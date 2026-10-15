export PinkNoise, noise_amplitudes, S_pink_noise, default_pink_noise
export calc_gs, calc_Γφ

# k_B T / h in GHz at T = 15 mK (1.38e-23 J/K * 15e-3 K / 6.63e-34 J s = 0.312 GHz),
# times 2π so that it is an angular frequency (rad/ns), like the frequencies passed to the PSD
const KT_15mK = 2π * 1.38 * 1.5 / 6.63
const DEFAULT_TNDC = 1.1e-6 # dielectric loss tangent
const DEFAULT_DF = 1.8e-6   # 1/f flux noise amplitude

"""
    noise_amplitudes(El, Ec, φ_ge; tndC=$DEFAULT_TNDC, df=$DEFAULT_DF)

Amplitudes `(Ad, Af)` of the dielectric-loss and 1/f flux noise terms of [`S_pink_noise`](@ref)
for a qubit with inductive energy `El`, charging energy `Ec` and flux matrix element
`φ_ge = ⟨g|φ|e⟩`:

    Af = df * El * |φ_ge|
    Ad = π² * tndC * |φ_ge|² / Ec
"""
function noise_amplitudes(El::Real, Ec::Real, φ_ge::Number; tndC::Real=DEFAULT_TNDC, df::Real=DEFAULT_DF)
    Af = df * El * abs(φ_ge)
    Ad = π^2 * tndC * abs2(φ_ge) / Ec
    return Ad, Af
end

"""
    S_pink_noise(w, Ad, Af; kT=$(round(KT_15mK, digits=4)))

Power spectral density of the environmental noise at angular frequency `w`:

    S(w) = Af² * T + α * Ad / T²,    T = |2π/w|,    α = |coth(w/2kT) + 1| / 2

The first term is 1/f flux noise (`Af² · 2π/|w|`), the second a thermally weighted
dielectric-loss term, which is larger for `w > 0` than for `w < 0` through `α`.
`w` may be a scalar; use broadcasting (`S_pink_noise.(ws, Ad, Af)`) for arrays.

`w = 0` is not regularized (`T = Inf`, so the result is `Inf` or `NaN`), the DC component
must be treated separately by the caller.
"""
function S_pink_noise(w::Real, Ad::Real, Af::Real; kT::Real=KT_15mK)
    α = abs(coth(w / (2kT)) + 1) / 2
    T = abs(2π / w)
    return Af^2 * T + α * Ad / T^2
end

"""
    PinkNoise(Ad, Af; kT=$(round(KT_15mK, digits=4)))
    PinkNoise(; El, Ec, φ_ge, tndC, df, kT)

Callable noise model: `S = PinkNoise(...)`, then `S(w)` is `S_pink_noise(w, S.Ad, S.Af; kT=S.kT)`,
and `S.(ws)` works on arrays. The second constructor computes `(Ad, Af)` with
[`noise_amplitudes`](@ref).
"""
struct PinkNoise{T<:Real}
    Ad::T
    Af::T
    kT::T
end

PinkNoise(Ad::Real, Af::Real; kT::Real=KT_15mK) = PinkNoise(promote(Ad, Af, kT)...)

function PinkNoise(; El::Real, Ec::Real, φ_ge::Number, kT::Real=KT_15mK, kwargs...)
    Ad, Af = noise_amplitudes(El, Ec, φ_ge; kwargs...)
    return PinkNoise(Ad, Af; kT=kT)
end

(S::PinkNoise)(w::Real) = S_pink_noise(w, S.Ad, S.Af; kT=S.kT)

"""
    default_pink_noise(; EJ=4.5, EL=1.0, EC=1.0, kwargs...)

[`PinkNoise`](@ref) for a fluxonium at its flux sweet spot with parameters in GHz
(default: `EJ=4.5, EL=1.0, EC=1.0`), the noise model of the bichromatic-Floquet study.
`El` and `Ec` are converted to angular units (×2π) and `φ_ge` is taken between the two lowest
eigenstates of the fluxonium. `kwargs` are forwarded to [`noise_amplitudes`](@ref).
"""
function default_pink_noise(; EJ::Real=4.5, EL::Real=1.0, EC::Real=1.0, kT::Real=KT_15mK, kwargs...)
    p = NamedTuple(tounit_ħ((; EJ, EL, EC)))
    fl = fluxonium(; p...)
    _, g = state(fl, :g)
    _, e = state(fl, :e)
    φ_ge = g' * fl.φ * e
    return PinkNoise(; El=p.EL, Ec=p.EC, φ_ge=φ_ge, kT=kT, kwargs...)
end


## Dephasing rate #######################

# Harmonic indices of an `n`-sample DFT, in the order of `numpy.fft.fftfreq(n) * n`
# (0, 1, …, -n÷2, …, -1), or of its `fftshift` (-n÷2, …, 0, …) if `centered`.
_harmonics(n::Int, centered::Bool) = centered ? collect(-(n ÷ 2):((n - 1) ÷ 2)) :
                                                vcat(0:((n - 1) ÷ 2), -(n ÷ 2):-1)

"""
    calc_gs(φ_upper, φ_lower)

Harmonic weights `gₖ` of the flux coupling of a Floquet doublet. `φ_upper` and `φ_lower` are the
expectation values `⟨φ⟩(tₙ)` of the flux operator in the upper and lower Floquet state (or, in the
adiabatic limit, the instantaneous eigenstates they follow), sampled at `N` equally spaced times
`tₙ = nT/N`, `n = 0, …, N-1`, over one drive period (the endpoint `T` must **not** be included).

Returns the `N` weights of `g(t) = (⟨φ⟩_upper - ⟨φ⟩_lower) / 2`, with the convention of
`numpy.fft.ifft`: `gₖ = (1/N) Σₙ g(tₙ) e^{+2πikn/N}`, ordered by harmonic index `0, 1, …, -1`.
`g₀` is the period average of `g`.
"""
function calc_gs(φ_upper::AbstractVector{<:Number}, φ_lower::AbstractVector{<:Number})
    length(φ_upper) == length(φ_lower) ||
        throw(DimensionMismatch("`φ_upper` and `φ_lower` must have the same length"))
    N = length(φ_upper)
    g = (φ_upper .- φ_lower) ./ 2
    return [sum(g[n + 1] * cis(2π * k * n / N) for n in 0:(N - 1)) / N for k in _harmonics(N, false)]
end

"""
    calc_Γφ(gs, S, Af, w_base; g0_centered=false)
    calc_Γφ(gs, S::PinkNoise, w_base; g0_centered=false)
    calc_Γφ(φ_upper, φ_lower, S::PinkNoise, w_base)

Dephasing rate of a driven qubit from the environmental noise `S`:

    Γφ = 8 Af |g₀| + Σ_{k≠0} |gₖ|² S(k w_base)

The first term is the dephasing from the zero-frequency (1/f) noise and `Af` its amplitude. The sum
runs over all non-zero harmonics `k` of the drive and uses `S(w)` at the angular frequency `k w_base`,
with `w_base = 2π/T` the angular frequency of the drive period.

# Arguments
- `gs`: harmonic weights, as returned by [`calc_gs`](@ref), ordered `0, 1, …, -1` (`fftfreq` order), or
        `-N÷2, …, N÷2-1` if `g0_centered`.
- `S`: noise PSD, a function `S(w)` of angular frequency that accepts a scalar or an array, e.g.
        a [`PinkNoise`](@ref).
- `Af`: amplitude of the zero-frequency noise. Taken from `S.Af` if `S` is a `PinkNoise`.
- `w_base`: angular frequency of the drive period.

The third method computes `gs` from the samples `⟨φ⟩(tₙ)` of the upper and lower Floquet state with
[`calc_gs`](@ref).
"""
function calc_Γφ(gs::AbstractVector{<:Number}, S, Af::Real, w_base::Real; g0_centered::Bool=false)
    harmonics = _harmonics(length(gs), g0_centered)
    zero_idx = findfirst(iszero, harmonics)
    g_0 = gs[zero_idx]
    nonzero = eachindex(gs) .!= zero_idx
    g_k = gs[nonzero]
    S_k = S.(harmonics[nonzero] .* w_base)
    rate_0 = 8 * Af * abs(g_0)
    rate_k = sum(abs2.(g_k) .* S_k)
    return rate_0 + rate_k
end

calc_Γφ(gs::AbstractVector{<:Number}, S::PinkNoise, w_base::Real; kwargs...) =
    calc_Γφ(gs, S, S.Af, w_base; kwargs...)

calc_Γφ(φ_upper::AbstractVector{<:Number}, φ_lower::AbstractVector{<:Number}, S::PinkNoise, w_base::Real) =
    calc_Γφ(calc_gs(φ_upper, φ_lower), S, w_base)
