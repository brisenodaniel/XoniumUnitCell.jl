export cosine_ramp_envelope, tanh_ramp_envelope, full_pulse

"""
Construct a flat-top pulse envelope with raised-cosine ramps.

Returns a function `f(t)` that
  - ramps 0 → 1 on [0, t_ramp],
  - is 1 on [t_ramp, t_ramp + t_plateau],
  - ramps 1 → 0 on [t_ramp + t_plateau, 2*t_ramp + t_plateau],
  - is 0 outside [0, 2*t_ramp + t_plateau].

Each ramp is a half-period raised cosine (Hann window): g(τ) = (1 - cos(πτ/t_ramp))/2
for τ ∈ [0, t_ramp]. Since g'(τ) = (π/2t_ramp)*sin(πτ/t_ramp) and sin(0) = sin(π) = 0,
the envelope has zero first derivative at every ramp/plateau/zero boundary, i.e. it is C¹.
"""
function cosine_ramp_envelope(t_ramp::Real, t_plateau::Real)
    t_down_start = t_ramp + t_plateau
    t_end = 2 * t_ramp + t_plateau

    return function(t::Real)
        if t < 0 || t > t_end
            return 0.0
        elseif t < t_ramp
            return (1 - cos(π * t / t_ramp)) / 2
        elseif t < t_down_start
            return 1.0
        else
            return (1 + cos(π * (t - t_down_start) / t_ramp)) / 2
        end
    end
end

"""
Flat-top envelope with tanh edges, as used in Fig. 2(c) of the paper
(characteristic rise time t_r = 2 ns there).

Same footprint as `cosine_ramp_envelope`: total length T = 2*t_ramp + t_plateau,
edges centered at t_ramp/2 and T - t_ramp/2 (so the half-amplitude points, and
hence the pulse area, match the cosine envelope). Each edge is
tanh((t - t_ramp/2)/t_rise); the default t_rise = t_ramp/6 puts the ends of each
ramp at ±3 t_rise, i.e. 99.5% of full swing. The raw sum of the two tanh edges is
shifted and rescaled so the envelope is exactly 0 at t = 0, T and exactly 1 at T/2,
avoiding a small flux step at the pulse boundaries.
"""
function tanh_ramp_envelope(t_ramp::Real, t_plateau::Real; t_rise::Real=t_ramp/6)
    t_end = 2 * t_ramp + t_plateau
    t_on, t_off = t_ramp / 2, t_end - t_ramp / 2
    raw(t) = (tanh((t - t_on) / t_rise) - tanh((t - t_off) / t_rise)) / 2
    r0, r_mid = raw(0.0), raw(t_end / 2)
    return function(t::Real)
        (t < 0 || t > t_end) && return 0.0
        return (raw(t) - r0) / (r_mid - r0)
    end
end

"""
Assemble a full carrier pulse: a flat-top envelope modulating a cosine carrier.

Returns a function `f(t) = envelope(t) * cos(ω * t + phase)`, where the envelope is
`cosine_ramp_envelope(t_ramp, t_plateau)` for `ramp_profile = :cosine` (default) or
`tanh_ramp_envelope(t_ramp, t_plateau; t_rise)` for `ramp_profile = :tanh`. `t_rise` is the
tanh characteristic rise time (default t_ramp/6) and is ignored for `:cosine`.
To match the paper's t_r = 2 ns, use `ramp_profile=:tanh, t_ramp=12.0` (or pass t_rise=2.0).

`ω` is an angular frequency (rad per unit time), matching the convention
used elsewhere in this project (EC, EL, EJ, J are all defined with an
explicit `*2π` at their point of definition, not inside the functions
that consume them). `phase` is in radians, default 0.
"""
function full_pulse(
    ω::Real,
    t_ramp::Real,
    t_plateau::Real;
    phase::Real=0.0,
    ramp_profile::Symbol=:cosine,
    t_rise::Maybe{Real}=nothing,
)
    envelope = if ramp_profile === :cosine
        cosine_ramp_envelope(t_ramp, t_plateau)
    elseif ramp_profile === :tanh
        isnothing(t_rise) ? tanh_ramp_envelope(t_ramp, t_plateau) :
                            tanh_ramp_envelope(t_ramp, t_plateau; t_rise=t_rise)
    else
        throw(ArgumentError("`ramp_profile` must be :cosine or :tanh, got :$ramp_profile"))
    end
    return function(t::Real)
        return envelope(t) * cos(ω * t + phase)
    end
end
