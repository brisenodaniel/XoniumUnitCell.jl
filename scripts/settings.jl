# Read and write julia_settings.TOML for the Makefile.
#
#   julia scripts/settings.jl <file> init             create <file> with defaults if missing
#   julia scripts/settings.jl <file> get <key>        print a value (default if file/key missing)
#   julia scripts/settings.jl <file> set <key> <val>  validate and store a value
#   julia scripts/settings.jl <file> resolve <val>    print <val> with "auto" replaced by Sys.CPU_THREADS
#
# Keys: threads (positive integer or "auto"), kernel_name (non-empty string).

using TOML

const DEFAULTS = Dict{String,Any}("threads" => "auto", "kernel_name" => "QuantumSims")

function parse_threads(s::AbstractString)
    s == "auto" && return "auto"
    n = tryparse(Int, s)
    (n === nothing || n < 1) && error("threads must be a positive integer or \"auto\", got \"$s\"")
    return n
end

function validate(key, val::AbstractString)
    haskey(DEFAULTS, key) || error("unknown setting \"$key\" (expected one of: $(join(sort!(collect(keys(DEFAULTS))), ", ")))")
    key == "threads" && return parse_threads(val)
    isempty(strip(val)) && error("kernel_name must not be empty")
    return String(val)
end

load(file) = merge(DEFAULTS, isfile(file) ? TOML.parsefile(file) : Dict{String,Any}())

function save(file, settings)
    open(file, "w") do io
        println(io, "# Defaults used by `make kernel` and `make session`.")
        println(io, "# threads: a positive integer, or \"auto\" (= Sys.CPU_THREADS, the logical CPU count).")
        println(io, "# `make kernel THREADS=... KERNEL_NAME=...` rewrites this file; you can also edit it by hand.")
        TOML.print(io, Dict(k => settings[k] for k in keys(DEFAULTS)); sorted=true)
    end
end

resolve(v) = v == "auto" ? Sys.CPU_THREADS : parse_threads(string(v))

function main(args)
    length(args) >= 2 || error("usage: settings.jl <file> init|get|set|resolve ...")
    file, cmd = args[1], args[2]
    if cmd == "init"
        isfile(file) || save(file, DEFAULTS)
    elseif cmd == "get"
        print(load(file)[validate_key(args[3])])
    elseif cmd == "set"
        key = args[3]
        settings = load(file)
        settings[key] = validate(key, args[4])
        save(file, settings)
    elseif cmd == "resolve"
        print(resolve(args[3]))
    else
        error("unknown command \"$cmd\"")
    end
end

validate_key(key) = (haskey(DEFAULTS, key) || error("unknown setting \"$key\""); key)

main(ARGS)
