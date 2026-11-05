import YAML

"""
    g2_entry_modes(path::AbstractString)
    g2_entry_modes(config::AbstractDict)

Inspect a simulation YAML file or its loaded configuration for scheduled G1 → G2
openings. Return a tuple containing enabled modes in the order `:regular`,
`:aligned`, `:anti_aligned`; `()` means no scheduled entry. Unrestricted
`OpenColocated` enables both colocated modes. Ignore zero-weight requests and
source sectors unreachable from Z through positive-weight requests.

Return `missing` with a warning if the file or required schedule fields cannot
be interpreted. Check schedule structure, not full simulation validity or physical
proposal feasibility/acceptance. The input and manifests are not modified.
"""
function g2_entry_modes(path::AbstractString)
    config = try
        YAML.load_file(path)
    catch err
        err isa InterruptException && rethrow()
        @warn "Cannot read configuration for G2 entry inspection" path exception=err
        return missing
    end
    if !(config isa AbstractDict)
        @warn "Cannot interpret configuration for G2 entry inspection: YAML root must be a mapping" path
        return missing
    end
    return g2_entry_modes(config)
end

function g2_entry_modes(config)
    try
        return inspect_g2_entry_modes(config)
    catch err
        err isa ArgumentError || rethrow()
        @warn "Cannot interpret configuration for G2 entry inspection" exception=err
        return missing
    end
end

"""Validate the required schedule fields and classify openings reachable from Z."""
function inspect_g2_entry_modes(config)
    config isa AbstractDict || throw(ArgumentError("config must be a mapping"))
    sampler = get(config, "sampler", nothing)
    sampler isa AbstractDict || throw(ArgumentError("sampler must be a mapping"))
    order = get(sampler, "multiworm_order", nothing)
    order isa Integer && !(order isa Bool) && order > 0 ||
        throw(ArgumentError("sampler.multiworm_order must be a positive integer"))
    moves = get(sampler, "moves", nothing)
    moves isa AbstractVector && !isempty(moves) ||
        throw(ArgumentError("sampler.moves must be a nonempty sequence"))
    tables = Dict{Tuple,Vector{Tuple{String,Tuple}}}()
    dimension = nothing
    function sector(value)
        value isa AbstractVector{<:Integer} && !(eltype(value) <: Bool) &&
            !isempty(value) && all(>=(0), value) ||
            throw(ArgumentError("source and target sectors must be nonnegative integer vectors"))
        dimension === nothing && (dimension = length(value))
        length(value) == dimension && sum(value) <= order ||
            throw(ArgumentError("sector dimensions or worm counts conflict with the sampler"))
        return Tuple(value)
    end
    opens = ("Open", "OpenColocated", "OpenColocatedAligned", "OpenColocatedAntiAligned")
    closes = ("Close", "CloseColocated", "CloseColocatedAligned", "CloseColocatedAntiAligned")
    preserving = ("ShiftTau", "ShiftTauColocated", "InsertKink", "InsertKinkColocated",
                  "DeleteKink", "DeleteKinkColocated")
    for entry in moves
        entry isa AbstractDict || throw(ArgumentError("each move table must be a mapping"))
        source = sector(get(entry, "source", nothing))
        haskey(tables, source) && throw(ArgumentError("duplicate source sector"))
        requests = get(entry, "requests", nothing)
        requests isa AbstractVector && !isempty(requests) ||
            throw(ArgumentError("requests must be a nonempty sequence"))
        active = Tuple{String,Tuple}[]
        for request in requests
            request isa AbstractDict || throw(ArgumentError("each request must be a mapping"))
            kind = get(request, "kind", nothing)
            kind isa AbstractString && kind in (opens..., closes..., preserving...) ||
                throw(ArgumentError("unknown or missing move kind"))
            target = sector(get(request, "target", nothing))
            weight = get(request, "weight", nothing)
            weight isa Real && !(weight isa Bool) && isfinite(weight) && weight >= 0 ||
                throw(ArgumentError("move weights must be finite and nonnegative"))
            delta = target .- source
            valid = kind in opens ? sum(delta) == 1 && all(>=(0), delta) :
                    kind in closes ? sum(delta) == -1 && all(<=(0), delta) : target == source
            valid || throw(ArgumentError("move kind conflicts with source and target sectors"))
            weight > 0 && push!(active, (String(kind), target))
        end
        isempty(active) && throw(ArgumentError("each source must have a positive-weight request"))
        tables[source] = active
    end
    z = ntuple(_ -> 0, dimension)
    haskey(tables, z) || throw(ArgumentError("schedule must contain the Z source sector"))
    all(haskey(tables, target) for requests in values(tables) for (_, target) in requests) ||
        throw(ArgumentError("schedule is missing a positive-weight target's source table"))
    pending, seen = [z], Set([z])
    enabled = Set{Symbol}()
    for source in pending
        for (kind, target) in tables[source]
            if sum(source) == 1 && sum(target) == 2 && kind in opens
                kind == "Open" && push!(enabled, :regular)
                kind in ("OpenColocated", "OpenColocatedAligned") && push!(enabled, :aligned)
                kind in ("OpenColocated", "OpenColocatedAntiAligned") && push!(enabled, :anti_aligned)
            end
            if !(target in seen)
                push!(seen, target)
                push!(pending, target)
            end
        end
    end
    return Tuple(mode for mode in (:regular, :aligned, :anti_aligned) if mode in enabled)
end
