"""
    observable_table(summaries; fields=nothing) -> DataFrame

Expand an observable column such as `table.Windings` into a DataFrame with one
row per input entry and one column per dictionary key. Inputs must be dictionaries
with string keys or `missing`. Missing entries and absent keys produce missing
cells. Measurements, arrays, and metadata are retained as stored; arrays remain
individual cells rather than being expanded into component columns.

By default, use the sorted union of keys across all dictionaries. Supply a
nonempty vector or tuple of unique string `fields` to select columns and their
order. Nonempty inputs without any dictionary keys (including all-missing
columns) require explicit fields so their rows can be preserved. Empty inputs
return an empty DataFrame, with the requested columns if fields are supplied.

The input dictionaries are not modified. Array-valued cells reference the
original arrays, as in a normal DataFrame column.

# Examples
```julia
windings = observable_table(table.Windings)
totals = observable_table(table.Windings; fields=["totalWinding", "squaredTotalWinding"])
```
"""
function observable_table(summaries::AbstractVector; fields=nothing)
    field_names = Set{String}()
    for (i, summary) in enumerate(summaries)
        ismissing(summary) && continue
        summary isa AbstractDict || throw(ArgumentError("entry $i must be a dictionary or missing"))
        all(key -> key isa AbstractString, keys(summary)) ||
            throw(ArgumentError("entry $i must have string keys"))
        union!(field_names, keys(summary))
    end

    if fields === nothing
        selected = sort!(collect(field_names))
        isempty(selected) && !isempty(summaries) && throw(ArgumentError(
            "cannot infer columns from summaries without keys; supply fields"))
    else
        fields isa Union{AbstractVector,Tuple} ||
            throw(ArgumentError("fields must be a vector or tuple of string names"))
        all(field -> field isa AbstractString, fields) ||
            throw(ArgumentError("fields must contain string names"))
        selected = String.(collect(fields))
        isempty(selected) && throw(ArgumentError("fields must not be empty"))
        allunique(selected) || throw(ArgumentError("fields must not contain duplicates"))
    end

    result = DataFrames.DataFrame()
    for field in selected
        result[!, field] = [ismissing(summary) ? missing : get(summary, field, missing)
                           for summary in summaries]
    end
    return result
end

"""
    describe_density(summary; io=stdout) -> nothing

Print one simulation's aggregated Density values, units, sample/block counts,
and saved uncertainty method. Pass a Density dictionary such as `table.Density[1]`,
or `missing` to report that no measurements are available. The input is not modified.

Species IDs label the rows when available; otherwise use component numbers and
state that species IDs are unavailable. Nonfinite uncertainties are displayed as
"uncertainty unavailable". Values must be Measurements, with matching species
vectors for density and particle number. No estimates are recomputed.
"""
function describe_density(summary::Union{Missing,AbstractDict}; io::IO=stdout)
    if ismissing(summary)
        println(io, "No Density measurements available.")
        return nothing
    end
    density, particles = summary["densityBySpecies"], summary["particleNumberBySpecies"]
    density isa AbstractVector && particles isa AbstractVector || throw(ArgumentError(
        "Density and particle number must be species vectors"))
    length(density) == length(particles) || throw(DimensionMismatch(
        "Density and particle number must have matching species counts"))
    species = get(summary, "species", nothing)
    has_species = species !== nothing && !ismissing(species)
    if has_species
        species isa AbstractVector{<:Integer} && !(eltype(species) <: Bool) &&
            all(>(0), species) && allunique(species) || throw(ArgumentError(
            "species must contain positive unique integer IDs"))
        length(species) == length(density) || throw(DimensionMismatch(
            "species IDs must match the Density species count"))
    end
    function format_estimate(value)
        value isa Measurements.Measurement || throw(ArgumentError(
            "Density estimates must be Measurements"))
        return isfinite(Measurements.uncertainty(value)) ? string(value) :
            "$(Measurements.value(value)) (uncertainty unavailable)"
    end
    headers = ["Species", "Density (particles/site)", "Particle number (particles)"]
    rows = [[has_species ? string(species[i]) : "Component $i",
             format_estimate(density[i]), format_estimate(particles[i])]
            for i in eachindex(density)]
    total = format_estimate(summary["totalParticleNumber"])
    widths = [maximum(textwidth, [headers[j]; [row[j] for row in rows]]) for j in 1:3]

    println(io, "Density\n")
    for row in [headers, rows...]
        println(io, rstrip(join(rpad.(row, widths), "  ")))
    end
    has_species || println(io, "Species IDs unavailable; rows are labeled by component.")
    println(io, "\nTotal particle number: ", total, " particles")
    println(io, "\nSaved production blocks: ", summary["nblocks"])
    println(io, "Density samples: ", summary["nsamples"])
    println(io, "Uncertainty: ", summary["error_method"])
    return nothing
end
