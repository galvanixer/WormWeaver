# Winding observable definition and summary construction. Keep winding-specific
# analysis here while using the shared reader and blocking statistics.

const WINDINGS_DEFINITION = (
    fields=("windingBySpecies", "squaredWindingBySpecies", "totalWinding", "squaredTotalWinding"),
    optional_fields=("windingCorrelationBySpecies",),
    metadata=("periodic_axes", "species_pairs"), reduction="mean",
)

function windings_summary(histories, counts, metadata; min_blocks, blocking_diagnostics::Bool=false)
    squared = histories["squaredWindingBySpecies"]
    ndims(squared) == 3 || throw(DimensionMismatch(
        "Windings/squaredWindingBySpecies must have species × spatial axis × block dimensions"))
    nspecies, naxes, nblocks = size(squared)
    nspecies >= 1 && naxes >= 1 || throw(DimensionMismatch(
        "Windings must have at least one species and spatial axis"))
    nblocks == length(counts) || throw(DimensionMismatch("Windings block counts do not match"))
    size(histories["windingBySpecies"]) == size(squared) || throw(DimensionMismatch(
        "Windings species moments must have matching dimensions"))
    for field in ("totalWinding", "squaredTotalWinding")
        size(histories[field]) == (naxes, nblocks) || throw(DimensionMismatch(
            "Windings/$field must have spatial axis × block dimensions matching the species moments"))
    end
    if haskey(metadata, "periodic_axes")
        mask = metadata["periodic_axes"]
        mask isa AbstractVector{Bool} && length(mask) == naxes || throw(DimensionMismatch(
            "Windings/periodic_axes must be a Boolean vector matching the spatial axes"))
    end

    products = get(histories, "windingCorrelationBySpecies", nothing)
    pairs = get(metadata, "species_pairs", nothing)
    if products !== nothing
        ndims(products) == 3 && size(products)[2:3] == (naxes, nblocks) ||
            throw(DimensionMismatch(
                "Windings/windingCorrelationBySpecies must have species pair × spatial axis × block dimensions"))
        size(products, 1) <= nspecies * (nspecies - 1) ÷ 2 || throw(DimensionMismatch(
            "Windings/windingCorrelationBySpecies has too many species pairs"))
    end
    if pairs !== nothing
        validate_winding_pairs(pairs, nspecies)
        products === nothing || size(products, 1) == size(pairs, 1) ||
            throw(DimensionMismatch("Windings/species_pairs and windingCorrelationBySpecies have different pair counts"))
    end

    # Work on copies: the raw saved moments and metadata remain unchanged.
    derived = copy(histories)
    retained_metadata = copy(metadata)
    if nspecies == 1
        retained_metadata["species_pairs"] = zeros(Int, 0, 2)
        coflow = zeros(Float64, 0, naxes, nblocks)
        counterflow = similar(coflow)
    elseif products !== nothing
        if pairs === nothing
            @warn "Windings pair products have no species_pairs metadata; omitting pair coflow and counterflow"
            return mean_observable_summary(histories, counts, metadata; min_blocks, blocking_diagnostics)
        end
        coflow, counterflow = winding_pair_histories(squared, products, pairs)
    elseif nspecies == 2
        # The total square determines the sole missing pair product exactly,
        # even when the inputs are block means rather than individual samples.
        pairs = reshape([1, 2], 1, 2)
        products = reshape((histories["squaredTotalWinding"] .-
                            squared[1, :, :] .- squared[2, :, :]) ./ 2,
                           1, naxes, nblocks)
        retained_metadata["species_pairs"] = pairs
        coflow, counterflow = winding_pair_histories(squared, products, pairs)
    else
        @warn "Windings with three or more species require pair products to resolve individual pairs; omitting pair coflow and counterflow"
        return mean_observable_summary(histories, counts, metadata; min_blocks, blocking_diagnostics)
    end
    derived["squaredCoflowBySpeciesPair"] = coflow
    derived["squaredCounterflowBySpeciesPair"] = counterflow
    return mean_observable_summary(derived, counts, retained_metadata; min_blocks, blocking_diagnostics)
end

function validate_winding_pairs(pairs, nspecies)
    pairs isa AbstractMatrix{<:Integer} && !(eltype(pairs) <: Bool) ||
        throw(ArgumentError("Windings/species_pairs must be an integer matrix"))
    size(pairs, 2) == 2 || throw(DimensionMismatch("Windings/species_pairs must have two columns"))
    seen = Set{Tuple{Int,Int}}()
    for row in eachrow(pairs)
        a, b = row
        1 <= a < b <= nspecies || throw(ArgumentError(
            "Windings/species_pairs must contain species indices 1 <= a < b <= $nspecies"))
        (a, b) in seen && throw(ArgumentError("Windings/species_pairs contains a duplicate pair"))
        push!(seen, (a, b))
    end
    return nothing
end

"""Build pair-channel block means before blocking, preserving joint fluctuations."""
function winding_pair_histories(squared, products, pairs)
    coflow = zeros(Float64, size(products))
    counterflow = similar(coflow)
    for (p, row) in enumerate(eachrow(pairs))
        a, b = row
        @views begin
            diagonal = squared[a, :, :] .+ squared[b, :, :]
            coflow[p, :, :] .= diagonal .+ 2 .* products[p, :, :]
            counterflow[p, :, :] .= diagonal .- 2 .* products[p, :, :]
        end
    end
    return coflow, counterflow
end
