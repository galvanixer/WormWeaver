# Read saved HDF5 observable data for aggregation. This file loads block histories,
# sample counts, and metadata, checking dataset shapes, reductions, and numeric
# contents. It selects required metadata based on the available measurements.
# Observable-specific validation, physical normalization, and statistical analysis
# are performed by the summary functions after reading.

import HDF5

"""Read histories, sample counts, and metadata; return missing for an absent group."""
function read_observable_data(file, observable, definition)
    group_path = "observables/$observable"
    if !haskey(file, group_path)
        @warn "$observable was not found; leaving the $observable cell missing" file=HDF5.filename(file)
        return missing
    end
    group = file[group_path]
    haskey(group, "sampleCount") || throw(ArgumentError("Missing $observable/sampleCount"))
    counts = read(group["sampleCount"])
    counts isa AbstractVector{<:Real} || throw(ArgumentError("$observable/sampleCount must be a numeric vector"))
    all(c -> isfinite(c) && c > 0 && isinteger(c), counts) ||
        throw(ArgumentError("$observable/sampleCount must contain positive finite integers"))

    optional_fields = filter(field -> haskey(group, field), definition.optional_fields)
    fields = (definition.fields..., optional_fields...)
    histories = Dict{String,Any}()
    for field in fields
        histories[field] = read_observable_history(file, group, observable, field, definition, length(counts))
    end
    required_metadata = required_observable_metadata(observable, definition, histories)
    metadata = read_observable_metadata(file, observable, definition.metadata, required_metadata)
    return (; histories, counts, metadata)
end

"""Read one field and check block shape, reduction, and finite numeric contents."""
function read_observable_history(file, group, observable, field, definition, block_count)
    haskey(group, field) || throw(ArgumentError("Missing $observable/$field"))
    dataset = group[field]
    ndims(dataset) >= 1 && size(dataset, ndims(dataset)) == block_count ||
        throw(DimensionMismatch("$observable/$field and sampleCount have different block counts"))
    attributes = HDF5.attributes(dataset)
    reduction_path = "observables/METADATA/$observable/reduction_$field"
    expected_reduction = definition.reduction
    if expected_reduction == "sum" && !haskey(attributes, "reduction") && !haskey(file, reduction_path)
        throw(ArgumentError("$observable/$field requires explicit sum reduction metadata"))
    end
    if haskey(attributes, "reduction")
        read(attributes["reduction"]) == expected_reduction ||
            throw(ArgumentError("$observable/$field requires $expected_reduction reduction"))
    end
    if haskey(file, reduction_path)
        read(file[reduction_path]) == expected_reduction ||
            throw(ArgumentError("$observable/$field requires $expected_reduction reduction"))
    end

    history = read(dataset)
    history isa AbstractArray{<:Real} || throw(ArgumentError("$observable/$field must be a numeric history"))
    all(isfinite, history) || throw(ArgumentError("$observable/$field contains nonfinite values"))
    if expected_reduction == "sum"
        if get(definition, :integer_sums, true)
            history isa AbstractArray{<:Integer} && !(eltype(history) <: Bool) && all(>=(0), history) ||
                throw(ArgumentError("$observable/$field must contain nonnegative integer counts"))
        else
            !(eltype(history) <: Bool) && all(>=(0), history) ||
                throw(ArgumentError("$observable/$field must contain nonnegative numeric counts"))
        end
    end
    return history
end

"""Select required metadata explicitly for the observable and available histories."""
function required_observable_metadata(observable, definition, histories)
    if observable == "SingleParticleGreenFunction"
        if haskey(histories, "gREqualTimeHits")
            return SINGLE_PARTICLE_GREEN_EQUAL_TIME_METADATA
        end
        return ()
    elseif definition.reduction == "sum"
        return definition.metadata
    end
    return ()
end

"""Retain available metadata and reject missing required fields."""
function read_observable_metadata(file, observable, fields, required_fields)
    metadata = Dict{String,Any}()
    for field in fields
        metadata_path = "observables/METADATA/$observable/$field"
        if haskey(file, metadata_path)
            metadata[field] = read(file[metadata_path])
        elseif field in required_fields
            throw(ArgumentError("Missing $observable/$field metadata"))
        end
    end
    return metadata
end
