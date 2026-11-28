# Read saved HDF5 observable data for aggregation. This file loads block histories,
# sample counts, and metadata, checking dataset shapes, reductions, and numeric
# contents. It selects required metadata based on the available measurements.
# Observable-specific validation, physical normalization, and statistical analysis
# are performed by the summary functions after reading.

import HDF5

"""Read histories, sample counts, metadata, and the G(k,τ) reduction when present."""
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
    gktau_reduction = nothing
    for field in fields
        history, reduction = read_observable_history(file, group, observable, field, definition, length(counts))
        history === nothing && continue
        histories[field] = history
        if observable == "SingleParticleGreenFunction" && field == "gKTau"
            gktau_reduction = reduction
        end
    end
    if observable == "CounterflowGreenFunction"
        omitted = filter(field -> haskey(group, field) && !haskey(histories, field),
                         ("gRTau", "gKTau", "gKOmega", "gROmega"))
        isempty(omitted) || @warn "CounterflowGreenFunction representations omitted" fields=omitted file=HDF5.filename(file)
    end
    required_metadata = required_observable_metadata(observable, definition, histories)
    metadata = read_observable_metadata(file, observable, definition.metadata, required_metadata)
    return (; histories, counts, metadata, gktau_reduction)
end

"""
Read a history and its resolved reduction, checking block shape and values.
Legacy counterflow Fourier means return `nothing` as the history after validating
the reduction metadata.
"""
function read_observable_history(file, group, observable, field, definition, block_count)
    haskey(group, field) || throw(ArgumentError("Missing $observable/$field"))
    dataset = group[field]
    attributes = HDF5.attributes(dataset)
    reduction_path = "observables/METADATA/$observable/reduction_$field"
    has_group_reduction = haskey(file, reduction_path)
    has_attribute = haskey(attributes, "reduction")
    group_reduction = has_group_reduction ? read(file[reduction_path]) : nothing
    reduction = definition.reduction
    is_gktau = field == "gKTau"
    if is_gktau
        # The writer uses group metadata to distinguish Fourier sums from legacy means.
        reduction = has_group_reduction ? group_reduction : "mean"
        reduction isa AbstractString && reduction in ("mean", "sum") ||
            throw(ArgumentError("$observable/$field requires mean or sum reduction"))
    end
    if reduction == "sum" && !has_attribute && !has_group_reduction
        throw(ArgumentError("$observable/$field requires explicit sum reduction metadata"))
    end
    if has_attribute
        attribute_reduction = read(attributes["reduction"])
        if is_gktau && !has_group_reduction && attribute_reduction == "sum"
            throw(ArgumentError("$observable/$field Fourier sums require $reduction_path"))
        end
        attribute_reduction == reduction ||
            throw(ArgumentError("$observable/$field requires $reduction reduction"))
    end
    if has_group_reduction
        group_reduction == reduction ||
            throw(ArgumentError("$observable/$field requires $reduction reduction"))
    end
    if is_gktau && observable == "CounterflowGreenFunction" && reduction == "mean"
        return nothing, reduction
    end
    ndims(dataset) >= 1 && size(dataset, ndims(dataset)) == block_count ||
        throw(DimensionMismatch("$observable/$field and sampleCount have different block counts"))

    history = read(dataset)
    numeric_type = field in get(definition, :complex_fields, ()) ? Number : Real
    history isa AbstractArray{<:numeric_type} || throw(ArgumentError("$observable/$field must be a numeric history"))
    all(isfinite, history) || throw(ArgumentError("$observable/$field contains nonfinite values"))
    if is_gktau && reduction == "sum"
        history isa AbstractArray{<:Complex} ||
            throw(ArgumentError("$observable/$field Fourier sums must be complex"))
    elseif reduction == "sum"
        if get(definition, :integer_sums, true)
            history isa AbstractArray{<:Integer} && !(eltype(history) <: Bool) && all(>=(0), history) ||
                throw(ArgumentError("$observable/$field must contain nonnegative integer counts"))
        else
            !(eltype(history) <: Bool) && all(>=(0), history) ||
                throw(ArgumentError("$observable/$field must contain nonnegative numeric counts"))
        end
    end
    return history, reduction
end

"""Select required metadata explicitly for the observable and available histories."""
function required_observable_metadata(observable, definition, histories)
    if observable == "CounterflowGreenFunction"
        required = ("species_pair", "channel")
        if haskey(histories, "gKTau")
            required = (required..., "momenta", "tau_edges")
        end
        if haskey(histories, "gREqualTimeHits")
            required = (required..., COUNTERFLOW_EQUAL_TIME_METADATA...)
        end
        return required
    elseif observable == "SingleParticleGreenFunction"
        required = String[]
        if haskey(histories, "gREqualTimeHits")
            append!(required, SINGLE_PARTICLE_GREEN_EQUAL_TIME_METADATA)
        end
        if haskey(histories, "gKTau")
            append!(required, ("species", "momenta", "tau_edges"))
        end
        return unique(required)
    elseif observable in ("DensityCorrelation", "StaticStructureFactor") || definition.reduction == "sum"
        return definition.metadata
    end
    return ()
end

"""Check the momentum and imaginary-time metadata against their history axes."""
function validate_fourier_grid(metadata, grid_shape, observable)
    momenta, edges = metadata["momenta"], metadata["tau_edges"]
    momenta isa AbstractMatrix{<:Real} && size(momenta, 2) > 0 && all(isfinite, momenta) ||
        throw(ArgumentError("$observable/momenta must be a finite numeric matrix"))
    edges isa AbstractVector{<:Real} && all(isfinite, edges) && all(>(0), diff(edges)) ||
        throw(ArgumentError("$observable/tau_edges must be finite and strictly increasing"))
    grid_shape == (size(momenta, 1), length(edges) - 1) ||
        throw(DimensionMismatch("$observable/gKTau axes do not match metadata"))
    if haskey(metadata, "momentum_indices")
        size(metadata["momentum_indices"]) == size(momenta) || throw(DimensionMismatch(
            "$observable/momentum_indices and momenta must have matching shapes"))
    end
    return nothing
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
