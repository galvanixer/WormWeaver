"""
    green_ratio_statistics(numerator_history, z_counts; reduction=:sum,
        sample_counts=nothing, transform=identity, min_blocks=10,
        blocking_diagnostics=false)

Analyze paired numerator/Z histories with blocks on the final axis. `:mean`
reconstructs numerator sums using `sample_counts`; `:sum` leaves sums unweighted.
Apply `transform` to the full ratio and each paired jackknife replicate. Return
`measurements`, `error_status`, and optional `blocking_diagnostics`.

The value uses every original block without bias correction. Errors are the
maximum valid jackknife errors over repeated blocking levels. Complex results
carry separate real/imaginary Measurement errors; use `transform` for joint
derived quantities because component covariance is not retained. Zero total Z
counts throw; unavailable errors are NaN. Diagnostics do not certify convergence
or include time-window bias. See the README for normalization and examples.
"""
function green_ratio_statistics(numerator_history::AbstractArray{<:Number},
                                z_counts::AbstractVector{<:Real};
                                reduction::Symbol=:sum, sample_counts=nothing,
                                transform=identity, min_blocks::Integer=10,
                                blocking_diagnostics::Bool=false)
    min_blocks >= 2 || throw(ArgumentError("min_blocks must be at least two"))
    histories = prepare_ratio_histories(numerator_history, z_counts; reduction, sample_counts)
    numerator_totals = vec(sum(histories.numerator_sums; dims=2))
    total_z = sum(histories.z_counts)
    central = transformed_ratio(numerator_totals, total_z, histories.input_shape, transform)
    errors = blocked_jackknife_errors(histories, central, transform; min_blocks, blocking_diagnostics)
    measurements = ratio_measurements(central, errors.real_error, errors.imaginary_error)
    return (; measurements, error_status=errors.error_status,
            blocking_diagnostics=errors.blocking_diagnostics)
end

"""Validate paired histories and convert saved numerator means to block sums."""
function prepare_ratio_histories(numerator_history, z_counts; reduction, sample_counts)
    reduction in (:sum, :mean) || throw(ArgumentError("reduction must be :sum or :mean"))
    ndims(numerator_history) >= 1 || throw(DimensionMismatch("numerator must have a block axis"))
    block_count = length(z_counts)
    block_count > 0 || throw(ArgumentError("Green history is empty"))
    size(numerator_history, ndims(numerator_history)) == block_count ||
        throw(DimensionMismatch("numerator and Z histories must have matching blocks"))
    !(eltype(numerator_history) <: Bool) && all(isfinite, numerator_history) ||
        throw(ArgumentError("numerator history must contain finite numbers"))
    !(eltype(z_counts) <: Bool) && all(z -> isfinite(z) && z >= 0, z_counts) ||
        throw(ArgumentError("Z counts must be finite and nonnegative"))

    samples = nothing
    if sample_counts !== nothing
        sample_counts isa AbstractVector{<:Real} && !(eltype(sample_counts) <: Bool) ||
            throw(ArgumentError("sample counts must be a numeric vector"))
        length(sample_counts) == block_count || throw(DimensionMismatch("sample counts must match blocks"))
        all(c -> isfinite(c) && c > 0 && isinteger(c), sample_counts) ||
            throw(ArgumentError("sample counts must be positive finite integers"))
        samples = Int64.(sample_counts)
    elseif reduction == :mean
        throw(ArgumentError("mean numerators require sample_counts"))
    end

    input_shape = size(numerator_history)[1:end-1]
    number_type = eltype(numerator_history) <: Real ? Float64 : ComplexF64
    numerator_sums = reshape(number_type.(numerator_history), prod(input_shape), block_count)
    if reduction == :mean
        numerator_sums = numerator_sums .* reshape(samples, 1, block_count)
    end
    all(isfinite, numerator_sums) || throw(ArgumentError("numerator sum reconstruction overflowed"))
    z = Float64.(z_counts)
    total_z = sum(z)
    isfinite(total_z) && total_z > 0 || throw(ArgumentError("require finite positive total Z visits"))
    return (; numerator_sums, z_counts=z, sample_counts=samples, input_shape)
end

"""Take the maximum usable jackknife error while merging adjacent blocks."""
function blocked_jackknife_errors(histories, central, transform; min_blocks, blocking_diagnostics)
    component_count = length(central.values)
    real_error = fill(NaN, component_count)
    imaginary_error = central.iscomplex ? fill(NaN, component_count) : nothing
    diagnostics = nothing
    if blocking_diagnostics
        diagnostics = create_ratio_diagnostics(central, histories.sample_counts)
    end
    block_width, valid_levels, invalid_levels = 1, 0, 0
    last_status = "fewer than two blocks"

    while true
        level = jackknife_standard_error(histories, central, transform)
        if level.status == "ok"
            if valid_levels == 0
                real_error = copy(level.real_error)
                if central.iscomplex
                    imaginary_error = copy(level.imaginary_error)
                end
            else
                real_error = max.(real_error, level.real_error)
                if central.iscomplex
                    imaginary_error = max.(imaginary_error, level.imaginary_error)
                end
            end
            valid_levels += 1
        else
            invalid_levels += 1
        end
        last_status = level.status
        if blocking_diagnostics
            record_ratio_level!(diagnostics, histories, level, block_width)
        end
        length(histories.z_counts) ÷ 2 >= min_blocks || break
        histories = merge_adjacent_blocks(histories)
        block_width *= 2
    end

    if valid_levels == 0
        error_status = last_status
    elseif invalid_levels == 0
        error_status = "ok"
    else
        error_status = "invalid blocking levels omitted"
    end
    diagnostic_results = nothing
    if blocking_diagnostics
        diagnostic_results = finish_ratio_diagnostics(diagnostics, central.shape)
    end
    return (; real_error, imaginary_error, error_status, blocking_diagnostics=diagnostic_results)
end

"""Compute one paired delete-one-block jackknife error using online variance."""
function jackknife_standard_error(histories, central, transform)
    z_counts = histories.z_counts
    block_count = length(z_counts)
    total_z = sum(z_counts)
    component_count = length(central.values)
    real_error = fill(NaN, component_count)
    imaginary_error = central.iscomplex ? fill(NaN, component_count) : nothing
    if block_count < 2
        return (; real_error, imaginary_error, status="fewer than two blocks")
    end
    if !(total_z > 0 && all(z -> total_z - z > 0, z_counts))
        return (; real_error, imaginary_error, status="zero leave-one-out Z count")
    end

    numerator_totals = vec(sum(histories.numerator_sums; dims=2))
    real_center, real_m2 = zeros(component_count), zeros(component_count)
    imaginary_center = central.iscomplex ? zeros(component_count) : nothing
    imaginary_m2 = central.iscomplex ? zeros(component_count) : nothing
    for block in 1:block_count
        remaining_numerator = numerator_totals .- histories.numerator_sums[:, block]
        remaining_z = total_z - z_counts[block]
        output = transformed_ratio(remaining_numerator, remaining_z, histories.input_shape, transform)
        output.shape == central.shape && output.iscomplex == central.iscomplex ||
            throw(DimensionMismatch("transform must preserve its output shape and real/complex type"))
        delta = real.(output.values) .- real_center
        real_center .+= delta ./ block
        real_m2 .+= delta .* (real.(output.values) .- real_center)
        if central.iscomplex
            imaginary_delta = imag.(output.values) .- imaginary_center
            imaginary_center .+= imaginary_delta ./ block
            imaginary_m2 .+= imaginary_delta .* (imag.(output.values) .- imaginary_center)
        end
    end
    jackknife_factor = (block_count - 1) / block_count
    real_error = sqrt.(jackknife_factor .* max.(real_m2, 0.0))
    if central.iscomplex
        imaginary_error = sqrt.(jackknife_factor .* max.(imaginary_m2, 0.0))
    end
    all(isfinite, real_error) && (!central.iscomplex || all(isfinite, imaginary_error)) ||
        throw(ArgumentError("jackknife variance overflowed"))
    return (; real_error, imaginary_error, status="ok")
end

"""Sum neighboring paired histories, omitting an unmatched trailing block."""
function merge_adjacent_blocks(histories)
    pair_count = length(histories.z_counts) ÷ 2
    left, right = 1:2:(2pair_count), 2:2:(2pair_count)
    numerator_sums = histories.numerator_sums[:, left] .+ histories.numerator_sums[:, right]
    z_counts = histories.z_counts[left] .+ histories.z_counts[right]
    sample_counts = nothing
    if histories.sample_counts !== nothing
        samples = histories.sample_counts
        sample_counts = [Base.Checked.checked_add(samples[2i-1], samples[2i]) for i in 1:pair_count]
    end
    return (; numerator_sums, z_counts, sample_counts, input_shape=histories.input_shape)
end

"""Shape a numerator/Z ratio and validate the result of the supplied transform."""
function transformed_ratio(numerator_totals, total_z, input_shape, transform)
    if isempty(input_shape)
        ratio = only(numerator_totals) / total_z
    else
        ratio = reshape(numerator_totals ./ total_z, input_shape)
    end
    return green_transform_output(transform(ratio))
end

"""Flatten finite transformed values while retaining their shape and complex type."""
function green_transform_output(output)
    output isa Union{Real,Complex,AbstractArray{<:Union{Real,Complex}}} ||
        throw(ArgumentError("transform must return a real/complex scalar or array"))
    scalar = output isa Number
    shape = scalar ? () : size(output)
    values = scalar ? [output] : vec(collect(output))
    all(isfinite, values) || throw(ArgumentError("transform returned nonfinite values"))
    iscomplex = scalar ? output isa Complex : eltype(output) <: Complex
    converted = iscomplex ? ComplexF64.(values) : Float64.(values)
    all(isfinite, converted) || throw(ArgumentError("transform output overflowed Float64"))
    return (; shape, values=converted, iscomplex)
end

"""Attach component errors to the full-data estimate and restore its output shape."""
function ratio_measurements(central, real_error, imaginary_error)
    values = Measurements.measurement.(real.(central.values), real_error)
    if central.iscomplex
        imaginary_values = Measurements.measurement.(imag.(central.values), imaginary_error)
        values = complex.(values, imaginary_values)
    end
    return isempty(central.shape) ? only(values) : reshape(values, central.shape)
end

"""Allocate optional level records separately from the statistical calculation."""
function create_ratio_diagnostics(central, sample_counts)
    data = Dict{String,Any}(
        "block_widths" => Int[], "block_counts" => Int[], "z_counts" => Float64[],
        "level_statuses" => String[],
    )
    if sample_counts !== nothing
        data["sample_counts"] = Int64[]
    end
    real_curves = Vector{Vector{Float64}}()
    imaginary_curves = central.iscomplex ? Vector{Vector{Float64}}() : nothing
    return (; data, real_curves, imaginary_curves)
end

"""Record one level's counts, status, and real/imaginary errors."""
function record_ratio_level!(diagnostics, histories, level, block_width)
    data = diagnostics.data
    push!(data["block_widths"], block_width)
    push!(data["block_counts"], length(histories.z_counts))
    push!(data["z_counts"], sum(histories.z_counts))
    push!(data["level_statuses"], level.status)
    if histories.sample_counts !== nothing
        push!(data["sample_counts"], checked_count_sum(histories.sample_counts))
    end
    push!(diagnostics.real_curves, level.real_error)
    if diagnostics.imaginary_curves !== nothing
        push!(diagnostics.imaginary_curves, level.imaginary_error)
    end
    return nothing
end

"""Restore component shapes with blocking level as the final diagnostic axis."""
function finish_ratio_diagnostics(diagnostics, output_shape)
    data = diagnostics.data
    curves = diagnostics.real_curves
    data["real_standard_errors"] = reshape(hcat(curves...), output_shape..., length(curves))
    if diagnostics.imaginary_curves !== nothing
        curves = diagnostics.imaginary_curves
        data["imaginary_standard_errors"] = reshape(hcat(curves...), output_shape..., length(curves))
    end
    return data
end
