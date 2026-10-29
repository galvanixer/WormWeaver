import Measurements

weighted_mean(values, counts) = sum(values .* counts) / sum(counts)

checked_count_sum(values) = foldl(Base.Checked.checked_add, (Int64(value) for value in values); init=Int64(0))

function blocking_sample_counts(counts, block_widths, block_counts)
    return [eltype(counts) <: Integer ? checked_count_sum(@view counts[1:(width * nblocks)]) :
            sum(@view counts[1:(width * nblocks)])
            for (width, nblocks) in zip(block_widths, block_counts)]
end

# Standard error across blocks. With equal counts this is std(values)/sqrt(B).
# With unequal counts, use the ratio-of-sums delta-method error across blocks.
function block_standard_error(values, counts)
    nb = length(values)
    nb >= 2 || return NaN
    avg = weighted_mean(values, counts)
    return sqrt(nb / (nb - 1) * sum(abs2, counts .* (values .- avg))) / sum(counts)
end

"""
    blocking_measurement(values, counts; min_blocks=10)

Combine saved block means into a sample-count-weighted mean ± repeated-blocking
standard error. Applicable to any scalar block-mean history, including individual
components of array observables. Raw block sums require their own normalization.
"""
function blocking_measurement(values::AbstractVector{<:Real}, counts::AbstractVector{<:Real};
                              min_blocks::Integer=10)
    return blocking_analysis(values, counts; min_blocks).measurement
end

"""Return the blocking measurement; construct per-level curves only with diagnostics=true."""
function blocking_analysis(values::AbstractVector{<:Real}, counts::AbstractVector{<:Real};
                           min_blocks::Integer=10, diagnostics::Bool=false)
    min_blocks >= 2 || throw(ArgumentError("min_blocks must be at least two"))
    length(values) == length(counts) || throw(DimensionMismatch("block values and counts must match"))
    isempty(values) && throw(ArgumentError("block history is empty"))
    all(isfinite, values) || throw(ArgumentError("block history contains nonfinite values"))
    all(c -> isfinite(c) && c > 0 && isinteger(c), counts) ||
        throw(ArgumentError("sample counts must be positive finite integers"))

    blocked = Float64.(values)
    weights = Float64.(counts)
    avg = weighted_mean(blocked, weights)
    stderr = block_standard_error(blocked, weights)
    standard_errors = diagnostics ? [stderr] : nothing
    block_widths = diagnostics ? [1] : nothing
    block_counts = diagnostics ? [length(blocked)] : nothing
    while length(blocked) ÷ 2 >= min_blocks
        # An odd trailing block is omitted only at this error-estimation level.
        # The reported mean above always includes all original saved blocks.
        npairs = length(blocked) ÷ 2
        paired_counts = [weights[2i - 1] + weights[2i] for i in 1:npairs]
        blocked = [(blocked[2i - 1] * weights[2i - 1] + blocked[2i] * weights[2i]) /
                   paired_counts[i] for i in 1:npairs]
        weights = paired_counts
        level_error = block_standard_error(blocked, weights)
        stderr = max(stderr, level_error)
        if diagnostics
            push!(standard_errors, level_error)
            push!(block_widths, 2 * last(block_widths))
            push!(block_counts, length(blocked))
        end
    end
    return (; measurement=Measurements.measurement(avg, stderr),
            standard_errors, block_widths, block_counts)
end

"""Aggregate scalar or array block means while retaining per-component blocking curves."""
function observable_blocking_statistics(history::AbstractArray{<:Real}, counts; min_blocks)
    shape = size(history)[1:end-1]
    size(history, ndims(history)) == length(counts) ||
        throw(DimensionMismatch("observable history must have one entry per saved block"))
    matrix = reshape(history, prod(shape), length(counts))
    analyses = [blocking_analysis(collect(row), counts; min_blocks, diagnostics=true)
                for row in eachrow(matrix)]
    # Empty component axes still have a valid blocking grid determined by counts.
    grid = isempty(analyses) ? blocking_analysis(zeros(length(counts)), counts; min_blocks, diagnostics=true) : first(analyses)
    measurements = Measurements.Measurement{Float64}[analysis.measurement for analysis in analyses]
    levels = length(grid.standard_errors)
    errors = Float64[analysis.standard_errors[level] for analysis in analyses, level in 1:levels]
    return (; measurements=isempty(shape) ? only(measurements) : reshape(measurements, shape),
            standard_errors=reshape(errors, shape..., levels),
            block_widths=grid.block_widths, block_counts=grid.block_counts)
end
