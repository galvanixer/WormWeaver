const COUNTERFLOW_EQUAL_TIME_METADATA = (
    "displacements", "displacement_vectors", "displacement_counts",
    "beta", "equal_time_window", "equal_time_convention",
)
const COUNTERFLOW_GREEN_DEFINITION = (
    fields=("zSectorVisitCount",), optional_fields=("gREqualTimeHits", "gKTau"),
    reduction="sum", integer_sums=false, complex_fields=("gKTau",),
    metadata=("species_pair", "channel", COUNTERFLOW_EQUAL_TIME_METADATA...,
              "production_log_eta", "endpoint_convention", "normalization",
              "equal_time_normalization", "momentum_indices", "momenta", "tau_edges",
              "frequency_indices", "frequencies"),
)

"""Validate counterflow counts, ordered species identity, and available measurement axes."""
function validate_counterflow_green(histories, counts, metadata)
    metadata["channel"] == "counterflow" || throw(ArgumentError(
        "CounterflowGreenFunction/channel must be counterflow"))
    pair = metadata["species_pair"]
    pair isa AbstractVector{<:Integer} && !(eltype(pair) <: Bool) &&
        length(pair) == 2 && all(>(0), pair) && allunique(pair) || throw(ArgumentError(
        "CounterflowGreenFunction/species_pair must contain two distinct positive integer IDs"))
    zcounts = histories["zSectorVisitCount"]
    size(zcounts) == (length(counts),) || throw(DimensionMismatch(
        "CounterflowGreenFunction Z counts must have one entry per block"))
    all(>(0), zcounts) || throw(ArgumentError(
        "CounterflowGreenFunction requires positive Z visits in every saved production block"))
    if haskey(metadata, "production_log_eta")
        weight = metadata["production_log_eta"]
        weight isa Real && !(weight isa Bool) && isfinite(weight) || throw(ArgumentError(
            "CounterflowGreenFunction/production_log_eta must be a finite scalar"))
    end
    if haskey(histories, "gKTau")
        history = histories["gKTau"]
        ndims(history) == 3 && all(>(0), size(history)[1:2]) || throw(DimensionMismatch(
            "CounterflowGreenFunction/gKTau must have momentum × time-bin × block dimensions"))
        validate_fourier_grid(metadata, size(history)[1:2], "CounterflowGreenFunction")
    end
    haskey(histories, "gREqualTimeHits") || return nothing
    hits = histories["gREqualTimeHits"]
    ndims(hits) == 2 && size(hits, 1) > 0 || throw(DimensionMismatch(
        "CounterflowGreenFunction/gREqualTimeHits must have displacement × block dimensions"))
    displacements, vectors = metadata["displacements"], metadata["displacement_vectors"]
    displacements isa AbstractMatrix{<:Integer} && !(eltype(displacements) <: Bool) ||
        throw(ArgumentError("CounterflowGreenFunction/displacements must be an integer matrix"))
    size(displacements, 1) == size(hits, 1) && size(displacements, 2) > 0 ||
        throw(DimensionMismatch("CounterflowGreenFunction displacement coordinates must match the hits axis"))
    vectors isa AbstractMatrix{<:Real} && all(isfinite, vectors) || throw(ArgumentError(
        "CounterflowGreenFunction/displacement_vectors must be a finite real matrix"))
    size(vectors) == size(displacements) || throw(DimensionMismatch(
        "CounterflowGreenFunction displacement coordinates and vectors must have matching shapes"))
    multiplicities = metadata["displacement_counts"]
    multiplicities isa AbstractVector{<:Integer} && !(eltype(multiplicities) <: Bool) &&
        all(>(0), multiplicities) || throw(ArgumentError(
        "CounterflowGreenFunction/displacement_counts must contain positive integers"))
    length(multiplicities) == size(hits, 1) || throw(DimensionMismatch(
        "CounterflowGreenFunction displacement counts must match the hits axis"))
    beta, window = metadata["beta"], metadata["equal_time_window"]
    beta isa Real && window isa Real && !(beta isa Bool || window isa Bool) &&
        isfinite(beta) && isfinite(window) && 0 < window < beta || throw(ArgumentError(
        "CounterflowGreenFunction requires finite 0 < equal_time_window < beta"))
    convention = metadata["equal_time_convention"]
    convention isa AbstractString && startswith(convention, "one-sided") || throw(ArgumentError(
        "CounterflowGreenFunction requires a one-sided equal-time convention"))
    return nothing
end

"""
    counterflow_green_summary(histories, counts, metadata)

Sum raw equal-time hits, Fourier blocks, and Z visits across saved blocks.
Retain `gKTauSum` as a momentum × time-bin matrix and preserve the ordered
species pair, grids, and sample/block counts. Sums carry no sample weighting,
normalization, or statistical uncertainties.
"""
function counterflow_green_summary(histories, counts, metadata)
    summary = Dict{String,Any}("nblocks" => length(counts),
        "sampleCount" => checked_count_sum(counts), "aggregation_method" => "sum")
    summary["nsamples"] = summary["sampleCount"]
    for (field, history) in histories
        is_fourier_sum = field == "gKTau"
        number_type = is_fourier_sum ? ComplexF64 : Float64
        values = eltype(history) <: Integer ? observable_counts(history) :
            dropdims(sum(number_type, history; dims=ndims(history)); dims=ndims(history))
        all(isfinite, values) || throw(ArgumentError("CounterflowGreenFunction/$field sums overflowed"))
        key = is_fourier_sum ? "gKTauSum" : field
        summary[key] = ndims(values) == 0 ? only(values) : values
    end
    merge!(summary, metadata)
    return summary
end
