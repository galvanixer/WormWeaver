# Equal-time normalization is an adapter over the shared paired-ratio engine.
const SINGLE_PARTICLE_GREEN_EQUAL_TIME_METADATA = (
    "species", "displacements", "displacement_vectors", "displacement_counts",
    "beta", "equal_time_window", "equal_time_convention", "production_log_eta",
)
const SINGLE_PARTICLE_GREEN_DEFINITION = (
    fields=("zSectorVisitCount",), optional_fields=("gREqualTimeHits",),
    reduction="sum", integer_sums=false,
    metadata=(SINGLE_PARTICLE_GREEN_EQUAL_TIME_METADATA...,
              "momentum_indices", "momenta", "tau_edges", "frequency_indices", "frequencies",
              "endpoint_convention", "normalization", "equal_time_normalization"),
)

"""
    validate_single_particle_green(histories, counts, metadata)

Check Z block counts and, when equal-time hits exist, their species/displacement
axes and normalization metadata. Return `nothing`; throw on invalid data.
"""
function validate_single_particle_green(histories, counts, metadata)
    zcounts = histories["zSectorVisitCount"]
    size(zcounts) == (length(counts),) || throw(DimensionMismatch(
        "SingleParticleGreenFunction Z counts must have one entry per block"))
    haskey(histories, "gREqualTimeHits") || return nothing
    hits = histories["gREqualTimeHits"]
    validate_equal_time_axes(hits, counts, metadata)
    validate_green_geometry(metadata, size(hits, 2))
    validate_equal_time_metadata(metadata)
    return nothing
end

"""Check equal-time histogram dimensions and selected-species labels."""
function validate_equal_time_axes(hits, counts, metadata)
    ndims(hits) == 3 || throw(DimensionMismatch(
        "SingleParticleGreenFunction/gREqualTimeHits must have species × displacement × block dimensions"))
    nspecies, ndisplacements, nblocks = size(hits)
    nspecies > 0 && ndisplacements > 0 || throw(DimensionMismatch(
        "SingleParticleGreenFunction species and displacement axes must not be empty"))
    nblocks == length(counts) || throw(DimensionMismatch(
        "SingleParticleGreenFunction counts must have one entry per block"))

    species = metadata["species"]
    species isa AbstractVector{<:Integer} && !(eltype(species) <: Bool) ||
        throw(ArgumentError("SingleParticleGreenFunction/species must be an integer vector"))
    length(species) == nspecies || throw(DimensionMismatch(
        "SingleParticleGreenFunction/species must match the hits species axis"))
    all(>(0), species) && allunique(species) || throw(ArgumentError(
        "SingleParticleGreenFunction/species must contain positive unique IDs"))
    weights = metadata["production_log_eta"]
    weights isa AbstractVector{<:Real} && all(isfinite, weights) || throw(ArgumentError(
        "SingleParticleGreenFunction/production_log_eta must be a finite numeric vector"))
    length(weights) == nspecies || throw(DimensionMismatch(
        "SingleParticleGreenFunction/production_log_eta must match the species axis"))
    return nothing
end

"""Check displacement coordinates, Cartesian vectors, and pair multiplicities."""
function validate_green_geometry(metadata, ndisplacements)
    displacements = metadata["displacements"]
    displacements isa AbstractMatrix{<:Integer} && !(eltype(displacements) <: Bool) ||
        throw(ArgumentError("SingleParticleGreenFunction/displacements must be an integer matrix"))
    size(displacements, 1) == ndisplacements && size(displacements, 2) > 0 ||
        throw(DimensionMismatch("SingleParticleGreenFunction/displacements must match the hits displacement axis"))
    vectors = metadata["displacement_vectors"]
    vectors isa AbstractMatrix{<:Real} && all(isfinite, vectors) || throw(ArgumentError(
        "SingleParticleGreenFunction/displacement_vectors must be a finite numeric matrix"))
    size(vectors) == size(displacements) || throw(DimensionMismatch(
        "SingleParticleGreenFunction displacement coordinates and vectors must have matching shapes"))
    multiplicities = metadata["displacement_counts"]
    multiplicities isa AbstractVector{<:Integer} && !(eltype(multiplicities) <: Bool) &&
        all(>(0), multiplicities) || throw(ArgumentError(
        "SingleParticleGreenFunction/displacement_counts must contain positive integers"))
    length(multiplicities) == ndisplacements || throw(DimensionMismatch(
        "SingleParticleGreenFunction/displacement_counts must match the hits displacement axis"))
    equal_time_site_count(multiplicities)
    return nothing
end

"""Infer site count from a checked total of all ordered displacement pairs."""
function equal_time_site_count(multiplicities)
    npairs = checked_count_sum(multiplicities)
    nsites = isqrt(npairs)
    nsites * nsites == npairs || throw(ArgumentError(
        "SingleParticleGreenFunction displacement counts must cover nsites^2 ordered pairs"))
    return nsites
end

"""Require a finite physical time window and one-sided equal-time convention."""
function validate_equal_time_metadata(metadata)
    beta, window = metadata["beta"], metadata["equal_time_window"]
    beta isa Real && window isa Real && !(beta isa Bool || window isa Bool) &&
        isfinite(beta) && isfinite(window) && 0 < window < beta ||
        throw(ArgumentError("SingleParticleGreenFunction requires finite 0 < equal_time_window < beta"))
    convention = metadata["equal_time_convention"]
    convention isa AbstractString && startswith(convention, "one-sided") ||
        throw(ArgumentError("SingleParticleGreenFunction requires one-sided equal-time counts"))
    return nothing
end

"""Return site count and species × displacement physical normalization factors."""
function equal_time_normalization(metadata)
    nsites = equal_time_site_count(metadata["displacement_counts"])
    beta, window = metadata["beta"], metadata["equal_time_window"]
    # Logarithms avoid forming a potentially overflowing sector weight first.
    factors = [exp(log(nsites) + log(beta) - log(count) - log(window) - log_weight)
               for log_weight in metadata["production_log_eta"], count in metadata["displacement_counts"]]
    all(isfinite, factors) || throw(ArgumentError("SingleParticleGreenFunction normalization overflowed"))
    return (; nsites, factors)
end

"""
    single_particle_green_summary(histories, counts, metadata;
        min_blocks=10, blocking_diagnostics=false)

Summarize validated histories as a dictionary retaining raw sums and metadata.
When equal-time hits exist, add normalized `grEqualTime` Measurements with paired
jackknife errors across blocking levels, plus optional diagnostic curves.
"""
function single_particle_green_summary(histories, counts, metadata;
                                       min_blocks=10, blocking_diagnostics::Bool=false)
    summary = summarize_raw_green_counts(histories, counts, metadata)
    haskey(histories, "gREqualTimeHits") || return summary

    normalization = equal_time_normalization(metadata)
    statistics = green_ratio_statistics(histories["gREqualTimeHits"], histories["zSectorVisitCount"];
        sample_counts=counts, transform=ratio -> normalization.factors .* ratio,
        min_blocks, blocking_diagnostics)
    summary["grEqualTime"] = statistics.measurements
    summary["nsites"] = normalization.nsites
    summary["normalization"] = "nsites*beta/(C_r*equal_time_window*exp(production_log_eta_s)) * sum(hits_s_r)/sum(z_counts)"
    summary["error_method"] = "repeated blocking: maximum paired jackknife standard error"
    summary["error_status"] = statistics.error_status
    summary["min_blocks"] = min_blocks
    if blocking_diagnostics
        diagnostics = statistics.blocking_diagnostics
        diagnostics["grEqualTimeStandardError"] = pop!(diagnostics, "real_standard_errors")
        summary["blocking_diagnostics"] = diagnostics
    end
    return summary
end

"""Sum raw Z, sample, and optional equal-time hit counts while retaining metadata."""
function summarize_raw_green_counts(histories, counts, metadata)
    zcounts = histories["zSectorVisitCount"]
    summed_z = eltype(zcounts) <: Integer ? checked_count_sum(zcounts) : sum(Float64, zcounts)
    isfinite(summed_z) || throw(ArgumentError(
        "SingleParticleGreenFunction count sums overflowed"))
    nsamples = checked_count_sum(counts)
    summary = Dict{String,Any}(
        "zSectorVisitCount" => summed_z,
        "sampleCount" => nsamples,
        "nsamples" => nsamples,
        "nblocks" => length(counts),
        "aggregation_method" => "sum",
    )
    merge!(summary, metadata)
    if haskey(histories, "gREqualTimeHits")
        hits = histories["gREqualTimeHits"]
        summed_hits = eltype(hits) <: Integer ? observable_counts(hits) :
            dropdims(sum(Float64, hits; dims=3); dims=3)
        all(isfinite, summed_hits) || throw(ArgumentError("SingleParticleGreenFunction hit sums overflowed"))
        summary["gREqualTimeHits"] = summed_hits
    end
    return summary
end
