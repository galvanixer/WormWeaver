# Particle-histogram counts, probabilities, moments, blocking diagnostics, and validation.

function observable_counts(history::AbstractArray{<:Integer})
    shape = size(history)[1:end-1]
    matrix = reshape(history, prod(shape), size(history, ndims(history)))
    totals = [checked_count_sum(row) for row in eachrow(matrix)]
    return isempty(shape) ? only(totals) : reshape(totals, shape)
end

function histogram_probability_statistics(history, totals, counts, nsamples;
                                         min_blocks, blocking_diagnostics::Bool=false)
    normalized = Float64.(history) ./ reshape(counts, ntuple(_ -> 1, ndims(history) - 1)..., length(counts))
    statistics = blocking_diagnostics ? observable_blocking_statistics(normalized, counts; min_blocks) : nothing
    measurements = blocking_diagnostics ? statistics.measurements : observable_measurements(normalized, counts; min_blocks)
    # Means come directly from the complete integer totals, independently of
    # truncation and floating-point pairing at the error-estimation levels.
    probabilities = Measurements.measurement.(totals ./ nsamples, Measurements.uncertainty.(measurements))
    return (; probabilities, statistics)
end

function histogram_moment_statistics(history, totals, bins, counts, nsamples, order;
                                    min_blocks, blocking_diagnostics::Bool=false)
    shape = size(history)[1:end-2]
    ncomponents = prod(shape)
    nbins = length(bins)
    matrix = reshape(history, ncomponents, nbins, length(counts))
    total_matrix = reshape(totals, ncomponents, nbins)
    # Floating-point bin weights make both n^2 and count-weight products safe
    # from integer overflow. Moments use raw counts, not probability Measurements.
    weights = Float64.(bins) .^ order
    block_moments = [sum(weights[n] * matrix[i, n, b] for n in 1:nbins) / counts[b]
                     for i in 1:ncomponents, b in eachindex(counts)]
    moment_history = reshape(block_moments, shape..., length(counts))
    statistics = blocking_diagnostics ? observable_blocking_statistics(moment_history, counts; min_blocks) : nothing
    estimates = blocking_diagnostics ? statistics.measurements : observable_measurements(moment_history, counts; min_blocks)
    means = [sum(weights[n] * total_matrix[i, n] for n in 1:nbins) / nsamples
             for i in 1:ncomponents]
    full_mean = isempty(shape) ? only(means) : reshape(means, shape)
    measurements = Measurements.measurement.(full_mean, Measurements.uncertainty.(estimates))
    return (; measurements, statistics)
end

function particle_histogram_summary(histories, counts, metadata; min_blocks, blocking_diagnostics::Bool=false)
    summary = Dict{String,Any}(field => observable_counts(history) for (field, history) in histories)
    nsamples = checked_count_sum(counts)
    species = histogram_probability_statistics(histories["particleNumberBySpeciesCounts"],
        summary["particleNumberBySpeciesCounts"], counts, nsamples; min_blocks, blocking_diagnostics)
    total = histogram_probability_statistics(histories["totalParticleNumberCounts"],
        summary["totalParticleNumberCounts"], counts, nsamples; min_blocks, blocking_diagnostics)
    summary["particleNumberBySpeciesProbability"] = species.probabilities
    summary["totalParticleNumberProbability"] = total.probabilities
    summary["nsamples"] = nsamples
    summary["nblocks"] = length(counts)
    summary["aggregation_method"] = "sum"
    summary["error_method"] = "repeated blocking: maximum standard error"
    summary["min_blocks"] = min_blocks
    if blocking_diagnostics
        summary["blocking_diagnostics"] = Dict{String,Any}(
            "block_widths" => species.statistics.block_widths,
            "block_counts" => species.statistics.block_counts,
            "sample_counts" => blocking_sample_counts(counts, species.statistics.block_widths, species.statistics.block_counts),
            "particleNumberBySpeciesStandardError" => species.statistics.standard_errors,
            "totalParticleNumberStandardError" => total.statistics.standard_errors,
        )
    end
    for (field, counts_field, bins_field, order) in (
        ("meanParticleNumberBySpecies", "particleNumberBySpeciesCounts", "species_particle_number_bins", 1),
        ("meanTotalParticleNumber", "totalParticleNumberCounts", "total_particle_number_bins", 1),
        ("meanSquaredParticleNumberBySpecies", "particleNumberBySpeciesCounts", "species_particle_number_bins", 2),
        ("meanSquaredTotalParticleNumber", "totalParticleNumberCounts", "total_particle_number_bins", 2),
    )
        moments = histogram_moment_statistics(histories[counts_field], summary[counts_field],
            metadata[bins_field], counts, nsamples, order; min_blocks, blocking_diagnostics)
        summary[field] = moments.measurements
        if blocking_diagnostics
            summary["blocking_diagnostics"][field * "StandardError"] = moments.statistics.standard_errors
        end
    end
    merge!(summary, metadata)
    return summary
end

function validate_particle_histogram(histories, counts, metadata)
    counts isa AbstractVector{<:Integer} && !(eltype(counts) <: Bool) ||
        throw(ArgumentError("ParticlesHistogram/sampleCount must contain integers"))
    species = histories["particleNumberBySpeciesCounts"]
    total = histories["totalParticleNumberCounts"]
    zcounts = histories["zSectorSampleCount"]
    ndims(species) == 3 && ndims(total) == 2 && ndims(zcounts) == 1 ||
        throw(DimensionMismatch("particle histograms require species × bin × block, bin × block, and block histories"))
    nspecies, nbins, nblocks = size(species)
    nspecies > 0 && nbins > 0 && size(total, 1) > 0 ||
        throw(ArgumentError("particle histogram species and bin axes must not be empty"))
    for (field, length_expected) in (
        ("species", nspecies), ("species_particle_number_bins", nbins),
        ("total_particle_number_bins", size(total, 1)), ("max_particle_number_by_species", nspecies),
    )
        value = metadata[field]
        value isa AbstractVector{<:Integer} && !(eltype(value) <: Bool) ||
            throw(ArgumentError("ParticlesHistogram/$field metadata must be an integer vector"))
        length(value) == length_expected || throw(DimensionMismatch("ParticlesHistogram/$field metadata has the wrong length"))
    end
    labels = metadata["species"]
    all(>(0), labels) && allunique(labels) ||
        throw(ArgumentError("particle histogram species labels must be positive and unique"))
    for field in ("species_particle_number_bins", "total_particle_number_bins")
        bins = metadata[field]
        bins == collect(0:(length(bins) - 1)) ||
            throw(ArgumentError("ParticlesHistogram/$field must enumerate consecutive bins starting at zero"))
    end
    capacities = metadata["max_particle_number_by_species"]
    all(capacity -> 0 <= capacity < nbins, capacities) ||
        throw(ArgumentError("species capacities must fit the histogram bins"))
    zcounts == counts || throw(ArgumentError("zSectorSampleCount must equal sampleCount in every block"))
    for b in 1:nblocks
        checked_count_sum(@view total[:, b]) == zcounts[b] ||
            throw(ArgumentError("total particle histogram counts do not match zSectorSampleCount in block $b"))
        for s in 1:nspecies
            checked_count_sum(@view species[s, :, b]) == zcounts[b] ||
                throw(ArgumentError("species $s histogram counts do not match zSectorSampleCount in block $b"))
            all(iszero, @view species[s, (capacities[s] + 2):nbins, b]) ||
                throw(ArgumentError("species $s has counts above its particle-number capacity"))
        end
    end
    return nothing
end
