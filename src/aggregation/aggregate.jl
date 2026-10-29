import DataFrames
import HDF5
import JLD2

const OBSERVABLE_DEFINITIONS = Dict(
    "Density" => (
        fields=("densityBySpecies", "particleNumberBySpecies", "totalParticleNumber"),
        optional_fields=(), metadata=("species",), reduction="mean",
    ),
    "Energy" => (
        fields=("diagonalEnergy", "diagonalEnergyPerSite", "kineticEnergy",
                "kineticEnergyPerSite", "totalEnergy", "totalEnergyPerSite"),
        optional_fields=(), metadata=(), reduction="mean",
    ),
    "Windings" => WINDINGS_DEFINITION,
    "SingleParticleGreenFunction" => SINGLE_PARTICLE_GREEN_DEFINITION,
    "ParticlesHistogram" => (
        fields=("particleNumberBySpeciesCounts", "totalParticleNumberCounts", "zSectorSampleCount"),
        optional_fields=(), reduction="sum",
        metadata=("species", "species_particle_number_bins", "total_particle_number_bins",
                  "max_particle_number_by_species"),
    ),
)

function observable_measurements(history::AbstractArray{<:Real}, counts; min_blocks)
    shape = size(history)[1:end-1]
    size(history, ndims(history)) == length(counts) ||
        throw(DimensionMismatch("observable history must have one entry per saved block"))
    matrix = reshape(history, prod(shape), length(counts))
    measurements = Measurements.Measurement{Float64}[
        blocking_measurement(collect(row), counts; min_blocks) for row in eachrow(matrix)]
    return isempty(shape) ? only(measurements) : reshape(measurements, shape)
end

"""Read and validate one observable, then dispatch to its summary calculation."""
function read_observable_summary(file, observable, expected_blocks; min_blocks, blocking_diagnostics::Bool=false)
    definition = OBSERVABLE_DEFINITIONS[observable]
    data = read_observable_data(file, observable, definition)
    ismissing(data) && return missing
    histories, counts, metadata = data.histories, data.counts, data.metadata
    if observable == "ParticlesHistogram"
        validate_particle_histogram(histories, counts, metadata)
    elseif observable == "SingleParticleGreenFunction"
        validate_single_particle_green(histories, counts, metadata)
    end
    if length(counts) < expected_blocks
        @warn "$observable run is incomplete; leaving the $observable cell missing" file=HDF5.filename(file) saved_blocks=length(counts) expected_blocks
        return missing
    end
    length(counts) == expected_blocks || throw(ArgumentError(
        "$observable has $(length(counts)) blocks, but the manifest expects $expected_blocks"))

    if observable == "ParticlesHistogram"
        return particle_histogram_summary(histories, counts, metadata; min_blocks, blocking_diagnostics)
    elseif observable == "SingleParticleGreenFunction"
        return single_particle_green_summary(histories, counts, metadata; min_blocks, blocking_diagnostics)
    elseif observable == "Windings"
        return windings_summary(histories, counts, metadata; min_blocks, blocking_diagnostics)
    end
    return mean_observable_summary(histories, counts, metadata; min_blocks, blocking_diagnostics)
end

"""Summarize validated block-mean histories with shared statistics and optional diagnostics."""
function mean_observable_summary(histories, counts, metadata; min_blocks, blocking_diagnostics::Bool=false)
    summary = Dict{String,Any}()
    diagnostics = blocking_diagnostics ? Dict{String,Any}() : nothing
    for (field, history) in histories
        if blocking_diagnostics
            statistics = observable_blocking_statistics(history, counts; min_blocks)
            summary[field] = statistics.measurements
            diagnostics[field * "StandardError"] = statistics.standard_errors
            if !haskey(diagnostics, "block_widths")
                diagnostics["block_widths"] = statistics.block_widths
                diagnostics["block_counts"] = statistics.block_counts
                diagnostics["sample_counts"] = blocking_sample_counts(counts, statistics.block_widths, statistics.block_counts)
            end
        else
            summary[field] = observable_measurements(history, counts; min_blocks)
        end
    end
    blocking_diagnostics && (summary["blocking_diagnostics"] = diagnostics)
    summary["nblocks"] = length(counts)
    summary["nsamples"] = sum(counts)
    summary["error_method"] = "repeated blocking: maximum standard error"
    summary["min_blocks"] = min_blocks
    merge!(summary, metadata)
    return summary
end

"""
    aggregate(study_dir; observables=["Density", "Energy", "Windings", "ParticlesHistogram", "SingleParticleGreenFunction"],
        min_blocks=10, blocking_diagnostics=false,
        output="aggregated.jld2", overwrite=false) -> DataFrame

Read `manifest.jld2` and each simulation's HDF5 block histories. Preserve manifest
rows and add a dictionary column for each selected observable. Density, energy,
and winding fields contain native Measurements.jl values (mean ± standard error).
Particle histograms retain checked integer sums across blocks and add normalized
probabilities and first/second particle-number moments with repeated-blocking
errors. Histogram species moments are vectors; total moments are scalars.
Set `blocking_diagnostics=true` to
calculate and retain per-level diagnostic curves for observables with error estimates.
Diagnostics are not constructed by default; blocking uncertainties are still
calculated in either mode.
Summaries retain block/sample counts and
observable metadata. Scalars and array shapes are preserved.
`observables` must be a nonempty collection of unique supported names: `Density`,
`Energy`, `Windings`, `ParticlesHistogram`, and `SingleParticleGreenFunction`.
All five are selected by default.
Each results file is opened once; each observable uses its own sample counts.

SingleParticleGreenFunction retains raw gREqualTimeHits summed over blocks as a
species × displacement numeric matrix, plus scalar zSectorVisitCount and
sampleCount sums. nsamples equals sampleCount. Species, signed displacements,
multiplicities, and normalization metadata are preserved. Counts are not weighted
by sampleCount. When equal-time hits exist, grEqualTime contains physically
normalized Measurements with repeated-blocking paired jackknife errors; nsites,
normalization, error_method, error_status, and min_blocks describe the analysis.
Optional diagnostics include grEqualTimeStandardError, retained Z/sample totals,
and level statuses. The central value uses all raw counts, without bias correction.
Time/frequency-only groups retain raw Z/sample summaries and available metadata.
Equal-time analysis requires one-sided sums and normalization metadata; zero total
Z visits throw. No derived HDF5 datasets are read and no time/frequency physical
normalization is performed.

Winding summaries add squaredCoflowBySpeciesPair and squaredCounterflowBySpeciesPair
as species-pair × spatial-axis Measurement matrices indexed by species_pairs.
They measure unnormalized (W_a + W_b)^2 and (W_a - W_b)^2, constructed from saved
block moments before blocking. squaredTotalWinding remains global coflow.
Single-species pair arrays are empty. Two-species files without pair products use
total and species squares to reconstruct the pair; larger systems without products,
or products without pair metadata, omit pair channels with a warning.

For mean observables and histogram probabilities, errors are the maximum across
successive pairwise blocking levels, retaining at least `min_blocks` blocks at
each additional level. Means use
all original blocks and sample-count weights. A single saved block gives an
undefined (NaN) error.
Missing files, absent groups, and incomplete observable histories leave missing
cells and emit warnings; complete observables in the same run still contribute.
Invalid data throw an error. Particle histograms require explicit sum reduction
metadata, as do SPGF sum fields; mean observables accept legacy files without reduction metadata as
block means. Histogram summaries contain counts with aggregation_method="sum",
probability Measurement arrays, meanParticleNumberBySpecies, meanTotalParticleNumber,
meanSquaredParticleNumberBySpecies, meanSquaredTotalParticleNumber, error_method,
and min_blocks. Moment histories come from raw block histograms using bin weights
n and n^2; reported means use full aggregated counts. Total moments use the saved
total histogram, preserving inter-species correlations.
Histogram means are ratios of full raw counts to the total sample count; errors
use each block's normalized probabilities weighted by its sample count. Diagnostic
curves record errors, block widths/counts, and exact retained sample totals at
every blocking level, without an automatic convergence assessment. With
`blocking_diagnostics=true`, each summary gains a `blocking_diagnostics` dictionary.
It contains block_widths, block_counts, sample_counts, and a StandardError array
for each estimate field, with blocking level appended as the final dimension
(a level vector for scalar fields). Histogram error fields refer to probabilities
and particle-number moments.
Burn-in is already excluded by the simulation and is not discarded again.
Component uncertainties do not reconstruct covariance; aggregate saved totals and
squares directly. Keep original block histograms for joint variance or
compressibility uncertainties.

Save the table under the `aggregate` key in `output`, relative to the study folder
unless absolute. Set `output=nothing` to only return the table. Existing outputs
require `overwrite=true`; source results and manifests are never overwritten.
"""
function aggregate(study_dir::AbstractString; observables=["Density", "Energy", "Windings", "ParticlesHistogram", "SingleParticleGreenFunction"],
                   min_blocks::Integer=10,
                   blocking_diagnostics::Bool=false,
                   output::Union{Nothing,AbstractString}="aggregated.jld2",
                   overwrite::Bool=false)
    min_blocks >= 2 || throw(ArgumentError("min_blocks must be at least two"))
    observables isa Union{AbstractVector,Tuple,AbstractSet} ||
        throw(ArgumentError("observables must be a collection of names"))
    selected = collect(observables)
    isempty(selected) && throw(ArgumentError("observables must not be empty"))
    all(name -> name isa AbstractString && haskey(OBSERVABLE_DEFINITIONS, name), selected) ||
        throw(ArgumentError("supported observables are Density, Energy, Windings, ParticlesHistogram, and SingleParticleGreenFunction"))
    allunique(selected) || throw(ArgumentError("observables must not contain duplicates"))
    root = abspath(study_dir)
    manifest_path = joinpath(root, "manifest.jld2")
    manifest = JLD2.load(manifest_path, "manifest")
    manifest isa DataFrames.AbstractDataFrame || throw(ArgumentError("manifest must be a DataFrame"))
    all(name -> name in names(manifest), ("output_file", "run.blocks")) ||
        throw(ArgumentError("manifest must contain output_file and run.blocks"))
    any(name -> name in names(manifest), selected) &&
        throw(ArgumentError("manifest already contains a selected observable column"))
    destination = output === nothing ? nothing : normpath(joinpath(root, output))
    if destination !== nothing
        isempty(strip(output)) && throw(ArgumentError("output must be a nonempty path"))
        protected = [manifest_path, joinpath(root, "manifest.csv")]
        append!(protected, [normpath(joinpath(root, file)) for file in manifest.output_file
                            if file isa AbstractString])
        for row in eachrow(manifest)
            "config" in names(manifest) && row.config isa AbstractString &&
                push!(protected, normpath(joinpath(root, row.config)))
        end
        real_destination = ispath(destination) ? realpath(destination) : destination
        any(path -> real_destination == (ispath(path) ? realpath(path) : path), protected) &&
            throw(ArgumentError("output cannot replace a source manifest, config, or results file"))
        islink(destination) && throw(ArgumentError("output must not be a symlink"))
        if ispath(destination)
            isfile(destination) || throw(ArgumentError("output must be a regular file"))
            overwrite || throw(ArgumentError("Output already exists: $destination; use overwrite=true"))
        end
    end

    table = DataFrames.DataFrame(manifest; copycols=true)
    for observable in selected
        table[!, observable] = Union{Missing,Dict{String,Any}}[missing for _ in 1:DataFrames.nrow(table)]
    end
    for (i, row) in enumerate(eachrow(manifest))
        expected_blocks = row["run.blocks"]
        expected_blocks isa Integer && !(expected_blocks isa Bool) && expected_blocks > 0 ||
            throw(ArgumentError("run.blocks must be a positive integer in manifest row $i"))
        if !(row.output_file isa AbstractString) || isempty(strip(row.output_file))
            @warn "No results path; leaving selected observable cells missing" manifest_row=i
            continue
        end
        path = normpath(joinpath(root, row.output_file))
        if !isfile(path)
            @warn "Results file is missing; leaving selected observable cells missing" path
            continue
        end
        HDF5.h5open(path, "r") do file
            for observable in selected
                table[i, observable] = read_observable_summary(file, observable, expected_blocks; min_blocks, blocking_diagnostics)
            end
        end
    end
    if destination !== nothing
        mkpath(dirname(destination))
        if all(==("SingleParticleGreenFunction"), selected)
            JLD2.jldsave(destination; aggregate=table)
        else
            JLD2.jldsave(destination; aggregate=table, error_method="repeated blocking", min_blocks)
        end
    end
    return table
end
