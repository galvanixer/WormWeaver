import DataFrames
import HDF5
import JLD2
import SHA

const OBSERVABLE_DEFINITIONS = Dict(
    "Density" => (
        fields=("densityBySpecies", "particleNumberBySpecies", "totalParticleNumber"),
        optional_fields=(), metadata=("species",), reduction="mean",
    ),
    "DensityCorrelation" => (
        fields=("full", "total"), optional_fields=(), reduction="mean",
        metadata=("displacements", "displacement_counts"),
    ),
    "Energy" => (
        fields=("diagonalEnergy", "diagonalEnergyPerSite", "kineticEnergy",
                "kineticEnergyPerSite", "totalEnergy", "totalEnergyPerSite"),
        optional_fields=(), metadata=(), reduction="mean",
    ),
    "Windings" => WINDINGS_DEFINITION,
    "StaticStructureFactor" => STATIC_STRUCTURE_FACTOR_DEFINITION,
    "CounterflowGreenFunction" => COUNTERFLOW_GREEN_DEFINITION,
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

"""Validate ordered component-pair histories and their saved displacement grid."""
function validate_density_correlation(histories, metadata)
    full, total = histories["full"], histories["total"]
    ndims(full) == 4 && ndims(total) == 2 || throw(DimensionMismatch(
        "DensityCorrelation requires component × component × displacement × block and displacement × block histories"))
    size(full, 1) == size(full, 2) && size(full, 1) > 0 && size(full, 3) > 0 ||
        throw(DimensionMismatch("DensityCorrelation requires nonempty, matching component axes and a nonempty displacement axis"))
    size(total) == (size(full, 3), size(full, 4)) || throw(DimensionMismatch(
        "DensityCorrelation full and total must have matching displacement and block axes"))
    displacements, multiplicities = metadata["displacements"], metadata["displacement_counts"]
    displacements isa AbstractMatrix{<:Integer} && !(eltype(displacements) <: Bool) &&
        size(displacements, 2) > 0 || throw(ArgumentError(
        "DensityCorrelation displacements must be an integer coordinate matrix"))
    multiplicities isa AbstractVector{<:Integer} && !(eltype(multiplicities) <: Bool) &&
        all(>(0), multiplicities) || throw(ArgumentError(
        "DensityCorrelation displacement_counts must be positive integer pair counts"))
    size(displacements, 1) == length(multiplicities) == size(full, 3) ||
        throw(DimensionMismatch("DensityCorrelation displacement metadata must match the history displacement axis"))
    return nothing
end

"""Read and validate one observable, then dispatch to its summary calculation."""
function read_observable_summary(file, observable, expected_blocks; min_blocks, blocking_diagnostics::Bool=false)
    definition = OBSERVABLE_DEFINITIONS[observable]
    data = read_observable_data(file, observable, definition)
    ismissing(data) && return missing
    histories, counts, metadata = data.histories, data.counts, data.metadata
    if isempty(counts)
        @warn "$observable has no saved blocks; leaving the $observable cell missing" file=HDF5.filename(file)
        return missing
    end
    length(counts) <= expected_blocks || throw(ArgumentError(
        "$observable has $(length(counts)) blocks, but the manifest expects $expected_blocks"))
    if observable == "ParticlesHistogram"
        validate_particle_histogram(histories, counts, metadata)
    elseif observable == "SingleParticleGreenFunction"
        validate_single_particle_green(histories, counts, metadata)
    elseif observable == "DensityCorrelation"
        validate_density_correlation(histories, metadata)
    elseif observable == "StaticStructureFactor"
        validate_static_structure_factor(histories, metadata)
    elseif observable == "CounterflowGreenFunction"
        validate_counterflow_green(histories, counts, metadata)
    end
    if observable == "ParticlesHistogram"
        return particle_histogram_summary(histories, counts, metadata; min_blocks, blocking_diagnostics)
    elseif observable == "SingleParticleGreenFunction"
        return single_particle_green_summary(histories, counts, metadata; min_blocks, blocking_diagnostics)
    elseif observable == "Windings"
        return windings_summary(histories, counts, metadata; min_blocks, blocking_diagnostics)
    elseif observable == "StaticStructureFactor"
        return static_structure_factor_summary(histories, counts, metadata; min_blocks, blocking_diagnostics)
    elseif observable == "CounterflowGreenFunction"
        return counterflow_green_summary(histories, counts, metadata)
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

"""Derive G2 entry modes from manifest configs, verifying recorded config hashes."""
function aggregate_g2_entry_modes(root, manifest)
    verify_hash = "config_sha256" in names(manifest)
    modes = Union{Missing,Tuple{Vararg{Symbol}}}[missing for _ in 1:DataFrames.nrow(manifest)]
    for (i, row) in enumerate(eachrow(manifest))
        path = normpath(joinpath(root, row.config))
        expected_hash = verify_hash ? row.config_sha256 : missing
        if !ismissing(expected_hash)
            actual_hash = try
                bytes2hex(SHA.sha256(read(path)))
            catch err
                err isa InterruptException && rethrow()
                @warn "Cannot verify configuration for G2 entry inspection; leaving g2_entry_modes missing" path manifest_row=i exception=err
                continue
            end
            if !(expected_hash isa AbstractString) ||
               lowercase(strip(expected_hash)) != actual_hash
                @warn "Configuration checksum does not match manifest; leaving g2_entry_modes missing" path manifest_row=i expected_hash actual_hash
                continue
            end
        end
        modes[i] = g2_entry_modes(path)
    end
    return modes
end

"""
    aggregate(study_dir; observables=["Density", "Energy", "Windings", "ParticlesHistogram", "SingleParticleGreenFunction", "DensityCorrelation", "StaticStructureFactor", "CounterflowGreenFunction"],
        min_blocks=10, blocking_diagnostics=false) -> DataFrame

Obtain the manifest with refreshed local-log progress from `status_manifest`,
then read each simulation's HDF5 block histories. Preserve its rows and all columns,
including `status`, `completed_blocks`, `failed`, and `log_source`, and add a
`g2_entry_modes` column derived from every simulation config, followed by a
dictionary column for each selected observable. Derive G2 entry modes even for
simulations without completed blocks. When `config_sha256` is recorded, require
the current config bytes to match it; otherwise leave the value missing with a
warning. Manifests without recorded hashes are inspected without checksum
verification. The source manifest must include
`study_id`, `simulation`, `config`, `output_file`, and `run.blocks`, with nonblank
study IDs and unique `(study_id, simulation)` keys. Simulations with zero
logged completed blocks are skipped, retaining rows with missing observable
summaries. Missing logs count as zero completed blocks.
Failed runs with saved production blocks remain eligible.
For eligible runs, summarize all valid saved blocks up to the planned count;
each observable's `nblocks` records its actual history length, which may differ
from the logged `completed_blocks`. Density, energy,
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
`Energy`, `Windings`, `ParticlesHistogram`, `SingleParticleGreenFunction`,
`DensityCorrelation`, `StaticStructureFactor`, and `CounterflowGreenFunction`. All eight are selected by default.
Each results file is opened once; each observable uses its own sample counts.

DensityCorrelation aggregates raw equal-time occupation correlations: full has
component × component × displacement shape and total is a displacement vector.
Both saved block-mean histories use sample-count weights and repeated-blocking
errors; total is analyzed independently to preserve its correlated fluctuations.
Required displacements and displacement_counts metadata retain their signed
coordinates, pair multiplicities, and saved order. Component indices are not
species IDs. Spatial pair averaging and imaginary-time averaging are already
applied; no further normalization, normal ordering, or disconnected subtraction
is performed. Same-site diagonal entries retain occupation squares.

StaticStructureFactor retains `structureFactorBySpeciesPair` as a component ×
component × momentum array of Complex{Measurement{Float64}} and the independently
analyzed real `structureFactorByMomentum` total. Each complex entry has separate
real and imaginary maximum repeated-blocking errors; no covariances are retained.
Both histories use their saved sample-count weights. Required `momentum_indices`
and `momenta` matrices retain the saved Fourier-grid order. Component indices are
not species IDs. Raw equal-time occupation products are already imaginary-time
averaged and site-normalized; no extra β/site factor, disconnected subtraction,
normal ordering, or symmetrization is applied. Only the current complex dataset
layout is supported. Optional diagnostics append blocking level to each array,
using `structureFactorBySpeciesPairRealStandardError`,
`structureFactorBySpeciesPairImaginaryStandardError`, and
`structureFactorByMomentumStandardError`.

CounterflowGreenFunction retains raw `gREqualTimeHits` summed over saved blocks
as a displacement vector, plus scalar sums of `zSectorVisitCount` and
`sampleCount`, `nsamples`, `nblocks`, and `aggregation_method="sum"`. Hits are
optional; absent hits do not imply a zero histogram. Require positive Z visits
in each saved block, preserve the ordered `species_pair` and signed displacement
metadata, and validate the one-sided finite-time window when hits exist.
`production_log_eta` is optional metadata and is not applied. No physical
normalization, uncertainty estimates, or blocking curves are computed, even with
`blocking_diagnostics=true`. Warn when other time/frequency representations are
present but omitted. These finite-window counts are not a normalized correlator.

SingleParticleGreenFunction retains raw gREqualTimeHits summed over blocks as a
species × displacement numeric matrix, plus scalar zSectorVisitCount and
sampleCount sums. nsamples equals sampleCount. Species, signed displacements,
multiplicities, and normalization metadata are preserved. Counts are not weighted
by sampleCount. When equal-time hits exist, grEqualTime contains physically
normalized Measurements with repeated-blocking paired jackknife errors; nsites,
normalization, error_method, error_status, and min_blocks describe the analysis.
Optional diagnostics include grEqualTimeStandardError, retained Z/sample totals,
and level statuses. The central value uses all raw counts, without bias correction.
When gKTau exists, its sample-count-weighted mean is retained as a species ×
momentum × time-bin array of Complex{Measurement{Float64}}, with independent
real/imaginary repeated-blocking errors. gKTauAggregationMethod and gKTauErrorMethod
describe this raw estimate. Optional diagnostics include gKTauRealStandardError
and gKTauImaginaryStandardError. No Z division or physical normalization is applied
to gKTau. Other time/frequency measurements are not summarized.
Equal-time analysis requires one-sided sums and normalization metadata; zero total
Z visits throw for equal-time analysis. No derived HDF5 datasets are read and no time/frequency physical
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
`min_blocks` controls additional blocking levels, not aggregation eligibility.
Missing files, absent groups, and empty observable histories leave missing
cells and emit warnings; other observables in the same run still contribute.
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

Return the DataFrame without writing files. Use `save_aggregate` to aggregate and
save the table as JLD2.
"""
function aggregate(study_dir::AbstractString; observables=["Density", "Energy", "Windings", "ParticlesHistogram", "SingleParticleGreenFunction", "DensityCorrelation", "StaticStructureFactor", "CounterflowGreenFunction"],
                   min_blocks::Integer=10,
                   blocking_diagnostics::Bool=false)
    min_blocks >= 2 || throw(ArgumentError("min_blocks must be at least two"))
    observables isa Union{AbstractVector,Tuple,AbstractSet} ||
        throw(ArgumentError("observables must be a collection of names"))
    selected = collect(observables)
    isempty(selected) && throw(ArgumentError("observables must not be empty"))
    all(name -> name isa AbstractString && haskey(OBSERVABLE_DEFINITIONS, name), selected) ||
        throw(ArgumentError("supported observables are Density, Energy, Windings, ParticlesHistogram, SingleParticleGreenFunction, DensityCorrelation, StaticStructureFactor, and CounterflowGreenFunction"))
    allunique(selected) || throw(ArgumentError("observables must not contain duplicates"))
    root = abspath(study_dir)
    manifest = status_manifest(root)
    all(name -> name in names(manifest), ("output_file", "run.blocks")) ||
        throw(ArgumentError("manifest must contain output_file and run.blocks"))
    any(name -> name in names(manifest), selected) &&
        throw(ArgumentError("manifest already contains a selected observable column"))
    "g2_entry_modes" in names(manifest) &&
        throw(ArgumentError("manifest must not contain the derived g2_entry_modes column"))

    table = DataFrames.DataFrame(manifest; copycols=true)
    table[!, "g2_entry_modes"] = aggregate_g2_entry_modes(root, manifest)
    for observable in selected
        table[!, observable] = Union{Missing,Dict{String,Any}}[missing for _ in 1:DataFrames.nrow(table)]
    end
    for (i, row) in enumerate(eachrow(manifest))
        expected_blocks = row["run.blocks"]
        expected_blocks isa Integer && !(expected_blocks isa Bool) && expected_blocks > 0 ||
            throw(ArgumentError("run.blocks must be a positive integer in manifest row $i"))
        row.completed_blocks >= 1 || continue
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
    return table
end

"""
    save_aggregate(study_dir; output="aggregated.jld2", overwrite=false,
        observables=["Density", "Energy", "Windings", "ParticlesHistogram", "SingleParticleGreenFunction", "DensityCorrelation", "StaticStructureFactor", "CounterflowGreenFunction"],
        min_blocks=10, blocking_diagnostics=false) -> DataFrame

Call `aggregate` with the selected observables and calculation options, save its
DataFrame under the `aggregate` key in a JLD2 file, and return that DataFrame.
The output path is relative to the study directory unless absolute. Create parent
directories as needed. Existing outputs require `overwrite=true`; source configs,
logs, and results are protected from replacement, and symlink outputs are rejected.

For selections containing observables other than `SingleParticleGreenFunction`
and `CounterflowGreenFunction`, also save
`error_method="repeated blocking"` and `min_blocks` as top-level metadata.
"""
function save_aggregate(study_dir::AbstractString;
                        output::AbstractString="aggregated.jld2",
                        overwrite::Bool=false,
                        observables=["Density", "Energy", "Windings", "ParticlesHistogram", "SingleParticleGreenFunction", "DensityCorrelation", "StaticStructureFactor", "CounterflowGreenFunction"],
                        min_blocks::Integer=10,
                        blocking_diagnostics::Bool=false)
    root = abspath(study_dir)
    isempty(strip(output)) && throw(ArgumentError("output must be a nonempty path"))
    destination = normpath(joinpath(root, output))
    islink(destination) && throw(ArgumentError("output must not be a symlink"))
    if ispath(destination)
        isfile(destination) || throw(ArgumentError("output must be a regular file"))
        overwrite || throw(ArgumentError("Output already exists: $destination; use overwrite=true"))
    end

    table = aggregate(root; observables, min_blocks, blocking_diagnostics)
    protected = [normpath(joinpath(root, file)) for file in table.output_file
                 if file isa AbstractString]
    for row in eachrow(table)
        config_path = normpath(joinpath(root, row.config))
        push!(protected, config_path, joinpath(dirname(config_path), "simulation.log"))
    end
    real_destination = ispath(destination) ? realpath(destination) : destination
    any(path -> real_destination == (ispath(path) ? realpath(path) : path), protected) &&
        throw(ArgumentError("output cannot replace a source config, log, or results file"))

    mkpath(dirname(destination))
    if all(name -> name in ("SingleParticleGreenFunction", "CounterflowGreenFunction"), observables)
        JLD2.jldsave(destination; aggregate=table)
    else
        JLD2.jldsave(destination; aggregate=table, error_method="repeated blocking", min_blocks)
    end
    return table
end
