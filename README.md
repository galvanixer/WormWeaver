# WormWeaver

WormWeaver is a Julia library for **job generation, aggregation, and analysis** of simulations run using **multiwormqmc.jl**. It provides the workflow layer around the simulation engine, helping turn a collection of simulation runs into organized, reproducible results.

## Purpose

Simulation studies often involve many parameter combinations and independent runs. Preparing jobs, keeping track of their outputs, and bringing those outputs together for analysis are essential parts of the work.

WormWeaver is being developed to handle these tasks in a consistent workflow, with multiwormqmc.jl responsible for executing the simulations.

## Scope

### Job generation

Define simulation studies and prepare jobs for multiwormqmc.jl:

- Describe simulation parameters and parameter sweeps.
- Prepare inputs for individual runs and independent repetitions.
- Organize job names, output locations, and metadata so results can be traced back to their inputs.

### Aggregation

Collect simulation outputs into datasets suitable for analysis:

- Discover and load results from multiple runs.
- Group outputs by simulation parameters and observables.
- Combine compatible results while retaining the metadata needed to interpret them.

### Analysis

Extract physical quantities and interpret simulation results:

- Summarize observables and estimate statistical uncertainties.
- Compare results across parameters and system sizes.
- Support fitting, visualization, and export for further scientific work.

These describe the intended scope of the library; they are not a statement of currently implemented features.

## Relationship to multiwormqmc.jl

The two libraries have complementary responsibilities:

| Library | Responsibility |
| --- | --- |
| **multiwormqmc.jl** | Simulation algorithms, simulation execution, and measurement of observables. |
| **WormWeaver** | Preparing simulation jobs, aggregating their outputs, and analyzing the resulting data. |

WormWeaver builds on multiwormqmc.jl. Simulation algorithms belong in the simulation library; the surrounding study workflow belongs here.

## Intended workflow

1. **Define a study:** choose model parameters, system sizes, simulation settings, and independent runs.
2. **Generate jobs:** prepare the inputs and job artifacts needed to run the study with multiwormqmc.jl.
3. **Run simulations:** execute the generated jobs using multiwormqmc.jl.
4. **Aggregate outputs:** collect results and organize them by their parameters and observables.
5. **Analyze results:** estimate quantities and uncertainties, compare runs, and produce figures or exported datasets.

Throughout this workflow, results should remain traceable to the simulation inputs and analysis choices that produced them.

## Design goals

- **Reproducibility:** make study definitions and analysis choices explicit and reusable.
- **Traceability:** preserve the relationship between jobs, parameters, raw outputs, and derived results.
- **Modularity:** allow job generation, aggregation, and analysis to be used independently or together.
- **Scientific reliability:** account for statistical uncertainty and the compatibility of results when combining simulation data.
- **Practical workflows:** support repeated studies without requiring a new set of ad hoc scripts for each one.

## Development status

WormWeaver is in early development. Study folder generation, density/energy/winding aggregation, and particle-histogram aggregation are available; the broader job execution and analysis workflow remains under development.

An editable [Unistra HPC study template](templates/studies/unistra-hpc.yaml) is available for scratch-mode runs. For installation and batch execution on Unistra, see [Running multiwormqmc.jl on Unistra HPC](docs/unistra-hpc.md).

## Generate a study

From the WormWeaver directory, install dependencies and start Julia:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=.
```

Create a commented study definition from the [general study template](templates/studies/study.yaml):

```julia
using WormWeaver
init_study("my_study.yaml")                   # General template (default).
# init_study("my_study.yaml"; template="unistra-hpc.yaml") # Unistra HPC template.

# Edit my_study.yaml: set output_dir, parameters, and job settings.
generate_study("my_study.yaml");
```

`template` is the filename inside `templates/studies/` and defaults to `"study.yaml"`. Adding a file there makes it available automatically. Unknown names report the available templates.

`init_study` copies the template verbatim, preserving comments, and creates missing parent directories. It returns the absolute YAML path. Existing files are protected unless `overwrite=true` is supplied; directories and symlinks are rejected. It only creates the study definition, without generating simulations or submitting jobs.

To copy a bundled simulation config into a new directory for editing:

```julia
init_config("my_run/config.yaml"; template="LRBHQM/config_LRBHQM_regular.yaml")
# Copy every LRBHQM config into a directory, preserving filenames:
init_config("my_configs"; template="LRBHQM")
```

`init_config` defaults to `template="config.yaml"`. Templates live in `templates/multiwormqmc/`, with the four LRBHQM variants grouped under `LRBHQM/`. A file template copies to an output file; a group template copies all its files into an output directory. Unique bare filenames such as `config_LRBHQM_regular.yaml` also work. Copies preserve comments and create missing directories. The function returns the absolute output file or directory. Existing files require `overwrite=true`; group copies check all destinations before writing and preserve unrelated directory contents. Symlink destinations are rejected. To use a copied config in a study definition in the same directory, set `template: ./config.yaml` in that study YAML.

After successful generation, WormWeaver prints the simulation count, parameter-combination and replica counts, output directory, job mode, and submission commands for the copied Slurm scripts. It still returns the absolute study directory and does not submit jobs.

The study file specifies constant overrides and a Cartesian parameter sweep:

```yaml
project: RFEBHM
output_dir: /path/to/my_study
overrides:
  lattice.extent: [32]
  run.blocks: 100
  run.samples_per_block: 10000
sweep:
  system.beta: [2.0, 4.0, 8.0]
  model.onsite.mu:
    - [1.0, 1.0]
    - [2.0, 2.0]
```

Relative custom paths are resolved from the study file. Omit `template` to use the bundled `config.yaml`. A unique bare filename or a path starting with a bundled group selects a bundled file, for example `template: LRBHQM/config_LRBHQM_regular.yaml`. Use `template: ./config.yaml`, `../configs/config.yaml`, or an absolute path for a custom file. Prefix local paths with `./` when their directory name matches a bundled group. Bare filenames never fall back to local files; unknown names list the bundled templates. Study generation requires one file, so select a variant within a group. Dotted field names address existing keys in the template. Arrays are replaced in full: each chemical-potential pair above is one sweep value.

The LRBHQM group contains regular, aligned, anti-aligned, and single-worm variants. Select `LRBHQM/config_LRBHQM_regular.yaml` for regular `Open`/`Close` moves between `[1, 0]` or `[0, 1]` and `[1, 1]`, or the aligned and anti-aligned variants for their corresponding colocated moves. The single-worm variant never enters `[1, 1]`. Keep `pair_green` and `counterflow_green` disabled with the regular and single-worm variants; those estimators require colocated moves.

This example generates six simulations:

```text
/path/to/my_study/
├── study.yaml
├── template.yaml
├── provenance.yaml
├── jobfile
├── grant_g2026a136c.slurm
├── grantgpu_g2026a136g.slurm
├── manifest.csv
├── manifest.jld2
├── sim1/config.yaml
├── sim2/config.yaml
├── sim3/config.yaml
├── sim4/config.yaml
├── sim5/config.yaml
└── sim6/config.yaml
```

By default, `append_date: true` appends the generation machine's local date to the output folder, for example `my_study_2026-09-30`. Set `append_date: false` to keep the folder name unchanged. An explicit `jobs.run_dir` receives the same suffix; otherwise it inherits the dated output directory. Supply undated base paths. Existing same-day folders still require `overwrite=true`. The saved `study.yaml` preserves the original definition, including its undated paths.

The required `project` field is a nonblank string identifying the broader research project, such as `RFEBHM` or `LRBHQM`. Multiple studies can share a project name. It is preserved in the study snapshot and repeated in the manifest's `project` column, without being added to simulation configs or changing the explicit `output_dir`.

The root YAML files preserve the original study definition and template verbatim. The copied study definition is a provenance record: its relative paths still refer to the original study location.

Both manifests have the same columns and one row per simulation:

| Columns | Contents |
| --- | --- |
| `project`, `study_id`, `simulation`, `replica` | Project name, a new UUID per generated study, local folder name, and replica number within each parameter combination. |
| `config`, `output_file` | Config and expected results paths relative to the study directory. An explicitly absolute results path stays absolute. |
| `config_sha256` | SHA-256 of the exact generated config file bytes, for detecting subsequent edits. |
| `system.*`, `lattice.*` | Temperature, ensemble, geometry, extents, and boundary conditions. |
| Model parameters | Model family, species count, hopping, onsite and offsite interactions, and disorder settings. |
| Initialization and seeds | Initial densities, occupation cutoffs, and initialization, simulation, and disorder seeds. |
| Sampling settings | Production blocks, samples per block, burn-in blocks, and multiworm order. |
| All overridden or swept fields | Additional selected parameters, taken from each final generated config. |

In `manifest.csv`, scalar values are written directly; arrays and mappings are YAML serialized into CSV cells. Missing or null fields are blank, allowing smaller custom templates. In `manifest.jld2`, the `manifest` key stores a DataFrame with native scalar, array, and dictionary values. Absent fields remain `missing`, while explicit YAML null values remain `nothing`. Large structures such as move schedules and tuning transitions are omitted unless explicitly overridden or swept. The core column list is defined in `CORE_FIELDS` in `src/study_generation/study.jl`.

Load the typed manifest for analysis with:

```julia
using JLD2, DataFrames
manifest = JLD2.load(joinpath(study_dir, "manifest.jld2"), "manifest")
mu = manifest[1, "model.onsite.mu"]  # A native array; no YAML parsing needed.
```

Here `study_dir` is the directory returned by `generate_study`. Both manifests are generated together from the same data; editing one afterward does not update the other.

Metadata column names (`project`, `study_id`, `simulation`, `replica`, `config`, `output_file`, and `config_sha256`) are reserved and cannot be override or sweep paths. Use `run.output_file` to change a simulation's output filename. Identical study definitions produce the same simulation ordering, but receive a fresh study ID when generated again. When replica randomization is enabled without `seed_generation_seed`, generated seeds and config hashes also change. The study ID is stored in the manifest and remains unchanged when the study directory is moved. Execution status and scheduler job IDs are not part of this generation manifest.

Ordinary `sweep` paths are sorted alphabetically, with the last path varying fastest and values following their listed order. Omitting both `sweep` and `correlated_sweeps` generates one parameter combination. Without `replicas`, it produces one simulation with `replica = 1`; seeds remain unchanged unless explicitly included in overrides or sweeps. Run each simulation from its own directory when using relative output filenames such as `results.h5`.

### Correlated sweeps

Use `correlated_sweeps` to vary fields together by matching list positions:

```yaml
project: RFEBHM
output_dir: /path/to/finite_size_study
sweep:
  model.onsite.mu:
    - [0.0, 0.0]
    - [1.0, 1.0]
correlated_sweeps:
  - lattice.extent: [[16], [32], [64]]
    system.beta: [16.0, 32.0, 64.0]
replicas: 5
```

This pairs extent `[16]` with beta `16.0`, `[32]` with `32.0`, and `[64]` with `64.0`. Combined with two chemical-potential vectors and five replicas, it generates 30 simulations.

Each group must be a nonempty mapping of existing field paths to nonempty lists of the same length. Arrays and mappings remain complete parameter values, just as in ordinary sweeps. Multiple correlated groups vary independently, forming a Cartesian product with each other and with ordinary `sweep` fields. Groups follow ordinary sweep fields in their listed order; the last group varies fastest, followed by replicas within each combination. Fields within a group are applied together before `validate_config` is called.

A field can appear only once across `overrides`, `sweep`, and all correlated groups. Parent/child path overlaps, reserved manifest names, and conflicts with randomized seed fields are rejected across all of these sections. Both manifests include correlated fields using their final config values. Omit `correlated_sweeps` or use `correlated_sweeps: []` when no correlated groups are needed.

### Replicas with random seeds

```yaml
project: RFEBHM
output_dir: /path/to/seed_study
replicas: 10
randomize_seeds:
  - run.seed
  - worldlines.init_seed
# Add model.disorder.seed above to vary the disorder realization too.
sweep:
  system.beta: [2.0, 4.0, 8.0]
```

This generates 30 simulations: ten replicas for each beta. Omit both sweep sections for ten simulations with identical parameters apart from the selected seeds. Replicas are numbered from 1 within each combination and vary fastest in simulation folder ordering.

`replicas` must be a positive integer. If `randomize_seeds` is omitted, it defaults to `[run.seed]`. An explicit list must be nonempty, contain no duplicates, and select only `run.seed`, `worldlines.init_seed`, or `model.disorder.seed`; it requires `replicas`. Unselected seeds retain their configured values. Overrides or sweeps that overlap selected seed fields are rejected, including replacement of their parent blocks.

Each selected seed gets a random positive Julia `Int`, distinct from every other generated seed in that study. Without `seed_generation_seed`, fresh seeds are drawn on each generation, including overwrite; keep the generated configs to rerun exactly the same simulations. Both manifests record `replica` and all three actual seed values. A private random number generator is used without reseeding the caller's default RNG. Unselected fixed seeds are not included in the uniqueness guarantee.

### Direct jobfile generation

Every study includes a study-root `jobfile` with one shell command per simulation in manifest order, plus every `.slurm` template in `templates/slurm/`. Add an optional `jobs` block to customize execution:

```yaml
jobs:
  command: [multiwormqmc]
  run_dir: /cluster/path/to/study
  mode: direct
```

`command` is an argument list, not a shell expression. It defaults to `[multiwormqmc]`; WormWeaver appends `config.yaml` and redirects stdout and stderr to each simulation's `simulation.log`. In direct mode, paths and arguments are shell quoted. Each direct-mode command changes into its simulation directory using `&&`, runs in a subshell, and returns the simulation command's exit status.

`run_dir` must be the absolute location of the study on the execution machine; it defaults to the generated study directory. It does not transfer files. Modes are `direct` (default) and `scratch`. Omit `jobs`, or use `jobs: {}`, for defaults. For a Julia script launcher, use an argument list such as `[julia, /cluster/path/multiwormqmc.jl/bin/multiwormqmc.jl]`.

Two account-specific templates are included, based on the HSCAnalysis launcher:

- `grant_g2026a136c.slurm`: `#SBATCH -p grant -A g2026a136c`
- `grantgpu_g2026a136g.slurm`: `#SBATCH -p grantgpu -A g2026a136g`

Both load `gcc/15.2.0` before invoking `hpc_multilauncher jobfile` and retain the original node, memory, and wall-time settings. The GPU-partition script does not add a GPU resource request or change the simulation executable. Choose the script appropriate to your allocation and submit from the generated study directory, for example `sbatch grant_g2026a136c.slurm`. The cluster must provide `hpc_multilauncher` and the configured command. Adding or removing `.slurm` files in `templates/slurm/` controls the scripts copied into future studies. WormWeaver does not submit jobs. Regenerating with `overwrite=true` replaces script edits with the bundled templates.

### Node-local scratch execution

```yaml
jobs:
  command: [multiwormqmc]
  run_dir: /cluster/path/to/study
  mode: scratch
  scratch_root: /scratch
```

Scratch mode emits the original HSCAnalysis command format exactly:

```bash
cp -r /path/study/sim1 /scratch/job.$SLURM_JOB_ID/sim1; cd /scratch/job.$SLURM_JOB_ID/sim1 ;multiwormqmc config.yaml &> simulation.log ; cp -r /scratch/job.$SLURM_JOB_ID/sim1/* /path/study/sim1 ;
```

The scratch job directory must already exist. Commands and paths are emitted as supplied, without automatic quoting, directory creation, or added failure handling. The shell expands `~` and `$SLURM_JOB_ID`. Copy-back uses `*`, as in HSCAnalysis. This mode requires Bash-compatible redirection. `scratch_root` defaults to `/scratch`.

### Reproducible seeds and software provenance

Set `seed_generation_seed` alongside `replicas` to regenerate the same replica seeds:

```yaml
replicas: 10
seed_generation_seed: 12345
randomize_seeds:
  - run.seed
  - worldlines.init_seed
# Optional: identify the simulation package checkout you plan to use.
multiwormqmc_path: /path/to/multiwormqmc.jl
```

The generation seed must be a nonnegative integer fitting Julia's `Int`. With the same study parameters, seed-field order, Julia version, and generator implementation, it reproduces seeds and configs, including when overwriting. Study UUIDs still change. Julia does not guarantee identical random sequences across versions; preserve generated configs for long-term reproducibility.

Every study includes `provenance.yaml` recording its study ID, Julia version, seed-generation algorithm and supplied seed, and WormWeaver and MultiwormQMC package versions, source paths, Git commits, and dirty states where available. This describes software identified at generation time, not proof of the software eventually used to execute jobs.

`multiwormqmc_path` is resolved relative to the study definition and must point to a MultiwormQMC package directory. If omitted, WormWeaver discovers MultiwormQMC from the Julia environment without loading it. If unavailable, the provenance marks it unavailable with null version/revision fields. Git commits are null for packages without their own repository or without a commit; dirty states include untracked files. No files are committed automatically.

Unknown fields, incompatible top-level value types, empty sweep lists, and overlapping override/sweep paths are rejected. Existing output directories are rejected by default. Use `generate_study("/path/to/study.yaml"; overwrite=true)` to replace the entire output directory, including old simulation outputs and any extra files. Replacement generates a fresh study ID. Files and symlinks are not accepted as replacement destinations. All configurations are prepared before output is written, and files are staged before publishing the study directory. During replacement, the old directory is retained until publication succeeds and restored if publication fails. Generated configs use compact formatting: short scalar arrays are inline, small records in lists use flow style, and larger sections remain indented blocks. Strings are quoted to preserve their meaning. Values are preserved, but original formatting and comments are not; the template snapshot retains both.

Validation here checks study structure and field types, not the full physical constraints enforced by multiwormqmc.jl. For additional checks, supply `generate_study(path; validate_config=my_validator)`, where `my_validator(config)` receives each generated configuration dictionary and throws if it is invalid. All callbacks run before writing output.

Run the tests with:

```sh
julia --project=. test/runtests.jl
```

## Aggregate observables

Aggregate saved production-block density, energy, winding, particle histogram,
and available single-particle Green measurements into one row per simulation, retaining every column of
the typed study manifest:

```julia
using WormWeaver, Measurements, JLD2

study_dir = "scratch/LRBHQM2_alpha1p5_2026-10-05"
table = aggregate(study_dir)  # Selects all five supported observables.
density = table.Density[1]["densityBySpecies"]  # Vector of mean ± uncertainty.
Measurements.value.(density)        # Species means.
Measurements.uncertainty.(density)  # Repeated-blocking standard errors.
energy = table.Energy[1]["totalEnergy"]         # Scalar mean ± uncertainty.
windings = table.Windings[1]["squaredWindingBySpecies"]  # Species × spatial axis.

# Select a subset; output=nothing only returns the table:
density_table = aggregate(study_dir; observables=["Density"], output=nothing)

# Read the saved table later:
table = JLD2.load(joinpath(study_dir, "aggregated.jld2"), "aggregate")
```

`observables` must be a nonempty vector, tuple, or set of unique supported names:
`"Density"`, `"Energy"`, `"Windings"`, `"ParticlesHistogram"`, and
`"SingleParticleGreenFunction"`. All five are
selected by default. Each selected observable adds a dictionary column.
Existing manifest columns with selected observable names are rejected.

- `Density` contains `densityBySpecies`, `particleNumberBySpecies`, and
  `totalParticleNumber`, plus `species` labels when present.
- `Energy` contains `diagonalEnergy`, `kineticEnergy`, and `totalEnergy`, and
  their `PerSite` counterparts.
- `Windings` contains `windingBySpecies`, `squaredWindingBySpecies`,
  `totalWinding`, and `squaredTotalWinding`. `windingCorrelationBySpecies` is
  included when saved; single-species runs may omit it. `periodic_axes` and
  `species_pairs` metadata are retained when present. Pair coflow and counterflow
  are derived as described below.
- `ParticlesHistogram` contains `particleNumberBySpeciesCounts`,
  `totalParticleNumberCounts`, and `zSectorSampleCount` as raw integer sums
  across blocks, plus normalized probabilities, particle-number moments, repeated-blocking errors,
  optional diagnostic curves, species labels, particle-number bins, and capacities.
- `SingleParticleGreenFunction` contains available raw `gREqualTimeHits` summed over the
  block axis as a species × displacement numeric matrix, plus scalar sums of
  `zSectorVisitCount` and `sampleCount`. `nsamples` equals `sampleCount`.
  Species IDs, signed displacement coordinates and vectors, multiplicities,
  beta, equal-time window and convention, and production sector weights are
  retained. Counts are not weighted by measurement-attempt counts. When equal-time
  hits exist, `grEqualTime` contains normalized species × displacement `Measurement`
  values with repeated-blocking paired jackknife errors. Analysis is computed
  from original histories, without reading `/derived/SingleParticleGreenFunction`.
  Time/frequency-only groups return Z/sample summaries and available metadata;
  equal-time results are omitted when that histogram was not measured.

Select only Green-function aggregation:

```julia
table = aggregate(study_dir; observables=["SingleParticleGreenFunction"], output=nothing)
green = observable_table(table.SingleParticleGreenFunction)
green.gREqualTimeHits[1]      # Species × signed displacement, summed over blocks.
green.zSectorVisitCount[1]   # Total closed-sector visits.
green.sampleCount[1]         # Total measurement attempts; no extra weighting.
green.grEqualTime[1]         # Normalized finite-window correlation ± error, if measured.
Measurements.value.(green.grEqualTime[1])
Measurements.uncertainty.(green.grEqualTime[1])
```

For equal-time sampling, normalization is
`nsites*beta/(displacement_counts[r]*equal_time_window*exp(production_log_eta[s]))`
times the ratio of total hits to total Z visits. `nsites` is inferred from the
checked ordered-pair total, which must equal `nsites^2`. The metadata must describe
one-sided sampling and satisfy `0 < equal_time_window < beta`. Normalization is
independent of spatial dimension and supports nonorthogonal lattice vectors.

The reported mean always includes every original block. Paired jackknife removes
the same block from hits and Z visits. Repeated blocking sums adjacent pairs of
both histories, taking the maximum valid jackknife standard error across levels.
Neither operation changes the reported mean or applies a bias correction. Zero
total Z visits throw; insufficient blocks or zero leave-one-out denominators
produce NaN errors when no valid level remains. Otherwise invalid levels are
omitted and reported by `error_status`. Error bars exclude finite-window bias.

The exported array helper supports future time/frequency adapters and custom
derived quantities without depending on MultiwormQMC:

```julia
stats = green_ratio_statistics(numerator_history, z_counts;
    reduction=:mean, sample_counts=sample_counts,
    transform=identity, min_blocks=10, blocking_diagnostics=true)
stats.measurements
stats.error_status
stats.blocking_diagnostics
```

The final numerator axis is block; a vector gives a scalar result. `:mean` requires
sample counts and reconstructs block sums once; `:sum` does not weight by sample
counts. `transform` acts on the full ratio and each jackknife replica and must
return a finite scalar or array with consistent shape and real/complex type.
Complex outputs are combined `Complex{Measurement{Float64}}` values with separate
real and imaginary uncertainties. Joint derived quantities should be computed
inside `transform`; individual output errors do not retain covariance.

The helper accepts the existing SPGF numerator representations `gRTauNumerator`,
the real/imaginary `gKTau`, `gROmega`, and `gKOmega` numerator pairs after joining
complex components. These are sector-reweighted block means, unlike equal-time
raw sums. Their physical normalization adapters are not implemented by `aggregate`
yet; the helper's default `identity` produces the numerator/Z ratio alone.

The aggregation code follows three steps: `aggregate.jl` coordinates runs and
selects summaries; `observable_io.jl` reads saved histories and metadata; each
observable file validates and analyzes its data. For SPGF,
`single_particle_green.jl` supplies equal-time normalization and
`ratio_statistics.jl` handles paired jackknife and repeated blocking. The latter
separates history preparation, one-level jackknife errors, adjacent-block merging,
and optional diagnostics. Mean-based statistics remain in `statistics.jl`.

Density, energy, and winding estimate fields contain native `Measurement`
values, preserving scalar, vector, and tensor shapes after removing the final
block dimension. Their summaries contain `error_method` and `min_blocks`.
All summaries contain `nblocks` and `nsamples` (the contributing sample count).
Each simulation's results file is opened once.

Winding summaries additionally contain `squaredCoflowBySpeciesPair` and
`squaredCounterflowBySpeciesPair` whenever the saved moments identify the pairs.
Both are `Measurement` matrices with axes `(species_pair, spatial_axis)`, indexed
by the rows of `species_pairs`. For pair `(a, b)` and spatial axis `d`, they report
`⟨(W[a,d] + W[b,d])²⟩` and `⟨(W[a,d] - W[b,d])²⟩`, respectively. Weights are
unnormalized: there is no division by two or by species count. These are
dimensionless winding moments; stiffness prefactors are applied separately.
`squaredTotalWinding` remains global coflow, `⟨(sum(W[:,d]))²⟩`. With two species,
pair and global coflow coincide; with more species, they generally differ. There
is no single global counterflow field for three or more species.

```julia
w = table.Windings[1]
w["species_pairs"]                       # Pair rows contain one-based a < b.
w["squaredCoflowBySpeciesPair"]           # Pair × axis, mean ± uncertainty.
w["squaredCounterflowBySpeciesPair"]      # Same pair order and shape.
```

Each pair's block histories are constructed as `A + B ± 2C`, where `A` and `B`
are the saved block means of squared species windings and `C` is their saved
pair product. These histories are blocked together with sample-count weights,
preserving the correlations needed for channel uncertainties. Squaring a block's
mean winding or combining separately aggregated errors gives a different result.
With `blocking_diagnostics=true`, the corresponding
`squaredCoflowBySpeciesPairStandardError` and
`squaredCounterflowBySpeciesPairStandardError` arrays have axes
`(species_pair, spatial_axis, blocking_level)`.

Single-species summaries contain empty `0 × spatial_dimension` pair-channel
matrices and empty `species_pairs` metadata. For two species without saved pair
products, the sole product is reconstructed from total and species squares and
`species_pairs` is set to `[1 2]`. For three or more species without products, or
when products lack pair metadata, pair-channel fields are omitted with a warning;
existing winding estimates remain available. Valid incomplete pair lists produce
only the listed channels, in their supplied order. Pair indices must be unique
integer rows satisfying `1 <= a < b <= nspecies`; malformed metadata and
incompatible dimensions raise errors. Open axes retain zero winding moments.

To aggregate particle histograms and expand the summaries into a table:

```julia
table = aggregate(study_dir; observables=["ParticlesHistogram"],
    blocking_diagnostics=true, output=nothing)
histograms = observable_table(table.ParticlesHistogram)
species_counts = histograms.particleNumberBySpeciesCounts[1]  # Species × bin.
total_counts = histograms.totalParticleNumberCounts[1]         # Bin vector.
sample_count = histograms.zSectorSampleCount[1]
probabilities = histograms.particleNumberBySpeciesProbability[1]  # Species × bin.
Measurements.value.(probabilities)
Measurements.uncertainty.(probabilities)
total_probabilities = histograms.totalParticleNumberProbability[1]
histograms.meanParticleNumberBySpecies[1]         # Species vector of mean ± error.
histograms.meanTotalParticleNumber[1]             # Scalar mean ± error.
histograms.meanSquaredParticleNumberBySpecies[1]  # Species second moments.
histograms.meanSquaredTotalParticleNumber[1]      # Total second moment.
diagnostics = histograms.blocking_diagnostics[1]
diagnostics["block_widths"]
diagnostics["particleNumberBySpeciesStandardError"]  # Species × bin × level.
```

For each simulation, particle histograms sum every bin across all saved production
blocks using checked integer addition. Species remain separate, zero-count bins
and padded bins remain intact, and the saved total-particle histogram is summed
directly. Each `ParticlesHistogram` summary contains `aggregation_method="sum"`,
`species`, `species_particle_number_bins`, `total_particle_number_bins`, and
`max_particle_number_by_species`. Bin index `n + 1` represents particle number
`n`. `zSectorSampleCount` and `nsamples` both contain the summed configuration
count, counted once per configuration rather than once per species.

Alongside the raw counts, every histogram summary contains
`particleNumberBySpeciesProbability` (species × bin) and
`totalParticleNumberProbability` (bin vector) as native `Measurement` arrays.
Their values are the full raw counts divided by `nsamples`. For uncertainties,
each original block is normalized by its own `zSectorSampleCount`, then repeatedly
blocked with sample-count weights. The largest standard error at any available
level is reported separately for each bin. All original blocks contribute to the
probability value, even if unmatched trailing blocks are omitted at coarser
error-estimation levels. Raw counts remain integers and carry no error bars.

Every histogram summary also contains four particle-number moments:

| Field | Meaning | Shape |
| --- | --- | --- |
| `meanParticleNumberBySpecies` | ⟨Nₛ⟩ | Species vector. |
| `meanTotalParticleNumber` | ⟨N_total⟩ | Scalar. |
| `meanSquaredParticleNumberBySpecies` | ⟨Nₛ²⟩ | Species vector. |
| `meanSquaredTotalParticleNumber` | ⟨N_total²⟩ | Scalar. |

Each is a native `Measurement` with a repeated-blocking standard error. For each
original block, first and second moments are calculated from its raw histogram
using the bin labels `n` and `n²`, divided by that block's sample count. The
moment histories are blocked independently with sample-count weights. Reported
means use the complete aggregated counts, including unmatched trailing blocks.
Second moments average squared particle numbers, rather than squaring average
particle numbers. Total-particle moments use the saved total histogram directly,
preserving inter-species correlations. These fields report individual moment
uncertainties; covariance-aware variance or compressibility errors require joint
analysis of the original blocks.

Histogram summaries also contain `error_method` and `min_blocks`.
Set `blocking_diagnostics=true` to calculate and retain `blocking_diagnostics`
for density, energy, windings, histograms, and measured equal-time Green functions:

```julia
table = aggregate(study_dir; blocking_diagnostics=true)
density_diagnostics = table.Density[1]["blocking_diagnostics"]
density_diagnostics["densityBySpeciesStandardError"]  # Species × level.
energy_diagnostics = table.Energy[1]["blocking_diagnostics"]
energy_diagnostics["totalEnergyStandardError"]         # Level vector.
winding_diagnostics = table.Windings[1]["blocking_diagnostics"]
winding_diagnostics["windingBySpeciesStandardError"]   # Species × axis × level.
# For a run with equal-time sampling:
green_diagnostics = table.SingleParticleGreenFunction[1]["blocking_diagnostics"]
green_diagnostics["grEqualTimeStandardError"]          # Species × displacement × level.
```

The switch defaults to `false`. In that mode, diagnostic curves and sampling
grids are not constructed or retained, and summaries have no `blocking_diagnostics`
key. Means, probabilities, and their repeated-blocking uncertainties are calculated
in both modes and give the same results.

When enabled, diagnostics include the initial level and every additional
level permitted by `min_blocks` (default 10), retaining at least that many blocks
at each additional level:

| Diagnostic key | Contents |
| --- | --- |
| `block_widths` | Number of original blocks combined at each level: 1, 2, 4, … |
| `block_counts` | Number of blocks remaining at each level. |
| `sample_counts` | Sample totals retained at each level; integer-count histories remain exact. |
| `particleNumberBySpeciesStandardError` | Species × bin × level error array. |
| `totalParticleNumberStandardError` | Bin × level error array. |
| `meanParticleNumberBySpeciesStandardError` | Species × level moment-error array. |
| `meanTotalParticleNumberStandardError` | Level vector of total first-moment errors. |
| `meanSquaredParticleNumberBySpeciesStandardError` | Species × level second-moment error array. |
| `meanSquaredTotalParticleNumberStandardError` | Level vector of total second-moment errors. |

For density, energy, and windings, each saved estimate field gets a corresponding
`<field>StandardError` key. The error array preserves the field's shape and appends
the blocking level as its final axis; scalar fields have a level vector. Optional
estimate fields are included only when present in the source results. Each
observable's diagnostics use its own sample counts.

Equal-time Green diagnostics additionally include `z_counts` (retained Z totals)
and `level_statuses`. The array helper uses `real_standard_errors` and, for complex
outputs, `imaginary_standard_errors` instead of a representation-specific name.
Its `sample_counts` curve is included only when sample counts are supplied.

These curves support inspection without an automatic convergence flag. A
single-block history retains its probabilities and moments with NaN uncertainties. Constant
or unobserved bins have zero estimated error when at least two blocks are
available; zero observed error does not establish convergence or rule out unseen
events. Diagnostic curves and individual `Measurement` values do not retain
cross-bin or cross-moment covariance. Keep the original block histograms for
reliable joint errors on variance or compressibility.

Independent runs remain separate rows. Original HDF5 block histories remain
unchanged, and existing saved aggregates are refreshed only by explicitly calling
`aggregate` (with `overwrite=true` if the output already exists).

To turn an observable's dictionary column into its own DataFrame:

```julia
windings_table = observable_table(table.Windings)
windings_table.totalWinding

# Select specific fields and their column order:
totals = observable_table(table.Windings;
    fields=["totalWinding", "squaredTotalWinding"])
```

`observable_table` works with all observable summaries. It preserves
input row order and uses the sorted union of dictionary keys as columns by
default. Measurements, arrays, and metadata remain native values; vectors and
matrices remain individual cells. Missing summaries and absent keys yield missing
cells. Supply `fields` for an all-missing column so its rows can be retained.
Input dictionaries are not modified; array-valued cells reference the original
arrays.

Mean-reduced observables use all saved production blocks, weighted by each observable's own
`sampleCount`; equal-time Green functions use ratios of raw count totals. Uncertainty is the largest valid standard error encountered while
repeatedly combining consecutive pairs of blocks. The default `min_blocks=10`
retains at least ten blocks at each additional blocking level. An unmatched
trailing block is omitted only at that error-estimation level; the reported mean
still includes it. Burn-in has already been excluded by multiwormqmc. A
single-block history retains its mean with a NaN uncertainty. The maximum-error
rule does not certify an error plateau. Component uncertainties do not reconstruct
covariance; use directly aggregated totals and saved squares for their errors.

The function returns the DataFrame and saves `aggregated.jld2` in the study
directory. Use `output=nothing` to only return the table, or specify another
output path (relative to the study directory unless absolute). Existing output
files require `overwrite=true`; source manifests, configs, and results are
protected from replacement.

Missing results, absent selected groups, and observable histories with fewer
blocks than `run.blocks` retain their manifest rows with missing cells and
warnings. Complete observables in the same simulation still contribute. Malformed
measurements, invalid sample counts, excess blocks, and missing required fields
raise errors. Density, energy, and windings accept files without reduction metadata
as legacy block means. Particle histograms require explicit `sum` reduction
metadata and all bin/species metadata. Counts must be nonnegative integers;
each block's total histogram and each species histogram must independently sum
to its positive `zSectorSampleCount`, which must equal `sampleCount`. Incompatible
dimensions, metadata lengths, populated bins above a species capacity, and integer
overflow raise errors.

## Author and contact

Tanul Gupta — [tanulgupta123@gmail.com](mailto:tanulgupta123@gmail.com)
