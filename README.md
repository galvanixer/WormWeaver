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

WormWeaver is in early development. Study folder generation is available; the broader job execution, aggregation, and analysis workflow remains under development.

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

Sweep paths are sorted alphabetically, with the last path varying fastest and values following their listed order. Omitting `sweep` generates one parameter combination. Without `replicas`, it produces one simulation with `replica = 1`; seeds remain unchanged unless explicitly included in overrides or sweeps. Run each simulation from its own directory when using relative output filenames such as `results.h5`.

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

This generates 30 simulations: ten replicas for each beta. Omit `sweep` for ten simulations with identical parameters apart from the selected seeds. Replicas are numbered from 1 within each combination and vary fastest in simulation folder ordering.

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

## Author and contact

Tanul Gupta — [tanulgupta123@gmail.com](mailto:tanulgupta123@gmail.com)
