# Running multiwormqmc.jl on Unistra HPC

This guide uses the multiwormqmc.jl checkout launcher and WormWeaver's generated jobfile. Replace the example username and study paths with your own.

## 1. Prepare the environment

Log in to the HPC and load GCC:

```bash
module load gcc/15.2.0
julia --version
```

MultiwormQMC requires Julia 1.12 or a compatible later Julia 1.x release. Ensure your Julia installation is on `PATH`, including inside batch jobs.

## 2. Install multiwormqmc.jl

Your HPC account needs GitLab SSH access to the repository.

```bash
mkdir -p ~/Code
cd ~/Code
git clone git@gitlab.com:galvanixer/multiwormqmc.jl.git
cd multiwormqmc.jl
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
julia bin/multiwormqmc.jl --help
```

If you already have a checkout, use that directory instead of cloning again. The checkout and Julia dependencies must be accessible from compute nodes. Keep its `Project.toml`, `Manifest.toml`, and Git revision to identify the simulation environment.

## 3. Install WormWeaver on the HPC

Clone WormWeaver into your code directory and install its dependencies:

```bash
module load gcc/15.2.0
mkdir -p ~/Code
cd ~/Code
git clone https://github.com/galvanixer/WormWeaver.git
cd WormWeaver
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
julia --project=. -e 'using WormWeaver; println("WormWeaver loaded")'
```

If the repository requires authentication, use your GitHub credentials or its SSH clone URL. If you already have a checkout, use it instead of cloning again. This setup runs WormWeaver from its checkout environment; it does not require registering the package or installing a separate executable.

Optionally verify the checkout:

```bash
julia --project=. test/runtests.jl
```

Keep study definitions and generated simulation data outside the WormWeaver repository. For example:

```bash
mkdir -p ~/Studies ~/Run
cp ~/Code/WormWeaver/templates/studies/unistra-hpc.yaml ~/Studies/my_study.yaml
```

Edit the [Unistra HPC study template](../templates/studies/unistra-hpc.yaml) for your project and paths. It uses WormWeaver's bundled simulation config automatically. For a custom config, place it beside your study YAML and uncomment `template: ./config.yaml`; the relative path is resolved from the study file's directory. The template generates 30 simulations (three beta values and ten replicas); the minimal example below generates ten.

## 4. Generate a study on the HPC

Save the following as `~/Studies/my_study.yaml`, adjusting all absolute paths for your account:

```yaml
project: RFEBHM
output_dir: /home2020/home/isis/tgupta/Run/my_study
multiwormqmc_path: /home2020/home/isis/tgupta/Code/multiwormqmc.jl
replicas: 10
seed_generation_seed: 12345

jobs:
  command:
    - julia
    - /home2020/home/isis/tgupta/Code/multiwormqmc.jl/bin/multiwormqmc.jl
  run_dir: /home2020/home/isis/tgupta/Run/my_study
  mode: scratch
  scratch_root: /scratch
```

Generate it from the shell with WormWeaver's environment explicitly selected:

```bash
module load gcc/15.2.0
julia --project="$HOME/Code/WormWeaver" -e 'using WormWeaver; generate_study(ARGS[1])' "$HOME/Studies/my_study.yaml"
```

Or use an interactive Julia session:

```bash
julia --project="$HOME/Code/WormWeaver"
```

```julia
using WormWeaver
study_dir = generate_study(joinpath(homedir(), "Studies", "my_study.yaml"))
```

This example creates ten simulation folders, both manifests, provenance, a jobfile, and the two account scripts under `~/Run/my_study`. Generation prepares files; it does not run or submit simulations. `multiwormqmc_path` records the HPC checkout's software provenance.

An existing output folder is rejected. To deliberately replace the entire study, including any previous simulation results:

```julia
generate_study(joinpath(homedir(), "Studies", "my_study.yaml"); overwrite=true)
```

If generating on your own computer instead, use a local `output_dir`, retain the HPC `jobs.run_dir` and command paths, and transfer the complete generated study to that HPC directory. For provenance, `multiwormqmc_path` must identify a checkout accessible on the generation machine.

Each simulation has its own `config.yaml`. The launcher selects the multiwormqmc.jl environment automatically and resolves relative output paths from the simulation's working directory.

## 5. Prepare the batch scripts

WormWeaver copies both account scripts into the study:

| Script | Partition | Account |
| --- | --- | --- |
| `grant_g2026a136c.slurm` | `grant` | `g2026a136c` |
| `grantgpu_g2026a136g.slurm` | `grantgpu` | `g2026a136g` |

Both bundled scripts already load GCC after all `#SBATCH` directives and before the launcher:

```bash
module load gcc/15.2.0
# Add any setup needed to make your Julia installation available here.
commandes=jobfile
hpc_multilauncher $commandes
```

Check the requested nodes, memory, and wall time against your allocation. The GPU-partition template only selects the supplied partition/account; it does not request a specific GPU resource or enable GPU computation.

Newly generated studies include this module command automatically. For an existing study, add it to its scripts manually or copy the updated templates; there is no need to overwrite simulation data. Regenerating a study with `overwrite=true` replaces the entire study, including script edits and results.

## 6. Submit and inspect results

From the study directory on the HPC, submit **one** account script:

```bash
cd /home2020/home/isis/tgupta/Run/my_study
sbatch grant_g2026a136c.slurm
```

Alternatively, choose `grantgpu_g2026a136g.slurm` when appropriate for that allocation. Both scripts use the same jobfile; submitting both would run the same simulations twice.

In scratch mode, each jobfile line copies a simulation to `/scratch/job.$SLURM_JOB_ID/simN`, runs it there, and copies its non-hidden contents back. The `/scratch/job.$SLURM_JOB_ID` directory must already exist in the job environment. Use relative output filenames in `config.yaml` to write results to scratch.

Inspect `simN/simulation.log` and the configured output files, normally `simN/results.h5`, after copy-back. The original-style jobfile has no added failure handling: interrupted jobs may leave results only on scratch, and the final shell status may reflect copy-back rather than simulation success. Use `mode: direct` to run in the study folders without scratch transfers.
