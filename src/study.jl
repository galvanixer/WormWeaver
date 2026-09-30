module Study

import YAML
import CSV
import JLD2
import DataFrames
import SHA
import UUIDs
import Random

include("config_writer.jl")
include("provenance.jl")
include("jobs.jl")

export generate_study

const DEFAULT_TEMPLATE = normpath(joinpath(@__DIR__, "..", "templates", "multiwormqmc", "config.yaml"))

const MANIFEST_METADATA = ("project", "study_id", "simulation", "replica", "config", "output_file", "config_sha256")
const CORE_FIELDS = [
    "system.beta", "system.ensemble", "lattice.geometry", "lattice.extent", "lattice.pbc",
    "model.family", "model.nspecies", "model.hopping.kind", "model.hopping.J",
    "model.hopping.alpha", "model.hopping.distance_cutoff", "model.onsite.U", "model.onsite.mu",
    "model.density_density.kind", "model.density_density.V", "model.density_density.alpha",
    "model.density_density.distance_cutoff", "model.disorder.kind", "model.disorder.strength",
    "model.disorder.distribution", "model.disorder.seed", "worldlines.initial_density",
    "worldlines.maxocc", "worldlines.init_seed", "run.seed", "run.blocks",
    "run.samples_per_block", "run.burnin_blocks", "sampler.multiworm_order",
]

function get_field(config, path)
    node = config
    for part in split(path, '.')
        node isa AbstractDict && haskey(node, part) || return missing
        node = node[part]
    end
    node
end

# Keep scalar columns directly usable; encode structured values as YAML cells.
manifest_value(value) = value === nothing || ismissing(value) ? missing :
    value isa AbstractDict || value isa AbstractVector ? strip(YAML.write(value)) : value

function result_path(config, simulation)
    value = get_field(config, "run.output_file")
    (ismissing(value) || value === nothing) && return missing
    value isa AbstractString && !isempty(strip(value)) ||
        throw(ArgumentError("run.output_file must be a nonempty path"))
    isabspath(value) ? normpath(value) : normpath(joinpath(simulation, value))
end

function check_destination(destination, overwrite)
    if ispath(destination) || islink(destination)
        overwrite || throw(ArgumentError("Output already exists: $destination; use overwrite=true to replace it"))
        !islink(destination) && isdir(destination) ||
            throw(ArgumentError("overwrite requires a directory, not a file or symlink: $destination"))
    end
end

function mapping(value, label)
    value isa AbstractDict || throw(ArgumentError("$label must be a mapping"))
    all(k -> k isa AbstractString, keys(value)) ||
        throw(ArgumentError("$label must have string keys"))
    value
end

function set_field!(config, path, value)
    parts = split(path, '.')
    any(isempty, parts) && throw(ArgumentError("Invalid field path: $path"))
    node = config
    for part in parts[1:end-1]
        node isa AbstractDict && haskey(node, part) ||
            throw(ArgumentError("Unknown field path: $path"))
        node = node[part]
    end
    node isa AbstractDict && haskey(node, last(parts)) ||
        throw(ArgumentError("Unknown field path: $path"))
    old = node[last(parts)]
    # Null-valued template fields deliberately allow a user-specified type.
    compatible = old === nothing || value === nothing ||
        (old isa Bool ? value isa Bool :
         old isa Integer ? value isa Integer && !(value isa Bool) :
         old isa Number ? value isa Number && !(value isa Bool) :
         old isa AbstractString ? value isa AbstractString :
         old isa AbstractVector ? value isa AbstractVector :
         old isa AbstractDict ? value isa AbstractDict : false)
    compatible || throw(ArgumentError("Incompatible value type for $path"))
    node[last(parts)] = deepcopy(value)
end

"""
    generate_study(path; overwrite=false, validate_config=nothing) -> String

Generate `sim1/config.yaml`, etc. from a study YAML file. Relative paths are
resolved against the study file. A nonblank `project` string is required and
stored in the study snapshot and manifest, not in simulation configs.
Writes a native DataFrame under the `manifest` key in `manifest.jld2`, plus
`manifest.csv` with structured values serialized as YAML cells. Each study also
includes `jobfile` and all bundled `.slurm` launchers.
Sweep paths are sorted alphabetically; the
last path varies fastest. `replicas` repeats each combination, randomizing
`run.seed` by default; `randomize_seeds` selects which of the three seed fields
to randomize. Without these options, seeds are preserved unless explicitly overridden.
An optional `validate_config(config)` callback validates each configuration
dictionary before any output is written; it must throw on invalid input.
Returns the absolute output directory. Set `overwrite=true` to replace an existing
directory and all its contents after generation succeeds. Files and symlinks are rejected.
"""
function generate_study(path::AbstractString; overwrite::Bool=false, validate_config=nothing)
    source = abspath(path)
    study_text = read(source, String)
    study = mapping(YAML.load(study_text), "Study")
    allowed = ("project", "template", "output_dir", "overrides", "sweep", "replicas", "randomize_seeds", "seed_generation_seed", "multiwormqmc_path", "jobs")
    for key in keys(study)
        key in allowed || throw(ArgumentError("Unknown study option: $key"))
    end
    project = get(study, "project", nothing)
    project isa AbstractString && !isempty(strip(project)) ||
        throw(ArgumentError("project must be a nonempty string"))
    output = get(study, "output_dir", nothing)
    output isa AbstractString && !isempty(strip(output)) ||
        throw(ArgumentError("output_dir must be a nonempty path"))
    resolve(p) = abspath(joinpath(dirname(source), p))
    destination = resolve(output)
    check_destination(destination, overwrite)
    template_path = get(study, "template", DEFAULT_TEMPLATE)
    template_path isa AbstractString && !isempty(strip(template_path)) ||
        throw(ArgumentError("template must be a nonempty filename or path"))
    if !isabspath(template_path) && !occursin('/', template_path) && !occursin('\\', template_path)
        template_dir = dirname(DEFAULT_TEMPLATE)
        available = sort(filter(name -> isfile(joinpath(template_dir, name)), readdir(template_dir)))
        template_path in available || throw(ArgumentError(
            "Unknown bundled template: $template_path; available templates: $(join(available, ", ")). Use ./filename.yaml for a local file."))
        template_path = joinpath(template_dir, template_path)
    else
        template_path = resolve(template_path)
    end
    template_text = read(template_path, String)
    template = mapping(YAML.load(template_text), "Template")
    overrides = mapping(get(study, "overrides", Dict()), "overrides")
    sweep = mapping(get(study, "sweep", Dict()), "sweep")
    fields = sort!(collect(union(keys(overrides), keys(sweep))))
    for field in fields
        field in MANIFEST_METADATA &&
            throw(ArgumentError("Field path conflicts with manifest metadata: $field"))
    end
    for (i, field) in enumerate(fields), other in fields[i+1:end]
        startswith(other, field * ".") &&
            throw(ArgumentError("Overlapping field paths: $field and $other"))
    end
    isempty(intersect(keys(overrides), keys(sweep))) ||
        throw(ArgumentError("A field cannot appear in both overrides and sweep"))
    replicas = get(study, "replicas", 1)
    replicas isa Integer && !(replicas isa Bool) && 1 <= replicas <= typemax(Int) ||
        throw(ArgumentError("replicas must be a positive integer"))
    haskey(study, "randomize_seeds") && !haskey(study, "replicas") &&
        throw(ArgumentError("randomize_seeds requires replicas"))
    seed_fields = get(study, "randomize_seeds", haskey(study, "replicas") ? ["run.seed"] : String[])
    valid_seeds = ("run.seed", "worldlines.init_seed", "model.disorder.seed")
    seed_fields isa AbstractVector && all(f -> f isa AbstractString && f in valid_seeds, seed_fields) ||
        throw(ArgumentError("randomize_seeds must be a list of supported seed paths"))
    haskey(study, "replicas") && isempty(seed_fields) &&
        throw(ArgumentError("randomize_seeds must select at least one seed field"))
    length(unique(seed_fields)) == length(seed_fields) ||
        throw(ArgumentError("randomize_seeds contains duplicate fields"))
    for seed_field in seed_fields
        for field in fields
            (field == seed_field || startswith(seed_field, field * ".") || startswith(field, seed_field * ".")) &&
                throw(ArgumentError("Randomized seed $seed_field conflicts with override or sweep $field"))
        end
        set_field!(deepcopy(template), seed_field, 1)
    end
    # A private RNG avoids changing the caller's default RNG state.
    generation_seed = get(study, "seed_generation_seed", nothing)
    if haskey(study, "seed_generation_seed")
        haskey(study, "replicas") || throw(ArgumentError("seed_generation_seed requires replicas"))
        generation_seed isa Integer && !(generation_seed isa Bool) && 0 <= generation_seed <= typemax(Int) ||
            throw(ArgumentError("seed_generation_seed must be a nonnegative integer fitting Int"))
    end
    seed_rng = generation_seed === nothing ? Random.Xoshiro() : Random.Xoshiro(generation_seed)
    used_seeds = Set{Int}()
    function fresh_seed()
        while true
            seed = rand(seed_rng, 1:typemax(Int))
            if !(seed in used_seeds)
                push!(used_seeds, seed)
                return seed
            end
        end
    end
    base = deepcopy(template)
    for field in sort!(collect(keys(overrides)))
        set_field!(base, field, overrides[field])
    end
    sweep_fields = sort!(collect(keys(sweep)))
    for field in sweep_fields
        values = sweep[field]
        values isa AbstractVector && !isempty(values) ||
            throw(ArgumentError("Sweep $field must be a nonempty list of values"))
        for value in values
            set_field!(deepcopy(base), field, value)
        end
    end

    configs = Any[]
    replica_numbers = Int[]
    function expand(config, index)
        if index > length(sweep_fields)
            for replica in 1:replicas
                instance = deepcopy(config)
                for field in seed_fields
                    set_field!(instance, field, fresh_seed())
                end
                validate_config === nothing || validate_config(deepcopy(instance))
                push!(configs, instance)
                push!(replica_numbers, replica)
            end
            return
        end
        field = sweep_fields[index]
        for value in sweep[field]
            next = deepcopy(config)
            set_field!(next, field, value)
            expand(next, index + 1)
        end
    end
    expand(base, 1)
    study_id = string(UUIDs.uuid4())
    simulations = ["sim$i" for i in eachindex(configs)]
    config_texts = [config_yaml(config) for config in configs]
    manifest = DataFrames.DataFrame(
        project=fill(project, length(configs)),
        study_id=fill(study_id, length(configs)),
        simulation=simulations,
        replica=replica_numbers,
        config=["$simulation/config.yaml" for simulation in simulations],
        output_file=[result_path(config, simulation) for (config, simulation) in zip(configs, simulations)],
        config_sha256=[bytes2hex(SHA.sha256(text)) for text in config_texts],
    )
    for field in unique(vcat(CORE_FIELDS, fields))
        manifest[!, field] = [deepcopy(get_field(config, field)) for config in configs]
    end
    csv_manifest = DataFrames.DataFrame()
    for field in names(manifest)
        csv_manifest[!, field] = manifest_value.(manifest[!, field])
    end

    jobfile = jobfile_text(study, destination, simulations)
    provenance = software_provenance(study, source)
    provenance["study_id"] = study_id

    # Stage the complete study before publishing it to the requested location.
    mkpath(dirname(destination))
    staging = mktempdir(dirname(destination); prefix=".wormweaver-")
    try
        write(joinpath(staging, "study.yaml"), study_text)
        write(joinpath(staging, "jobfile"), jobfile)
        slurm_dir = joinpath(@__DIR__, "..", "templates", "slurm")
        for name in sort(readdir(slurm_dir))
            template = joinpath(slurm_dir, name)
            endswith(name, ".slurm") && isfile(template) || continue
            cp(template, joinpath(staging, name))
        end
        YAML.write_file(joinpath(staging, "provenance.yaml"), provenance)
        write(joinpath(staging, "template.yaml"), template_text)
        for (i, config) in enumerate(configs)
            folder = joinpath(staging, "sim$i")
            mkdir(folder)
            write(joinpath(folder, "config.yaml"), config_texts[i])
        end
        JLD2.jldsave(joinpath(staging, "manifest.jld2"); manifest)
        CSV.write(joinpath(staging, "manifest.csv"), csv_manifest)
        check_destination(destination, overwrite)
        if isdir(destination)
            # Retain the old directory until publication succeeds; restore on error.
            backup_root = mktempdir(dirname(destination); prefix=".wormweaver-backup-")
            backup = joinpath(backup_root, "previous")
            try
                mv(destination, backup; force=false)
                try
                    mv(staging, destination; force=false)
                catch
                    mv(backup, destination; force=false)
                    rethrow()
                end
                rm(backup; recursive=true)
            finally
                # If recovery failed, leave the backup available for recovery.
                isempty(readdir(backup_root)) && rm(backup_root)
            end
        else
            mv(staging, destination; force=false)
        end
    finally
        isdir(staging) && rm(staging; recursive=true)
    end
    combinations = div(length(configs), replicas)
    jobs = get(study, "jobs", Dict())
    mode = get(jobs, "mode", "direct")
    run_dir = get(jobs, "run_dir", destination)
    println("Created $(length(configs)) simulations: $combinations parameter combinations × $replicas replicas")
    println("Directory: $destination")
    println("Job mode: $mode")
    println("\nSubmit from the study directory on the execution machine:")
    println("  cd ", shell_quote(run_dir))
    for name in sort(readdir(destination))
        endswith(name, ".slurm") && isfile(joinpath(destination, name)) || continue
        println("  sbatch ", name)
    end
    println("Choose one account script. No jobs have been submitted.")
    destination
end

end
