# Each line is an independent shell command consumed by hpc_multilauncher.
shell_quote(value::AbstractString) = "'" * replace(value, "'" => "'\"'\"'") * "'"

function jobfile_text(study, destination, simulations)
    jobs = mapping(get(study, "jobs", Dict()), "jobs")
    for key in keys(jobs)
        key in ("command", "run_dir", "mode", "scratch_root") || throw(ArgumentError("Unknown jobs option: $key"))
    end
    mode = get(jobs, "mode", "direct")
    mode in ("direct", "scratch") || throw(ArgumentError("jobs.mode must be direct or scratch"))
    scratch_root = get(jobs, "scratch_root", "/scratch")
    scratch_root isa AbstractString && isabspath(scratch_root) ||
        throw(ArgumentError("jobs.scratch_root must be an absolute directory"))
    haskey(jobs, "scratch_root") && mode != "scratch" &&
        throw(ArgumentError("jobs.scratch_root requires scratch mode"))
    command = get(jobs, "command", ["multiwormqmc"])
    command isa AbstractVector && !isempty(command) &&
        all(arg -> arg isa AbstractString, command) && !isempty(strip(first(command))) ||
        throw(ArgumentError("jobs.command must be a nonempty list of strings with a nonblank executable"))
    run_dir = get(jobs, "run_dir", destination)
    run_dir isa AbstractString && isabspath(run_dir) ||
        throw(ArgumentError("jobs.run_dir must be an absolute execution directory"))
    for value in [command; run_dir; scratch_root]
        any(c -> c in ('\n', '\r', '\0'), value) &&
            throw(ArgumentError("Job paths and arguments cannot contain newlines or NUL bytes"))
    end
    invocation = join(shell_quote.(command), " ") * " config.yaml > simulation.log 2>&1"
    if mode == "direct"
        return join(["(cd " * shell_quote(joinpath(run_dir, simulation)) * " && " * invocation * ")\n"
                     for simulation in simulations])
    end
    scratch_invocation = join(command, " ") * " config.yaml &> simulation.log"
    join([scratch_jobline(scratch_invocation, joinpath(run_dir, simulation), scratch_root, simulation)
          for simulation in simulations])
end

function scratch_jobline(invocation, source, scratch_root, simulation)
    work = rstrip(scratch_root, '/') * "/job.\$SLURM_JOB_ID/" * simulation
    "cp -r " * source * " " * work * "; cd " * work * " ;" * invocation *
        " ; cp -r " * work * "/* " * source * " ;\n"
end
