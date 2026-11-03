import DataFrames
import JLD2

"""
    slurm_array_tasks(id; scratch_root="/scratch") -> DataFrame

Query `squeue` for all currently queued tasks in a Slurm array. Accept an array
ID such as `18523211` or a task ID such as `"18523211_2"`; either queries the
entire array. Run on an HPC host with `squeue` available and access to Slurm.

Return `task`, `index`, numeric allocation `jobid`, `state`, allocated `nodes`,
and expected `scratch_dir`, sorted by array index. Node lists retain Slurm's
compact notation; an empty node list is `missing`. Scratch directories follow
`scratch_root/job.<jobid>` and are not checked for existence.

Completed tasks no longer visible in `squeue` are excluded. Pending tasks may
have no nodes, and their allocation details can change. Slurm command failures
raise an error. No simulation files are read or modified.
"""
function slurm_array_tasks(id; scratch_root::AbstractString="/scratch")
    m = match(r"^(\d+)(?:_\d+)?$", strip(string(id)))
    isnothing(m) &&
        throw(ArgumentError("Expected an array ID or array_task ID"))
    arrayid = m.captures[1]

    output = read(`squeue -h -r -j $arrayid -o "%i|%A|%T|%N"`, String)
    table = DataFrames.DataFrame(
        task=String[], index=Int[], jobid=Int[], state=String[],
        nodes=Union{Missing,String}[], scratch_dir=String[])

    for line in split(output, '\n'; keepempty=false)
        task, jobid_text, state, nodes = strip.(split(line, '|'; limit=4))
        task_match = match(r"^(\d+)_(\d+)$", task)
        isnothing(task_match) && continue
        task_match.captures[1] == arrayid || continue

        index = parse(Int, task_match.captures[2])
        jobid = parse(Int, jobid_text)
        push!(table, (
            task=String(task), index=index, jobid=jobid, state=String(state),
            nodes=isempty(nodes) ? missing : String(nodes),
            scratch_dir=joinpath(scratch_root, "job.$jobid")))
    end

    return sort!(table, :index)
end

"""
    find_scratch_logs(study_dir, arrayid; scratch_root="/scratch", manifest=nothing) -> DataFrame

Read `manifest.jld2` and locate its simulations'
`simulation.log` files on running tasks of the given Slurm array. A task ID such
as `"18523211_2"` also inspects the whole array. Run on the HPC login node.

Return `simulation`, `task`, numeric `jobid`, `node`, and `log_path`, in manifest
simulation order. Keep every matching location when a simulation appears in
multiple allocations or nodes. Simulation names must be single directory
components; repeated manifest names are checked once.

Use one read-only `srun` check per running allocation and node. Pending tasks
are skipped; Slurm/inspection failures raise errors. Missing logs are omitted.
The function checks file existence only: it does not read log contents, parse
status, copy files, or change simulation execution. Completed allocations no
longer returned by Slurm are not inspected.
An already-loaded JLD2 manifest can be supplied through `manifest` to reuse
the same snapshot as `study_status`. CSV manifests are not read.
"""
function find_scratch_logs(study_dir::AbstractString, arrayid;
                           scratch_root::AbstractString="/scratch", manifest=nothing)
    isabspath(scratch_root) && !occursin('\0', scratch_root) ||
        throw(ArgumentError("scratch_root must be an absolute path without NUL bytes"))
    if manifest === nothing
        path = joinpath(abspath(study_dir), "manifest.jld2")
        isfile(path) || throw(ArgumentError("no manifest.jld2 found in $(abspath(study_dir))"))
        manifest = JLD2.load(path, "manifest")
    end
    manifest isa DataFrames.AbstractDataFrame && hasproperty(manifest, :simulation) ||
        throw(ArgumentError("manifest must be a DataFrame containing simulation"))
    all(name -> name isa AbstractString && !isempty(strip(name)) &&
                name ∉ (".", "..") && !occursin(r"[/\x00]", name), manifest.simulation) ||
        throw(ArgumentError("manifest simulation names must be single directory components"))
    simulations = unique(String.(manifest.simulation))
    tasks = slurm_array_tasks(arrayid; scratch_root)
    table = DataFrames.DataFrame(simulation=String[], task=String[], jobid=Int[],
                                 node=String[], log_path=String[])
    isempty(simulations) && return table

    # Return manifest indices rather than paths, so whitespace and shell
    # characters in names cannot affect the output format or shell source.
    script = raw"""
    scratch_root=$1
    job_dir=$2
    shift 2
    if [[ ! -d "$scratch_root" || ! -x "$scratch_root" ]]; then
        echo "Cannot inspect scratch root: $scratch_root" >&2
        exit 1
    fi
    if [[ -d "$job_dir" && ! -x "$job_dir" ]]; then
        echo "Cannot inspect scratch directory: $job_dir" >&2
        exit 1
    fi
    i=0
    for simulation do
        i=$((i + 1))
        folder="$job_dir/$simulation"
        if [[ -d "$folder" && ! -x "$folder" ]]; then
            echo "Cannot inspect simulation directory: $folder" >&2
            exit 1
        fi
        if [[ -f "$folder/simulation.log" ]]; then
            printf '%s\n' "$i"
        fi
    done
    """
    for task in eachrow(tasks)
        task.state == "RUNNING" || continue
        ismissing(task.nodes) && error("Running task $(task.task) has no allocated nodes")
        nodes = occursin(r"[\[\],]", task.nodes) ?
            String.(split(read(`scontrol show hostnames $(task.nodes)`, String))) : [task.nodes]
        isempty(nodes) && error("No nodes resolved for task $(task.task)")
        for node in unique(nodes)
            command = `srun --jobid=$(task.jobid) --nodelist=$node --overlap --nodes=1 --ntasks=1 --cpus-per-task=1 --quiet --chdir=/ bash -c $script wormweaver-log-check $scratch_root $(task.scratch_dir) $simulations`
            indices = split(read(command, String))
            for value in indices
                i = tryparse(Int, value)
                i !== nothing && 1 <= i <= length(simulations) ||
                    error("Invalid scratch-log response from $(task.task) on $node")
                simulation = simulations[i]
                push!(table, (; simulation, task=task.task, jobid=task.jobid, node,
                                log_path=joinpath(task.scratch_dir, simulation, "simulation.log")))
            end
        end
    end
    order = Dict(name => i for (i, name) in enumerate(simulations))
    return table[sortperm([order[name] for name in table.simulation]), :]
end

"""
    read_scratch_log(row) -> String

Read the full log identified by a row from `find_scratch_logs`, using its
`jobid`, `node`, and `log_path`. Run on the HPC login node. A named tuple with
these fields is also accepted.

Use a read-only `srun` step inside the allocation and return the log text
unchanged. No files are copied or modified and no status is parsed. Slurm or
file-read failures raise an error, including when a log disappears or its
allocation ends after discovery.
"""
function read_scratch_log(row)
    all(name -> hasproperty(row, name), (:jobid, :node, :log_path)) ||
        throw(ArgumentError("row must contain jobid, node, and log_path"))
    jobid, node, path = row.jobid, row.node, row.log_path
    jobid isa Integer && !(jobid isa Bool) && jobid > 0 ||
        throw(ArgumentError("jobid must be a positive integer"))
    node isa AbstractString && !isempty(strip(node)) && !occursin('\0', node) ||
        throw(ArgumentError("node must be a nonblank string without NUL bytes"))
    path isa AbstractString && isabspath(path) && !occursin('\0', path) ||
        throw(ArgumentError("log_path must be an absolute path without NUL bytes"))
    return read(`srun --jobid=$jobid --nodelist=$node --overlap --nodes=1 --ntasks=1 --cpus-per-task=1 --quiet --chdir=/ cat -- $path`, String)
end
