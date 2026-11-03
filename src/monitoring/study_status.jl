import Dates
import DataFrames
import JLD2
import YAML

"""
    study_log_status(path::AbstractString)
    study_log_status(io::IO)

Parse a simulation log and return `(status, completed_blocks, failed)` as a named
tuple. Accept a file path or an IO stream such as `IOBuffer(text)` containing a
scratch log. Streams are read from their current position and left open.

A missing file returns `status="No Log File"`, zero blocks, and `failed=false`.
An empty log is `Pre Tuning`. For appended logs, report only the latest attempt.
Checkpoint stages require saved/restored checkpoint evidence; `Finished`
requires a success summary. A failure retains its stage and sets `failed=true`.
Unreadable files raise an error. No files are modified.
"""
function study_log_status(path::AbstractString)
    try
        (ispath(path) || islink(path)) && !isfile(path) &&
            throw(ArgumentError("log is not a regular file: $path"))
        return open(study_log_status, path)
    catch err
        # Opening distinguishes a missing log from an inaccessible parent.
        if err isa SystemError && err.errnum == Base.Libc.ENOENT && !islink(path)
            return (; status="No Log File", completed_blocks=0, failed=false)
        end
        rethrow()
    end
end

# Parse from the stream's current position; ownership stays with the caller.
function study_log_status(io::IO)
    status, completed_blocks, failed = "Pre Tuning", 0, false
    production, target = false, 0

    for line in eachline(io)
        line = strip(line)
        if startswith(line, "MultiWormQMC v")
            # Appended logs may contain earlier successful or failed attempts.
            status, completed_blocks, failed = "Pre Tuning", 0, false
            production, target = false, 0
        elseif line == "TUNING"
            status, production = "Tuning", false
        elseif line == "THERMALIZATION"
            status, production = "Thermalizing", false
        elseif line == "PRODUCTION"
            status, production = "Production Started", true
        elseif line in ("FINALIZATION", "SUMMARY")
            production = false
        elseif occursin(r"^(Simulation failed:|ERROR:)", line)
            failed = true
        elseif occursin(r"^Status\s+SUCCESS$", line)
            status = "Finished"
        elseif production
            checkpoint = false
            if (m = match(r"^Blocks\s+([\d,]+)$", line)) !== nothing
                target = parse(Int, replace(m[1], "," => ""))
            elseif (m = match(r"^Restored blocks\s+([\d,]+)$", line)) !== nothing
                completed_blocks = parse(Int, replace(m[1], "," => ""))
            elseif startswith(line, "Checkpoints") &&
                   !occursin("Disabled", line) && target > 0
                # A resumed log can start after previously saved quarters.
                checkpoint = true
            elseif (m = match(r"^(\d+)/(\d+)\s+\d+%\s", line)) !== nothing
                completed_blocks, target = parse.(Int, m.captures)
                checkpoint = occursin('◆', line)
            end
            if checkpoint
                quarter = count(q -> cld(widemul(target, q), 4) <= completed_blocks, 1:4)
                status = quarter == 0 ? "Production Started" : "Checkpoint $quarter"
            end
        end
    end
    return (; status, completed_blocks, failed)
end

"""
    study_status(root; arrayid=nothing, scratch_root="/scratch", manifest=nothing) -> DataFrame

Read `manifest.jld2` and scan each simulation's
`simulation.log`, located beside its manifest `config` path. Return columns
`study_id`, `simulation`, `lattice_size`, `beta`, `status`, `completed_blocks`,
`failed`, and `log_source`, in manifest order. The manifest's `study_id` and `simulation`
together identify a simulation across generated studies. Lattice sizes remain
arrays, including for a 1D lattice.

Statuses are `No Log File`, `Pre Tuning`, `Tuning`, `Thermalizing`, `Production Started`,
`Checkpoint 1` through `Checkpoint 4`, and `Finished`. Checkpoints require a
saved-checkpoint marker or restored checkpoint evidence; `Finished` requires a
success summary. A failure retains its last stage and sets `failed=true`.
For appended logs, only the latest run is reported.

Missing logs are `No Log File`. `Pre Tuning` means the simulation has started
logging but tuning has not started. Absent or null study ID, lattice size, or
beta values are `missing`.
Unreadable logs raise an error. This read-only snapshot does not determine
whether a process is currently running; `failed=false` only means no explicit
failure was found in the log. Results and checkpoint files are not opened.

With `arrayid`, use `find_scratch_logs` to locate logs on running allocations
of that array, preferring a unique scratch log for each simulation. Otherwise
read the study-directory log. Multiple scratch matches and read failures raise
errors. A task ID such as `"18523211_2"` inspects the whole array. Without
`arrayid`, only study-directory logs are read and Slurm is not called.

`log_source` is `"scratch"` for a log read from a compute node, `"local"` for
a study-directory log, or `"none"` when no log was found. A local fallback
describes that log's progress and is not proof of the current scratch attempt.
Callers that already loaded the JLD2 manifest can pass it through `manifest`
to reuse that snapshot. CSV manifests are not read.
"""
function study_status(root::AbstractString; arrayid=nothing,
                      scratch_root::AbstractString="/scratch", manifest=nothing)
    root = abspath(root)
    if manifest === nothing
        path = joinpath(root, "manifest.jld2")
        isfile(path) || throw(ArgumentError("no manifest.jld2 found in $root"))
        manifest = JLD2.load(path, "manifest")
    end
    manifest isa DataFrames.AbstractDataFrame ||
        throw(ArgumentError("manifest must be a DataFrame"))
    all(name -> name in names(manifest), ("simulation", "config")) ||
        throw(ArgumentError("manifest must contain simulation and config"))

    for row in eachrow(manifest)
        row.simulation isa AbstractString && !isempty(strip(row.simulation)) &&
            row.config isa AbstractString &&
            !isempty(strip(row.config)) ||
            throw(ArgumentError("each manifest row must have a simulation name and config path"))
    end

    scratch_logs = Dict{String,NamedTuple}()
    if arrayid !== nothing
        for log in eachrow(find_scratch_logs(root, arrayid; scratch_root, manifest))
            if haskey(scratch_logs, log.simulation)
                previous = scratch_logs[log.simulation]
                throw(ArgumentError("multiple scratch logs for $(log.simulation): " *
                                    "$(previous.node):$(previous.log_path), $(log.node):$(log.log_path)"))
            end
            scratch_logs[log.simulation] = NamedTuple(log)
        end
    end

    table = DataFrames.DataFrame(study_id=Union{Missing,String}[], simulation=String[],
                                lattice_size=Any[], beta=Any[], status=String[],
                                completed_blocks=Int[], failed=Bool[], log_source=String[])
    for row in eachrow(manifest)
        folder = dirname(joinpath(root, row.config))
        if haskey(scratch_logs, row.simulation)
            text = read_scratch_log(scratch_logs[row.simulation])
            progress = study_log_status(IOBuffer(text))
            log_source = "scratch"
        else
            progress = study_log_status(joinpath(folder, "simulation.log"))
            log_source = progress.status == "No Log File" ? "none" : "local"
        end
        study_id = hasproperty(row, :study_id) ? something(row.study_id, missing) : missing
        lattice_size = hasproperty(row, Symbol("lattice.extent")) ?
            something(row["lattice.extent"], missing) : missing
        beta = hasproperty(row, Symbol("system.beta")) ?
            something(row["system.beta"], missing) : missing
        if lattice_size isa AbstractString
            lattice_size = YAML.load(String(lattice_size))
        end
        push!(table, (; study_id, simulation=String(row.simulation), lattice_size, beta,
                      progress..., log_source))
    end
    return table
end

"""
    status_manifest(root; save=false) -> DataFrame

Read the study manifest and current simulation logs, then return all manifest
columns with refreshed `status`, `completed_blocks`, and `failed` columns.
Join on `(study_id, simulation)`, requiring nonblank study IDs and unique keys,
and preserve manifest row order. Read only `manifest.jld2`.

By default, no files are written. With `save=true`, also write
`manifest_status_yyyymmdd_HHMMSS_sssZ.jld2` in the study directory. The timestamp
is UTC. Store the joined DataFrame under `manifest` and the capture-start time
under `captured_at_utc`. The return value is the DataFrame in either mode.
Existing files are never overwritten; a timestamp collision raises an error.
Source manifests, configs, logs, and simulation outputs are not modified.
This is a log-progress snapshot, not a running-process or failure-cause check.
"""
function status_manifest(root::AbstractString; save::Bool=false)
    root = abspath(root)
    captured_at_utc = Dates.now(Dates.UTC)
    path = joinpath(root, "manifest.jld2")
    isfile(path) || throw(ArgumentError("no manifest.jld2 found in $root"))
    manifest = JLD2.load(path, "manifest")
    statuses = study_status(root; manifest)
    hasproperty(manifest, :study_id) &&
        all(id -> id isa AbstractString && !isempty(strip(id)), manifest.study_id) ||
        throw(ArgumentError("manifest must contain a nonblank study_id for every simulation"))

    keys = [:study_id, :simulation]
    progress = [:status, :completed_blocks, :failed]
    # Refresh any previous status columns without duplicating parameter aliases.
    base = DataFrames.select(manifest, DataFrames.Not(intersect(names(manifest), string.(progress))))
    table = DataFrames.leftjoin(base, DataFrames.select(statuses, [keys; progress]);
                                on=keys, validate=(true, true), order=:left)
    if save
        timestamp = Dates.format(captured_at_utc, "yyyymmdd_HHMMSS_sss")
        destination = joinpath(root, "manifest_status_$(timestamp)Z.jld2")
        (ispath(destination) || islink(destination)) &&
            throw(ArgumentError("status snapshot already exists: $destination; retry to use a new timestamp"))
        JLD2.jldsave(destination; manifest=table, captured_at_utc)
    end
    return table
end

"""
    status_summary(df::AbstractDataFrame; io=stdout) -> DataFrame

Summarize a `study_status` table by `lattice_size`, `beta`, `status`, and `failed`,
with a count column named `n_simulations`. Failed runs remain separate from
other runs at the same stage. Only observed combinations are returned, in order
of first appearance; missing lattice sizes or beta values are retained.
Print the total number of simulations, successful finishes, failures, and runs
left before returning the grouped table. Successful runs have `status="Finished"`
and `failed=false`; failures take precedence. Runs left are neither successful
nor failed, regardless of whether their processes are currently running.
Use `io` to redirect this output. The input table is not modified.
"""
function status_summary(df::DataFrames.AbstractDataFrame; io::IO=stdout)
    columns = [:lattice_size, :beta, :status, :failed]
    all(column -> hasproperty(df, column), columns) ||
        throw(ArgumentError("status table must contain lattice_size, beta, status, and failed"))
    summary = DataFrames.combine(DataFrames.groupby(df, columns; sort=false),
                                 DataFrames.nrow => :n_simulations)
    total = DataFrames.nrow(df)
    failed = count(identity, df.failed)
    finished = count(row -> row.status == "Finished" && !row.failed, eachrow(df))
    left = total - finished - failed
    println()
    println(io, "Total simulations: ", total)
    println(io, "Finished successfully: ", finished)
    println(io, "Failed: ", failed)
    println(io, "Left: ", left)
    println(io)
    return summary
end
