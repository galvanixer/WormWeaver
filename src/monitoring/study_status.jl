import CSV
import Dates
import DataFrames
import JLD2
import YAML

_status_count(value) = parse(Int, replace(value, "," => ""))

function _checkpoint_status(blocks, target)
    quarter = count(q -> cld(widemul(target, q), 4) <= blocks, 1:4)
    return quarter == 0 ? "Production Started" : "Checkpoint $quarter"
end

function _study_log_status(path::AbstractString)
    status, completed_blocks, failed = "Pre Tuning", 0, false
    production, target = false, 0
    isfile(path) || return (; status="No Log File", completed_blocks, failed)

    for line in eachline(path)
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
            if (m = match(r"^Blocks\s+([\d,]+)$", line)) !== nothing
                target = _status_count(m[1])
            elseif (m = match(r"^Restored blocks\s+([\d,]+)$", line)) !== nothing
                completed_blocks = _status_count(m[1])
            elseif startswith(line, "Checkpoints") &&
                   !occursin("Disabled", line) && target > 0
                # A resumed log can start after previously saved quarters.
                status = _checkpoint_status(completed_blocks, target)
            elseif (m = match(r"^(\d+)/(\d+)\s+\d+%\s", line)) !== nothing
                completed_blocks, target = parse.(Int, m.captures)
                if occursin('◆', line)
                    status = _checkpoint_status(completed_blocks, target)
                end
            end
        end
    end
    return (; status, completed_blocks, failed)
end

function _status_metadata(row, name; yaml=false)
    hasproperty(row, Symbol(name)) || return missing
    value = row[name]
    value === nothing && return missing
    return yaml && value isa AbstractString ? YAML.load(String(value)) : value
end

"""
    study_status(root) -> DataFrame

Read `manifest.jld2` (or `manifest.csv` if absent) and scan each simulation's
`simulation.log`, located beside its manifest `config` path. Return columns
`study_id`, `simulation`, `lattice_size`, `beta`, `status`, `completed_blocks`,
and `failed`, in manifest order. The manifest's `study_id` and `simulation`
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
"""
function study_status(root::AbstractString)
    root = abspath(root)
    return _study_status(root, _read_study_manifest(root))
end

function _read_study_manifest(root::AbstractString)
    native = joinpath(root, "manifest.jld2")
    csv = joinpath(root, "manifest.csv")
    return if isfile(native)
        JLD2.load(native, "manifest")
    elseif isfile(csv)
        CSV.read(csv, DataFrames.DataFrame)
    else
        throw(ArgumentError("no manifest.jld2 or manifest.csv found in $root"))
    end
end

# Reuse the manifest snapshot already loaded by aggregation.
function _study_status(root::AbstractString, manifest)
    manifest isa DataFrames.AbstractDataFrame ||
        throw(ArgumentError("manifest must be a DataFrame"))
    all(name -> name in names(manifest), ("simulation", "config")) ||
        throw(ArgumentError("manifest must contain simulation and config"))

    table = DataFrames.DataFrame(study_id=Union{Missing,String}[], simulation=String[],
                                lattice_size=Any[], beta=Any[],
                                status=String[], completed_blocks=Int[], failed=Bool[])
    for row in eachrow(manifest)
        row.simulation isa AbstractString && !isempty(strip(row.simulation)) &&
            row.config isa AbstractString &&
            !isempty(strip(row.config)) ||
            throw(ArgumentError("each manifest row must have a simulation name and config path"))
        folder = dirname(joinpath(root, row.config))
        progress = _study_log_status(joinpath(folder, "simulation.log"))
        push!(table, (; study_id=_status_metadata(row, "study_id"),
                      simulation=String(row.simulation),
                      lattice_size=_status_metadata(row, "lattice.extent"; yaml=true),
                      beta=_status_metadata(row, "system.beta"), progress...))
    end
    return table
end

"""
    save_status_manifest(root) -> String

Read the study manifest and current simulation logs, then save all manifest
columns with refreshed `status`, `completed_blocks`, and `failed` columns.
Join on `(study_id, simulation)`, requiring nonblank study IDs and unique keys,
and preserve manifest row order. Prefer `manifest.jld2`, falling back to CSV.

Write `manifest_status_yyyymmdd_HHMMSS_sssZ.jld2` in the study directory and
return its absolute path. The timestamp is UTC. Store the joined DataFrame
under `manifest` and the capture-start time under `captured_at_utc`.
Existing files are never overwritten; a timestamp collision raises an error.
Source manifests, configs, logs, and simulation outputs are not modified.
This is a log-progress snapshot, not a running-process or failure-cause check.
"""
function save_status_manifest(root::AbstractString)
    root = abspath(root)
    captured_at_utc = Dates.now(Dates.UTC)
    manifest = _read_study_manifest(root)
    statuses = _study_status(root, manifest)
    hasproperty(manifest, :study_id) &&
        all(id -> id isa AbstractString && !isempty(strip(id)), manifest.study_id) ||
        throw(ArgumentError("manifest must contain a nonblank study_id for every simulation"))

    keys = [:study_id, :simulation]
    progress = [:status, :completed_blocks, :failed]
    # Refresh any previous status columns without duplicating parameter aliases.
    base = DataFrames.select(manifest, DataFrames.Not(intersect(names(manifest), string.(progress))))
    table = DataFrames.leftjoin(base, DataFrames.select(statuses, [keys; progress]);
                                on=keys, validate=(true, true), order=:left)
    timestamp = Dates.format(captured_at_utc, "yyyymmdd_HHMMSS_sss")
    destination = joinpath(root, "manifest_status_$(timestamp)Z.jld2")
    (ispath(destination) || islink(destination)) &&
        throw(ArgumentError("status snapshot already exists: $destination; retry to use a new timestamp"))
    JLD2.jldsave(destination; manifest=table, captured_at_utc)
    return destination
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
