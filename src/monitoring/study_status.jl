import CSV
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
    isfile(path) || return (; status, completed_blocks, failed)

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
`simulation`, `lattice_size`, `beta`, `status`, `completed_blocks`, and `failed`,
in manifest order. Lattice sizes remain arrays, including for a 1D lattice.

Statuses are `Pre Tuning`, `Tuning`, `Thermalizing`, `Production Started`,
`Checkpoint 1` through `Checkpoint 4`, and `Finished`. Checkpoints require a
saved-checkpoint marker or restored checkpoint evidence; `Finished` requires a
success summary. A failure retains its last stage and sets `failed=true`.
For appended logs, only the latest run is reported.

Missing logs are `Pre Tuning`; missing lattice size or beta is `missing`.
Unreadable logs raise an error. This read-only snapshot does not determine
whether a process is currently running; `failed=false` only means no explicit
failure was found in the log. Results and checkpoint files are not opened.
"""
function study_status(root::AbstractString)
    root = abspath(root)
    native = joinpath(root, "manifest.jld2")
    csv = joinpath(root, "manifest.csv")
    manifest = if isfile(native)
        JLD2.load(native, "manifest")
    elseif isfile(csv)
        CSV.read(csv, DataFrames.DataFrame)
    else
        throw(ArgumentError("no manifest.jld2 or manifest.csv found in $root"))
    end
    return _study_status(root, manifest)
end

# Reuse the manifest snapshot already loaded by aggregation.
function _study_status(root::AbstractString, manifest)
    manifest isa DataFrames.AbstractDataFrame ||
        throw(ArgumentError("manifest must be a DataFrame"))
    all(name -> name in names(manifest), ("simulation", "config")) ||
        throw(ArgumentError("manifest must contain simulation and config"))

    table = DataFrames.DataFrame(simulation=String[], lattice_size=Any[], beta=Any[],
                                status=String[], completed_blocks=Int[], failed=Bool[])
    for row in eachrow(manifest)
        row.simulation isa AbstractString && !isempty(strip(row.simulation)) &&
            row.config isa AbstractString &&
            !isempty(strip(row.config)) ||
            throw(ArgumentError("each manifest row must have a simulation name and config path"))
        folder = dirname(joinpath(root, row.config))
        progress = _study_log_status(joinpath(folder, "simulation.log"))
        push!(table, (; simulation=String(row.simulation),
                      lattice_size=_status_metadata(row, "lattice.extent"; yaml=true),
                      beta=_status_metadata(row, "system.beta"), progress...))
    end
    return table
end

"""
    status_summary(df::AbstractDataFrame) -> DataFrame

Summarize a `study_status` table by `lattice_size`, `beta`, `status`, and `failed`,
with a count column named `n_simulations`. Failed runs remain separate from
other runs at the same stage. Only observed combinations are returned, in order
of first appearance; missing lattice sizes or beta values are retained.
The input table is not modified.
"""
function status_summary(df::DataFrames.AbstractDataFrame)
    columns = [:lattice_size, :beta, :status, :failed]
    all(column -> hasproperty(df, column), columns) ||
        throw(ArgumentError("status table must contain lattice_size, beta, status, and failed"))
    return DataFrames.combine(DataFrames.groupby(df, columns; sort=false),
                              DataFrames.nrow => :n_simulations)
end
