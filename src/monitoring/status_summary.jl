import DataFrames

"""
    status_summary_data(df::AbstractDataFrame) -> NamedTuple

Calculate log-derived study progress without printing, file access, or Slurm
queries. Require `lattice_size`, `beta`, `status`, and `failed`; `log_source`
is optional. Failure takes precedence over `status="Finished"`.

Return `summary` (the original grouped status table), `overall` (counts named
`total`, `finished`, `failed`, `remaining`), `parameters` (those counts per
lattice/beta pair), `stages` (counts per lattice/beta/status/failure), and
`has_source`. Stage counts are `n_local`, `n_scratch`, `n_none`,
`n_unspecified`, and `n_simulations`. Missing or unrecognized sources are
unspecified; without a source column all stage counts are unspecified.

The grouped `summary` retains the original column types and first-appearance
order, including `log_source` when present. Input data and status strings are
not modified. Lattice arrays and missing metadata are retained in all tables.
"""
function status_summary_data(df::DataFrames.AbstractDataFrame)
    columns = [:lattice_size, :beta, :status, :failed]
    all(column -> hasproperty(df, column), columns) ||
        throw(ArgumentError("status table must contain lattice_size, beta, status, and failed"))
    has_source = hasproperty(df, :log_source)
    has_source && push!(columns, :log_source)
    summary = DataFrames.combine(DataFrames.groupby(df, columns; sort=false),
                                 DataFrames.nrow => :n_simulations)

    parameters = DataFrames.DataFrame(lattice_size=Any[], beta=Any[], total=Int[],
                                      finished=Int[], failed=Int[], remaining=Int[])
    for group in DataFrames.groupby(df, [:lattice_size, :beta]; sort=false)
        total = DataFrames.nrow(group)
        failed = count(identity, group.failed)
        finished = count(row -> row.status == "Finished" && !row.failed, eachrow(group))
        push!(parameters, (; lattice_size=group.lattice_size[1], beta=group.beta[1],
                            total, finished, failed, remaining=total - finished - failed))
    end
    overall = (; total=sum(parameters.total), finished=sum(parameters.finished),
                failed=sum(parameters.failed), remaining=sum(parameters.remaining))

    stages = DataFrames.DataFrame(lattice_size=Any[], beta=Any[], status=String[],
                                  failed=Bool[], n_local=Int[], n_scratch=Int[],
                                  n_none=Int[], n_unspecified=Int[], n_simulations=Int[])
    for group in DataFrames.groupby(df, [:lattice_size, :beta, :status, :failed]; sort=false)
        total = DataFrames.nrow(group)
        n_local = has_source ? count(isequal("local"), group.log_source) : 0
        n_scratch = has_source ? count(isequal("scratch"), group.log_source) : 0
        n_none = has_source ? count(isequal("none"), group.log_source) : 0
        push!(stages, (; lattice_size=group.lattice_size[1], beta=group.beta[1],
                        status=group.status[1], failed=group.failed[1], n_local, n_scratch,
                        n_none, n_unspecified=total - n_local - n_scratch - n_none,
                        n_simulations=total))
    end
    return (; summary, overall, parameters, stages, has_source)
end

"""
    status_summary(df::AbstractDataFrame; io=stdout, detailed=true) -> DataFrame

Print an adaptive dashboard of log-derived progress and return the original
grouped status DataFrame from `status_summary_data`. Always show overall counts
and percentages; show a parameter table only for multiple lattice/beta pairs.
With `detailed=true`, also print stages in logical order, combining local and
scratch counts and labeling failures with their last stage. `No log` and
`Unspecified` source columns appear only when needed. Legacy inputs without
`log_source` show stage totals alone.
Overall counts end with a divided Total row. Stage tables also end with source
totals; multiple parameter groups have a subtotal each and a final grand total.

Full lattice dimensions and missing metadata are preserved in the display.
Colors follow the IO context's `:color` setting, using Dark2_8 green for success
and orange for failure. All output goes to `io`; `io=devnull` suppresses it.
Use `status_summary(df);` to suppress automatic REPL display of the returned
table. Neither log progress nor log source establishes Slurm process state.
"""
function status_summary(df::DataFrames.AbstractDataFrame; io::IO=stdout, detailed::Bool=true)
    data = status_summary_data(df)
    overall, parameters, stages = data.overall, data.parameters, data.stages
    multiple = DataFrames.nrow(parameters) > 1
    colored = get(io, :color, false)
    green = colored ? "\e[38;2;27;158;119m" : ""
    orange = colored ? "\e[38;2;217;95;2m" : ""
    reset = colored ? "\e[0m" : ""

    # Prepare plain cells before padding/coloring. Both dashboard tables use
    # this same rendering loop, without a separate private printer layer.
    sections = NamedTuple[]
    if multiple
        headers = ["Lattice", "β", "Total", "Finished", "Failed", "Left"]
        rows = Vector{String}[]
        for row in eachrow(parameters)
            lattice = ismissing(row.lattice_size) ? "—" : string(row.lattice_size)
            beta = ismissing(row.beta) ? "—" : string(row.beta)
            push!(rows, [lattice, beta, string(row.total), string(row.finished),
                         string(row.failed), string(row.remaining)])
        end
        push!(rows, ["TOTAL", "", string(overall.total), string(overall.finished),
                     string(overall.failed), string(overall.remaining)])
        push!(sections, (; title="PROGRESS BY PARAMETERS", headers, rows, numeric_start=3,
                          row_colors=fill("", length(rows)),
                          column_colors=["", "", "", green, orange, ""],
                          breaks=Set{Int}(), dividers=Set([length(rows)])))
    end

    if detailed && overall.total > 0
        headers = multiple ? ["Lattice", "β", "Stage"] : ["Stage"]
        source_columns = Symbol[]
        if data.has_source
            append!(headers, ["Local", "Scratch"])
            append!(source_columns, [:n_local, :n_scratch])
            if any(>(0), stages.n_none)
                push!(headers, "No log")
                push!(source_columns, :n_none)
            end
            if any(>(0), stages.n_unspecified)
                push!(headers, "Unspecified")
                push!(source_columns, :n_unspecified)
            end
        end
        push!(headers, "Total")
        push!(source_columns, :n_simulations)
        rows, row_colors, breaks = Vector{String}[], String[], Set{Int}()
        dividers = Set{Int}()
        known = ["No Log File", "Pre Tuning", "Tuning", "Thermalizing", "Production Started",
                 "Checkpoint 1", "Checkpoint 2", "Checkpoint 3", "Checkpoint 4", "Finished"]
        for group in DataFrames.groupby(stages, [:lattice_size, :beta]; sort=false)
            order = Dict(stage => i for (i, stage) in enumerate(unique([known; group.status])))
            indices = sortperm(1:DataFrames.nrow(group); by=i -> (order[group.status[i]], group.failed[i]))
            !isempty(rows) && push!(breaks, length(rows) + 1)
            for (position, i) in enumerate(indices)
                row = group[i, :]
                label = row.failed ? "Failed ($(row.status))" : row.status
                cells = String[]
                if multiple
                    lattice = ismissing(row.lattice_size) ? "—" : string(row.lattice_size)
                    beta = ismissing(row.beta) ? "—" : string(row.beta)
                    append!(cells, position == 1 ? [lattice, beta] : ["", ""])
                end
                push!(cells, label)
                append!(cells, [string(row[column]) for column in source_columns])
                push!(rows, cells)
                push!(row_colors, row.failed ? orange : row.status == "Finished" ? green : "")
            end
            if multiple
                cells = ["", "", "Subtotal"]
                append!(cells, [string(sum(group[!, column])) for column in source_columns])
                push!(rows, cells)
                push!(row_colors, "")
                push!(dividers, length(rows))
            end
        end
        totals = multiple ? ["", "", "Grand total"] : ["Total"]
        append!(totals, [string(sum(stages[!, column])) for column in source_columns])
        push!(rows, totals)
        push!(row_colors, "")
        push!(dividers, length(rows))
        push!(sections, (; title="PROGRESS BY STAGE", headers, rows,
                          numeric_start=multiple ? 4 : 2, row_colors,
                          column_colors=fill("", length(headers)), breaks, dividers))
    end

    # Widths use complete, uncolored text, including full lattice dimensions.
    widths = [[maximum(textwidth, [section.headers[j]; [row[j] for row in section.rows]])
               for j in eachindex(section.headers)] for section in sections]
    rule_width = max(52, maximum([sum(w) + 2 * (length(w) - 1) for w in widths]; init=0))
    rule = repeat("─", rule_width)
    println(io, "STUDY STATUS")
    if DataFrames.nrow(parameters) == 1
        row = parameters[1, :]
        lattice = ismissing(row.lattice_size) ? "—" : string(row.lattice_size)
        beta = ismissing(row.beta) ? "—" : string(row.beta)
        println(io, "Lattice ", lattice, "  β ", beta)
    end
    println(io, rule)
    count_width = max(5, textwidth(string(overall.total)))
    for (label, value, color) in (("Finished", overall.finished, green),
                                  ("Failed", overall.failed, orange),
                                  ("Remaining", overall.remaining, ""),
                                  ("Total", overall.total, ""))
        label == "Total" && println(io, rule)
        print(io, rpad(label, 13), color, lpad(string(value), count_width))
        percentage = overall.total == 0 ? 0.0 : round(100 * value / overall.total; digits=1)
        print(io, "  ", lpad(string(percentage, "%"), 6))
        isempty(color) || print(io, reset)
        println(io)
    end

    for (section, column_widths) in zip(sections, widths)
        println(io, "\n", section.title)
        println(io, rule)
        # Row zero holds the header; all other rows are numeric summary data.
        for i in 0:length(section.rows)
            if i > 0 && i in section.breaks
                println(io)
            end
            if i in section.dividers
                println(io, rule)
            end
            cells = i == 0 ? section.headers : section.rows[i]
            for j in eachindex(cells)
                j == 1 || print(io, "  ")
                padding = repeat(" ", column_widths[j] - textwidth(cells[j]))
                right = j >= section.numeric_start
                color = i == 0 ? "" : isempty(section.row_colors[i]) ?
                    section.column_colors[j] : section.row_colors[i]
                right && print(io, padding)
                print(io, color, cells[j])
                isempty(color) || print(io, reset)
                !right && j != length(cells) && print(io, padding)
            end
            println(io)
        end
    end
    println(io, "\nNOTE: Progress is based on logs, not Slurm job state.")
    println(io)
    return data.summary
end
