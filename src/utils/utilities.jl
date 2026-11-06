import DataFrames

"""
    describe_column_sizes(table; io=stdout, unit=:auto, top=nothing) -> nothing

Print the estimated in-memory size of DataFrame columns in descending order.

Sizes are computed independently using `Base.summarysize` and may not equal
the serialized JLD2 file size or sum to the DataFrame's total memory usage.

# Keyword Arguments
- `io::IO=stdout`: Output stream.
- `unit::Symbol=:auto`: Display unit (`:B`, `:KiB`, `:MiB`, `:GiB`, `:TiB`).
  With `:auto`, the most appropriate unit is selected for each column.
- `top::Union{Int,Nothing}=nothing`: Number of largest columns to display.
  `nothing` displays all columns.

Units use binary multiples (1 KiB = 1024 bytes).
"""
function describe_column_sizes(table::DataFrames.AbstractDataFrame;
                               io::IO=stdout,
                               unit::Symbol=:auto,
                               top::Union{Int,Nothing}=nothing)
    units = (:B, :KiB, :MiB, :GiB, :TiB)
    factors = (1, 1024, 1024^2, 1024^3, 1024^4)

    unit == :auto || unit in units ||
        throw(ArgumentError("Unsupported unit: $unit"))

    isnothing(top) || top >= 0 ||
        throw(ArgumentError("top must be non-negative"))

    sizes = [(name=String(name), bytes=Base.summarysize(table[!, name]))
             for name in DataFrames.names(table)]

    sort!(sizes; by=entry -> entry.bytes, rev=true)
    width = max(32, maximum(length(entry.name) for entry in sizes; init=0))

    selected = isnothing(top) ? sizes : @view sizes[1:min(top, length(sizes))]

    for entry in selected
        idx = unit == :auto ?
              something(findlast(f -> entry.bytes >= f, factors), 1) :
              findfirst(==(unit), units)

        value = entry.bytes / factors[idx]
        formatted = idx == 1 ? string(entry.bytes) :
                    string(round(value; digits=1))

        println(io, rpad(entry.name, width), lpad(formatted, 10), " ", units[idx])
    end

    return nothing
end
