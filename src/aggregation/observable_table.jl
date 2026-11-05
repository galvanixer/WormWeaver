"""
    observable_table(summaries; fields=nothing) -> DataFrame

Expand an observable column such as `table.Windings` into a DataFrame with one
row per input entry and one column per dictionary key. Inputs must be dictionaries
with string keys or `missing`. Missing entries and absent keys produce missing
cells. Measurements, arrays, and metadata are retained as stored; arrays remain
individual cells rather than being expanded into component columns.

By default, use the sorted union of keys across all dictionaries. Supply a
nonempty vector or tuple of unique string `fields` to select columns and their
order. Nonempty inputs without any dictionary keys (including all-missing
columns) require explicit fields so their rows can be preserved. Empty inputs
return an empty DataFrame, with the requested columns if fields are supplied.

The input dictionaries are not modified. Array-valued cells reference the
original arrays, as in a normal DataFrame column.

# Examples
```julia
windings = observable_table(table.Windings)
totals = observable_table(table.Windings; fields=["totalWinding", "squaredTotalWinding"])
```
"""
function observable_table(summaries::AbstractVector; fields=nothing)
    field_names = Set{String}()
    for (i, summary) in enumerate(summaries)
        ismissing(summary) && continue
        summary isa AbstractDict || throw(ArgumentError("entry $i must be a dictionary or missing"))
        all(key -> key isa AbstractString, keys(summary)) ||
            throw(ArgumentError("entry $i must have string keys"))
        union!(field_names, keys(summary))
    end

    if fields === nothing
        selected = sort!(collect(field_names))
        isempty(selected) && !isempty(summaries) && throw(ArgumentError(
            "cannot infer columns from summaries without keys; supply fields"))
    else
        fields isa Union{AbstractVector,Tuple} ||
            throw(ArgumentError("fields must be a vector or tuple of string names"))
        all(field -> field isa AbstractString, fields) ||
            throw(ArgumentError("fields must contain string names"))
        selected = String.(collect(fields))
        isempty(selected) && throw(ArgumentError("fields must not be empty"))
        allunique(selected) || throw(ArgumentError("fields must not contain duplicates"))
    end

    result = DataFrames.DataFrame()
    for field in selected
        result[!, field] = [ismissing(summary) ? missing : get(summary, field, missing)
                           for summary in summaries]
    end
    return result
end
