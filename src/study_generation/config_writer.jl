# Compact formatting for YAML-compatible simulation configuration values.
# Quote strings explicitly so YAML never interprets them as numbers or booleans.
function yaml_string(value::AbstractString)
    io = IOBuffer()
    print(io, '"')
    for c in value
        if c == '"' || c == '\\'
            print(io, '\\', c)
        elseif Int(c) < 0x20 || 0x7f <= Int(c) <= 0x9f || c in ('\u2028', '\u2029')
            print(io, "\\u", string(Int(c); base=16, pad=4))
        else
            print(io, c)
        end
    end
    print(io, '"')
    String(take!(io))
end

function yaml_scalar(value)
    value === nothing && return "null"
    value isa AbstractString && return yaml_string(value)
    value isa Bool && return value ? "true" : "false"
    value isa Integer && return string(value)
    if value isa AbstractFloat
        isnan(value) && return ".nan"
        isinf(value) && return value > 0 ? ".inf" : "-.inf"
        return string(value)
    end
    throw(ArgumentError("Unsupported configuration value type: $(typeof(value))"))
end

yaml_key(key) = key isa AbstractString && occursin(r"^[A-Za-z_][A-Za-z_0-9]*$", key) &&
    !(lowercase(key) in ("null", "true", "false", "yes", "no", "on", "off", "y", "n")) ?
    String(key) : yaml_scalar(key)

function inline_yaml(value; record=false)
    if value isa AbstractVector
        isempty(value) && return "[]"
        length(value) <= 8 || return nothing
        any(v -> v isa AbstractVector || v isa AbstractDict, value) && return nothing
        text = "[" * join(yaml_scalar.(value), ", ") * "]"
        return length(text) <= 100 ? text : nothing
    elseif value isa AbstractDict
        isempty(value) && return "{}"
        record && length(value) <= 4 || return nothing
        parts = String[]
        for (key, item) in value
            item isa AbstractDict && return nothing
            text = inline_yaml(item)
            text === nothing && return nothing
            push!(parts, yaml_key(key) * ": " * text)
        end
        text = "{" * join(parts, ", ") * "}"
        return length(text) <= 120 ? text : nothing
    end
    yaml_scalar(value)
end

function write_config_node(io, value, indent=0)
    padding = repeat(" ", indent)
    if value isa AbstractDict && !isempty(value)
        for (key, item) in value
            text = inline_yaml(item)
            print(io, padding, yaml_key(key), ":")
            if text === nothing
                println(io)
                write_config_node(io, item, indent + 2)
            else
                println(io, " ", text)
            end
        end
    elseif value isa AbstractVector && !isempty(value)
        for item in value
            text = inline_yaml(item; record=true)
            if text === nothing
                println(io, padding, "-")
                write_config_node(io, item, indent + 2)
            else
                println(io, padding, "- ", text)
            end
        end
    else
        println(io, padding, inline_yaml(value))
    end
end

function config_yaml(config)
    io = IOBuffer()
    write_config_node(io, config)
    String(take!(io))
end
