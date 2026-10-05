"""
    init_config(path; template="config.yaml", overwrite=false) -> String

Copy a bundled simulation config or template group, preserving comments.
For a file such as `"LRBHQM/config_LRBHQM_regular.yaml"`, `path` is the output
file. For a group such as `"LRBHQM"`, `path` is the output directory and all
group files are copied with their relative filenames. Unique bare filenames
are also accepted. Relative destinations resolve from the current directory;
missing directories are created. Returns the absolute output file or directory.
Existing files require `overwrite=true`; unrelated directory contents are kept.
All destinations are checked before copying a group. Symlink destinations are
rejected. No simulations are run or submitted.
"""
function init_config(path::AbstractString; template::AbstractString="config.yaml", overwrite::Bool=false)
    source = Study.resolve_bundled_template(template; allow_group=true)
    isempty(strip(path)) && throw(ArgumentError("Config path must not be empty"))
    destination = abspath(path)
    group = isdir(source)
    copies = Pair{String,String}[]
    if group
        if ispath(destination) || islink(destination)
            isdir(destination) && !islink(destination) ||
                throw(ArgumentError("Config group destination must be a directory: $destination"))
        end
        files, _ = Study.bundled_template_catalog()
        for file in files
            input = joinpath(Study.CONFIG_TEMPLATE_DIR, file)
            relative = relpath(input, source)
            first(splitpath(relative)) == ".." && continue
            push!(copies, input => joinpath(destination, relative))
        end
        isempty(copies) && throw(ArgumentError("Template group contains no files: $template"))
    else
        push!(copies, source => destination)
    end
    for (_, output) in copies
        if ispath(output) || islink(output)
            isfile(output) && !islink(output) ||
                throw(ArgumentError("Config destination must be a regular file: $output"))
            overwrite || throw(ArgumentError("Config file already exists: $output; use overwrite=true to replace it"))
        end
        parent = dirname(output)
        while group && parent != dirname(destination)
            islink(parent) && throw(ArgumentError("Config destination directory must not be a symlink: $parent"))
            ispath(parent) && !isdir(parent) &&
                throw(ArgumentError("Config destination parent must be a directory: $parent"))
            parent = dirname(parent)
        end
    end
    contents = [read(input, String) for (input, _) in copies]
    for ((_, output), text) in zip(copies, contents)
        mkpath(dirname(output))
        open(output, "w") do io
            write(io, text)
        end
    end
    if group
        println("Copied $(length(copies)) simulation configs to: $destination (template: $template)")
    else
        println("Created simulation config: $destination (template: $template)")
    end
    destination
end
