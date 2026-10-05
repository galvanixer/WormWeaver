const CONFIG_TEMPLATE_DIR = normpath(joinpath(@__DIR__, "..", "..", "templates", "multiwormqmc"))

function bundled_template_catalog()
    files = String[]
    groups = String[]
    for (root, dirs, names) in walkdir(CONFIG_TEMPLATE_DIR)
        filter!(name -> !islink(joinpath(root, name)), dirs)
        root == CONFIG_TEMPLATE_DIR || push!(groups, relpath(root, CONFIG_TEMPLATE_DIR))
        for name in names
            path = joinpath(root, name)
            !islink(path) && isfile(path) && push!(files, relpath(path, CONFIG_TEMPLATE_DIR))
        end
    end
    sort!(files), sort!(groups)
end

function resolve_bundled_template(template::AbstractString; allow_group::Bool=false)
    files, groups = bundled_template_catalog()
    available = allow_group ? sort(vcat(files, groups)) : files
    valid = !isabspath(template) && !isempty(strip(template)) &&
        !any(part -> part in (".", ".."), splitpath(template))
    if valid
        template in files && return joinpath(CONFIG_TEMPLATE_DIR, template)
        allow_group && template in groups && return joinpath(CONFIG_TEMPLATE_DIR, template)
        if basename(template) == template
            matches = filter(path -> basename(path) == template, files)
            length(matches) == 1 && return joinpath(CONFIG_TEMPLATE_DIR, only(matches))
            length(matches) > 1 && throw(ArgumentError(
                "Ambiguous bundled template: $template; use one of: $(join(matches, ", "))"))
        end
    end
    throw(ArgumentError(
        "Unknown bundled template: $template; available templates: $(join(available, ", ")). Use ./filename.yaml for a local file in a study definition."))
end

function is_bundled_template_reference(path::AbstractString)
    isabspath(path) && return false
    !occursin('/', path) && !occursin('\\', path) && return true
    _, groups = bundled_template_catalog()
    first(splitpath(path)) in groups
end
