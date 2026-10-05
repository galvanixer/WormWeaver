"""
    init_study(path; template="study.yaml", overwrite=false) -> String

Copy a commented study template to `path` for editing. `template` is a filename
in `templates/studies/`, such as `"study.yaml"` or `"unistra-hpc.yaml"`.
Relative destinations resolve from the current directory;
missing parent directories are created. Returns the absolute destination path.
Existing files require `overwrite=true`; directories and symlinks are rejected.
No simulations are generated or submitted.
"""
function init_study(path::AbstractString; template::AbstractString="study.yaml", overwrite::Bool=false)
    template_dir = joinpath(@__DIR__, "..", "..", "templates", "studies")
    available = sort(filter(name -> isfile(joinpath(template_dir, name)), readdir(template_dir)))
    template in available || throw(ArgumentError(
        "Unknown template: $template; available templates: $(join(available, ", "))"))
    isempty(strip(path)) && throw(ArgumentError("Study path must not be empty"))
    destination = abspath(path)
    if ispath(destination) || islink(destination)
        isfile(destination) && !islink(destination) ||
            throw(ArgumentError("Study destination must be a regular file: $destination"))
        overwrite || throw(ArgumentError("Study file already exists: $destination; use overwrite=true to replace it"))
    end
    template_path = joinpath(template_dir, template)
    contents = read(template_path, String)
    mkpath(dirname(destination))
    open(destination, "w") do io
        write(io, contents)
    end
    println("Created study definition: $destination (template: $template)")
    println("Edit its paths and parameters, then run generate_study(", repr(destination), ").")
    destination
end
