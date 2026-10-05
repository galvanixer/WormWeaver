import TOML

function git_output(root, args...)
    try
        strip(read(pipeline(`git -C $root $args`; stderr=devnull), String))
    catch
        nothing
    end
end

function package_provenance(root, expected_name)
    metadata = TOML.parsefile(joinpath(root, "Project.toml"))
    get(metadata, "name", nothing) == expected_name ||
        throw(ArgumentError("Expected $expected_name project at $root"))
    git_root = git_output(root, "rev-parse", "--show-toplevel")
    # Avoid attributing an enclosing repository's revision to an installed package.
    own_repo = git_root !== nothing && realpath(git_root) == realpath(root)
    commit = own_repo ? git_output(root, "rev-parse", "HEAD") : nothing
    status = own_repo ? git_output(root, "status", "--porcelain", "--untracked-files=normal") : nothing
    Dict{String,Any}(
        "name" => expected_name, "version" => get(metadata, "version", nothing),
        "source_path" => root, "git_commit" => commit,
        "git_dirty" => status === nothing ? nothing : !isempty(status),
    )
end

function software_provenance(study, source)
    wormweaver = package_provenance(normpath(joinpath(@__DIR__, "..", "..")), "WormWeaver")
    if haskey(study, "multiwormqmc_path")
        path = study["multiwormqmc_path"]
        path isa AbstractString && !isempty(strip(path)) ||
            throw(ArgumentError("multiwormqmc_path must be a nonempty package directory path"))
        root = abspath(joinpath(dirname(source), path))
        isfile(joinpath(root, "Project.toml")) ||
            throw(ArgumentError("No Project.toml at multiwormqmc_path: $root"))
        multiworm = package_provenance(root, "MultiwormQMC")
        multiworm["discovery"] = "explicit_path"
    else
        entry = Base.find_package("MultiwormQMC")
        if entry === nothing
            multiworm = Dict{String,Any}("name" => "MultiwormQMC", "version" => nothing,
                "source_path" => nothing, "git_commit" => nothing, "git_dirty" => nothing,
                "discovery" => "unavailable")
        else
            multiworm = package_provenance(dirname(dirname(entry)), "MultiwormQMC")
            multiworm["discovery"] = "julia_environment"
        end
    end
    Dict{String,Any}("julia_version" => string(VERSION), "wormweaver" => wormweaver,
        "multiwormqmc" => multiworm,
        "seed_generation" => Dict("algorithm" => "Random.Xoshiro",
            "seed" => get(study, "seed_generation_seed", nothing)))
end
