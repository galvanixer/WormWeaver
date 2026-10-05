@testset "Simulation config initialization" begin
    template_dir = joinpath(@__DIR__, "..", "templates", "multiwormqmc")
    mktempdir() do root
        path = joinpath(root, "nested", "config.yaml")
        template = "LRBHQM/config_LRBHQM_regular.yaml"
        @test init_config(path; template) == path
        @test read(path) == read(joinpath(template_dir, template))
        bare_path = joinpath(root, "bare.yaml")
        @test init_config(bare_path; template=basename(template)) == bare_path
        @test read(bare_path) == read(path)

        write(path, "edited config\n")
        @test_throws ArgumentError init_config(path; template)
        @test read(path, String) == "edited config\n"
        @test init_config(path; template, overwrite=true) == path
        @test read(path) == read(joinpath(template_dir, template))

        cd(root) do
            destination = init_config("default/config.yaml")
            @test destination == abspath("default/config.yaml")
            @test read(destination) == read(joinpath(template_dir, "config.yaml"))
        end

        invalid_path = joinpath(root, "invalid", "config.yaml")
        for invalid_template in ("unknown.yaml", "../multiwormqmc/config.yaml", "")
            @test_throws ArgumentError init_config(invalid_path; template=invalid_template)
            @test !ispath(dirname(invalid_path))
        end
        @test_throws ArgumentError init_config("")
        @test_throws ArgumentError init_config("   ")
        @test_throws ArgumentError init_config(root; overwrite=true)

        link = joinpath(root, "link.yaml")
        symlink(path, link)
        @test_throws ArgumentError init_config(link; overwrite=true)
        @test islink(link)
        @test read(path) == read(joinpath(template_dir, template))
        dangling_link = joinpath(root, "dangling.yaml")
        symlink(joinpath(root, "missing.yaml"), dangling_link)
        @test_throws ArgumentError init_config(dangling_link; overwrite=true)
        @test !ispath(joinpath(root, "missing.yaml"))
    end
end

@testset "Simulation config group initialization" begin
    group_dir = joinpath(@__DIR__, "..", "templates", "multiwormqmc", "LRBHQM")
    names = sort(readdir(group_dir))
    mktempdir() do root
        output = joinpath(root, "nested", "configs")
        @test init_config(output; template="LRBHQM") == output
        @test readdir(output) == names
        @test all(read(joinpath(output, name)) == read(joinpath(group_dir, name)) for name in names)

        sentinel = joinpath(output, "notes.txt")
        write(sentinel, "keep me")
        edited = joinpath(output, first(names))
        write(edited, "user edits")
        @test_throws ArgumentError init_config(output; template="LRBHQM")
        @test read(edited, String) == "user edits"
        @test init_config(output; template="LRBHQM", overwrite=true) == output
        @test read(sentinel, String) == "keep me"
        @test all(read(joinpath(output, name)) == read(joinpath(group_dir, name)) for name in names)

        # A conflict in the last filename must leave earlier files uncreated.
        conflict_dir = joinpath(root, "conflict")
        mkpath(conflict_dir)
        conflict = joinpath(conflict_dir, last(names))
        write(conflict, "existing config")
        @test_throws ArgumentError init_config(conflict_dir; template="LRBHQM")
        @test readdir(conflict_dir) == [last(names)]
        @test read(conflict, String) == "existing config"

        rm(conflict)
        symlink(sentinel, conflict)
        @test_throws ArgumentError init_config(conflict_dir; template="LRBHQM", overwrite=true)
        @test readdir(conflict_dir) == [last(names)]
        @test read(sentinel, String) == "keep me"
        @test_throws ArgumentError init_config(sentinel; template="LRBHQM", overwrite=true)
        link = joinpath(root, "linked_directory")
        symlink(output, link)
        @test_throws ArgumentError init_config(link; template="LRBHQM", overwrite=true)

        cd(root) do
            @test init_config("relative"; template="LRBHQM") == abspath("relative")
            @test readdir("relative") == names
        end
        @test_throws ArgumentError init_config(joinpath(root, "invalid"); template="LRBHQM/../config.yaml")
        @test !ispath(joinpath(root, "invalid"))
    end
end
