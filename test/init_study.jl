@testset "Initialize study definition" begin
    mktempdir() do root
        for template in ("study.yaml", "unistra-hpc.yaml")
            path = joinpath(root, template, "study.yaml")
            @test init_study(path; template) == path
            original = read(joinpath(@__DIR__, "..", "templates", "studies", template))
            @test read(path) == original
            write(path, "user edits")
            @test_throws ArgumentError init_study(path; template)
            @test read(path, String) == "user edits"
            @test init_study(path; template, overwrite=true) == path
            @test read(path) == original
            @test readdir(dirname(path)) == ["study.yaml"]
        end
        cd(root) do
            @test realpath(init_study("default.yaml")) == realpath(joinpath(root, "default.yaml"))
            @test read("default.yaml") == read(joinpath(@__DIR__, "..", "templates", "studies", "study.yaml"))
        end
        @test_throws ArgumentError init_study(joinpath(root, "bad", "study.yaml"); template="unknown")
        @test !ispath(joinpath(root, "bad"))
        for template in ("../multiwormqmc/config.yaml", "/tmp/study.yaml", "general", "unistra")
            @test_throws ArgumentError init_study(joinpath(root, "bad.yaml"); template)
        end
        err = try
            init_study(joinpath(root, "bad.yaml"); template="unknown.yaml")
        catch e
            e
        end
        @test occursin("study.yaml", sprint(showerror, err))
        @test occursin("unistra-hpc.yaml", sprint(showerror, err))
        @test_throws ArgumentError init_study("")
        @test_throws ArgumentError init_study(root; overwrite=true)
        link = joinpath(root, "link.yaml")
        symlink(joinpath(root, "default.yaml"), link)
        @test_throws ArgumentError init_study(link; overwrite=true)
        @test islink(link)
    end
end
