import SHA
import YAML
import CSV
import JLD2
import DataFrames

@testset "Study generation" begin
    mktempdir() do root
        template = joinpath(root, "template.yaml")
        original = Dict("system" => Dict("beta" => 4.0),
                        "model" => Dict("mu" => [0.0, 0.0]),
                        "run" => Dict("seed" => 42, "blocks" => 10))
        YAML.write_file(template, original)
        definition = Dict{String,Any}(
            "project" => "RFEBHM",
            "template" => "./template.yaml", "output_dir" => "generated",
            "overrides" => Dict("run.blocks" => 100),
            "sweep" => Dict("system.beta" => [2.0, 4.0, 8.0],
                            "model.mu" => [[1.0, 1.0], [2.0, 2.0]]))
        source = joinpath(root, "study.yaml")
        YAML.write_file(source, definition)
        count = Ref(0)
        output = generate_study(source; validate_config=c -> (count[] += 1))
        @test output == joinpath(root, "generated")
        @test count[] == 6
        configs = [YAML.load_file(joinpath(output, "sim$i", "config.yaml")) for i in 1:6]
        @test [c["system"]["beta"] for c in configs] == [2, 4, 8, 2, 4, 8]
        @test [c["model"]["mu"] for c in configs] == [[1, 1], [1, 1], [1, 1], [2, 2], [2, 2], [2, 2]]
        @test all(c -> c["run"] == Dict("seed" => 42, "blocks" => 100), configs)
        @test read(joinpath(output, "template.yaml")) == read(template)
        @test read(joinpath(output, "study.yaml")) == read(source)
        @test YAML.load_file(template) == original
        table = CSV.read(joinpath(output, "manifest.csv"), DataFrames.DataFrame)
        native = JLD2.load(joinpath(output, "manifest.jld2"), "manifest")
        @test native isa DataFrames.DataFrame
        @test names(native) == names(table)
        @test !("simulation_id" in names(native))
        @test !("simulation_id" in names(table))
        @test native.config_sha256 == table.config_sha256
        @test native[1, "run.seed"] isa Integer
        @test native[1, "system.beta"] isa AbstractFloat
        @test all(ismissing, native[!, "lattice.extent"])
        @test all(ismissing, native.output_file)
        for i in 1:6
            @test native[i, "model.mu"] == configs[i]["model"]["mu"]
            @test typeof(native[i, "model.mu"]) == typeof(configs[i]["model"]["mu"])
        end
        @test table.simulation == ["sim$i" for i in 1:6]
        @test table.project == fill("RFEBHM", 6)
        @test YAML.load_file(joinpath(output, "study.yaml"))["project"] == "RFEBHM"
        @test all(c -> !haskey(c, "project"), configs)
        @test YAML.load(String(table[1, "model.mu"])) == [1, 1]
        @test_throws ArgumentError generate_study(source)
        @test YAML.load_file(joinpath(output, "sim1", "config.yaml")) == configs[1]

        definition["output_dir"] = "repeat"
        YAML.write_file(source, definition)
        second = generate_study(source)
        repeated = CSV.read(joinpath(second, "manifest.csv"), DataFrames.DataFrame)
        @test repeated.study_id[1] != table.study_id[1]
        @test length(unique(table.study_id)) == 1
        stable_columns = filter(n -> n != "study_id", names(table))
        @test isequal(repeated[:, stable_columns], table[:, stable_columns])
        @test table[!, "run.seed"] == fill(42, 6)
        @test table[!, "run.blocks"] == fill(100, 6)
        @test table[!, "system.beta"] == [2, 4, 8, 2, 4, 8]
        @test all(ismissing, table.output_file)
        for i in 1:6
            @test table.config_sha256[i] == bytes2hex(SHA.sha256(read(joinpath(output, table.config[i]))))
        end
        definition["output_dir"] = "invalid"
        for bad in (Dict("missing.field" => [1]), Dict("system.beta" => []),
                    Dict("system.beta" => 2), Dict("system.beta" => ["bad"]))
            definition["sweep"] = bad
            YAML.write_file(source, definition)
            @test_throws ArgumentError generate_study(source)
            @test !ispath(joinpath(root, "invalid"))
        end
        definition["sweep"] = Dict("run.blocks" => [20])
        YAML.write_file(source, definition)
        @test_throws ArgumentError generate_study(source)
        definition["sweep"] = Dict("run" => [Dict("blocks" => 20)])
        YAML.write_file(source, definition)
        @test_throws ArgumentError generate_study(source)
        definition["sweep"] = Dict()
        YAML.write_file(source, definition)
        @test_throws ErrorException generate_study(source; validate_config=c -> error("Invalid physics"))
        @test !ispath(joinpath(root, "invalid"))
        single = generate_study(source)
        @test isfile(joinpath(single, "sim1", "config.yaml"))
        @test !ispath(joinpath(single, "sim2"))

        # Exercise the actual bundled simulation template, including nested arrays.
        YAML.write_file(source, Dict("project" => "LRBHQM", "output_dir" => "bundled"))
        bundled = generate_study(source)
        @test YAML.load_file(joinpath(bundled, "sim1", "config.yaml")) ==
              YAML.load_file(joinpath(@__DIR__, "..", "templates", "multiwormqmc", "config.yaml"))
        bundled_manifest = CSV.read(joinpath(bundled, "manifest.csv"), DataFrames.DataFrame)
        @test bundled_manifest.project == ["LRBHQM"]
        @test bundled_manifest.output_file == ["sim1/results.h5"]
        @test bundled_manifest[1, "model.family"] == "bose_hubbard"
        @test YAML.load(String(bundled_manifest[1, "model.onsite.U"])) == [[4.0, 1.5], [1.5, 5.0]]
        @test YAML.load(String(bundled_manifest[1, "worldlines.initial_density"])) ==
              [Dict("num" => 2, "den" => 16), Dict("num" => 2, "den" => 16)]
        @test !("sampler.moves" in names(bundled_manifest))
        native_bundled = JLD2.load(joinpath(bundled, "manifest.jld2"), "manifest")
        bundled_config = YAML.load_file(joinpath(bundled, "sim1", "config.yaml"))
        for (field, value) in (
            "model.onsite.U" => bundled_config["model"]["onsite"]["U"],
            "worldlines.initial_density" => bundled_config["worldlines"]["initial_density"],
            "lattice.pbc" => bundled_config["lattice"]["pbc"],
        )
            @test isequal(native_bundled[1, field], value)
            @test typeof(native_bundled[1, field]) == typeof(value)
        end
        @test native_bundled[1, "model.density_density.alpha"] === nothing
        @test native_bundled[1, "lattice.pbc"][1] isa Bool
        for file in ("../shared/results.h5", joinpath(root, "absolute.h5"))
            YAML.write_file(source, Dict("project" => "paths", "output_dir" => "path$(basename(file))",
                "overrides" => Dict("run.output_file" => file)))
            generated = generate_study(source)
            expected = isabspath(file) ? file : normpath(joinpath("sim1", file))
            @test CSV.read(joinpath(generated, "manifest.csv"), DataFrames.DataFrame).output_file == [expected]
        end

        for invalid in (nothing, "", "   ", 123, ["RFEBHM"])
            bad = Dict{String,Any}("output_dir" => "bad_project", "project" => invalid)
            YAML.write_file(source, bad)
            @test_throws ArgumentError generate_study(source)
            @test !ispath(joinpath(root, "bad_project"))
        end
        YAML.write_file(source, Dict("output_dir" => "missing_project"))
        @test_throws ArgumentError generate_study(source)
        @test !ispath(joinpath(root, "missing_project"))
    end
end

@testset "Study replacement" begin
    mktempdir() do root
        source = joinpath(root, "study.yaml")
        definition = Dict{String,Any}("project" => "rewrite", "output_dir" => "output",
            "sweep" => Dict("system.beta" => [2.0, 4.0]))
        YAML.write_file(source, definition)
        output = generate_study(source; overwrite=true)
        old_manifest = read(joinpath(output, "manifest.csv"))
        write(joinpath(output, "sim1", "results.h5"), "old result")
        @test_throws ArgumentError generate_study(source)
        @test_throws ErrorException generate_study(source; overwrite=true,
            validate_config=c -> error("validation failure"))
        @test read(joinpath(output, "manifest.csv")) == old_manifest
        @test read(joinpath(output, "sim1", "results.h5"), String) == "old result"
        definition["sweep"] = Dict("system.beta" => [8.0])
        YAML.write_file(source, definition)
        @test generate_study(source; overwrite=true) == output
        @test !ispath(joinpath(output, "sim2"))
        @test !ispath(joinpath(output, "sim1", "results.h5"))
        @test YAML.load_file(joinpath(output, "sim1", "config.yaml"))["system"]["beta"] == 8.0
        @test read(joinpath(output, "manifest.csv")) != old_manifest
        native = JLD2.load(joinpath(output, "manifest.jld2"), "manifest")
        @test size(native, 1) == 1
        @test native.config_sha256[1] == bytes2hex(SHA.sha256(read(joinpath(output, "sim1", "config.yaml"))))
        @test isempty(filter(n -> startswith(n, ".wormweaver-"), readdir(root)))
        for target in ("file", "link")
            path = joinpath(root, target)
            target == "file" ? write(path, "keep") : symlink(output, path)
            definition["output_dir"] = target
            YAML.write_file(source, definition)
            @test_throws ArgumentError generate_study(source; overwrite=true)
        end
        @test read(joinpath(root, "file"), String) == "keep"
        @test islink(joinpath(root, "link"))
        @test isdir(output)
    end
end

@testset "Random seed replicas" begin
    mktempdir() do root
        source = joinpath(root, "study.yaml")
        definition = Dict{String,Any}("project" => "seeds", "output_dir" => "default", "replicas" => 10)
        YAML.write_file(source, definition)
        output = generate_study(source)
        native = JLD2.load(joinpath(output, "manifest.jld2"), "manifest")
        @test native.replica == collect(1:10)
        @test length(unique(native[!, "run.seed"])) == 10
        @test all(>(0), native[!, "run.seed"])
        template = YAML.load_file(joinpath(output, "template.yaml"))
        @test all(==(template["worldlines"]["init_seed"]), native[!, "worldlines.init_seed"])
        @test all(==(template["model"]["disorder"]["seed"]), native[!, "model.disorder.seed"])
        paths = ["run.seed", "worldlines.init_seed", "model.disorder.seed"]
        definition["output_dir"] = "all"
        definition["replicas"] = 3
        definition["randomize_seeds"] = paths
        definition["sweep"] = Dict("system.beta" => [2.0, 4.0])
        YAML.write_file(source, definition)
        output = generate_study(source)
        native = JLD2.load(joinpath(output, "manifest.jld2"), "manifest")
        csv = CSV.read(joinpath(output, "manifest.csv"), DataFrames.DataFrame)
        @test native.replica == [1, 2, 3, 1, 2, 3]
        @test native[!, "system.beta"] == [2, 2, 2, 4, 4, 4]
        @test length(unique(vcat([native[!, field] for field in paths]...))) == 18
        @test csv.replica == native.replica
        for i in 1:6
            config = YAML.load_file(joinpath(output, "sim$i", "config.yaml"))
            for field in paths
                parts = split(field, '.')
                value = config
                for part in parts
                    value = value[part]
                end
                @test value == native[i, field] == csv[i, field]
            end
        end
        # Selected subset leaves an explicitly fixed sampling seed untouched.
        definition["output_dir"] = "subset"
        definition["randomize_seeds"] = ["worldlines.init_seed", "model.disorder.seed"]
        definition["overrides"] = Dict("run.seed" => 7)
        YAML.write_file(source, definition)
        output = generate_study(source)
        @test all(==(7), JLD2.load(joinpath(output, "manifest.jld2"), "manifest")[!, "run.seed"])
        for changes in (
            Dict("replicas" => 0), Dict("replicas" => -1), Dict("replicas" => true), Dict("replicas" => 2.5),
            Dict("randomize_seeds" => []), Dict("randomize_seeds" => "run.seed"),
            Dict("randomize_seeds" => ["unknown.seed"]), Dict("randomize_seeds" => ["run.seed", "run.seed"]),
            Dict("randomize_seeds" => ["run.seed"], "overrides" => Dict("run.seed" => 4)),
            Dict("randomize_seeds" => ["run.seed"], "sweep" => Dict("run.seed" => [4, 5])),
            Dict("randomize_seeds" => ["run.seed"], "overrides" => Dict("run" => Dict("seed" => 4))),
        )
            bad = merge(Dict{String,Any}("project" => "bad", "output_dir" => "bad", "replicas" => 2), changes)
            YAML.write_file(source, bad)
            @test_throws ArgumentError generate_study(source)
            @test !ispath(joinpath(root, "bad"))
        end
        YAML.write_file(source, Dict("project" => "bad", "output_dir" => "bad", "randomize_seeds" => ["run.seed"]))
        @test_throws ArgumentError generate_study(source)
    end
end

@testset "Compact configuration YAML" begin
    encode = WormWeaver.Study.config_yaml
    template = YAML.load_file(joinpath(@__DIR__, "..", "templates", "multiwormqmc", "config.yaml"))
    text = encode(template)
    @test isequal(YAML.load(text), template)
    @test occursin("target: [0, 1]", text)
    @test occursin(r"- \{[^\n]*kind: \"Open\"[^\n]*\}", text)
    @test occursin("- [4.0, 1.5]", text)
    @test count(==('\n'), text) < count(==('\n'), YAML.write(template))
    tricky = Dict{String,Any}(
        "strings" => ["true", "null", "01", "2026-01-01", "a: b", "#comment", "", "yes"],
        "special" => ["quote\"", "slash\\", "line\nnext\n", "\t\r\0\u001b", "α🐛", "\u0085\u2028", "\$var"],
        "null" => nothing, "yes" => "no", "a,b: []" => "x",
        "empty" => Any[Dict(), Any[]], "long" => collect(1:20),
        "numbers" => [1.0, -0.0, Inf, -Inf, NaN],
        "records" => [Dict("a" => Dict("nested" => true)), Dict("list" => [1, 2], "flag" => false)],
    )
    @test isequal(YAML.load(encode(tricky)), tricky)
    @test_throws ArgumentError encode(Dict("bad" => missing))
end

@testset "Reproducible seeds and provenance" begin
    mktempdir() do root
        package = joinpath(root, "MultiwormQMC")
        mkdir(package)
        write(joinpath(package, "Project.toml"), "name = \"MultiwormQMC\"\nversion = \"0.2.3\"\n")
        source = joinpath(root, "study.yaml")
        definition = Dict{String,Any}("project" => "reproducible", "output_dir" => "one",
            "replicas" => 3, "seed_generation_seed" => 42, "multiwormqmc_path" => "MultiwormQMC",
            "randomize_seeds" => ["run.seed", "worldlines.init_seed", "model.disorder.seed"],
            "sweep" => Dict("system.beta" => [2.0, 4.0]))
        YAML.write_file(source, definition)
        first = generate_study(source)
        a = JLD2.load(joinpath(first, "manifest.jld2"), "manifest")
        definition["output_dir"] = "two"
        YAML.write_file(source, definition)
        second = generate_study(source)
        b = JLD2.load(joinpath(second, "manifest.jld2"), "manifest")
        @test a.config_sha256 == b.config_sha256
        @test a.study_id != b.study_id
        for field in definition["randomize_seeds"]
            @test a[!, field] == b[!, field]
        end
        generate_study(source; overwrite=true)
        @test JLD2.load(joinpath(second, "manifest.jld2"), "manifest").config_sha256 == a.config_sha256
        provenance = YAML.load_file(joinpath(first, "provenance.yaml"))
        @test provenance["multiwormqmc"]["version"] == "0.2.3"
        @test provenance["multiwormqmc"]["git_commit"] === nothing
        @test provenance["multiwormqmc"]["discovery"] == "explicit_path"
        @test provenance["wormweaver"]["version"] == "0.1.0"
        @test provenance["julia_version"] == string(VERSION)
        @test provenance["seed_generation"]["seed"] == 42
        @test provenance["study_id"] == a.study_id[1]
        definition["seed_generation_seed"] = 43
        generate_path = joinpath(root, "different")
        definition["output_dir"] = generate_path
        YAML.write_file(source, definition)
        generate_study(source)
        @test JLD2.load(joinpath(generate_path, "manifest.jld2"), "manifest")[!, "run.seed"] != a[!, "run.seed"]
        for bad in (-1, true, 1.5, "42", nothing)
            definition["seed_generation_seed"] = bad
            definition["output_dir"] = "invalid_seed"
            YAML.write_file(source, definition)
            @test_throws ArgumentError generate_study(source)
            @test !ispath(joinpath(root, "invalid_seed"))
        end
        definition["seed_generation_seed"] = 42
        delete!(definition, "replicas")
        delete!(definition, "randomize_seeds")
        YAML.write_file(source, definition)
        @test_throws ArgumentError generate_study(source)
    end
end

@testset "Direct jobfiles" begin
    mktempdir() do root
        source = joinpath(root, "study.yaml")
        script = joinpath(root, "fake runner.sh")
        write(script, "printf '%s\\n' \"\$PWD\" \"\$1\" \"\$2\"; echo stderr >&2; test -f \"\$2\"\n")
        literal = "a 'quoted' \$(touch SHOULD_NOT_EXIST) argument"
        definition = Dict{String,Any}("project" => "jobs", "output_dir" => "study 'quoted' space",
            "sweep" => Dict("system.beta" => [2.0, 4.0]),
            "jobs" => Dict("command" => ["bash", script, literal]))
        YAML.write_file(source, definition)
        output = generate_study(source)
        lines = readlines(joinpath(output, "jobfile"))
        @test length(lines) == 2
        for (i, line) in enumerate(lines)
            run(`bash -c $line`)
            log = readlines(joinpath(output, "sim$i", "simulation.log"))
            @test log == [joinpath(output, "sim$i"), literal, "config.yaml", "stderr"]
            @test !ispath(joinpath(output, "sim$i", "SHOULD_NOT_EXIST"))
        end
        @test !success(pipeline(`bash -c $(replace(lines[1], "sim1" => "absent"))`; stderr=devnull))
        definition["output_dir"] = "remote"
        definition["jobs"] = Dict("run_dir" => "/cluster/study", "mode" => "direct")
        YAML.write_file(source, definition)
        remote = generate_study(source)
        @test readlines(joinpath(remote, "jobfile"))[1] ==
            "(cd '/cluster/study/sim1' && 'multiwormqmc' config.yaml > simulation.log 2>&1)"
        definition["output_dir"] = "failure"
        definition["jobs"] = Dict("command" => ["bash", "-c", "exit 7"])
        YAML.write_file(source, definition)
        failed = generate_study(source)
        proc = run(ignorestatus(`bash -c $(readlines(joinpath(failed, "jobfile"))[1])`))
        @test proc.exitcode == 7
        for jobs in (Dict("mode" => "unsupported"), Dict("command" => "multiwormqmc"),
                     Dict("command" => []), Dict("command" => [""]),
                     Dict("command" => ["cmd\nother"]), Dict("run_dir" => "relative"),
                     Dict("run_dir" => "/bad\npath"), Dict("unknown" => true))
            definition["output_dir"] = "invalid_jobs"
            definition["jobs"] = jobs
            YAML.write_file(source, definition)
            @test_throws ArgumentError generate_study(source)
            @test !ispath(joinpath(root, "invalid_jobs"))
        end
        delete!(definition, "jobs")
        definition["output_dir"] = "no_jobs"
        YAML.write_file(source, definition)
        default_output = generate_study(source)
        @test length(readlines(joinpath(default_output, "jobfile"))) == 2
        @test occursin("'multiwormqmc' config.yaml", read(joinpath(default_output, "jobfile"), String))
        grant = joinpath(@__DIR__, "..", "templates", "slurm", "grant_g2026a136c.slurm")
        @test read(joinpath(default_output, "grant_g2026a136c.slurm")) == read(grant)
        @test read(joinpath(remote, "grant_g2026a136c.slurm")) == read(grant)
        templates = joinpath(@__DIR__, "..", "templates", "slurm")
        for name in ("grant_g2026a136c.slurm", "grantgpu_g2026a136g.slurm")
            @test read(joinpath(default_output, name)) == read(joinpath(templates, name))
        end
        @test occursin("#SBATCH -p grant -A g2026a136c", read(grant, String))
        @test occursin("#SBATCH -p grantgpu -A g2026a136g", read(joinpath(default_output, "grantgpu_g2026a136g.slurm"), String))
        @test !isfile(joinpath(default_output, "my_script_grant.slurm"))
        write(joinpath(default_output, "grant_g2026a136c.slurm"), "old script")
        generate_study(source; overwrite=true)
        @test read(joinpath(default_output, "grant_g2026a136c.slurm")) == read(grant)
        @test isfile(joinpath(default_output, "jobfile"))
    end
end

@testset "Scratch jobfile format" begin
    source = "/home/user/Run/study/sim600"
    command = "~/Code/multiwormqmc config.yaml &> simulation.log"
    expected = "cp -r /home/user/Run/study/sim600 /scratch/job.\$SLURM_JOB_ID/sim600; cd /scratch/job.\$SLURM_JOB_ID/sim600 ;~/Code/multiwormqmc config.yaml &> simulation.log ; cp -r /scratch/job.\$SLURM_JOB_ID/sim600/* /home/user/Run/study/sim600 ;\n"
    @test WormWeaver.Study.scratch_jobline(command, source, "/scratch", "sim600") == expected
    study = Dict("jobs" => Dict("mode" => "scratch", "command" => ["~/Code/multiwormqmc"]))
    @test WormWeaver.Study.jobfile_text(study, "/home/user/Run/study", ["sim600"]) == expected
    @test WormWeaver.Study.scratch_jobline(command, source, "/scratch/", "sim600") == expected
end

@testset "Generation summary" begin
    mktempdir() do root
        source = joinpath(root, "study.yaml")
        YAML.write_file(source, Dict("project" => "summary", "output_dir" => "generated",
            "replicas" => 2, "sweep" => Dict("system.beta" => [2.0, 4.0, 8.0]),
            "jobs" => Dict("mode" => "scratch", "run_dir" => "/cluster/study")))
        capture = joinpath(root, "stdout.txt")
        result = open(capture, "w") do io
            redirect_stdout(io) do
                generate_study(source)
            end
        end
        text = read(capture, String)
        @test result == joinpath(root, "generated")
        @test occursin("Created 6 simulations: 3 parameter combinations × 2 replicas", text)
        @test occursin("Directory: $result", text)
        @test occursin("Job mode: scratch", text)
        @test occursin("cd '/cluster/study'", text)
        @test occursin("sbatch grant_g2026a136c.slurm", text)
        @test occursin("sbatch grantgpu_g2026a136g.slurm", text)
        @test occursin("Choose one account script. No jobs have been submitted.", text)
        open(capture, "w") do io
            redirect_stdout(io) do
                @test_throws ArgumentError generate_study(source)
            end
        end
        @test isempty(read(capture, String))
    end
end

@testset "Bundled simulation template selection" begin
    mktempdir() do root
        source = joinpath(root, "study.yaml")
        local_config = Dict("system" => Dict("beta" => 123.0))
        YAML.write_file(joinpath(root, "config_LRBHQM.yaml"), local_config)
        for name in ("config.yaml", "config_LRBHQM.yaml")
            YAML.write_file(source, Dict("project" => "templates", "output_dir" => name * "_output", "template" => name))
            output = generate_study(source)
            bundled = joinpath(@__DIR__, "..", "templates", "multiwormqmc", name)
            @test read(joinpath(output, "template.yaml")) == read(bundled)
            @test YAML.load_file(joinpath(output, "sim1", "config.yaml")) == YAML.load_file(bundled)
        end
        for (i, template) in enumerate(("./config_LRBHQM.yaml", joinpath(root, "config_LRBHQM.yaml")))
            YAML.write_file(source, Dict("project" => "local", "output_dir" => "local$i", "template" => template))
            output = generate_study(source)
            @test YAML.load_file(joinpath(output, "sim1", "config.yaml")) == local_config
        end
        YAML.write_file(joinpath(root, "custom.yaml"), local_config)
        YAML.write_file(source, Dict("project" => "bad", "output_dir" => "bad", "template" => "custom.yaml"))
        err = try
            generate_study(source)
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("config_LRBHQM.yaml", sprint(showerror, err))
        @test occursin("./filename.yaml", sprint(showerror, err))
        @test !ispath(joinpath(root, "bad"))
    end
end
