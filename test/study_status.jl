import CSV
import Dates
import DataFrames
import JLD2

@testset "Status summary" begin
    df = DataFrames.DataFrame(
        lattice_size=Any[[64], [64], [64], [64], [32], [64], [16, 16], missing, missing],
        beta=[64.0, 64.0, 64.0, 32.0, 64.0, 64.0, 16.0, missing, missing],
        status=["Finished", "Finished", "Tuning", "Finished", "Finished", "Tuning",
                "Finished", "Pre Tuning", "Pre Tuning"],
        failed=[false, false, true, false, false, false, false, false, false])
    original = deepcopy(df)
    io = IOBuffer()
    summary = status_summary(df; io)
    printed = String(take!(io))
    @test startswith(printed, "STUDY STATUS\n")
    @test occursin(r"Total\s+9\s+100\.0%", printed)
    @test occursin(r"Finished\s+5\s+55\.6%", printed)
    @test occursin(r"Failed\s+1\s+11\.1%", printed)
    @test occursin(r"Remaining\s+3\s+33\.3%", printed)
    @test names(summary) == ["lattice_size", "beta", "status", "failed", "n_simulations"]
    @test summary.n_simulations == [2, 1, 1, 1, 1, 1, 2]
    @test isequal(summary[:, 1:4], df[[1, 3, 4, 5, 6, 7, 8], :])
    @test sum(summary.n_simulations) == DataFrames.nrow(df)
    @test isequal(df, original)

    empty = status_summary(df[1:0, :]; io)
    @test occursin(r"Remaining\s+0\s+0\.0%", String(take!(io)))
    @test DataFrames.nrow(empty) == 0
    @test names(empty) == names(summary)
    @test eltype(empty.n_simulations) == Int
    @test_throws ArgumentError status_summary(DataFrames.select(df, DataFrames.Not(:failed)))

    # A failure after finishing is not counted as a success; the final
    # checkpoint alone also does not establish successful finalization.
    terminal = DataFrames.DataFrame(lattice_size=[[64], [64], [64]], beta=fill(64.0, 3),
                                   status=["Finished", "Finished", "Checkpoint 4"],
                                   failed=[false, true, false])
    status_summary(terminal; io)
    printed = String(take!(io))
    @test occursin(r"Finished\s+1\s+33\.3%", printed)
    @test occursin(r"Failed\s+1\s+33\.3%", printed)
    @test occursin(r"Remaining\s+1\s+33\.3%", printed)
    @test occursin("Failed (Finished)", printed)

    sourced = copy(df)
    sourced.log_source = ["scratch", "local", "scratch", "scratch", "local",
                          "local", "scratch", "none", "none"]
    original = deepcopy(sourced)
    summary = status_summary(sourced; io)
    printed = String(take!(io))
    @test occursin(r"Local\s+Scratch\s+No log\s+Total", printed)
    @test names(summary) == ["lattice_size", "beta", "status", "failed", "log_source", "n_simulations"]
    @test summary.n_simulations == [1, 1, 1, 1, 1, 1, 1, 2]
    @test isequal(summary[:, 1:5], sourced[1:8, :])
    @test sum(summary.n_simulations) == DataFrames.nrow(sourced)
    @test isequal(sourced, original)

    empty = status_summary(sourced[1:0, :]; io)
    @test occursin(r"Total\s+0\s+0\.0%", String(take!(io)))
    @test DataFrames.nrow(empty) == 0
    @test names(empty) == names(summary)
    @test eltype(empty.n_simulations) == Int
end

@testset "Status manifest and optional snapshots" begin
    mktempdir() do root
        # Reused simulation labels require both key columns, not simulation alone.
        manifest = DataFrames.DataFrame(
            study_id=["study-b", "study-a", "study-a"],
            simulation=["sim1", "sim1", "sim2"],
            config=["job-b/config.yaml", "job-a/config.yaml", "pending/config.yaml"],
            status=fill("stale", 3), failed=trues(3), completed_blocks=fill(-1, 3),
            log_source=fill("scratch", 3))
        manifest[!, "lattice.extent"] = [[16, 16], [32], [64]]
        manifest[!, "system.beta"] = [8.0, 16.0, 32.0]
        manifest[!, "custom.settings"] = [Dict("weights" => [1, 2]) for _ in 1:3]
        native = joinpath(root, "manifest.jld2")
        JLD2.jldsave(native; manifest)
        # A conflicting CSV must not replace the native manifest.
        csv = joinpath(root, "manifest.csv")
        CSV.write(csv, DataFrames.DataFrame(study_id=["ignored"], simulation=["ignored"],
                                          config=["ignored/config.yaml"]))
        for folder in ("job-b", "job-a")
            mkpath(joinpath(root, folder))
            write(joinpath(root, folder, "config.yaml"), "preserved config\n")
        end
        failed_log = joinpath(root, "job-b", "simulation.log")
        success_log = joinpath(root, "job-a", "simulation.log")
        write(failed_log, "TUNING\nSimulation failed: tuning exhausted\n")
        write(success_log, "PRODUCTION\n  100/100 100% 1.000s ◆\nSUMMARY\n  Status SUCCESS\n")
        sources = [native, csv, failed_log, success_log,
                   joinpath(root, "job-b", "config.yaml"), joinpath(root, "job-a", "config.yaml")]
        original = read.(sources)

        existing = readdir(root)
        table = status_manifest(root)
        @test table isa DataFrames.DataFrame
        @test readdir(root) == existing
        @test isequal(status_manifest(root; save=false), table)
        @test readdir(root) == existing

        before = Dates.now(Dates.UTC)
        saved_table = status_manifest(root; save=true)
        path = joinpath(root, only(setdiff(readdir(root), existing)))
        saved = JLD2.load(path)
        @test isequal(saved_table, table)
        @test isequal(saved["manifest"], saved_table)
        @test isabspath(path) && dirname(path) == root
        @test occursin(r"^manifest_status_\d{8}_\d{6}_\d{3}Z\.jld2$", basename(path))
        @test before <= saved["captured_at_utc"] <= Dates.now(Dates.UTC)
        @test basename(path) == "manifest_status_" *
              Dates.format(saved["captured_at_utc"], "yyyymmdd_HHMMSS_sss") * "Z.jld2"
        parameters = DataFrames.Not([:status, :completed_blocks, :failed, :log_source])
        @test isequal(table[:, parameters], manifest[:, parameters])
        @test table.status == ["Tuning", "Finished", "No Log File"]
        @test table.completed_blocks == [0, 100, 0]
        @test table.failed == [true, false, false]
        @test table.log_source == ["local", "local", "none"]
        @test read.(sources) == original

        # A later snapshot reflects new progress and preserves the old snapshot.
        snapshot_bytes = read(path)
        write(failed_log, "MultiWormQMC v0.2.0\nPRODUCTION\n  25/100 25% 1.000s ◆\n")
        existing = readdir(root)
        refreshed = status_manifest(root; save=true)
        next_path = joinpath(root, only(setdiff(readdir(root), existing)))
        @test next_path != path
        @test read(path) == snapshot_bytes
        @test isequal(JLD2.load(next_path, "manifest"), refreshed)
        @test refreshed.status[1] == "Checkpoint 1"
        @test refreshed.completed_blocks[1] == 25
        @test !refreshed.failed[1]
        @test refreshed.log_source == ["local", "local", "none"]

        # Missing IDs and duplicate composite keys must not produce snapshots.
        existing = sort(readdir(root))
        for invalid in (DataFrames.select(manifest, DataFrames.Not(:study_id)),
                        DataFrames.transform(manifest, :study_id => (_ -> ["", "study-a", "study-a"]) => :study_id),
                        DataFrames.transform(manifest, :study_id => (_ -> [missing, "study-a", "study-a"]) => :study_id),
                        vcat(manifest, manifest[1:1, :]))
            JLD2.jldsave(native; manifest=invalid)
            @test_throws ArgumentError status_manifest(root)
            @test_throws ArgumentError status_manifest(root; save=true)
            @test sort(readdir(root)) == existing
        end
    end

    mktempdir() do root
        manifest = DataFrames.DataFrame(study_id=["csv-study"], simulation=["sim1"],
                                        config=["sim1/config.yaml"])
        CSV.write(joinpath(root, "manifest.csv"), manifest)
        existing = readdir(root)
        @test_throws ArgumentError status_manifest(root)
        @test_throws ArgumentError status_manifest(root; save=true)
        @test readdir(root) == existing
        JLD2.jldsave(joinpath(root, "manifest.jld2"); manifest)
        existing = readdir(root)
        table = status_manifest(root)
        @test readdir(root) == existing
        @test isequal(table[:, names(manifest)], manifest)
        @test table.status == ["No Log File"]
        @test table.completed_blocks == [0]
        @test table.failed == [false]
        @test table.log_source == ["none"]
        saved_table = status_manifest(root; save=true)
        path = joinpath(root, only(setdiff(readdir(root), existing)))
        @test isequal(saved_table, table)
        @test isequal(JLD2.load(path, "manifest"), table)
    end
end

@testset "Study status" begin
    mktempdir() do root
        # Deliberately use folders different from the simulation labels.
        logs = [
            nothing,
            "SETUP\n  ✓ Input validated\n",
            "TUNING\n  Tuning 1/10000 | stable 0/5\n",
            "TUNING\nTHERMALIZATION\n  ✓ Thermalization block 1/10 completed\n",
            "PRODUCTION\n  Blocks 100\n  25/100 25% 1.000s\n",
            "PRODUCTION\n  25/100 25% 1.000s ◆\n  48/100 48% 1.000s\n",
            "PRODUCTION\n  50/100 50% 1.000s ◆\n",
            "PRODUCTION\n  75/100 75% 1.000s ◆\n",
            "PRODUCTION\n  100/100 100% 1.000s ◆\nFINALIZATION\n",
            "PRODUCTION\n  100/100 100% 1.000s ◆\nSUMMARY\n  Status SUCCESS\n  Restored blocks 0\n",
            "TUNING\nSimulation failed:\nTuning exhausted max_attempts\n",
            "PRODUCTION\n  100/100 100% 1.000s ◆\nFINALIZATION\nERROR: finalization failed\n",
        ]
        manifest = DataFrames.DataFrame(
            study_id=fill("67db9eed-b978-47f7-a90f-e6604b195df5", length(logs)),
            simulation=["sim$i" for i in eachindex(logs)],
            config=["job$i/config.yaml" for i in eachindex(logs)])
        manifest[!, "lattice.extent"] = [[64] for _ in logs]
        manifest[!, "system.beta"] = fill(64.0, length(logs))
        manifest[!, "run.blocks"] = fill(100, length(logs))
        path = joinpath(root, "manifest.jld2")
        JLD2.jldsave(path; manifest)
        for (i, log) in enumerate(logs)
            log === nothing && continue
            mkpath(joinpath(root, "job$i"))
            write(joinpath(root, "job$i", "simulation.log"), log)
        end
        sources = [path; [joinpath(root, "job$i", "simulation.log") for i in 2:length(logs)]]
        original = read.(sources)
        table = study_status(root)
        @test names(table) == ["study_id", "simulation", "lattice_size", "beta", "status", "completed_blocks", "failed", "log_source"]
        @test table.study_id == manifest.study_id
        @test table.simulation == manifest.simulation
        @test table.lattice_size == fill([64], length(logs))
        @test table.beta == fill(64.0, length(logs))
        @test table.status == ["No Log File", "Pre Tuning", "Tuning", "Thermalizing",
                               "Production Started", "Checkpoint 1", "Checkpoint 2",
                               "Checkpoint 3", "Checkpoint 4", "Finished", "Tuning", "Checkpoint 4"]
        @test table.completed_blocks == [0, 0, 0, 0, 25, 48, 50, 75, 100, 100, 0, 100]
        @test table.failed == [falses(10); true; true]
        @test table.log_source == ["none"; fill("local", length(logs) - 1)]
        @test read.(sources) == original

        # Scratch log text uses the same rules as every local-file fixture.
        for (i, log) in enumerate(logs)
            log === nothing && continue
            @test study_log_status(IOBuffer(log)) ==
                (; status=table.status[i], completed_blocks=table.completed_blocks[i],
                   failed=table.failed[i])
        end
        empty_stream = IOBuffer()
        @test study_log_status(empty_stream) ==
            (; status="Pre Tuning", completed_blocks=0, failed=false)
        @test isopen(empty_stream)
        stream = IOBuffer("TUNING\nERROR: earlier attempt\nTHERMALIZATION\n")
        readline(stream)
        readline(stream)
        @test study_log_status(stream) ==
            (; status="Thermalizing", completed_blocks=0, failed=false)

        logpath = joinpath(root, "job2", "simulation.log")
        # A fresh attempt must not inherit an earlier attempt's state.
        write(logpath, "SUMMARY\n  Status SUCCESS\nSimulation failed:\nMultiWormQMC v0.2.0\nSETUP\nTUNING\n")
        fresh = study_status(root)[2, :]
        @test fresh.status == "Tuning"
        @test fresh.completed_blocks == 0
        @test !fresh.failed

        # Resume output lists only the remaining checkpoints.
        write(logpath, "PRODUCTION\n  Blocks 100\n  Restored blocks 50\n  Checkpoints 75, 100\n  51/100 51% 1.000s\n")
        resumed = study_status(root)[2, :]
        @test resumed.status == "Checkpoint 2"
        @test resumed.completed_blocks == 51

        write(logpath, "PRODUCTION\n  Blocks 100\n  Restored blocks 50\n  Checkpoints Disabled\n")
        @test study_status(root).status[2] == "Production Started"

        # Rounded-up quarters can share one checkpoint for short runs.
        write(logpath, "PRODUCTION\n  1/3 33% 1.000s ◆\n")
        @test study_status(root).status[2] == "Checkpoint 1"
        write(logpath, "PRODUCTION\n  3/3 100% 1.000s ◆\n")
        @test study_status(root).status[2] == "Checkpoint 4"
    end

    mktempdir() do root
        @test_throws ArgumentError study_status(root)
        manifest = DataFrames.DataFrame(study_id=["1cbbd8cf-b1aa-48ad-a822-ea3ac9df77b0"],
                                        simulation=["custom"], config=["nested/config.yaml"])
        manifest[!, "lattice.extent"] = ["- 16\n- 16"]
        manifest[!, "system.beta"] = [8.0]
        CSV.write(joinpath(root, "manifest.csv"), manifest)
        @test_throws ArgumentError study_status(root)
        JLD2.jldsave(joinpath(root, "manifest.jld2"); manifest)
        table = study_status(root)
        @test table.study_id == manifest.study_id
        @test table.lattice_size == [[16, 16]]
        @test table.beta == [8.0]
        @test table.status == ["No Log File"]

        # Optional metadata may be absent from the JLD2 manifest.
        DataFrames.select!(manifest, :simulation, :config)
        JLD2.jldsave(joinpath(root, "manifest.jld2"); manifest)
        @test ismissing(study_status(root).study_id[1])
        @test ismissing(study_status(root).lattice_size[1])
        @test ismissing(study_status(root).beta[1])
        manifest.study_id = [nothing]
        JLD2.jldsave(joinpath(root, "manifest.jld2"); manifest)
        @test ismissing(study_status(root).study_id[1])
        manifest.config = [missing]
        JLD2.jldsave(joinpath(root, "manifest.jld2"); manifest)
        @test_throws ArgumentError study_status(root)

        JLD2.jldsave(joinpath(root, "manifest.jld2"); manifest=DataFrames.DataFrame())
        @test_throws ArgumentError study_status(root)
        JLD2.jldsave(joinpath(root, "manifest.jld2"); manifest=[1, 2])
        @test_throws ArgumentError study_status(root)
    end
end
