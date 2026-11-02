import CSV
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
    @test String(take!(io)) == "Total simulations: 9\nFinished successfully: 5\nFailed: 1\nLeft: 3\n\n"
    @test names(summary) == ["lattice_size", "beta", "status", "failed", "n_simulations"]
    @test summary.n_simulations == [2, 1, 1, 1, 1, 1, 2]
    @test isequal(summary[:, 1:4], df[[1, 3, 4, 5, 6, 7, 8], :])
    @test sum(summary.n_simulations) == DataFrames.nrow(df)
    @test isequal(df, original)

    empty = status_summary(df[1:0, :]; io)
    @test String(take!(io)) == "Total simulations: 0\nFinished successfully: 0\nFailed: 0\nLeft: 0\n\n"
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
    @test String(take!(io)) == "Total simulations: 3\nFinished successfully: 1\nFailed: 1\nLeft: 1\n\n"
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
        @test names(table) == ["simulation", "lattice_size", "beta", "status", "completed_blocks", "failed"]
        @test table.simulation == manifest.simulation
        @test table.lattice_size == fill([64], length(logs))
        @test table.beta == fill(64.0, length(logs))
        @test table.status == ["Pre Tuning", "Pre Tuning", "Tuning", "Thermalizing",
                               "Production Started", "Checkpoint 1", "Checkpoint 2",
                               "Checkpoint 3", "Checkpoint 4", "Finished", "Tuning", "Checkpoint 4"]
        @test table.completed_blocks == [0, 0, 0, 0, 25, 48, 50, 75, 100, 100, 0, 100]
        @test table.failed == [falses(10); true; true]
        @test read.(sources) == original

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
        manifest = DataFrames.DataFrame(simulation=["custom"], config=["nested/config.yaml"])
        manifest[!, "lattice.extent"] = ["- 16\n- 16"]
        manifest[!, "system.beta"] = [8.0]
        CSV.write(joinpath(root, "manifest.csv"), manifest)
        table = study_status(root)
        @test table.lattice_size == [[16, 16]]
        @test table.beta == [8.0]
        @test table.status == ["Pre Tuning"]

        # The typed manifest takes precedence; optional metadata may be absent.
        DataFrames.select!(manifest, :simulation, :config)
        JLD2.jldsave(joinpath(root, "manifest.jld2"); manifest)
        @test ismissing(study_status(root).lattice_size[1])
        @test ismissing(study_status(root).beta[1])
        manifest.config = [missing]
        JLD2.jldsave(joinpath(root, "manifest.jld2"); manifest)
        @test_throws ArgumentError study_status(root)

        JLD2.jldsave(joinpath(root, "manifest.jld2"); manifest=DataFrames.DataFrame())
        @test_throws ArgumentError study_status(root)
        JLD2.jldsave(joinpath(root, "manifest.jld2"); manifest=[1, 2])
        @test_throws ArgumentError study_status(root)
    end
end
