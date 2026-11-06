@testset "Aggregation of runs with saved production blocks" begin
    mktempdir() do root
        labels = ["waiting", "tuning_failed", "one", "partial", "finished",
                  "production_failed", "lagging_log", "empty", "missing_log"]
        manifest = DataFrames.DataFrame(
            study_id=fill("partial-aggregation-fixture", length(labels)),
            simulation=labels,
            config=["job$i/config.yaml" for i in eachindex(labels)],
            output_file=["job$i/results.h5" for i in eachindex(labels)])
        manifest[!, "run.blocks"] = fill(4, length(labels))
        manifest[!, "lattice.extent"] = [[64] for _ in labels]
        manifest[!, "system.beta"] = fill(64.0, length(labels))
        manifest_path = joinpath(root, "manifest.jld2")
        JLD2.jldsave(manifest_path; manifest)
        logs = [
            "SETUP\n",
            "TUNING\nSimulation failed:\nTuning exhausted max_attempts\n",
            "PRODUCTION\n  1/4 25% 0.001s\n",
            "PRODUCTION\n  2/4 50% 0.001s ◆\n",
            "PRODUCTION\n  4/4 100% 0.001s ◆\nFINALIZATION\nSUMMARY\n  Status SUCCESS\n  Restored blocks 0\n",
            "PRODUCTION\n  1/4 25% 0.001s\nSimulation failed:\nInterrupted\n",
            "PRODUCTION\n  1/4 25% 0.001s\n",
            "PRODUCTION\n  1/4 25% 0.001s\n",
            nothing,
        ]
        saved_blocks = [0, 0, 1, 2, 4, 1, 3, 0, 0]
        sources = [manifest_path]
        for i in eachindex(labels)
            folder = joinpath(root, "job$i")
            mkpath(folder)
            config_path = joinpath(root, manifest.config[i])
            write(config_path, AGGREGATION_CONFIG)
            push!(sources, config_path)
            if logs[i] !== nothing
                logpath = joinpath(folder, "simulation.log")
                write(logpath, logs[i])
                push!(sources, logpath)
            end
            path = joinpath(root, manifest.output_file[i])
            if i in (1, 2, 9)
                # These files must never be opened: the logs report zero blocks.
                write(path, "not an HDF5 file")
            else
                n = saved_blocks[i]
                write_combined_fixture(path; density_blocks=n, energy_blocks=n,
                                       winding_blocks=n, histogram_blocks=n, correlation_blocks=n, structure_blocks=n, counterflow_blocks=n)
                write_single_particle_green_fixture(path; counts=fill(10, n),
                                                    zcounts=ones(n), append=true)
            end
            push!(sources, path)
        end
        original = read.(sources)
        table = @test_logs (:warn, r"Density has no saved blocks") (:warn, r"Energy has no saved blocks") (:warn, r"Windings has no saved blocks") (:warn, r"ParticlesHistogram has no saved blocks") (:warn, r"SingleParticleGreenFunction has no saved blocks") (:warn, r"DensityCorrelation has no saved blocks") (:warn, r"StaticStructureFactor has no saved blocks") (:warn, r"CounterflowGreenFunction has no saved blocks") save_aggregate(root; blocking_diagnostics=true)
        @test table.simulation == labels
        @test isequal(table[:, names(manifest)], manifest)
        @test table.status == ["Pre Tuning", "Tuning", "Production Started", "Checkpoint 2",
                               "Finished", "Production Started", "Production Started",
                               "Production Started", "No Log File"]
        @test table.failed == [false, true, false, false, false, true, false, false, false]
        @test table.completed_blocks == [0, 0, 1, 2, 4, 1, 1, 1, 0]
        @test table.log_source == [fill("local", 8); "none"]
        @test table.g2_entry_modes == fill((), length(labels))
        snapshot = status_manifest(root)
        @test isequal(table[:, names(snapshot)], snapshot)
        @test !any(startswith("manifest_status_"), readdir(root))
        for name in ("Density", "Energy", "Windings", "ParticlesHistogram", "SingleParticleGreenFunction", "DensityCorrelation", "StaticStructureFactor", "CounterflowGreenFunction")
            @test all(ismissing, table[[1, 2, 8, 9], name])
            @test [table[i, name]["nblocks"] for i in 3:7] == saved_blocks[3:7]
        end
        # One block remains eligible even with the default min_blocks=10.
        @test Measurements.value.(table.Density[3]["densityBySpecies"]) == [0, 1]
        @test all(isnan, Measurements.uncertainty.(table.Density[3]["densityBySpecies"]))
        @test isnan(Measurements.uncertainty(table.Energy[3]["totalEnergy"]))
        @test all(isnan, Measurements.uncertainty.(table.Windings[3]["squaredWindingBySpecies"]))
        @test all(isnan, Measurements.uncertainty.(table.ParticlesHistogram[3]["totalParticleNumberProbability"]))
        @test all(isnan, Measurements.uncertainty.(table.SingleParticleGreenFunction[3]["grEqualTime"]))
        @test all(isnan, Measurements.uncertainty.(table.DensityCorrelation[3]["full"]))
        # HDF5 histories, rather than a stale log count, determine the summary.
        @test table.Energy[7]["nblocks"] == 3
        @test Measurements.value(table.Energy[7]["totalEnergy"]) ≈ 1.5
        @test table.ParticlesHistogram[4]["zSectorSampleCount"] == 3
        saved = JLD2.load(default_aggregate_path(root), "aggregate")
        metadata_columns = names(snapshot)
        @test isequal(saved[:, metadata_columns], table[:, metadata_columns])
        # Measurements with NaN uncertainties compare unequal, so check the
        # stored values and uncertainties separately for the one-block run.
        for (name, field) in (("Density", "densityBySpecies"), ("Energy", "totalEnergy"),
                              ("Windings", "squaredWindingBySpecies"),
                              ("ParticlesHistogram", "totalParticleNumberProbability"),
                              ("SingleParticleGreenFunction", "grEqualTime"),
                              ("DensityCorrelation", "full"))
            @test isequal(Measurements.value.(saved[3, name][field]),
                          Measurements.value.(table[3, name][field]))
            @test isequal(Measurements.uncertainty.(saved[3, name][field]),
                          Measurements.uncertainty.(table[3, name][field]))
        end
        @test read.(sources) == original
        @test_throws ArgumentError save_aggregate(root; output="job3/simulation.log", overwrite=true)

        # Refresh stale monitoring columns before deciding which runs to read.
        manifest.status = fill("Finished", length(labels))
        manifest.completed_blocks = zeros(Int, length(labels))
        manifest.failed = trues(length(labels))
        manifest.log_source = fill("scratch", length(labels))
        JLD2.jldsave(manifest_path; manifest)
        stale_bytes = read(manifest_path)
        refreshed = @test_logs (:warn, r"Density has no saved blocks") aggregate(
            root; observables=["Density"])
        @test isequal(refreshed[:, names(snapshot)], snapshot)
        @test all(ismissing, refreshed.Density[[1, 2, 8, 9]])
        @test refreshed.Density[7]["nblocks"] == 3
        @test read(manifest_path) == stale_bytes
    end
end
