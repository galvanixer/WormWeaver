import HDF5
import Measurements

@testset "Repeated blocking statistics" begin
    # Long correlated stretches: larger blocks must increase the reported error.
    history = vcat(zeros(32), ones(32))
    result = WormWeaver.blocking_measurement(history, fill(1000, 64); min_blocks=4)
    @test Measurements.value(result) == 0.5
    @test Measurements.uncertainty(result) ≈ 1 / sqrt(12)
    analysis = WormWeaver.blocking_analysis(history, fill(1000, 64); min_blocks=4, diagnostics=true)
    @test analysis.block_widths == [1, 2, 4, 8, 16]
    @test analysis.block_counts == [64, 32, 16, 8, 4]
    @test analysis.standard_errors ≈ [0.5 / sqrt(n - 1) for n in analysis.block_counts]
    @test isequal(analysis.measurement, result)
    plain = WormWeaver.blocking_analysis(history, fill(1000, 64); min_blocks=4)
    @test isequal(plain.measurement, analysis.measurement)
    @test plain.standard_errors === nothing
    @test plain.block_widths === nothing
    @test plain.block_counts === nothing

    # Unequal counts travel with the block means when pairing them.
    result = WormWeaver.blocking_measurement([0.0, 1.0, 3.0, 5.0], [1, 3, 2, 4]; min_blocks=2)
    @test Measurements.value(result) ≈ 2.9
    @test Measurements.uncertainty(result) ≈ 1.72

    # The full-data mean retains an unmatched final block, and even histories
    # never lose their last entry before the first standard-error calculation.
    odd = WormWeaver.blocking_measurement([0.0, 0, 0, 0, 10], ones(Int, 5); min_blocks=2)
    @test Measurements.value(odd) == 2
    @test Measurements.uncertainty(odd) ≈ 2
    even = WormWeaver.blocking_measurement([0.0, 0, 0, 4], ones(Int, 4))
    @test Measurements.value(even) == 1
    @test Measurements.uncertainty(even) ≈ 1

    constant = WormWeaver.blocking_measurement(fill(0.5, 100), fill(1000, 100))
    @test Measurements.uncertainty(constant) == 0
    one = WormWeaver.blocking_measurement([0.5], [1000])
    @test Measurements.value(one) == 0.5
    @test isnan(Measurements.uncertainty(one))
    @test_throws ArgumentError WormWeaver.blocking_measurement([0.5, NaN], [1000, 1000])
    @test_throws ArgumentError WormWeaver.blocking_measurement([0.5], [0])
    @test_throws ArgumentError WormWeaver.blocking_measurement([0.5], [1.5])
    @test_throws DimensionMismatch WormWeaver.blocking_measurement([0.5], [1000, 1000])
    @test_throws ArgumentError WormWeaver.blocking_measurement([0.5], [1000]; min_blocks=1)
end

function write_observable_fixture(file, name, fields, counts; metadata=(;),
                                  reduction="mean")
    group = HDF5.create_group(file, "observables/$name")
    group["sampleCount"] = counts
    md = HDF5.create_group(file, "observables/METADATA/$name")
    for (field, history) in pairs(fields)
        field = string(field)
        group[field] = history
        if reduction !== nothing
            HDF5.attributes(group[field])["reduction"] = reduction
            md["reduction_$field"] = reduction
        end
    end
    for (field, value) in pairs(metadata)
        md[string(field)] = value
    end
end

function write_particle_histogram_fixture(file; nblocks=4, nspecies=2)
    species = nspecies == 1 ? reshape([1 0 1; 0 0 1]', 1, 3, 2) :
        cat([1 0 1 0 0; 2 0 0 0 0], [0 0 1 0 0; 0 0 0 0 1]; dims=3)
    total = nspecies == 1 ? [1 0; 0 0; 1 1] : [1 0; 0 0; 1 0; 0 0; 0 0; 0 0; 0 1]
    repetitions = cld(nblocks, 2)
    write_observable_fixture(file, "ParticlesHistogram", (
        particleNumberBySpeciesCounts=repeat(species, 1, 1, repetitions)[:, :, 1:nblocks],
        totalParticleNumberCounts=repeat(total, 1, repetitions)[:, 1:nblocks],
        zSectorSampleCount=repeat([2, 1], repetitions)[1:nblocks],
    ), repeat([2, 1], repetitions)[1:nblocks]; reduction="sum", metadata=(
        species=collect(1:nspecies),
        species_particle_number_bins=collect(0:(nspecies == 1 ? 2 : 4)),
        total_particle_number_bins=collect(0:(nspecies == 1 ? 2 : 6)),
        max_particle_number_by_species=nspecies == 1 ? [2] : [2, 4],
    ))
end

function write_combined_fixture(path; density_blocks=4, nspecies=2,
                                groups=("Density", "Energy", "Windings", "ParticlesHistogram"),
                                legacy=false, histogram_blocks=4)
    mkpath(dirname(path))
    HDF5.h5open(path, "w") do file
        HDF5.create_group(file, "observables")
        HDF5.create_group(file, "observables/METADATA")
        reduction = legacy ? nothing : "mean"
        if "Density" in groups
            history = [0.0 0 0 4; 1.0 1 1 1][:, 1:density_blocks]
            write_observable_fixture(file, "Density", (
                densityBySpecies=history, particleNumberBySpecies=64 .* history,
                totalParticleNumber=vec(sum(64 .* history; dims=1)),
            ), fill(1000, density_blocks); metadata=(species=[1, 2],), reduction)
        end
        if "Energy" in groups
            history = [0.0, 1, 3, 5]
            write_observable_fixture(file, "Energy", (
                diagonalEnergy=2 .* history, diagonalEnergyPerSite=2 .* history ./ 64,
                kineticEnergy=-history, kineticEnergyPerSite=-history ./ 64,
                totalEnergy=history, totalEnergyPerSite=history ./ 64,
            ), [1, 3, 2, 4]; reduction)
        end
        if "Windings" in groups
            history = zeros(nspecies, 2, 4)
            history[1, 1, 4] = 4
            nspecies == 2 && (history[2, 1, 4] = -4)
            total = dropdims(sum(history; dims=1); dims=1)
            fields = (
                windingBySpecies=history, squaredWindingBySpecies=history .^ 2,
                totalWinding=total, squaredTotalWinding=total .^ 2,
            )
            if nspecies == 2
                correlation = reshape(history[1, :, :] .* history[2, :, :], 1, 2, 4)
                fields = (; fields..., windingCorrelationBySpecies=correlation)
            end
            write_observable_fixture(file, "Windings", fields, fill(1000, 4);
                metadata=(periodic_axes=[true, false],
                          species_pairs=nspecies == 2 ? reshape([1, 2], 1, 2) : zeros(Int, 0, 2)),
                reduction)
        end
        if "ParticlesHistogram" in groups
            write_particle_histogram_fixture(file; nblocks=histogram_blocks, nspecies)
        end
    end
end

@testset "Shared observable aggregation" begin
    mktempdir() do root
        path = joinpath(root, "results.h5")
        write_combined_fixture(path)
        manifest = DataFrames.DataFrame(output_file=["results.h5"], config=["config.yaml"], replica=[1])
        manifest[!, "run.blocks"] = [4]
        manifest_path = joinpath(root, "manifest.jld2")
        JLD2.jldsave(manifest_path; manifest)
        write(joinpath(root, "config.yaml"), "original config")
        write(joinpath(root, "manifest.csv"), "original manifest")
        sources = [path, manifest_path, joinpath(root, "config.yaml"), joinpath(root, "manifest.csv")]
        original_bytes = read.(sources)

        table = aggregate(root; min_blocks=2)
        @test names(table) == [names(manifest); ["Density", "Energy", "Windings", "ParticlesHistogram", "SingleParticleGreenFunction"]]
        @test ismissing(table.SingleParticleGreenFunction[1])
        @test isequal(table[:, names(manifest)], manifest)
        @test table.Density[1]["species"] == [1, 2]
        @test Measurements.value.(table.Density[1]["densityBySpecies"]) == [1, 1]
        energy = table.Energy[1]
        @test energy["totalEnergy"] isa Measurements.Measurement{Float64}
        @test Measurements.value(energy["totalEnergy"]) ≈ 2.9
        @test Measurements.uncertainty(energy["totalEnergy"]) ≈ 1.72
        @test Measurements.value(energy["diagonalEnergy"]) ≈ 5.8
        @test Measurements.value(energy["kineticEnergy"]) ≈ -2.9
        @test Measurements.value(energy["totalEnergyPerSite"]) ≈ 2.9 / 64
        @test Measurements.uncertainty(energy["diagonalEnergyPerSite"]) ≈ 3.44 / 64
        @test Measurements.uncertainty(energy["kineticEnergyPerSite"]) ≈ 1.72 / 64
        @test energy["nsamples"] == 10
        @test energy["nblocks"] == 4
        @test energy["min_blocks"] == 2
        @test energy["error_method"] == "repeated blocking: maximum standard error"
        windings = table.Windings[1]
        @test windings["nsamples"] == 4000
        @test windings["windingBySpecies"] isa Matrix{Measurements.Measurement{Float64}}
        @test Measurements.value.(windings["windingBySpecies"]) == [1 0; -1 0]
        @test Measurements.uncertainty.(windings["windingBySpecies"]) ≈ [1 0; 1 0]
        @test Measurements.value.(windings["squaredWindingBySpecies"]) == [4 0; 4 0]
        @test Measurements.uncertainty.(windings["squaredWindingBySpecies"]) ≈ [4 0; 4 0]
        @test Measurements.value.(windings["totalWinding"]) == [0, 0]
        @test Measurements.uncertainty.(windings["squaredTotalWinding"]) == [0, 0]
        @test Measurements.value.(windings["windingCorrelationBySpecies"]) == [-4 0]
        @test Measurements.uncertainty.(windings["windingCorrelationBySpecies"]) ≈ [4 0]
        @test windings["periodic_axes"] == [true, false]
        @test windings["species_pairs"] == [1 2]

        expanded = observable_table(table.Windings)
        @test DataFrames.nrow(expanded) == DataFrames.nrow(table)
        @test expanded.windingBySpecies[1] === windings["windingBySpecies"]
        @test expanded.species_pairs[1] == [1 2]
        @test expanded.nsamples == [4000]
        @test table.ParticlesHistogram[1]["zSectorSampleCount"] == 6

        saved = JLD2.load(joinpath(root, "aggregated.jld2"), "aggregate")
        @test isequal(saved, table)
        @test read.(sources) == original_bytes
        @test_throws ArgumentError aggregate(root)
        @test isequal(aggregate(root; min_blocks=2, overwrite=true), table)
        subset = aggregate(root; observables=("Windings", "Energy"), output=nothing)
        @test names(subset) == [names(manifest); ["Windings", "Energy"]]
        custom_output = joinpath(root, "nested", "custom.jld2")
        aggregate(root; observables=["Energy"], output=custom_output)
        @test names(JLD2.load(custom_output, "aggregate")) == [names(manifest); ["Energy"]]
        for source in sources
            @test_throws ArgumentError aggregate(root; output=source, overwrite=true)
        end
        symlink(path, joinpath(root, "alias.h5"))
        @test_throws ArgumentError aggregate(root; output="alias.h5", overwrite=true)
        symlink(joinpath(root, "new.jld2"), joinpath(root, "dangling.jld2"))
        @test_throws ArgumentError aggregate(root; output="dangling.jld2", overwrite=true)
        @test_throws ArgumentError aggregate(root; output="nested", overwrite=true)
        @test_throws ArgumentError aggregate(root; output=" ")
        @test read.(sources) == original_bytes

        for selection in (String[], ["Density", "Density"], ["Unknown"], [1], "Density", nothing)
            @test_throws ArgumentError aggregate(root; observables=selection, output=nothing)
        end
        @test_throws ArgumentError aggregate(root; min_blocks=1, output=nothing)
        manifest[!, "Energy"] = [missing]
        JLD2.jldsave(manifest_path; manifest)
        @test_throws ArgumentError aggregate(root; output=nothing)
        @test "Density" in names(aggregate(root; observables=["Density"], output=nothing))
        DataFrames.select!(manifest, DataFrames.Not("Energy"))
        JLD2.jldsave(manifest_path; manifest)

        # Incompleteness is specific to each observable, not the entire file.
        write_combined_fixture(path; density_blocks=2)
        partial = @test_logs (:warn, r"Density run is incomplete") (:warn, r"SingleParticleGreenFunction was not found") aggregate(root; output=nothing)
        @test ismissing(partial.Density[1])
        @test !ismissing(partial.Energy[1])
        @test !ismissing(partial.Windings[1])
        @test !ismissing(partial.ParticlesHistogram[1])
        write_combined_fixture(path; groups=("Energy",))
        absent = @test_logs (:warn, r"Density was not found") (:warn, r"Windings was not found") (:warn, r"ParticlesHistogram was not found") (:warn, r"SingleParticleGreenFunction was not found") aggregate(root; output=nothing)
        @test ismissing(absent.Density[1]) && ismissing(absent.Windings[1])
        @test !ismissing(absent.Energy[1])
        @test ismissing(absent.ParticlesHistogram[1])

        missing_manifest = DataFrames.DataFrame(output_file=[missing, "missing.h5"])
        missing_manifest[!, "run.blocks"] = [4, 4]
        JLD2.jldsave(manifest_path; manifest=missing_manifest)
        absent_files = @test_logs (:warn, r"No results path") (:warn, r"Results file is missing") aggregate(root; output=nothing)
        @test DataFrames.nrow(absent_files) == 2
        @test all(ismissing, absent_files.Density)
        @test all(ismissing, absent_files.Energy)
        @test all(ismissing, absent_files.Windings)
        @test all(ismissing, absent_files.ParticlesHistogram)
        JLD2.jldsave(manifest_path; manifest)

        write_combined_fixture(path; nspecies=1, legacy=true)
        single = aggregate(root; output=nothing).Windings[1]
        @test size(single["windingBySpecies"]) == (1, 2)
        @test !haskey(single, "windingCorrelationBySpecies")
        @test size(single["species_pairs"]) == (0, 2)
        @test Measurements.value.(single["squaredTotalWinding"]) == [4, 0]
    end

    # The shared shape handler also preserves higher-dimensional observable arrays.
    history = repeat(reshape(collect(1.0:8), 2, 2, 2, 1), 1, 1, 1, 4)
    measurements = WormWeaver.observable_measurements(history, ones(Int, 4); min_blocks=2)
    @test size(measurements) == (2, 2, 2)
    @test Measurements.value.(measurements) == reshape(collect(1.0:8), 2, 2, 2)
    @test all(iszero, Measurements.uncertainty.(measurements))
end

@testset "Particle counts and probabilities" begin
    mktempdir() do root
        path = joinpath(root, "results.h5")
        single_path = joinpath(root, "single.h5")
        write_combined_fixture(path)
        write_combined_fixture(single_path; nspecies=1)
        manifest = DataFrames.DataFrame(output_file=["results.h5", "single.h5"], replica=[1, 2])
        manifest[!, "run.blocks"] = [4, 4]
        manifest_path = joinpath(root, "manifest.jld2")
        JLD2.jldsave(manifest_path; manifest)
        before = read.([path, single_path, manifest_path])

        table = aggregate(root; observables=["ParticlesHistogram"], blocking_diagnostics=true)
        @test names(table) == [names(manifest); ["ParticlesHistogram"]]
        @test isequal(table[:, names(manifest)], manifest)
        @test DataFrames.nrow(table) == 2  # Replicas remain separate rows.
        histogram = table.ParticlesHistogram[1]
        @test histogram["particleNumberBySpeciesCounts"] isa Matrix{Int64}
        @test histogram["totalParticleNumberCounts"] isa Vector{Int64}
        @test histogram["particleNumberBySpeciesCounts"] == [2 0 4 0 0; 4 0 0 0 2]
        @test histogram["totalParticleNumberCounts"] == [2, 0, 2, 0, 0, 0, 2]
        @test histogram["zSectorSampleCount"] === Int64(6)
        @test histogram["nsamples"] === Int64(6)
        @test histogram["nblocks"] == 4
        @test histogram["aggregation_method"] == "sum"
        @test histogram["species"] == [1, 2]
        @test histogram["species_particle_number_bins"] == collect(0:4)
        @test histogram["total_particle_number_bins"] == collect(0:6)
        @test histogram["max_particle_number_by_species"] == [2, 4]
        @test histogram["min_blocks"] == 10
        @test histogram["error_method"] == "repeated blocking: maximum standard error"
        @test histogram["particleNumberBySpeciesProbability"] isa Matrix{Measurements.Measurement{Float64}}
        @test histogram["totalParticleNumberProbability"] isa Vector{Measurements.Measurement{Float64}}
        @test Measurements.value.(histogram["particleNumberBySpeciesProbability"]) ≈ [1/3 0 2/3 0 0; 2/3 0 0 0 1/3]
        @test Measurements.value.(histogram["totalParticleNumberProbability"]) ≈ [1/3, 0, 1/3, 0, 0, 0, 1/3]
        @test histogram["meanParticleNumberBySpecies"] isa Vector{Measurements.Measurement{Float64}}
        @test histogram["meanSquaredParticleNumberBySpecies"] isa Vector{Measurements.Measurement{Float64}}
        @test histogram["meanTotalParticleNumber"] isa Measurements.Measurement{Float64}
        @test histogram["meanSquaredTotalParticleNumber"] isa Measurements.Measurement{Float64}
        @test Measurements.value.(histogram["meanParticleNumberBySpecies"]) ≈ [4/3, 4/3]
        @test Measurements.value.(histogram["meanSquaredParticleNumberBySpecies"]) ≈ [8/3, 16/3]
        @test Measurements.value(histogram["meanTotalParticleNumber"]) ≈ 8/3
        @test Measurements.value(histogram["meanSquaredTotalParticleNumber"]) ≈ 40/3
        # The first block samples N=0 and N=2, so its second moment is 2,
        # not the square of its mean (1). Unequal block counts also matter.
        moment_error = sqrt(4/3 * 16/9) / 6
        @test Measurements.uncertainty.(histogram["meanParticleNumberBySpecies"]) ≈ [moment_error, 4moment_error]
        @test Measurements.uncertainty.(histogram["meanSquaredParticleNumberBySpecies"]) ≈ [2moment_error, 16moment_error]
        @test Measurements.uncertainty(histogram["meanTotalParticleNumber"]) ≈ 5moment_error
        @test Measurements.uncertainty(histogram["meanSquaredTotalParticleNumber"]) ≈ 34moment_error
        @test Measurements.value.(histogram["meanSquaredParticleNumberBySpecies"]) !=
            Measurements.value.(histogram["meanParticleNumberBySpecies"]) .^ 2
        @test Measurements.value(histogram["meanSquaredTotalParticleNumber"]) !=
            sum(Measurements.value.(histogram["meanSquaredParticleNumberBySpecies"]))
        # Unequal block sizes require sample weighting, rather than averaging
        # [1/2, 0, 1/2, 0] to obtain 1/4 for total particle-number bin zero.
        @test Measurements.uncertainty(histogram["totalParticleNumberProbability"][1]) ≈ sqrt(4/3 * 4/9) / 6
        @test histogram["blocking_diagnostics"]["block_widths"] == [1]
        @test histogram["blocking_diagnostics"]["block_counts"] == [4]
        @test histogram["blocking_diagnostics"]["sample_counts"] == [6]
        @test size(histogram["blocking_diagnostics"]["particleNumberBySpeciesStandardError"]) == (2, 5, 1)
        @test size(histogram["blocking_diagnostics"]["totalParticleNumberStandardError"]) == (7, 1)
        @test size(histogram["blocking_diagnostics"]["meanParticleNumberBySpeciesStandardError"]) == (2, 1)
        @test size(histogram["blocking_diagnostics"]["meanSquaredParticleNumberBySpeciesStandardError"]) == (2, 1)
        @test size(histogram["blocking_diagnostics"]["meanTotalParticleNumberStandardError"]) == (1,)
        @test size(histogram["blocking_diagnostics"]["meanSquaredTotalParticleNumberStandardError"]) == (1,)
        @test all(iszero, histogram["particleNumberBySpeciesCounts"][1, 4:5])
        @test vec(sum(histogram["particleNumberBySpeciesCounts"]; dims=2)) == [6, 6]
        @test sum(histogram["totalParticleNumberCounts"]) == histogram["nsamples"]
        @test table.ParticlesHistogram[2]["particleNumberBySpeciesCounts"] == [2 0 4]
        @test table.ParticlesHistogram[2]["totalParticleNumberCounts"] == [2, 0, 4]
        single_species = table.ParticlesHistogram[2]
        @test size(single_species["meanParticleNumberBySpecies"]) == (1,)
        @test size(single_species["meanSquaredParticleNumberBySpecies"]) == (1,)
        @test isequal(single_species["meanParticleNumberBySpecies"][1], single_species["meanTotalParticleNumber"])
        @test isequal(single_species["meanSquaredParticleNumberBySpecies"][1], single_species["meanSquaredTotalParticleNumber"])

        saved = JLD2.load(joinpath(root, "aggregated.jld2"))
        @test isequal(saved["aggregate"], table)
        @test saved["error_method"] == "repeated blocking"
        @test saved["min_blocks"] == 10
        @test read.([path, single_path, manifest_path]) == before
        more_levels = aggregate(root; observables=["ParticlesHistogram"], blocking_diagnostics=true, min_blocks=2, output=nothing).ParticlesHistogram[1]
        @test more_levels["particleNumberBySpeciesCounts"] == histogram["particleNumberBySpeciesCounts"]
        @test isequal(more_levels["particleNumberBySpeciesProbability"], histogram["particleNumberBySpeciesProbability"])
        @test more_levels["blocking_diagnostics"]["block_widths"] == [1, 2]
        @test more_levels["blocking_diagnostics"]["block_counts"] == [4, 2]
        @test more_levels["blocking_diagnostics"]["sample_counts"] == [6, 6]
        @test all(iszero, more_levels["blocking_diagnostics"]["totalParticleNumberStandardError"][:, 2])
        flattened = observable_table(table.ParticlesHistogram)
        @test flattened.nsamples == [6, 6]
        @test flattened.particleNumberBySpeciesCounts[1] === histogram["particleNumberBySpeciesCounts"]
        @test flattened.totalParticleNumberProbability[1] === histogram["totalParticleNumberProbability"]
        @test flattened.blocking_diagnostics[1] === histogram["blocking_diagnostics"]
        @test flattened.meanParticleNumberBySpecies[1] === histogram["meanParticleNumberBySpecies"]
        @test flattened.meanSquaredTotalParticleNumber[1] === histogram["meanSquaredTotalParticleNumber"]

        write_combined_fixture(path; histogram_blocks=2)
        partial = @test_logs (:warn, r"ParticlesHistogram run is incomplete") (:warn, r"SingleParticleGreenFunction was not found") (:warn, r"SingleParticleGreenFunction was not found") aggregate(root; output=nothing)
        @test ismissing(partial.ParticlesHistogram[1])
        @test !ismissing(partial.Density[1]) && !ismissing(partial.Energy[1]) && !ismissing(partial.Windings[1])
        @test !ismissing(partial.ParticlesHistogram[2])
        manifest[!, "run.blocks"] = [1, 4]
        JLD2.jldsave(manifest_path; manifest)
        @test_throws ArgumentError aggregate(root; observables=["ParticlesHistogram"], blocking_diagnostics=true, output=nothing)
        write_combined_fixture(path; histogram_blocks=1)
        single_block = aggregate(root; observables=["ParticlesHistogram"], blocking_diagnostics=true, output=nothing).ParticlesHistogram[1]
        @test single_block["nblocks"] == 1
        @test single_block["nsamples"] == 2
        @test single_block["totalParticleNumberCounts"] == [1, 0, 1, 0, 0, 0, 0]
        @test Measurements.value.(single_block["totalParticleNumberProbability"]) == [0.5, 0, 0.5, 0, 0, 0, 0]
        @test all(isnan, Measurements.uncertainty.(single_block["totalParticleNumberProbability"]))
        @test all(isnan, Measurements.uncertainty.(single_block["particleNumberBySpeciesProbability"]))
        @test single_block["blocking_diagnostics"]["block_widths"] == [1]
        @test single_block["blocking_diagnostics"]["block_counts"] == [1]
        @test single_block["blocking_diagnostics"]["sample_counts"] == [2]
        @test all(isnan, single_block["blocking_diagnostics"]["totalParticleNumberStandardError"])
        @test Measurements.value.(single_block["meanParticleNumberBySpecies"]) == [1, 0]
        @test Measurements.value.(single_block["meanSquaredParticleNumberBySpecies"]) == [2, 0]
        @test Measurements.value(single_block["meanTotalParticleNumber"]) == 1
        @test Measurements.value(single_block["meanSquaredTotalParticleNumber"]) == 2
        for field in ("meanParticleNumberBySpecies", "meanSquaredParticleNumberBySpecies")
            @test all(isnan, Measurements.uncertainty.(single_block[field]))
            @test all(isnan, single_block["blocking_diagnostics"][field * "StandardError"])
        end
        for field in ("meanTotalParticleNumber", "meanSquaredTotalParticleNumber")
            @test isnan(Measurements.uncertainty(single_block[field]))
            @test all(isnan, single_block["blocking_diagnostics"][field * "StandardError"])
        end
    end
end

function write_histogram_trajectory_fixture(path, counts, occupations; capacities)
    nspecies, nblocks = size(occupations)
    species = zeros(Int64, nspecies, maximum(capacities) + 1, nblocks)
    total = zeros(Int64, sum(capacities) + 1, nblocks)
    for b in 1:nblocks
        for s in 1:nspecies
            species[s, occupations[s, b] + 1, b] = counts[b]
        end
        total[sum(occupations[:, b]) + 1, b] = counts[b]
    end
    HDF5.h5open(path, "w") do file
        HDF5.create_group(file, "observables")
        HDF5.create_group(file, "observables/METADATA")
        write_observable_fixture(file, "ParticlesHistogram", (
            particleNumberBySpeciesCounts=species, totalParticleNumberCounts=total,
            zSectorSampleCount=counts,
        ), counts; reduction="sum", metadata=(species=collect(1:nspecies),
            species_particle_number_bins=collect(0:maximum(capacities)),
            total_particle_number_bins=collect(0:sum(capacities)),
            max_particle_number_by_species=capacities))
    end
end

@testset "Histogram blocking curves" begin
    mktempdir() do root
        path = joinpath(root, "results.h5")
        manifest = DataFrames.DataFrame(output_file=["results.h5"])
        manifest[!, "run.blocks"] = [64]
        JLD2.jldsave(joinpath(root, "manifest.jld2"); manifest)
        first_species = vcat(zeros(Int, 32), ones(Int, 32))
        occupations = permutedims(hcat(first_species, 1 .- first_species))
        write_histogram_trajectory_fixture(path, fill(1000, 64), occupations; capacities=[1, 2])
        before = read(path)
        histogram = aggregate(root; observables=["ParticlesHistogram"], blocking_diagnostics=true, min_blocks=4).ParticlesHistogram[1]
        diagnostics = histogram["blocking_diagnostics"]
        @test diagnostics["block_widths"] == [1, 2, 4, 8, 16]
        @test diagnostics["block_counts"] == [64, 32, 16, 8, 4]
        @test diagnostics["sample_counts"] == fill(64000, 5)
        @test diagnostics["sample_counts"] isa Vector{Int64}
        @test size(diagnostics["particleNumberBySpeciesStandardError"]) == (2, 3, 5)
        @test size(diagnostics["totalParticleNumberStandardError"]) == (4, 5)
        @test Set(keys(diagnostics)) == Set(["block_widths", "block_counts", "sample_counts",
            "particleNumberBySpeciesStandardError", "totalParticleNumberStandardError",
            "meanParticleNumberBySpeciesStandardError", "meanTotalParticleNumberStandardError",
            "meanSquaredParticleNumberBySpeciesStandardError", "meanSquaredTotalParticleNumberStandardError"])
        expected_errors = [0.5 / sqrt(n - 1) for n in [64, 32, 16, 8, 4]]
        @test diagnostics["particleNumberBySpeciesStandardError"][1, 2, :] ≈ expected_errors
        @test all(>(0), diff(diagnostics["particleNumberBySpeciesStandardError"][1, 2, :]))
        @test Measurements.uncertainty(histogram["particleNumberBySpeciesProbability"][1, 2]) ≈ 1 / sqrt(12)
        @test Measurements.value.(histogram["particleNumberBySpeciesProbability"]) == [0.5 0.5 0; 0.5 0.5 0]
        @test vec(sum(Measurements.value.(histogram["particleNumberBySpeciesProbability"]); dims=2)) == [1, 1]
        @test histogram["particleNumberBySpeciesCounts"] == [32000 32000 0; 32000 32000 0]
        # Species populations change, but their total is constant; use its
        # saved histogram directly, preserving this correlation.
        @test Measurements.value.(histogram["totalParticleNumberProbability"]) == [0, 1, 0, 0]
        @test all(iszero, Measurements.uncertainty.(histogram["totalParticleNumberProbability"]))
        @test all(iszero, diagnostics["totalParticleNumberStandardError"])
        @test all(iszero, diagnostics["particleNumberBySpeciesStandardError"][:, 3, :])
        @test Measurements.value.(histogram["meanParticleNumberBySpecies"]) == [0.5, 0.5]
        @test Measurements.value.(histogram["meanSquaredParticleNumberBySpecies"]) == [0.5, 0.5]
        @test Measurements.uncertainty.(histogram["meanParticleNumberBySpecies"]) ≈ fill(1/sqrt(12), 2)
        @test Measurements.value(histogram["meanTotalParticleNumber"]) == 1
        @test Measurements.value(histogram["meanSquaredTotalParticleNumber"]) == 1
        @test Measurements.uncertainty(histogram["meanTotalParticleNumber"]) == 0
        @test Measurements.uncertainty(histogram["meanSquaredTotalParticleNumber"]) == 0
        @test diagnostics["meanParticleNumberBySpeciesStandardError"][1, :] ≈ expected_errors
        @test diagnostics["meanSquaredParticleNumberBySpeciesStandardError"][1, :] ≈ expected_errors
        @test all(iszero, diagnostics["meanTotalParticleNumberStandardError"])
        @test all(iszero, diagnostics["meanSquaredTotalParticleNumberStandardError"])
        @test read(path) == before
        @test isequal(JLD2.load(joinpath(root, "aggregated.jld2"), "aggregate").ParticlesHistogram[1], histogram)
        for (probability_field, error_field) in (
            ("particleNumberBySpeciesProbability", "particleNumberBySpeciesStandardError"),
            ("totalParticleNumberProbability", "totalParticleNumberStandardError"),
        )
            errors = diagnostics[error_field]
            @test Measurements.uncertainty.(histogram[probability_field]) ≈
                dropdims(maximum(errors; dims=ndims(errors)); dims=ndims(errors))
        end
        for field in ("meanParticleNumberBySpecies", "meanSquaredParticleNumberBySpecies",
                      "meanTotalParticleNumber", "meanSquaredTotalParticleNumber")
            curve = diagnostics[field * "StandardError"]
            maximum_error = dropdims(maximum(curve; dims=ndims(curve)); dims=ndims(curve))
            expected = ndims(curve) == 1 ? only(maximum_error) : maximum_error
            @test Measurements.uncertainty.(histogram[field]) ≈ expected
        end

        # The last block contains the only N=1 events; it must stay in the mean
        # even though both coarser error levels omit those events.
        counts = collect(Int64, 1:11)
        write_histogram_trajectory_fixture(path, counts, reshape([zeros(Int, 10); 1], 1, 11); capacities=[1])
        manifest[!, "run.blocks"] = [11]
        JLD2.jldsave(joinpath(root, "manifest.jld2"); manifest)
        odd = aggregate(root; observables=["ParticlesHistogram"], blocking_diagnostics=true, min_blocks=2, output=nothing).ParticlesHistogram[1]
        @test odd["totalParticleNumberCounts"] == [55, 11]
        @test Measurements.value.(odd["totalParticleNumberProbability"]) ≈ [5/6, 1/6]
        @test odd["blocking_diagnostics"]["block_widths"] == [1, 2, 4]
        @test odd["blocking_diagnostics"]["block_counts"] == [11, 5, 2]
        @test odd["blocking_diagnostics"]["sample_counts"] == [66, 55, 36]
        @test all(iszero, odd["blocking_diagnostics"]["totalParticleNumberStandardError"][:, 2:3])
        @test Measurements.uncertainty(odd["totalParticleNumberProbability"][2]) > 0
        @test Measurements.uncertainty(odd["totalParticleNumberProbability"][2]) ==
            maximum(odd["blocking_diagnostics"]["totalParticleNumberStandardError"][2, :])
        @test Measurements.value.(odd["meanParticleNumberBySpecies"]) ≈ [1/6]
        @test Measurements.value.(odd["meanSquaredParticleNumberBySpecies"]) ≈ [1/6]
        @test Measurements.value(odd["meanTotalParticleNumber"]) ≈ 1/6
        @test Measurements.value(odd["meanSquaredTotalParticleNumber"]) ≈ 1/6
        for field in ("meanParticleNumberBySpecies", "meanSquaredParticleNumberBySpecies",
                      "meanTotalParticleNumber", "meanSquaredTotalParticleNumber")
            curve = odd["blocking_diagnostics"][field * "StandardError"]
            @test all(iszero, selectdim(curve, ndims(curve), 2:3))
            @test all(>(0), selectdim(curve, ndims(curve), 1))
        end

        # Diagnostic sample totals must remain exact above Float64's integer precision.
        counts = [Int64(2)^53 + 1, 1, 1, 1, 1]
        write_histogram_trajectory_fixture(path, counts, zeros(Int, 1, 5); capacities=[1])
        manifest[!, "run.blocks"] = [5]
        JLD2.jldsave(joinpath(root, "manifest.jld2"); manifest)
        large = aggregate(root; observables=["ParticlesHistogram"], blocking_diagnostics=true, min_blocks=2, output=nothing).ParticlesHistogram[1]
        @test large["nsamples"] === Int64(2)^53 + 5
        @test large["blocking_diagnostics"]["sample_counts"] == [Int64(2)^53 + 5, Int64(2)^53 + 4]
        @test large["totalParticleNumberCounts"] == [Int64(2)^53 + 5, 0]
        @test Measurements.value.(large["totalParticleNumberProbability"]) == [1, 0]
        @test all(iszero, Measurements.uncertainty.(large["totalParticleNumberProbability"]))
        for field in ("meanParticleNumberBySpecies", "meanSquaredParticleNumberBySpecies",
                      "meanTotalParticleNumber", "meanSquaredTotalParticleNumber")
            @test all(iszero, Measurements.value.(large[field]))
            @test all(iszero, Measurements.uncertainty.(large[field]))
        end
    end

    # Exercise large bin weights without allocating billions of empty bins.
    large_bins = [0, Int64(2)^32]
    history = cat([0 2], [0 3]; dims=3)
    result = WormWeaver.histogram_moment_statistics(history, [0 5], large_bins, [2, 3], 5, 2; min_blocks=2)
    @test Measurements.value(result.measurements[1]) == 2.0^64
    @test Measurements.uncertainty(result.measurements[1]) == 0
    @test result.statistics === nothing
end

@testset "Particle histogram validation" begin
    mktempdir() do root
        path = joinpath(root, "results.h5")
        manifest = DataFrames.DataFrame(output_file=["results.h5"])
        manifest[!, "run.blocks"] = [4]
        JLD2.jldsave(joinpath(root, "manifest.jld2"); manifest)
        for (dataset, value, expected_error) in (
            ("sampleCount", [2.0, 1, 2, 1], ArgumentError),
            ("sampleCount", [2, 0, 2, 1], ArgumentError),
            ("zSectorSampleCount", [2, 2, 2, 1], ArgumentError),
            ("zSectorSampleCount", [2.0, 1, 2, 1], ArgumentError),
            ("zSectorSampleCount", reshape([2, 1, 2, 1], 1, 4), DimensionMismatch),
            ("totalParticleNumberCounts", zeros(Int, 7, 4), ArgumentError),
            ("totalParticleNumberCounts", fill(-1, 7, 4), ArgumentError),
            ("totalParticleNumberCounts", zeros(Float64, 7, 4), ArgumentError),
            ("totalParticleNumberCounts", fill(false, 7, 4), ArgumentError),
            ("totalParticleNumberCounts", zeros(Int, 7, 3), DimensionMismatch),
            ("totalParticleNumberCounts", zeros(Int, 1, 7, 4), DimensionMismatch),
            ("particleNumberBySpeciesCounts", zeros(Int, 2, 5, 4), ArgumentError),
            ("particleNumberBySpeciesCounts", zeros(Int, 2, 4, 4), DimensionMismatch),
        )
            write_combined_fixture(path)
            HDF5.h5open(path, "r+") do file
                HDF5.delete_object(file, "observables/ParticlesHistogram/$dataset")
                file["observables/ParticlesHistogram/$dataset"] = value
                dataset != "sampleCount" && (HDF5.attributes(file["observables/ParticlesHistogram/$dataset"])["reduction"] = "sum")
            end
            @test_throws expected_error aggregate(root; observables=["ParticlesHistogram"], output=nothing)
        end
        for field in ("particleNumberBySpeciesCounts", "totalParticleNumberCounts", "zSectorSampleCount",
                      "sampleCount")
            write_combined_fixture(path)
            HDF5.h5open(path, "r+") do file
                HDF5.delete_object(file, "observables/ParticlesHistogram/$field")
            end
            @test_throws ArgumentError aggregate(root; observables=["ParticlesHistogram"], output=nothing)
        end
        for field in ("species", "species_particle_number_bins", "total_particle_number_bins",
                      "max_particle_number_by_species")
            write_combined_fixture(path)
            HDF5.h5open(path, "r+") do file
                HDF5.delete_object(file, "observables/METADATA/ParticlesHistogram/$field")
            end
            @test_throws ArgumentError aggregate(root; observables=["ParticlesHistogram"], output=nothing)
        end
        for (field, value, expected_error) in (
            ("species", [1], DimensionMismatch),
            ("species", [1, 1], ArgumentError),
            ("species", [0, 2], ArgumentError),
            ("species_particle_number_bins", collect(0:3), DimensionMismatch),
            ("species_particle_number_bins", [0, 1, 2, 4, 5], ArgumentError),
            ("total_particle_number_bins", collect(0:5), DimensionMismatch),
            ("max_particle_number_by_species", [2], DimensionMismatch),
            ("max_particle_number_by_species", [2.0, 4], ArgumentError),
            ("max_particle_number_by_species", [1, 4], ArgumentError),
            ("max_particle_number_by_species", [2, 5], ArgumentError),
        )
            write_combined_fixture(path)
            HDF5.h5open(path, "r+") do file
                HDF5.delete_object(file, "observables/METADATA/ParticlesHistogram/$field")
                file["observables/METADATA/ParticlesHistogram/$field"] = value
            end
            @test_throws expected_error aggregate(root; observables=["ParticlesHistogram"], output=nothing)
        end
        for location in (:attribute, :metadata, :missing, :attribute_only, :metadata_only)
            write_combined_fixture(path)
            HDF5.h5open(path, "r+") do file
                dataset = file["observables/ParticlesHistogram/totalParticleNumberCounts"]
                mdpath = "observables/METADATA/ParticlesHistogram/reduction_totalParticleNumberCounts"
                if location in (:attribute, :missing, :metadata_only)
                    HDF5.delete_attribute(dataset, "reduction")
                    location == :attribute && (HDF5.attributes(dataset)["reduction"] = "mean")
                end
                if location in (:metadata, :missing, :attribute_only)
                    HDF5.delete_object(file, mdpath)
                    location == :metadata && (file[mdpath] = "mean")
                end
            end
            if location in (:attribute_only, :metadata_only)
                @test aggregate(root; observables=["ParticlesHistogram"], output=nothing).ParticlesHistogram[1]["nsamples"] == 6
            else
                @test_throws ArgumentError aggregate(root; observables=["ParticlesHistogram"], output=nothing)
            end
        end

        # Each block is individually valid, but its accumulated bin count overflows Int64.
        HDF5.h5open(path, "w") do file
            HDF5.create_group(file, "observables")
            HDF5.create_group(file, "observables/METADATA")
            write_observable_fixture(file, "ParticlesHistogram", (
                particleNumberBySpeciesCounts=reshape([typemax(Int64), 1], 1, 1, 2),
                totalParticleNumberCounts=reshape([typemax(Int64), 1], 1, 2),
                zSectorSampleCount=[typemax(Int64), 1],
            ), [typemax(Int64), 1]; reduction="sum", metadata=(species=[1],
                species_particle_number_bins=[0], total_particle_number_bins=[0],
                max_particle_number_by_species=[0]))
        end
        manifest[!, "run.blocks"] = [2]
        JLD2.jldsave(joinpath(root, "manifest.jld2"); manifest)
        @test_throws OverflowError aggregate(root; observables=["ParticlesHistogram"], output=nothing)
        # Individual bins fit, but summing bins inside a block also must be checked.
        HDF5.h5open(path, "r+") do file
            HDF5.delete_object(file, "observables/ParticlesHistogram/totalParticleNumberCounts")
            file["observables/ParticlesHistogram/totalParticleNumberCounts"] = [typemax(Int64) 1; 1 0]
            HDF5.attributes(file["observables/ParticlesHistogram/totalParticleNumberCounts"])["reduction"] = "sum"
            HDF5.delete_object(file, "observables/METADATA/ParticlesHistogram/total_particle_number_bins")
            file["observables/METADATA/ParticlesHistogram/total_particle_number_bins"] = [0, 1]
        end
        @test_throws OverflowError aggregate(root; observables=["ParticlesHistogram"], output=nothing)
    end
end

@testset "Observable dictionary columns to tables" begin
    winding = reshape(Measurements.measurement.([0.006, -0.033, -0.016], [0.0086, 0.03, 0.024]), 1, 3)
    by_species = reshape(Measurements.measurement.([-0.0031, 0.0091], [0.0028, 0.0088]), 2, 1)
    summaries = Union{Missing,Dict{String,Any}}[
        Dict("totalWinding" => winding[:, 1], "windingBySpecies" => by_species,
             "nsamples" => 100000, "species_pairs" => [1 2], "label" => nothing),
        missing,
        Dict("totalWinding" => winding[:, 2], "nsamples" => 100000),
        Dict("totalWinding" => winding[:, 3], "nsamples" => 100000,
             "windingCorrelationBySpecies" => reshape([Measurements.measurement(0.0019, 0.0022)], 1, 1)),
    ]
    original = deepcopy(summaries)
    table = observable_table(summaries)
    @test DataFrames.nrow(table) == 4
    @test names(table) == ["label", "nsamples", "species_pairs", "totalWinding",
                          "windingBySpecies", "windingCorrelationBySpecies"]
    @test table.totalWinding[1] === summaries[1]["totalWinding"]
    @test Measurements.value.(table.totalWinding[1]) == [0.006]
    @test Measurements.uncertainty.(table.totalWinding[1]) == [0.0086]
    @test table.windingBySpecies[1] === by_species
    @test size(table.windingBySpecies[1]) == (2, 1)
    @test table.species_pairs[1] == [1 2]
    @test all(ismissing, collect(table[2, :]))
    @test ismissing(table.windingBySpecies[3])
    @test ismissing(table.windingCorrelationBySpecies[1])
    @test Measurements.value.(table.windingCorrelationBySpecies[4]) == reshape([0.0019], 1, 1)
    @test table.label[1] === nothing
    @test isequal(summaries, original)

    selected = observable_table(summaries; fields=("totalWinding", "nsamples", "nblocks"))
    @test names(selected) == ["totalWinding", "nsamples", "nblocks"]
    @test DataFrames.nrow(selected) == 4
    @test all(ismissing, selected.nblocks)
    @test isequal(selected.nsamples, [100000, missing, 100000, 100000])
    all_missing = observable_table(fill(missing, 3); fields=["totalEnergy"])
    @test DataFrames.nrow(all_missing) == 3
    @test all(ismissing, all_missing.totalEnergy)
    @test size(observable_table(Dict{String,Any}[])) == (0, 0)
    @test size(observable_table(Missing[]; fields=["totalEnergy"])) == (0, 1)
    energy = Measurements.measurement(2.9, 1.72)
    scalar_table = observable_table([Dict("totalEnergy" => energy)])
    @test scalar_table.totalEnergy[1] === energy
    @test names(observable_table([Dict{String,Any}(), Dict("x" => 1)])) == ["x"]
    for input in (fill(missing, 3), [Dict{String,Any}()], [1], [Dict(:x => 1)])
        @test_throws ArgumentError observable_table(input)
    end
    for fields in (String[], ["x", "x"], [1], "x")
        @test_throws ArgumentError observable_table(summaries; fields)
    end
end

@testset "Observable history validation" begin
    mktempdir() do root
        path = joinpath(root, "results.h5")
        manifest = DataFrames.DataFrame(output_file=["results.h5"])
        manifest[!, "run.blocks"] = [4]
        JLD2.jldsave(joinpath(root, "manifest.jld2"); manifest)
        for (dataset, value, expected_error) in (
            ("sampleCount", [1, 0, 1, 1], ArgumentError),
            ("sampleCount", [1.0, 1.5, 1, 1], ArgumentError),
            ("sampleCount", [1.0, NaN, 1, 1], ArgumentError),
            ("sampleCount", ones(Int, 2, 2), ArgumentError),
            ("sampleCount", ["one", "two", "three", "four"], ArgumentError),
            ("totalEnergy", [0.0, NaN, 3, 5], ArgumentError),
            ("totalEnergy", ["a", "b", "c", "d"], ArgumentError),
            ("totalEnergy", [0.0, 1, 3], DimensionMismatch),
            ("totalEnergy", 1.0, DimensionMismatch),
        )
            write_combined_fixture(path)
            HDF5.h5open(path, "r+") do file
                HDF5.delete_object(file, "observables/Energy/$dataset")
                file["observables/Energy/$dataset"] = value
            end
            @test_throws expected_error aggregate(root; output=nothing)
        end
        for field in ("sampleCount", "diagonalEnergy")
            write_combined_fixture(path)
            HDF5.h5open(path, "r+") do file
                HDF5.delete_object(file, "observables/Energy/$field")
            end
            @test_throws ArgumentError aggregate(root; output=nothing)
        end
        for location in (:attribute, :metadata)
            write_combined_fixture(path)
            HDF5.h5open(path, "r+") do file
                if location == :attribute
                    HDF5.delete_attribute(file["observables/Energy/totalEnergy"], "reduction")
                    HDF5.attributes(file["observables/Energy/totalEnergy"])["reduction"] = "sum"
                else
                    HDF5.delete_object(file, "observables/METADATA/Energy/reduction_totalEnergy")
                    file["observables/METADATA/Energy/reduction_totalEnergy"] = "sum"
                end
            end
            @test_throws ArgumentError aggregate(root; output=nothing)
        end
        write_combined_fixture(path)
        HDF5.h5open(path, "r+") do file
            HDF5.delete_object(file, "observables/Windings/windingCorrelationBySpecies")
            file["observables/Windings/windingCorrelationBySpecies"] = zeros(1, 2, 3)
        end
        @test_throws DimensionMismatch aggregate(root; output=nothing)
        write_combined_fixture(path; density_blocks=2)
        HDF5.h5open(path, "r+") do file
            HDF5.delete_object(file, "observables/Density/sampleCount")
            file["observables/Density/sampleCount"] = [0, 1]
        end
        @test_throws ArgumentError aggregate(root; output=nothing)
        write_combined_fixture(path)
        manifest[!, "run.blocks"] = [3]
        JLD2.jldsave(joinpath(root, "manifest.jld2"); manifest)
        @test_throws ArgumentError aggregate(root; output=nothing)
        for expected_blocks in (0, true, 4.0, missing)
            manifest[!, "run.blocks"] = [expected_blocks]
            JLD2.jldsave(joinpath(root, "manifest.jld2"); manifest)
            @test_throws ArgumentError aggregate(root; output=nothing)
        end
    end
end

function write_density_fixture(path, history, counts; reduction="mean")
    mkpath(dirname(path))
    HDF5.h5open(path, "w") do file
        obs = HDF5.create_group(file, "observables")
        group = HDF5.create_group(obs, "Density")
        group["sampleCount"] = counts
        group["densityBySpecies"] = history
        group["particleNumberBySpecies"] = 64 .* history
        group["totalParticleNumber"] = vec(sum(64 .* history; dims=1))
        for field in ("densityBySpecies", "particleNumberBySpecies", "totalParticleNumber")
            HDF5.attributes(group[field])["reduction"] = reduction
        end
        metadata = HDF5.create_group(obs, "METADATA")
        density_metadata = HDF5.create_group(metadata, "Density")
        density_metadata["species"] = collect(1:size(history, 1))
        density_metadata["reduction_densityBySpecies"] = reduction
    end
end

@testset "Density study aggregation and Measurements storage" begin
    mktempdir() do root
        history = [0.0 0 0 4; 1.0 1 1 1]
        path = joinpath(root, "sim1", "results.h5")
        write_density_fixture(path, history, fill(1000, 4))
        incomplete = joinpath(root, "sim2", "results.h5")
        write_density_fixture(incomplete, history[:, 1:2], fill(1000, 2))
        absolute = joinpath(root, "absolute.h5")
        write_density_fixture(absolute, history, fill(1000, 4))
        manifest = DataFrames.DataFrame(
            project=fill("LRBHQM", 4), study_id=fill("fixture", 4),
            simulation=["sim1", "sim2", "sim3", "sim4"],
            output_file=["sim1/results.h5", "sim2/results.h5", "sim3/results.h5", absolute],
        )
        manifest[!, "run.blocks"] = fill(4, 4)
        manifest[!, "model.hopping.J"] = [[0.62, 0.62] for _ in 1:4]
        JLD2.jldsave(joinpath(root, "manifest.jld2"); manifest)
        source_bytes = read(path)
        aggregate = @test_logs (:warn, r"incomplete") (:warn, r"missing") WormWeaver.aggregate(root; observables=["Density"])
        @test DataFrames.nrow(aggregate) == 4
        @test isequal(aggregate[:, names(manifest)], manifest)
        @test ismissing(aggregate.Density[2])
        @test ismissing(aggregate.Density[3])
        summary = aggregate.Density[1]
        density = summary["densityBySpecies"]
        @test density isa Vector{Measurements.Measurement{Float64}}
        @test Measurements.value.(density) == [1.0, 1.0]
        @test Measurements.uncertainty.(density) ≈ [1.0, 0.0]
        @test Measurements.value.(summary["particleNumberBySpecies"]) == [64.0, 64.0]
        @test Measurements.value(summary["totalParticleNumber"]) == 128
        @test Measurements.uncertainty(summary["totalParticleNumber"]) ≈ 64
        @test summary["nblocks"] == 4
        @test summary["nsamples"] == 4000
        @test summary["species"] == [1, 2]
        @test Measurements.value.(aggregate.Density[4]["densityBySpecies"]) == [1.0, 1.0]

        saved = JLD2.load(joinpath(root, "aggregated.jld2"), "aggregate")
        @test Measurements.value.(saved.Density[1]["densityBySpecies"]) == [1.0, 1.0]
        @test Measurements.uncertainty.(saved.Density[1]["densityBySpecies"]) ≈ [1.0, 0.0]
        @test read(path) == source_bytes
        @test_throws ArgumentError WormWeaver.aggregate(root; observables=["Density"])
        @test_logs (:warn, r"incomplete") (:warn, r"missing") WormWeaver.aggregate(root; observables=["Density"], overwrite=true)
        @test_logs (:warn, r"incomplete") (:warn, r"missing") WormWeaver.aggregate(root; observables=["Density"], output=nothing)
        @test_throws ArgumentError WormWeaver.aggregate(root; observables=["Density"], output="manifest.jld2", overwrite=true)
        @test_throws ArgumentError WormWeaver.aggregate(root; observables=["Density"], output="sim1/results.h5", overwrite=true)

        # Fail on malformed histories instead of silently deleting measurements.
        write_density_fixture(path, history, [1000, 1000, 1000])
        @test_throws DimensionMismatch WormWeaver.aggregate(root; observables=["Density"], output=nothing)
        write_density_fixture(path, history, fill(1000, 4); reduction="sum")
        @test_throws ArgumentError WormWeaver.aggregate(root; observables=["Density"], output=nothing)
    end
end

@testset "Optional diagnostics for all observables" begin
    mktempdir() do root
        path = joinpath(root, "results.h5")
        write_combined_fixture(path)
        manifest = DataFrames.DataFrame(output_file=["results.h5"])
        manifest[!, "run.blocks"] = [4]
        manifest_path = joinpath(root, "manifest.jld2")
        JLD2.jldsave(manifest_path; manifest)
        before = read.([path, manifest_path])
        plain = aggregate(root; min_blocks=2, output=nothing)
        disabled = aggregate(root; min_blocks=2, blocking_diagnostics=false, output=nothing)
        @test isequal(plain, disabled)
        enabled = aggregate(root; min_blocks=2, blocking_diagnostics=true)
        for observable in ("Density", "Energy", "Windings", "ParticlesHistogram")
            without_diagnostics = plain[1, observable]
            with_diagnostics = enabled[1, observable]
            @test !haskey(without_diagnostics, "blocking_diagnostics")
            @test haskey(with_diagnostics, "blocking_diagnostics")
            retained = copy(with_diagnostics)
            delete!(retained, "blocking_diagnostics")
            @test isequal(retained, without_diagnostics)
            diagnostics = with_diagnostics["blocking_diagnostics"]
            @test diagnostics["block_widths"] == [1, 2]
            @test diagnostics["block_counts"] == [4, 2]
        end
        density = enabled.Density[1]["blocking_diagnostics"]
        @test density["sample_counts"] == [4000, 4000]
        @test size(density["densityBySpeciesStandardError"]) == (2, 2)
        @test size(density["totalParticleNumberStandardError"]) == (2,)
        @test density["densityBySpeciesStandardError"][1, :] ≈ [1, 1]
        energy = enabled.Energy[1]["blocking_diagnostics"]
        @test energy["sample_counts"] == [10, 10]
        @test size(energy["totalEnergyStandardError"]) == (2,)
        @test energy["totalEnergyStandardError"][2] ≈ 1.72
        @test energy["totalEnergyStandardError"][1] < energy["totalEnergyStandardError"][2]
        windings = enabled.Windings[1]["blocking_diagnostics"]
        @test windings["sample_counts"] == [4000, 4000]
        @test size(windings["windingBySpeciesStandardError"]) == (2, 2, 2)
        @test size(windings["totalWindingStandardError"]) == (2, 2)
        @test enabled.ParticlesHistogram[1]["blocking_diagnostics"]["sample_counts"] == [6, 6]
        for (observable, fields) in (
            ("Density", ("densityBySpecies", "particleNumberBySpecies", "totalParticleNumber")),
            ("Energy", ("diagonalEnergy", "diagonalEnergyPerSite", "kineticEnergy", "kineticEnergyPerSite", "totalEnergy", "totalEnergyPerSite")),
            ("Windings", ("windingBySpecies", "squaredWindingBySpecies", "totalWinding", "squaredTotalWinding", "windingCorrelationBySpecies",
                          "squaredCoflowBySpeciesPair", "squaredCounterflowBySpeciesPair")),
        )
            summary = enabled[1, observable]
            for field in fields
                curve = summary["blocking_diagnostics"][field * "StandardError"]
                max_error = dropdims(maximum(curve; dims=ndims(curve)); dims=ndims(curve))
                expected = ndims(curve) == 1 ? only(max_error) : max_error
                @test Measurements.uncertainty.(summary[field]) ≈ expected
            end
        end
        saved = JLD2.load(joinpath(root, "aggregated.jld2"), "aggregate")
        @test isequal(saved, enabled)
        @test observable_table(enabled.Energy).blocking_diagnostics[1] === energy
        @test read.([path, manifest_path]) == before

        write_combined_fixture(path; nspecies=1)
        single_species = aggregate(root; observables=["Windings"], blocking_diagnostics=true, output=nothing).Windings[1]
        @test !haskey(single_species["blocking_diagnostics"], "windingCorrelationBySpeciesStandardError")

        # Every observable supports a single saved block, including scalar
        # energy fields and matrix-valued winding fields.
        HDF5.h5open(path, "r+") do file
            for observable in ("Density", "Energy", "Windings", "ParticlesHistogram")
                group = file["observables/$observable"]
                for field in keys(group)
                    dataset = group[field]
                    history = read(dataset)
                    first_block = Array(selectdim(history, ndims(history), 1:1))
                    HDF5.delete_object(group, field)
                    group[field] = first_block
                    if field != "sampleCount"
                        HDF5.attributes(group[field])["reduction"] = observable == "ParticlesHistogram" ? "sum" : "mean"
                    end
                end
            end
        end
        manifest[!, "run.blocks"] = [1]
        JLD2.jldsave(manifest_path; manifest)
        single_block = aggregate(root; blocking_diagnostics=true, output=nothing)
        for observable in ("Density", "Energy", "Windings", "ParticlesHistogram")
            diagnostics = single_block[1, observable]["blocking_diagnostics"]
            @test diagnostics["block_widths"] == [1]
            @test diagnostics["block_counts"] == [1]
            for (key, values) in diagnostics
                if endswith(key, "StandardError")
                    @test all(isnan, values)
                end
            end
        end

        # Odd histories retain all samples in their mean and exact prefix
        # counts in diagnostics for density as well as histograms.
        counts = collect(Int64, 1:11)
        write_density_fixture(path, reshape([zeros(10); 4.0], 1, 11), counts)
        manifest[!, "run.blocks"] = [11]
        JLD2.jldsave(manifest_path; manifest)
        odd = aggregate(root; observables=["Density"], min_blocks=2, blocking_diagnostics=true, output=nothing).Density[1]
        @test Measurements.value(odd["densityBySpecies"][1]) ≈ 2/3
        @test odd["blocking_diagnostics"]["block_widths"] == [1, 2, 4]
        @test odd["blocking_diagnostics"]["block_counts"] == [11, 5, 2]
        @test odd["blocking_diagnostics"]["sample_counts"] == [66, 55, 36]
        @test all(iszero, odd["blocking_diagnostics"]["densityBySpeciesStandardError"][:, 2:3])
    end
end

# Each block represents equal numbers of samples W and -W. First moments vanish
# while the squared moments and products remain nonzero, unlike squared means.
function winding_channel_fixture(windings)
    nspecies, naxes, nblocks = size(windings)
    pairs = zeros(Int, nspecies * (nspecies - 1) ÷ 2, 2)
    products = zeros(Float64, size(pairs, 1), naxes, nblocks)
    for (p, (a, b)) in enumerate([(a, b) for a in 1:nspecies for b in (a + 1):nspecies])
        pairs[p, :] .= (a, b)
        products[p, :, :] .= windings[a, :, :] .* windings[b, :, :]
    end
    total = dropdims(sum(windings; dims=1); dims=1)
    histories = Dict{String,Any}(
        "windingBySpecies" => zeros(size(windings)),
        "squaredWindingBySpecies" => windings .^ 2,
        "totalWinding" => zeros(size(total)),
        "squaredTotalWinding" => total .^ 2,
    )
    nspecies > 1 && (histories["windingCorrelationBySpecies"] = products)
    return histories, Dict{String,Any}("species_pairs" => pairs,
                                      "periodic_axes" => fill(true, naxes))
end

@testset "Winding pair channels" begin
    fields = ("squaredCoflowBySpeciesPair", "squaredCounterflowBySpeciesPair")
    counts = [1, 3, 2, 4]
    @testset "$nspecies species" for nspecies in 1:4
        windings = [Float64((-1)^s * s * (b - 2)) for s in 1:nspecies, d in 1:2, b in 1:4]
        windings[:, 2, :] .= 0
        histories, metadata = winding_channel_fixture(windings)
        metadata["periodic_axes"][2] = false
        before = deepcopy((histories, metadata))
        plain = WormWeaver.windings_summary(histories, counts, metadata; min_blocks=2)
        enabled = WormWeaver.windings_summary(histories, counts, metadata;
                                             min_blocks=2, blocking_diagnostics=true)
        retained = copy(enabled)
        diagnostics = pop!(retained, "blocking_diagnostics")
        @test isequal(plain, retained)
        @test (histories, metadata) == before
        @test enabled["species_pairs"] == metadata["species_pairs"]
        @test enabled["nsamples"] == 10
        @test diagnostics["sample_counts"] == [10, 10]
        for field in fields
            @test size(enabled[field]) == (nspecies * (nspecies - 1) ÷ 2, 2)
            @test eltype(enabled[field]) == Measurements.Measurement{Float64}
            @test size(diagnostics[field * "StandardError"]) == (nspecies * (nspecies - 1) ÷ 2, 2, 2)
            @test all(iszero, Measurements.value.(enabled[field][:, 2]))
            @test all(iszero, Measurements.uncertainty.(enabled[field][:, 2]))
        end
        for (p, row) in enumerate(eachrow(metadata["species_pairs"]))
            a, b = row
            for (field, sign) in zip(fields, (1, -1)), d in 1:2
                expected_history = (windings[a, d, :] .+ sign .* windings[b, d, :]) .^ 2
                expected = WormWeaver.blocking_analysis(expected_history, counts;
                                                       min_blocks=2, diagnostics=true)
                @test isequal(enabled[field][p, d], expected.measurement)
                @test diagnostics[field * "StandardError"][p, d, :] ≈ expected.standard_errors
            end
        end

        # Exercise the public HDF5 reader, JLD2 persistence and table expansion.
        mktempdir() do root
            path = joinpath(root, "windings.h5")
            HDF5.h5open(path, "w") do file
                HDF5.create_group(file, "observables")
                HDF5.create_group(file, "observables/METADATA")
                write_observable_fixture(file, "Windings", histories, counts; metadata)
            end
            manifest = DataFrames.DataFrame(output_file=["windings.h5"])
            manifest[!, "run.blocks"] = [4]
            manifest_path = joinpath(root, "manifest.jld2")
            JLD2.jldsave(manifest_path; manifest)
            before_bytes = read.([path, manifest_path])
            table = aggregate(root; observables=["Windings"], min_blocks=2,
                              blocking_diagnostics=true)
            @test isequal(table.Windings[1], enabled)
            @test isequal(JLD2.load(joinpath(root, "aggregated.jld2"), "aggregate"), table)
            expanded = observable_table(table.Windings)
            for field in fields
                @test expanded[1, field] === table.Windings[1][field]
            end
            @test read.([path, manifest_path]) == before_bytes
        end
    end

    @testset "Joint fluctuations and cancellation" begin
        amplitude = vcat(zeros(32), ones(32))
        for sign in (1, -1)
            windings = zeros(2, 1, 64)
            windings[1, 1, :] .= amplitude
            windings[2, 1, :] .= sign .* amplitude
            histories, metadata = winding_channel_fixture(windings)
            summary = WormWeaver.windings_summary(histories, ones(Int, 64), metadata;
                                                 min_blocks=4, blocking_diagnostics=true)
            zero_field, nonzero_field = sign == 1 ? reverse(fields) : fields
            @test Measurements.value(summary[zero_field][1, 1]) == 0
            @test Measurements.uncertainty(summary[zero_field][1, 1]) == 0
            @test Measurements.value(summary[nonzero_field][1, 1]) == 2
            @test Measurements.uncertainty(summary[nonzero_field][1, 1]) ≈ 4 / sqrt(12)
            @test Measurements.uncertainty(summary["squaredWindingBySpecies"][1, 1]) > 0
            @test all(iszero, summary["blocking_diagnostics"][zero_field * "StandardError"])
            @test Measurements.value(summary["windingBySpecies"][1, 1]) == 0
            @test isequal(summary[fields[1]][1, 1], summary["squaredTotalWinding"][1])
        end
        histories, metadata = winding_channel_fixture(ones(2, 1, 1))
        one_block = WormWeaver.windings_summary(histories, [10], metadata;
                                               min_blocks=2, blocking_diagnostics=true)
        for field in fields
            @test isnan(Measurements.uncertainty(one_block[field][1, 1]))
            @test all(isnan, one_block["blocking_diagnostics"][field * "StandardError"])
        end
    end

    @testset "Pair ordering and partial coverage" begin
        windings = reshape([1.0, 2, -3], 3, 1, 1) .* reshape([1.0, 2, 3, 4], 1, 1, 4)
        histories, metadata = winding_channel_fixture(windings)
        baseline = WormWeaver.windings_summary(histories, counts, metadata; min_blocks=2)
        order = [3, 1, 2]
        metadata["species_pairs"] = metadata["species_pairs"][order, :]
        histories["windingCorrelationBySpecies"] = histories["windingCorrelationBySpecies"][order, :, :]
        reordered = WormWeaver.windings_summary(histories, counts, metadata; min_blocks=2)
        @test reordered["species_pairs"] == [2 3; 1 2; 1 3]
        for field in fields
            @test isequal(reordered[field], baseline[field][order, :])
        end
        metadata["species_pairs"] = metadata["species_pairs"][1:1, :]
        histories["windingCorrelationBySpecies"] = histories["windingCorrelationBySpecies"][1:1, :, :]
        partial = WormWeaver.windings_summary(histories, counts, metadata; min_blocks=2)
        for field in fields
            @test isequal(partial[field], reordered[field][1:1, :])
        end
        metadata["species_pairs"] = zeros(Int, 0, 2)
        histories["windingCorrelationBySpecies"] = zeros(0, 1, 4)
        empty_pairs = WormWeaver.windings_summary(histories, counts, metadata;
                                                 min_blocks=2, blocking_diagnostics=true)
        @test size(empty_pairs[fields[1]]) == (0, 1)
    end

    @testset "Legacy moments and missing metadata" begin
        histories, metadata = winding_channel_fixture(reshape([1.0, -2, 2, -3, 3, -4, 4, -5], 2, 1, 4))
        baseline = WormWeaver.windings_summary(histories, counts, metadata;
                                              min_blocks=2, blocking_diagnostics=true)
        without_pairs = copy(metadata)
        delete!(without_pairs, "species_pairs")
        unmapped = @test_logs (:warn, r"no species_pairs metadata") WormWeaver.windings_summary(
            histories, counts, without_pairs; min_blocks=2, blocking_diagnostics=true)
        @test !haskey(unmapped, fields[1]) && !haskey(unmapped, fields[2])
        @test haskey(unmapped, "windingCorrelationBySpecies")
        delete!(histories, "windingCorrelationBySpecies")
        for md in (metadata, without_pairs)
            legacy = WormWeaver.windings_summary(histories, counts, md;
                                                 min_blocks=2, blocking_diagnostics=true)
            @test legacy["species_pairs"] == [1 2]
            @test !haskey(legacy, "windingCorrelationBySpecies")
            for field in fields
                @test isequal(legacy[field], baseline[field])
                @test legacy["blocking_diagnostics"][field * "StandardError"] ≈
                      baseline["blocking_diagnostics"][field * "StandardError"]
            end
        end
        for nspecies in (3, 4)
            histories, metadata = winding_channel_fixture(ones(nspecies, 1, 4))
            delete!(histories, "windingCorrelationBySpecies")
            missing_products = @test_logs (:warn, r"require pair products") WormWeaver.windings_summary(
                histories, counts, metadata; min_blocks=2, blocking_diagnostics=true)
            @test !haskey(missing_products, fields[1]) && !haskey(missing_products, fields[2])
            @test !haskey(missing_products["blocking_diagnostics"], fields[1] * "StandardError")
            @test Measurements.value.(missing_products["squaredTotalWinding"]) == [nspecies^2]
        end
        histories, metadata = winding_channel_fixture(ones(1, 1, 4))
        delete!(metadata, "species_pairs")
        single = WormWeaver.windings_summary(histories, counts, metadata; min_blocks=2)
        @test size(single["species_pairs"]) == (0, 2)
        @test size(single[fields[1]]) == (0, 1)
    end

    @testset "Invalid dimensions and pair metadata" begin
        histories, metadata = winding_channel_fixture(ones(3, 2, 4))
        for (field, replacement) in (
            ("squaredWindingBySpecies", zeros(3, 8)),
            ("squaredWindingBySpecies", zeros(0, 2, 4)),
            ("squaredWindingBySpecies", zeros(3, 2, 3)),
            ("windingBySpecies", zeros(2, 2, 4)),
            ("totalWinding", zeros(1, 4)),
            ("squaredTotalWinding", zeros(2, 1, 4)),
            ("windingCorrelationBySpecies", zeros(3, 2, 1, 4)),
            ("windingCorrelationBySpecies", zeros(3, 1, 4)),
            ("windingCorrelationBySpecies", zeros(4, 2, 4)),
            ("windingCorrelationBySpecies", zeros(2, 2, 4)),
        )
            bad = copy(histories)
            bad[field] = replacement
            @test_throws DimensionMismatch WormWeaver.windings_summary(bad, counts, metadata; min_blocks=2)
        end
        for (pairs, error) in (
            ([1, 2], ArgumentError),
            (ones(Int, 3, 3), DimensionMismatch),
            ([1.0 2; 1 3; 2 3], ArgumentError),
            (trues(3, 2), ArgumentError),
            ([0 2; 1 3; 2 3], ArgumentError),
            ([1 4; 1 3; 2 3], ArgumentError),
            ([1 1; 1 3; 2 3], ArgumentError),
            ([2 1; 1 3; 2 3], ArgumentError),
            ([1 2; 1 2; 2 3], ArgumentError),
            (zeros(Int, 0, 2), DimensionMismatch),
        )
            bad = copy(metadata)
            bad["species_pairs"] = pairs
            @test_throws error WormWeaver.windings_summary(histories, counts, bad; min_blocks=2)
        end
        bad = copy(metadata)
        bad["periodic_axes"] = [true]
        @test_throws DimensionMismatch WormWeaver.windings_summary(histories, counts, bad; min_blocks=2)
        delete!(histories, "windingCorrelationBySpecies")
        bad = copy(metadata)
        bad["species_pairs"] = [1 2; 1 2]
        @test_throws ArgumentError WormWeaver.windings_summary(histories, counts, bad; min_blocks=2)
    end
end
