function write_single_particle_green_fixture(path; nspecies=2, counts=[11, 23, 5, 41],
        hits=reshape(collect(1.0:(nspecies * 3 * length(counts))), nspecies, 3, length(counts)),
        zcounts=Float64.(1:length(counts)), species=collect(1:nspecies),
        metadata=(;), reduction="sum", append=false, gktau=nothing)
    md = merge((;
        species, displacements=reshape([0, 1, -1], 3, 1),
        displacement_vectors=reshape([0.0, 1, -1], 3, 1),
        displacement_counts=[3, 3, 3], beta=2.0, equal_time_window=0.05,
        equal_time_convention="one-sided inclusive window",
        production_log_eta=collect(1.0:nspecies),
    ), gktau === nothing ? (;) : (;
        momenta=reshape(collect(0.0:(size(gktau, 2) - 1)), :, 1),
        momentum_indices=reshape(collect(0:(size(gktau, 2) - 1)), :, 1),
        tau_edges=collect(range(0.0, 2.0; length=size(gktau, 3) + 1)),
    ), metadata)
    HDF5.h5open(path, append ? "r+" : "w") do file
        if !append
            HDF5.create_group(file, "observables")
            HDF5.create_group(file, "observables/METADATA")
        end
        write_observable_fixture(file, "SingleParticleGreenFunction",
            (; gREqualTimeHits=hits, zSectorVisitCount=zcounts), counts; metadata=md, reduction)
        if gktau !== nothing
            group = file["observables/SingleParticleGreenFunction"]
            group["gKTau"] = gktau
            HDF5.attributes(group["gKTau"])["reduction"] = "mean"
            file["observables/METADATA/SingleParticleGreenFunction/reduction_gKTau"] = "mean"
        end
        # Deliberately unusable derived data must not be consulted by raw aggregation.
        derived = HDF5.create_group(file, "derived")
        derived["SingleParticleGreenFunction"] = "do not read this"
    end
    return (; hits, zcounts, counts, metadata=md)
end

@testset "Raw complex G(k,tau) aggregation" begin
    mktempdir() do root
        path = joinpath(root, "results.h5")
        manifest = DataFrames.DataFrame(output_file=["results.h5"])
        manifest[!, "run.blocks"] = [5]
        write_aggregation_manifest(joinpath(root, "manifest.jld2"); manifest)
        counts, zcounts = [1, 3, 2, 4], [1.0, 2, 1, 1]
        offsets = reshape(ComplexF64.(1:8), 2, 2, 2)
        history = offsets .+ reshape([0-1im, 1+0im, 3+2im, 5-3im], 1, 1, 1, 4)
        kwargs = (; counts, zcounts, species=[4, 2],
                  metadata=(tau_edges=[0.0, 0.5, 2.0],))
        write_single_particle_green_fixture(path; kwargs...)
        baseline = aggregate(root; observables=["SingleParticleGreenFunction"], min_blocks=2).SingleParticleGreenFunction[1]
        write_single_particle_green_fixture(path; kwargs..., gktau=history)
        source_bytes = read(path)
        table = save_aggregate(root; observables=["SingleParticleGreenFunction"],
                               min_blocks=2, blocking_diagnostics=true)
        summary = table.SingleParticleGreenFunction[1]
        result = summary["gKTau"]
        @test result isa Array{Complex{Measurements.Measurement{Float64}},3}
        @test size(result) == (2, 2, 2)
        @test Measurements.value.(real.(result)) ≈ real.(offsets) .+ 2.9
        @test Measurements.value.(imag.(result)) ≈ fill(-0.9, 2, 2, 2)
        # The maximum error occurs at different blocking levels for each component.
        @test Measurements.uncertainty.(real.(result)) ≈ fill(1.72, 2, 2, 2)
        @test Measurements.uncertainty.(imag.(result)) ≈ fill(sqrt(446/3)/10, 2, 2, 2)
        @test summary["sampleCount"] == 10
        @test summary["zSectorVisitCount"] == 5
        @test summary["nblocks"] == 4
        @test summary["species"] == [4, 2]
        @test summary["tau_edges"] == [0.0, 0.5, 2.0]
        @test isequal(summary["grEqualTime"], baseline["grEqualTime"])
        @test summary["gKTauAggregationMethod"] == "sample-count-weighted mean"
        @test summary["gKTauErrorMethod"] == "repeated blocking: maximum standard error"
        diagnostics = summary["blocking_diagnostics"]
        @test size(diagnostics["gKTauRealStandardError"]) == (2, 2, 2, 2)
        @test size(diagnostics["gKTauImaginaryStandardError"]) == (2, 2, 2, 2)
        @test haskey(diagnostics, "grEqualTimeStandardError")
        @test diagnostics["block_counts"] == [4, 2]
        @test diagnostics["sample_counts"] == [10, 10]
        @test isequal(JLD2.load(default_aggregate_path(root), "aggregate"), table)
        @test observable_table(table.SingleParticleGreenFunction).gKTau[1] === result
        @test read(path) == source_bytes
        plain = aggregate(root; observables=["SingleParticleGreenFunction"], min_blocks=2).SingleParticleGreenFunction[1]
        @test isequal(plain["gKTau"], result)
        @test !haskey(plain, "blocking_diagnostics")

        # Raw means need neither equal-time metadata nor positive total Z visits.
        HDF5.h5open(path, "r+") do file
            HDF5.delete_object(file, "observables/SingleParticleGreenFunction/gREqualTimeHits")
            file["observables/SingleParticleGreenFunction/zSectorVisitCount"][:] = zeros(4)
            for name in ("equal_time_window", "equal_time_convention", "production_log_eta")
                HDF5.delete_object(file, "observables/METADATA/SingleParticleGreenFunction/$name")
            end
        end
        raw_only = aggregate(root; observables=["SingleParticleGreenFunction"], min_blocks=2).SingleParticleGreenFunction[1]
        @test isequal(raw_only["gKTau"], result)
        @test raw_only["zSectorVisitCount"] == 0
        @test !haskey(raw_only, "grEqualTime")

        for equal_time in (false, true), species in ([1, 1], [-1, 2], [true, false], [1.0, 2.0])
            write_single_particle_green_fixture(path; merge(kwargs, (; species))..., gktau=history)
            if !equal_time
                HDF5.h5open(path, "r+") do file
                    HDF5.delete_object(file, "observables/SingleParticleGreenFunction/gREqualTimeHits")
                end
            end
            @test_throws ArgumentError aggregate(root; observables=["SingleParticleGreenFunction"])
        end

        odd = cat(history, offsets .+ (10+4im); dims=4)
        write_single_particle_green_fixture(path; counts=[counts; 5], gktau=odd)
        summary = aggregate(root; observables=["SingleParticleGreenFunction"], min_blocks=2).SingleParticleGreenFunction[1]
        @test Measurements.value.(real.(summary["gKTau"])) ≈ real.(offsets) .+ 79/15
        @test Measurements.value.(imag.(summary["gKTau"])) ≈ fill(11/15, 2, 2, 2)

        write_single_particle_green_fixture(path; counts=[9], gktau=reshape([1+2im], 1, 1, 1, 1), nspecies=1)
        single = aggregate(root; observables=["SingleParticleGreenFunction"]).SingleParticleGreenFunction[1]["gKTau"][1]
        @test Measurements.value(real(single)) == 1
        @test Measurements.value(imag(single)) == 2
        @test isnan(Measurements.uncertainty(real(single)))
        @test isnan(Measurements.uncertainty(imag(single)))

        for (data, expected_error) in ((history[:, :, :, 1:3], DimensionMismatch),
                                      (fill(ComplexF64(NaN), 2, 2, 2, 4), ArgumentError),
                                      (ones(2, 2, 4), DimensionMismatch))
            write_single_particle_green_fixture(path; kwargs..., gktau=data)
            @test_throws expected_error aggregate(root; observables=["SingleParticleGreenFunction"])
        end
        for name in ("species", "momenta", "tau_edges")
            write_single_particle_green_fixture(path; kwargs..., gktau=history)
            HDF5.h5open(path, "r+") do file
                HDF5.delete_object(file, "observables/METADATA/SingleParticleGreenFunction/$name")
            end
            @test_throws ArgumentError aggregate(root; observables=["SingleParticleGreenFunction"])
        end
        write_single_particle_green_fixture(path; kwargs..., gktau=history)
        HDF5.h5open(path, "r+") do file
            HDF5.attrs(file["observables/SingleParticleGreenFunction/gKTau"])["reduction"] = "sum"
        end
        @test_throws ArgumentError aggregate(root; observables=["SingleParticleGreenFunction"])
    end
end

@testset "Equal-time SPGF analysis and optional representations" begin
    mktempdir() do root
        path = joinpath(root, "results.h5")
        manifest = DataFrames.DataFrame(output_file=["results.h5"])
        manifest[!, "run.blocks"] = [3]
        write_aggregation_manifest(joinpath(root, "manifest.jld2"); manifest)
        write_single_particle_green_fixture(path; nspecies=1, counts=[7, 11, 9],
            hits=reshape([2.0, 8, 10], 1, 1, 3), zcounts=[1.0, 2, 1],
            metadata=(displacements=zeros(Int, 1, 1), displacement_vectors=zeros(1, 1),
                displacement_counts=[1], equal_time_window=0.5, production_log_eta=[log(2.0)]))
        sources = read.([path, joinpath(root, "manifest.jld2")])
        table = save_aggregate(root; observables=["SingleParticleGreenFunction"])
        summary = table.SingleParticleGreenFunction[1]
        @test summary["grEqualTime"] isa Matrix{Measurements.Measurement{Float64}}
        @test Measurements.value(summary["grEqualTime"][1]) == 10
        @test Measurements.uncertainty(summary["grEqualTime"][1]) ≈ 32/9
        diagnosed = aggregate(root; observables=["SingleParticleGreenFunction"], blocking_diagnostics=true)
        @test isequal(diagnosed.SingleParticleGreenFunction[1]["grEqualTime"], summary["grEqualTime"])
        @test diagnosed.SingleParticleGreenFunction[1]["blocking_diagnostics"]["grEqualTimeStandardError"] ≈ reshape([32/9], 1, 1, 1)
        @test read.([path, joinpath(root, "manifest.jld2")]) == sources
        @test isequal(JLD2.load(default_aggregate_path(root), "aggregate"), table)
        @test observable_table(table.SingleParticleGreenFunction).grEqualTime[1] === summary["grEqualTime"]

        # Additional numerator representations are optional and do not affect
        # the equal-time result. Their metadata is retained for future adapters.
        HDF5.h5open(path, "r+") do file
            obs = file["observables/SingleParticleGreenFunction"]
            obs["gKTauRealNumerator"] = ones(1, 2, 4, 3)
            obs["gKTauImaginaryNumerator"] = zeros(1, 2, 4, 3)
            for name in ("gKTauRealNumerator", "gKTauImaginaryNumerator")
                HDF5.attrs(obs[name])["reduction"] = "mean"
            end
            md = file["observables/METADATA/SingleParticleGreenFunction"]
            md["tau_edges"] = collect(0.0:0.5:2.0)
            md["momenta"] = reshape([0.0, pi], 2, 1)
        end
        mixed = aggregate(root; observables=["SingleParticleGreenFunction"]).SingleParticleGreenFunction[1]
        @test isequal(mixed["grEqualTime"], summary["grEqualTime"])
        @test mixed["tau_edges"] == collect(0.0:0.5:2.0)
        HDF5.h5open(path, "r+") do file
            HDF5.delete_object(file, "observables/SingleParticleGreenFunction/gREqualTimeHits")
            for name in ("equal_time_window", "equal_time_convention", "production_log_eta")
                HDF5.delete_object(file, "observables/METADATA/SingleParticleGreenFunction/$name")
            end
        end
        time_only = aggregate(root; observables=["SingleParticleGreenFunction"], blocking_diagnostics=true).SingleParticleGreenFunction[1]
        @test time_only["zSectorVisitCount"] == 4
        @test time_only["sampleCount"] == 27
        @test time_only["momenta"] == reshape([0.0, pi], 2, 1)
        @test !haskey(time_only, "grEqualTime")
        @test !haskey(time_only, "gREqualTimeHits")
        @test !haskey(time_only, "blocking_diagnostics")

        # Dimensions and primitive-vector geometry do not alter normalization.
        for dimension in 1:3
            displacements = zeros(Int, 3, dimension)
            displacements[2, 1], displacements[3, 1] = 1, -1
            vectors = Float64.(displacements)
            if dimension >= 2
                vectors[2, 2], vectors[3, 2] = 0.5, -0.5
            end
            write_single_particle_green_fixture(path; nspecies=1, counts=[1, 2, 3],
                hits=ones(1, 3, 3), zcounts=ones(3),
                metadata=(displacements=displacements, displacement_vectors=vectors, production_log_eta=[0.0]))
            result = aggregate(root; observables=["SingleParticleGreenFunction"]).SingleParticleGreenFunction[1]
            @test Measurements.value.(result["grEqualTime"]) ≈ fill(40.0, 1, 3)
            @test all(iszero, Measurements.uncertainty.(result["grEqualTime"]))
            @test size(result["displacement_vectors"]) == (3, dimension)
            @test result["displacement_vectors"] == vectors
        end
        write_single_particle_green_fixture(path; nspecies=1, counts=[1, 2, 3],
            hits=ones(1, 3, 3), zcounts=ones(3),
            metadata=(displacement_counts=[2, 1, 1], equal_time_window=0.5, production_log_eta=[0.0]))
        unequal = aggregate(root; observables=["SingleParticleGreenFunction"]).SingleParticleGreenFunction[1]
        @test unequal["nsites"] == 2
        @test Measurements.value.(unequal["grEqualTime"]) ≈ [4.0 8.0 8.0]

        manifest[!, "run.blocks"] = [1]
        write_aggregation_manifest(joinpath(root, "manifest.jld2"); manifest)
        write_single_particle_green_fixture(path; nspecies=1, counts=[2], zcounts=[1])
        one = aggregate(root; observables=["SingleParticleGreenFunction"]).SingleParticleGreenFunction[1]
        @test all(isnan, Measurements.uncertainty.(one["grEqualTime"]))
        @test one["error_status"] == "fewer than two blocks"
    end
end

@testset "Raw single-particle Green-function aggregation" begin
    mktempdir() do root
        path = joinpath(root, "results.h5")
        manifest = DataFrames.DataFrame(output_file=["results.h5"])
        manifest[!, "run.blocks"] = [4]
        manifest_path = joinpath(root, "manifest.jld2")
        write_aggregation_manifest(manifest_path; manifest)
        fixture = write_single_particle_green_fixture(path; species=[4, 2])
        source_bytes = read.([path, manifest_path])
        table = save_aggregate(root; observables=["SingleParticleGreenFunction"])
        summary = table.SingleParticleGreenFunction[1]
        @test names(table) == [names(manifest); ["status", "completed_blocks", "failed", "log_source", "g2_entry_modes"]; ["SingleParticleGreenFunction"]]
        @test summary["gREqualTimeHits"] == dropdims(sum(fixture.hits; dims=3); dims=3)
        @test size(summary["gREqualTimeHits"]) == (2, 3)
        @test eltype(summary["gREqualTimeHits"]) == Float64
        @test summary["zSectorVisitCount"] === 10.0
        @test summary["sampleCount"] == summary["nsamples"] == 80
        @test summary["sampleCount"] isa Int64
        @test summary["nblocks"] == 4
        @test summary["aggregation_method"] == "sum"
        for (field, value) in pairs(fixture.metadata)
            @test summary[string(field)] == value
        end
        @test !haskey(summary, "gr_equal_time")
        @test summary["error_method"] == "repeated blocking: maximum paired jackknife standard error"
        @test summary["min_blocks"] == 10
        @test summary["nsites"] == 3
        factors = exp.(-fixture.metadata.production_log_eta) .* (2.0 / 0.05)
        @test Measurements.value.(summary["grEqualTime"]) ≈ factors .* summary["gREqualTimeHits"] ./ 10
        @test !haskey(summary, "blocking_diagnostics")
        enabled = aggregate(root; observables=["SingleParticleGreenFunction"],
            blocking_diagnostics=true, min_blocks=2)
        @test Measurements.value.(enabled.SingleParticleGreenFunction[1]["grEqualTime"]) == Measurements.value.(summary["grEqualTime"])
        @test haskey(enabled.SingleParticleGreenFunction[1], "blocking_diagnostics")
        saved = JLD2.load(default_aggregate_path(root), "aggregate")
        @test isequal(saved, table)
        @test Set(keys(JLD2.load(default_aggregate_path(root)))) == Set(["aggregate"])
        expanded = observable_table(table.SingleParticleGreenFunction)
        @test expanded.gREqualTimeHits[1] === summary["gREqualTimeHits"]
        @test expanded.zSectorVisitCount == [10.0]
        @test expanded.sampleCount == [80]
        @test expanded.species[1] == [4, 2]
        @test read.([path, manifest_path]) == source_bytes

        # Integer-valued storage remains supported, including checked integer sums.
        fixture = write_single_particle_green_fixture(path;
            hits=reshape(collect(1:24), 2, 3, 4), zcounts=[3, 0, 1, 2])
        integer_summary = aggregate(root; observables=["SingleParticleGreenFunction"]).SingleParticleGreenFunction[1]
        @test integer_summary["gREqualTimeHits"] == dropdims(sum(fixture.hits; dims=3); dims=3)
        @test integer_summary["zSectorVisitCount"] == 6

        for nspecies in 1:3, nblocks in (1, 4)
            counts = collect(1:nblocks)
            manifest[!, "run.blocks"] = [nblocks]
            write_aggregation_manifest(manifest_path; manifest)
            fixture = write_single_particle_green_fixture(path; nspecies, counts,
                species=reverse(collect(1:nspecies)), zcounts=ones(nblocks))
            result = aggregate(root; observables=["SingleParticleGreenFunction"]).SingleParticleGreenFunction[1]
            @test size(result["gREqualTimeHits"]) == (nspecies, 3)
            @test result["gREqualTimeHits"] == dropdims(sum(fixture.hits; dims=3); dims=3)
            @test result["species"] == reverse(collect(1:nspecies))
            @test result["zSectorVisitCount"] == nblocks
            @test result["sampleCount"] == sum(counts)
            @test result["nblocks"] == nblocks
        end
        write_single_particle_green_fixture(path; hits=zeros(2, 3, 4), zcounts=zeros(4))
        manifest[!, "run.blocks"] = [4]
        write_aggregation_manifest(manifest_path; manifest)
        @test_throws ArgumentError aggregate(root; observables=["SingleParticleGreenFunction"])

        write_single_particle_green_fixture(path; counts=[Int64(2)^53 + 1, 1, 1, 1])
        large = aggregate(root; observables=["SingleParticleGreenFunction"]).SingleParticleGreenFunction[1]
        @test large["sampleCount"] == Int64(2)^53 + 4

        write_single_particle_green_fixture(path; counts=[11, 23])
        partial = aggregate(
            root; observables=["SingleParticleGreenFunction"])
        @test partial.SingleParticleGreenFunction[1]["nblocks"] == 2
        @test partial.SingleParticleGreenFunction[1]["sampleCount"] == 34
        write_combined_fixture(path)
        absent = @test_logs (:warn, r"SingleParticleGreenFunction was not found") aggregate(
            root; observables=["SingleParticleGreenFunction"])
        @test ismissing(absent.SingleParticleGreenFunction[1])
        @test ismissing(aggregate(root).SingleParticleGreenFunction[1])

        manifest.output_file = ["missing.h5"]
        write_aggregation_manifest(manifest_path; manifest)
        missing_file = @test_logs (:warn, r"Results file is missing") aggregate(
            root; observables=["SingleParticleGreenFunction"])
        @test ismissing(missing_file.SingleParticleGreenFunction[1])
    end
end

@testset "Raw Green-function validation" begin
    mktempdir() do root
        path = joinpath(root, "results.h5")
        manifest = DataFrames.DataFrame(output_file=["results.h5"])
        manifest[!, "run.blocks"] = [4]
        write_aggregation_manifest(joinpath(root, "manifest.jld2"); manifest)
        read_summary() = aggregate(root; observables=["SingleParticleGreenFunction"])
        cases = (
            ((; hits=ones(2, 3)), DimensionMismatch),
            ((; hits=ones(2, 3, 3)), DimensionMismatch),
            ((; hits=zeros(0, 3, 4)), DimensionMismatch),
            ((; zcounts=ones(1, 4)), DimensionMismatch),
            ((; zcounts=[0.0, -1, 0, 0]), ArgumentError),
            ((; hits=fill(NaN, 2, 3, 4)), ArgumentError),
            ((; hits=fill(-1.0, 2, 3, 4)), ArgumentError),
            ((; hits=fill(false, 2, 3, 4)), ArgumentError),
            ((; counts=[1, 0, 2, 3]), ArgumentError),
            ((; counts=[1.0, 1.5, 2, 3]), ArgumentError),
            ((; species=[1]), DimensionMismatch),
            ((; species=[1, 1]), ArgumentError),
            ((; metadata=(; displacements=zeros(Int, 2, 1))), DimensionMismatch),
            ((; metadata=(; displacement_vectors=zeros(3, 2))), DimensionMismatch),
            ((; metadata=(; displacement_counts=[8, 8])), DimensionMismatch),
            ((; metadata=(; displacement_counts=[8, 0, 8])), ArgumentError),
            ((; metadata=(; production_log_eta=[0.0])), DimensionMismatch),
            ((; metadata=(; displacement_counts=[8, 8, 8])), ArgumentError),
            ((; metadata=(; beta=Inf)), ArgumentError),
            ((; metadata=(; equal_time_window=2.0)), ArgumentError),
            ((; metadata=(; equal_time_window=0.0)), ArgumentError),
            ((; metadata=(; equal_time_convention="two-sided")), ArgumentError),
            ((; reduction="mean"), ArgumentError),
            ((; reduction=nothing), ArgumentError),
        )
        for (kwargs, exception) in cases
            write_single_particle_green_fixture(path; kwargs...)
            @test_throws exception read_summary()
        end
        write_single_particle_green_fixture(path)
        HDF5.h5open(path, "r+") do file
            HDF5.delete_object(file, "observables/METADATA/SingleParticleGreenFunction/species")
        end
        @test_throws ArgumentError read_summary()
        write_single_particle_green_fixture(path)
        HDF5.h5open(path, "r+") do file
            HDF5.attrs(file["observables/SingleParticleGreenFunction/zSectorVisitCount"])["reduction"] = "mean"
        end
        @test_throws ArgumentError read_summary()
        write_single_particle_green_fixture(path; counts=[1, 2, 3, 4, 5])
        @test_throws ArgumentError read_summary()
    end
end
