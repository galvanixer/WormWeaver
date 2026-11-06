@testset "Raw counterflow equal-time aggregation" begin
    mktempdir() do root
        path = joinpath(root, "results.h5")
        manifest = DataFrames.DataFrame(output_file=["results.h5", "results.h5"], seed=[11, 22])
        manifest[!, "run.blocks"] = [4, 4]
        write_aggregation_manifest(joinpath(root, "manifest.jld2"); manifest)
        counts, zcounts = [10, 100, 20, 200], [2.0, 2, 2, 2]
        hits = [1.0 10 0 20; 0 0 0 0; 2 3 4 5]
        metadata = (species_pair=[2, 1], channel="counterflow",
            displacements=reshape([0, 1, -1], 3, 1),
            displacement_vectors=reshape([0.0, 0.5, -0.5], 3, 1), displacement_counts=[4, 2, 2],
            beta=10.0, equal_time_window=0.05, equal_time_convention="one-sided inclusive window",
            endpoint_convention="first species sets direction", production_log_eta=-2.0)
        function write_fixture(; n=4, hit_history=hits[:, 1:n], z=zcounts[1:n], samples=counts[1:n],
                               md=metadata, reduction="sum", extras=false)
            HDF5.h5open(path, "w") do file
                HDF5.create_group(file, "observables")
                HDF5.create_group(file, "observables/METADATA")
                fields = hit_history === nothing ? (zSectorVisitCount=z,) :
                    (zSectorVisitCount=z, gREqualTimeHits=hit_history)
                write_observable_fixture(file, "CounterflowGreenFunction", fields, samples;
                                         metadata=md, reduction)
                if extras
                    for field in ("gRTau", "gKTau", "gKOmega", "gROmega")
                        file["observables/CounterflowGreenFunction/$field"] = zeros(3, 2, n)
                    end
                end
                derived = HDF5.create_group(file, "derived")
                derived["CounterflowGreenFunction"] = "must not be read"
            end
        end
        write_fixture()
        original = read(path)
        table = aggregate(root; observables=["CounterflowGreenFunction"])
        summary = table.CounterflowGreenFunction[1]
        @test table.seed == [11, 22]
        @test isequal(table[:, names(manifest)], manifest)
        @test summary["gREqualTimeHits"] == [31, 0, 14]
        @test summary["zSectorVisitCount"] == 8
        @test summary["sampleCount"] == summary["nsamples"] == 330
        @test summary["nblocks"] == 4
        @test summary["aggregation_method"] == "sum"
        for (key, value) in pairs(metadata)
            @test summary[string(key)] == value
        end
        for field in ("grEqualTime", "error_method", "error_status", "min_blocks", "blocking_diagnostics")
            @test !haskey(summary, field)
        end
        enabled = save_aggregate(root; observables=["CounterflowGreenFunction"], blocking_diagnostics=true)
        @test isequal(enabled, table)
        @test isequal(JLD2.load(default_aggregate_path(root), "aggregate"), table)
        @test Set(keys(JLD2.load(default_aggregate_path(root)))) == Set(["aggregate"])
        @test read(path) == original
        @test observable_table(table.CounterflowGreenFunction).species_pair[1] == [2, 1]

        write_fixture(; n=1)
        one = aggregate(root; observables=["CounterflowGreenFunction"]).CounterflowGreenFunction[1]
        @test one["gREqualTimeHits"] == hits[:, 1]
        @test one["nblocks"] == 1 && one["zSectorVisitCount"] == 2
        write_fixture(; n=2)
        partial = aggregate(root; observables=["CounterflowGreenFunction"]).CounterflowGreenFunction[1]
        @test partial["gREqualTimeHits"] == [11, 0, 5]
        @test partial["sampleCount"] == 110 && partial["nblocks"] == 2
        write_fixture(; hit_history=zeros(3, 4))
        @test all(iszero, aggregate(root; observables=["CounterflowGreenFunction"]).CounterflowGreenFunction[1]["gREqualTimeHits"])
        write_fixture(; hit_history=Int.(hits), z=Int.(zcounts))
        integer = aggregate(root; observables=["CounterflowGreenFunction"]).CounterflowGreenFunction[1]
        @test integer["gREqualTimeHits"] == [31, 0, 14] && integer["zSectorVisitCount"] == 8

        write_fixture(; hit_history=nothing, md=(species_pair=[2, 1], channel="counterflow"))
        only_counts = aggregate(root; observables=["CounterflowGreenFunction"]).CounterflowGreenFunction[1]
        @test !haskey(only_counts, "gREqualTimeHits")
        @test only_counts["zSectorVisitCount"] == 8
        no_weight = Base.structdiff(metadata, (; production_log_eta=metadata.production_log_eta))
        write_fixture(; md=no_weight)
        @test !haskey(aggregate(root; observables=["CounterflowGreenFunction"]).CounterflowGreenFunction[1], "production_log_eta")
        write_fixture(; extras=true)
        omitted = @test_logs (:warn, r"representations omitted") (:warn, r"representations omitted") aggregate(root; observables=["CounterflowGreenFunction"])
        @test omitted.CounterflowGreenFunction[1]["gREqualTimeHits"] == [31, 0, 14]
        @test !haskey(omitted.CounterflowGreenFunction[1], "gRTau")

        invalid_cases = (
            (ArgumentError, (; z=[2.0, 0, 2, 2])),
            (ArgumentError, (; z=[2.0, -1, 2, 2])),
            (ArgumentError, (; z=[2.0, NaN, 2, 2])),
            (DimensionMismatch, (; z=ones(1, 4))),
            (DimensionMismatch, (; z=ones(3))),
            (ArgumentError, (; hit_history=fill(-1.0, 3, 4))),
            (ArgumentError, (; hit_history=fill(Inf, 3, 4))),
            (DimensionMismatch, (; hit_history=zeros(1, 3, 4))),
            (DimensionMismatch, (; hit_history=zeros(0, 4))),
            (ArgumentError, (; md=merge(metadata, (channel="pair",)))),
            (ArgumentError, (; md=merge(metadata, (species_pair=[1, 1],)))),
            (ArgumentError, (; md=merge(metadata, (species_pair=[0, 2],)))),
            (ArgumentError, (; md=merge(metadata, (species_pair=[1.0, 2.0],)))),
            (ArgumentError, (; md=Base.structdiff(metadata, (; channel=metadata.channel)))),
            (ArgumentError, (; md=Base.structdiff(metadata, (; beta=metadata.beta)))),
            (ArgumentError, (; md=merge(metadata, (displacements=zeros(3, 1),)))),
            (DimensionMismatch, (; md=merge(metadata, (displacements=zeros(Int, 2, 1),)))),
            (DimensionMismatch, (; md=merge(metadata, (displacement_vectors=zeros(3, 2),)))),
            (ArgumentError, (; md=merge(metadata, (displacement_counts=[4, 0, 2],)))),
            (DimensionMismatch, (; md=merge(metadata, (displacement_counts=[4, 2],)))),
            (ArgumentError, (; md=merge(metadata, (equal_time_window=10.0,)))),
            (ArgumentError, (; md=merge(metadata, (equal_time_convention="two-sided",)))),
            (ArgumentError, (; md=merge(metadata, (production_log_eta=[-2.0],)))),
            (ArgumentError, (; samples=[1, 2, 0, 4])),
            (ArgumentError, (; reduction="mean")),
            (ArgumentError, (; reduction=nothing)),
        )
        for (error_type, kwargs) in invalid_cases
            write_fixture(; kwargs...)
            @test_throws error_type aggregate(root; observables=["CounterflowGreenFunction"])
        end
        write_fixture(; n=0)
        empty = @test_logs (:warn, r"CounterflowGreenFunction has no saved blocks") (:warn, r"CounterflowGreenFunction has no saved blocks") aggregate(root; observables=["CounterflowGreenFunction"])
        @test all(ismissing, empty.CounterflowGreenFunction)
        write_combined_fixture(path; groups=())
        absent = @test_logs (:warn, r"CounterflowGreenFunction was not found") (:warn, r"CounterflowGreenFunction was not found") aggregate(root; observables=["CounterflowGreenFunction"])
        @test all(ismissing, absent.CounterflowGreenFunction)
    end
end
