@testset "DensityCorrelation aggregation" begin
    mktempdir() do root
        path = joinpath(root, "results.h5")
        manifest = DataFrames.DataFrame(output_file=["results.h5", "results.h5"], seed=[11, 22])
        manifest[!, "run.blocks"] = [5, 5]
        write_aggregation_manifest(joinpath(root, "manifest.jld2"); manifest)
        counts = [1, 3, 2, 4, 5]
        trend = [0.0, 1, 3, 5, 10]
        history = zeros(2, 2, 3, 5)
        for r in 1:3
            history[1, 1, r, :] = 10r .+ trend
            history[2, 2, r, :] = 20r .- trend
            history[1, 2, r, :] .= r
            history[2, 1, r, :] .= r
        end
        total = dropdims(sum(history; dims=(1, 2)); dims=(1, 2))
        metadata = (displacements=reshape([0, 1, -1], 3, 1), displacement_counts=[4, 2, 2])
        function write_fixture(; nblocks=5, fields=(full=history[:, :, :, 1:nblocks], total=total[:, 1:nblocks]),
                               md=metadata, reduction="mean")
            HDF5.h5open(path, "w") do file
                HDF5.create_group(file, "observables")
                HDF5.create_group(file, "observables/METADATA")
                write_observable_fixture(file, "DensityCorrelation", fields, counts[1:nblocks]; metadata=md, reduction)
                derived = HDF5.create_group(file, "derived")
                derived["DensityCorrelation"] = "must not be read"
            end
        end
        write_fixture()
        original = read(path)
        table = aggregate(root; observables=["DensityCorrelation"], min_blocks=2)
        @test DataFrames.nrow(table) == 2
        @test table.seed == [11, 22]
        @test isequal(table[:, names(manifest)], manifest)
        @test names(table)[end:end] == ["DensityCorrelation"]
        summary = table.DensityCorrelation[1]
        @test summary["full"] isa Array{Measurements.Measurement{Float64},3}
        @test size(summary["full"]) == (2, 2, 3)
        @test summary["total"] isa Vector{Measurements.Measurement{Float64}}
        average = sum(counts .* trend) / sum(counts)
        original_error = sqrt(5 / 4 * sum(abs2, counts .* (trend .- average))) / sum(counts)
        expected_error = max(original_error, 1.72)
        for r in 1:3
            @test Measurements.value(summary["full"][1, 1, r]) ≈ 10r + average
            @test Measurements.value(summary["full"][2, 2, r]) ≈ 20r - average
            @test Measurements.value(summary["full"][1, 2, r]) == r
            @test Measurements.value(summary["full"][2, 1, r]) == r
            @test Measurements.uncertainty(summary["full"][1, 1, r]) ≈ expected_error
            @test Measurements.uncertainty(summary["full"][2, 2, r]) ≈ expected_error
            @test Measurements.value(summary["total"][r]) == 32r
            @test Measurements.uncertainty(summary["total"][r]) == 0
            @test Measurements.uncertainty(sum(summary["full"][:, :, r])) > 0
        end
        @test summary["displacements"] == metadata.displacements
        @test summary["displacement_counts"] == [4, 2, 2]
        @test !haskey(summary, "species")
        @test summary["nblocks"] == 5
        @test summary["nsamples"] == 15
        @test summary["min_blocks"] == 2
        @test summary["error_method"] == "repeated blocking: maximum standard error"
        @test !haskey(summary, "blocking_diagnostics")

        diagnosed = save_aggregate(root; observables=["DensityCorrelation"], min_blocks=2,
                                   blocking_diagnostics=true)
        diagnostic_summary = diagnosed.DensityCorrelation[1]
        for field in ("full", "total")
            @test isequal(Measurements.value.(diagnostic_summary[field]), Measurements.value.(summary[field]))
            @test Measurements.uncertainty.(diagnostic_summary[field]) ≈ Measurements.uncertainty.(summary[field])
        end
        diagnostics = diagnostic_summary["blocking_diagnostics"]
        @test diagnostics["block_widths"] == [1, 2]
        @test diagnostics["block_counts"] == [5, 2]
        @test diagnostics["sample_counts"] == [15, 10]
        @test size(diagnostics["fullStandardError"]) == (2, 2, 3, 2)
        @test size(diagnostics["totalStandardError"]) == (3, 2)
        @test diagnostics["fullStandardError"][1, 1, 1, :] ≈ [original_error, 1.72]
        @test all(iszero, diagnostics["totalStandardError"])
        saved = JLD2.load(joinpath(root, "aggregated.jld2"), "aggregate")
        @test isequal(saved, diagnosed)
        @test read(path) == original

        # Missing reduction markers on legacy mean histories remain supported.
        write_fixture(; nblocks=4, reduction=nothing)
        partial = aggregate(root; observables=["DensityCorrelation"], min_blocks=2).DensityCorrelation[1]
        @test partial["nblocks"] == 4
        @test partial["nsamples"] == 10
        @test Measurements.value(partial["full"][1, 1, 1]) ≈ 12.9
        @test Measurements.uncertainty(partial["full"][1, 1, 1]) ≈ 1.72
        write_fixture(; nblocks=1)
        one = aggregate(root; observables=["DensityCorrelation"]).DensityCorrelation[1]
        @test Measurements.value.(one["full"]) == history[:, :, :, 1]
        @test all(isnan, Measurements.uncertainty.(one["full"]))
        @test all(isnan, Measurements.uncertainty.(one["total"]))

        # One component and multiple coordinate dimensions retain their axes.
        write_fixture(; fields=(full=history[1:1, 1:1, :, :], total=reshape(history[1, 1, :, :], 3, 5)),
                      md=(displacements=[0 0; 1 0; -1 0], displacement_counts=[4, 2, 2]))
        single = aggregate(root; observables=["DensityCorrelation"]).DensityCorrelation[1]
        @test size(single["full"]) == (1, 1, 3)
        @test size(single["displacements"]) == (3, 2)

        # Shared reader checks and correlation-specific shape/metadata checks.
        invalid_cases = (
            (DimensionMismatch, (; fields=(full=zeros(2, 3, 3, 5), total))),
            (DimensionMismatch, (; fields=(full=zeros(2, 2, 5), total))),
            (DimensionMismatch, (; fields=(full=history, total=zeros(2, 5)))),
            (DimensionMismatch, (; fields=(full=zeros(0, 0, 3, 5), total))),
            (DimensionMismatch, (; fields=(full=zeros(2, 2, 3, 4), total))),
            (DimensionMismatch, (; md=(displacements=zeros(Int, 2, 1), displacement_counts=[4, 2, 2]))),
            (DimensionMismatch, (; md=(displacements=metadata.displacements, displacement_counts=[4, 2]))),
            (ArgumentError, (; md=(displacement_counts=[4, 2, 2],))),
            (ArgumentError, (; md=(displacements=metadata.displacements,))),
            (ArgumentError, (; md=(displacements=Float64.(metadata.displacements), displacement_counts=[4, 2, 2]))),
            (ArgumentError, (; md=(displacements=zeros(Int, 3, 0), displacement_counts=[4, 2, 2]))),
            (ArgumentError, (; md=(displacements=metadata.displacements, displacement_counts=[4, 0, 2]))),
            (ArgumentError, (; fields=(total,))),
            (ArgumentError, (; fields=(full=fill(NaN, size(history)), total))),
            (ArgumentError, (; reduction="sum")),
        )
        for (error, kwargs) in invalid_cases
            write_fixture(; kwargs...)
            @test_throws error aggregate(root; observables=["DensityCorrelation"])
        end
        write_fixture(; nblocks=0)
        empty = @test_logs (:warn, r"DensityCorrelation has no saved blocks") (:warn, r"DensityCorrelation has no saved blocks") aggregate(root; observables=["DensityCorrelation"])
        @test all(ismissing, empty.DensityCorrelation)
        write_combined_fixture(path; groups=())
        absent = @test_logs (:warn, r"DensityCorrelation was not found") (:warn, r"DensityCorrelation was not found") aggregate(root; observables=["DensityCorrelation"])
        @test all(ismissing, absent.DensityCorrelation)
    end
end
