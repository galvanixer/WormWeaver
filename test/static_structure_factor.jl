@testset "StaticStructureFactor aggregation" begin
    mktempdir() do root
        path = joinpath(root, "results.h5")
        manifest = DataFrames.DataFrame(output_file=["results.h5", "results.h5"], seed=[11, 22])
        manifest[!, "run.blocks"] = [5, 5]
        write_aggregation_manifest(joinpath(root, "manifest.jld2"); manifest)
        counts = [1, 3, 2, 4, 5]
        trend = [0.0, 1, 3, 5, 10]
        pairs = zeros(ComplexF64, 2, 2, 3, 5)
        for k in 1:3
            pairs[1, 1, k, :] = 10k .+ trend
            pairs[2, 2, k, :] = 20k .- trend
            pairs[1, 2, k, :] = k .+ 2im .* trend
            pairs[2, 1, k, :] = k .- 2im .* trend
        end
        total = real.(dropdims(sum(pairs; dims=(1, 2)); dims=(1, 2)))
        md = (momentum_indices=[0 0; 2 0; 1 1], momenta=[0.0 0.0; 4π/3 0.0; 2π/3 π])
        function write_fixture(; n=5, pair_history=pairs[:, :, :, 1:n], total_history=total[:, 1:n],
                               metadata=md, reduction="mean", weights=counts[1:n], split=false)
            HDF5.h5open(path, "w") do file
                HDF5.create_group(file, "observables")
                HDF5.create_group(file, "observables/METADATA")
                fields = split ? (structureFactorBySpeciesPairReal=real.(pair_history),
                                  structureFactorBySpeciesPairImaginary=imag.(pair_history),
                                  structureFactorByMomentum=total_history) :
                    (structureFactorBySpeciesPair=pair_history, structureFactorByMomentum=total_history)
                write_observable_fixture(file, "StaticStructureFactor", fields, weights; metadata, reduction)
                derived = HDF5.create_group(file, "derived")
                derived["StaticStructureFactor"] = "must not be read"
            end
        end
        write_fixture()
        original = read(path)
        table = aggregate(root; observables=["StaticStructureFactor"], min_blocks=2)
        summary = table.StaticStructureFactor[1]
        @test DataFrames.nrow(table) == 2
        @test table.seed == [11, 22]
        @test isequal(table[:, names(manifest)], manifest)
        @test summary["structureFactorBySpeciesPair"] isa Array{Complex{Measurements.Measurement{Float64}},3}
        @test summary["structureFactorByMomentum"] isa Vector{Measurements.Measurement{Float64}}
        avg = sum(counts .* trend) / sum(counts)
        base_error = sqrt(5/4 * sum(abs2, counts .* (trend .- avg))) / sum(counts)
        error = max(base_error, 1.72)
        result = summary["structureFactorBySpeciesPair"]
        for k in 1:3
            @test Measurements.value(real(result[1, 1, k])) ≈ 10k + avg
            @test Measurements.uncertainty(real(result[1, 1, k])) ≈ error
            @test Measurements.value(imag(result[1, 2, k])) ≈ 2avg
            @test Measurements.value(imag(result[2, 1, k])) ≈ -2avg
            @test Measurements.uncertainty(imag(result[1, 2, k])) ≈ 2error
            @test Measurements.value(summary["structureFactorByMomentum"][k]) == 32k
            @test Measurements.uncertainty(summary["structureFactorByMomentum"][k]) == 0
            @test Measurements.uncertainty(real(sum(result[:, :, k]))) > 0
        end
        @test summary["momentum_indices"] == md.momentum_indices
        @test summary["momenta"] == md.momenta
        @test summary["nblocks"] == 5 && summary["nsamples"] == 15
        @test summary["min_blocks"] == 2
        @test summary["error_method"] == "repeated blocking: maximum standard error"
        @test !haskey(summary, "species") && !haskey(summary, "blocking_diagnostics")
        diagnosed = save_aggregate(root; observables=["StaticStructureFactor"], min_blocks=2,
                                   blocking_diagnostics=true)
        retained = copy(diagnosed.StaticStructureFactor[1])
        diagnostics = pop!(retained, "blocking_diagnostics")
        @test isequal(retained, summary)
        @test diagnostics["block_widths"] == [1, 2]
        @test diagnostics["block_counts"] == [5, 2]
        @test diagnostics["sample_counts"] == [15, 10]
        @test size(diagnostics["structureFactorBySpeciesPairRealStandardError"]) == (2, 2, 3, 2)
        @test diagnostics["structureFactorBySpeciesPairRealStandardError"][1, 1, 1, :] ≈ [base_error, 1.72]
        @test diagnostics["structureFactorBySpeciesPairImaginaryStandardError"][1, 2, 1, :] ≈ 2 .* [base_error, 1.72]
        @test size(diagnostics["structureFactorByMomentumStandardError"]) == (3, 2)
        @test all(iszero, diagnostics["structureFactorByMomentumStandardError"])
        @test isequal(JLD2.load(default_aggregate_path(root), "aggregate"), diagnosed)
        @test read(path) == original
        @test observable_table(table.StaticStructureFactor).momenta[1] == md.momenta

        write_fixture(; n=4, reduction=nothing)
        partial = aggregate(root; observables=["StaticStructureFactor"], min_blocks=2).StaticStructureFactor[1]
        @test partial["nblocks"] == 4 && partial["nsamples"] == 10
        @test Measurements.value(imag(partial["structureFactorBySpeciesPair"][1, 2, 1])) ≈ 5.8
        write_fixture(; n=1)
        one = aggregate(root; observables=["StaticStructureFactor"]).StaticStructureFactor[1]
        @test all(isnan, Measurements.uncertainty.(real.(one["structureFactorBySpeciesPair"])))
        @test all(isnan, Measurements.uncertainty.(imag.(one["structureFactorBySpeciesPair"])))
        @test all(isnan, Measurements.uncertainty.(one["structureFactorByMomentum"]))
        @test Measurements.value.(real.(one["structureFactorBySpeciesPair"])) == real.(pairs[:, :, :, 1])

        # Stored entries survive independently, even when physical identities do not hold.
        asymmetric = copy(pairs)
        asymmetric[1, 1, :, :] .+= 3im
        asymmetric[2, 1, :, :] .= -7 + 4im
        write_fixture(; pair_history=asymmetric, total_history=fill(-9.0, 3, 5))
        raw = aggregate(root; observables=["StaticStructureFactor"]).StaticStructureFactor[1]
        @test Measurements.value(imag(raw["structureFactorBySpeciesPair"][1, 1, 1])) == 3
        @test Measurements.value(real(raw["structureFactorBySpeciesPair"][2, 1, 1])) == -7
        @test Measurements.value(raw["structureFactorByMomentum"][1]) == -9

        write_fixture(; pair_history=pairs[1:1, 1:1, :, :])
        @test size(aggregate(root; observables=["StaticStructureFactor"]).StaticStructureFactor[1]["structureFactorBySpeciesPair"]) == (1, 1, 3)
        invalid_cases = (
            (DimensionMismatch, (; pair_history=zeros(ComplexF64, 2, 3, 3, 5))),
            (DimensionMismatch, (; pair_history=zeros(ComplexF64, 2, 3, 5))),
            (DimensionMismatch, (; pair_history=zeros(ComplexF64, 0, 0, 3, 5))),
            (DimensionMismatch, (; pair_history=zeros(ComplexF64, 2, 2, 0, 5))),
            (DimensionMismatch, (; total_history=zeros(2, 5))),
            (DimensionMismatch, (; total_history=zeros(3, 4))),
            (DimensionMismatch, (; metadata=(momentum_indices=zeros(Int, 3, 1), momenta=md.momenta))),
            (DimensionMismatch, (; metadata=(momentum_indices=zeros(Int, 2, 2), momenta=zeros(2, 2)))),
            (DimensionMismatch, (; metadata=(momentum_indices=zeros(Int, 3, 0), momenta=zeros(3, 0)))),
            (ArgumentError, (; metadata=(momenta=md.momenta,))),
            (ArgumentError, (; metadata=(momentum_indices=md.momentum_indices,))),
            (ArgumentError, (; metadata=(momentum_indices=Float64.(md.momentum_indices), momenta=md.momenta))),
            (ArgumentError, (; metadata=(momentum_indices=md.momentum_indices, momenta=fill(NaN, 3, 2)))),
            (ArgumentError, (; pair_history=real.(pairs))),
            (ArgumentError, (; pair_history=fill(ComplexF64(NaN), size(pairs)))),
            (ArgumentError, (; total_history=fill(Inf, 3, 5))),
            (ArgumentError, (; weights=[1, 3, 0, 4, 5])),
            (ArgumentError, (; reduction="sum")),
            (ArgumentError, (; split=true)),
        )
        for (error_type, kwargs) in invalid_cases
            write_fixture(; kwargs...)
            @test_throws error_type aggregate(root; observables=["StaticStructureFactor"])
        end
        write_fixture(; n=0)
        empty = @test_logs (:warn, r"StaticStructureFactor has no saved blocks") (:warn, r"StaticStructureFactor has no saved blocks") aggregate(root; observables=["StaticStructureFactor"])
        @test all(ismissing, empty.StaticStructureFactor)
        write_combined_fixture(path; groups=())
        absent = @test_logs (:warn, r"StaticStructureFactor was not found") (:warn, r"StaticStructureFactor was not found") aggregate(root; observables=["StaticStructureFactor"])
        @test all(ismissing, absent.StaticStructureFactor)
    end
end
