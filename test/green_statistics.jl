function reference_green_jackknife(numerator, z; transform=identity)
    n = length(z)
    replicas = [transform((sum(numerator) - numerator[b]) / (sum(z) - z[b])) for b in 1:n]
    center = sum(replicas) / n
    real_error = sqrt((n - 1) / n * sum((real(x) - real(center))^2 for x in replicas))
    imag_error = sqrt((n - 1) / n * sum((imag(x) - imag(center))^2 for x in replicas))
    return real_error, imag_error
end

@testset "Reusable paired Green statistics" begin
    stats = green_ratio_statistics([2.0, 8, 10], [1, 2, 1]; transform=x -> 2x)
    @test Measurements.value(stats.measurements) == 10
    @test Measurements.uncertainty(stats.measurements) ≈ 32/9
    @test stats.error_status == "ok"
    @test stats.blocking_diagnostics === nothing

    proportional = green_ratio_statistics(3 .* [0, 1, 2, 3], [0, 1, 2, 3]; min_blocks=2)
    @test Measurements.value(proportional.measurements) == 3
    @test Measurements.uncertainty(proportional.measurements) == 0

    numerator = vcat(zeros(32), ones(32))
    stats = green_ratio_statistics(numerator, ones(64); min_blocks=4, blocking_diagnostics=true)
    diag = stats.blocking_diagnostics
    @test diag["block_widths"] == [1, 2, 4, 8, 16]
    @test diag["block_counts"] == [64, 32, 16, 8, 4]
    @test diag["real_standard_errors"] ≈ [0.5 / sqrt(n - 1) for n in diag["block_counts"]]
    @test Measurements.value(stats.measurements) == 0.5
    @test Measurements.uncertainty(stats.measurements) ≈ 1 / sqrt(12)
    @test !haskey(diag, "sample_counts")
    plain = green_ratio_statistics(numerator, ones(64); min_blocks=4)
    @test isequal(plain.measurements, stats.measurements)

    odd = green_ratio_statistics([0, 0, 0, 0, 10], ones(5);
        sample_counts=[1, 2, 3, 4, 5], min_blocks=2, blocking_diagnostics=true)
    @test Measurements.value(odd.measurements) == 2
    @test Measurements.uncertainty(odd.measurements) ≈ 2
    @test odd.blocking_diagnostics["sample_counts"] == [15, 10]
    @test odd.blocking_diagnostics["z_counts"] == [5, 4]
    @test Measurements.value(green_ratio_statistics([0, 0, 0, 0, 10], ones(5)).measurements) == 2

    counts = [1, 3, 2, 4]
    means, z = [-2.0, 1, 3, 5], [0, 3, 2, 1]
    summed = means .* counts
    mean_stats = green_ratio_statistics(means, z; reduction=:mean, sample_counts=counts, min_blocks=2)
    sum_stats = green_ratio_statistics(summed, z; sample_counts=counts, min_blocks=2)
    @test isequal(mean_stats.measurements, sum_stats.measurements)
    @test Measurements.value(mean_stats.measurements) == sum(summed) / sum(z)
    @test Measurements.value(green_ratio_statistics(summed, z; sample_counts=100 .* counts).measurements) == sum(summed) / sum(z)
    large = green_ratio_statistics([1, 1, 1, 1], ones(4);
        sample_counts=[Int64(2)^53 + 1, 1, 1, 1], min_blocks=2, blocking_diagnostics=true)
    @test large.blocking_diagnostics["sample_counts"] == [Int64(2)^53 + 4, Int64(2)^53 + 4]

    complex_values = ComplexF64[1 + 2im, -3 + 4im, 5 - 1im]
    stats = green_ratio_statistics(complex_values, [1, 2, 3]; blocking_diagnostics=true)
    @test stats.measurements isa Complex{Measurements.Measurement{Float64}}
    expected_real, expected_imag = reference_green_jackknife(complex_values, [1, 2, 3])
    @test Measurements.value(real(stats.measurements)) == real(sum(complex_values) / 6)
    @test Measurements.value(imag(stats.measurements)) == imag(sum(complex_values) / 6)
    @test Measurements.uncertainty(real(stats.measurements)) ≈ expected_real
    @test Measurements.uncertainty(imag(stats.measurements)) ≈ expected_imag
    @test stats.blocking_diagnostics["imaginary_standard_errors"] ≈ [expected_imag]

    # A joint, nonlinear transform is evaluated on paired replicas, not on
    # independently uncertain output components.
    history = reshape(collect(1.0:24), 2, 3, 4)
    stats = green_ratio_statistics(history, [1, 2, 1, 3]; transform=x -> sum(x)^2)
    totals = vec(sum(history; dims=(1, 2)))
    expected, _ = reference_green_jackknife(totals, [1, 2, 1, 3]; transform=x -> x^2)
    @test Measurements.value(stats.measurements) == (sum(history) / 7)^2
    @test Measurements.uncertainty(stats.measurements) ≈ expected
    array_stats = green_ratio_statistics(complex.(history, -history), [1, 2, 1, 3]; blocking_diagnostics=true)
    @test size(array_stats.measurements) == (2, 3)
    @test size(array_stats.blocking_diagnostics["real_standard_errors"]) == (2, 3, 1)
    mktempdir() do root
        path = joinpath(root, "complex.jld2")
        JLD2.jldsave(path; statistics=array_stats)
        @test isequal(JLD2.load(path, "statistics"), array_stats)
    end

    partial = green_ratio_statistics([1, 2, 3, 4], [1, 1, 0, 0]; min_blocks=2, blocking_diagnostics=true)
    @test isfinite(Measurements.uncertainty(partial.measurements))
    @test partial.error_status == "invalid blocking levels omitted"
    @test isnan(partial.blocking_diagnostics["real_standard_errors"][2])
    unavailable = green_ratio_statistics([1, 2, 3, 4, 5], [0, 0, 0, 0, 1]; min_blocks=2)
    @test Measurements.value(unavailable.measurements) == 15
    @test isnan(Measurements.uncertainty(unavailable.measurements))
    @test unavailable.error_status == "zero leave-one-out Z count"
    one = green_ratio_statistics([2], [1])
    @test Measurements.value(one.measurements) == 2
    @test isnan(Measurements.uncertainty(one.measurements))
    @test one.error_status == "fewer than two blocks"

    @test_throws ArgumentError green_ratio_statistics([1, 2], [0, 0])
    @test_throws ArgumentError green_ratio_statistics([1, 2], [-1, 2])
    @test_throws ArgumentError green_ratio_statistics([1, NaN], [1, 1])
    @test_throws ArgumentError green_ratio_statistics([1, 2], [1, Inf])
    @test_throws ArgumentError green_ratio_statistics([1, 2], [1, 1]; reduction=:other)
    @test_throws ArgumentError green_ratio_statistics([1, 2], [1, 1]; reduction=:mean)
    @test_throws ArgumentError green_ratio_statistics([1, 2], [1, 1]; sample_counts=[1, 1.5])
    @test_throws ArgumentError green_ratio_statistics([1, 2], [1, 1]; sample_counts=[1, 0])
    @test_throws ArgumentError green_ratio_statistics([1, 2], [1, 1]; transform=x -> Inf)
    @test_throws ArgumentError green_ratio_statistics([1], [1]; min_blocks=1)
    @test_throws DimensionMismatch green_ratio_statistics([1, 2], [1])
    @test_throws DimensionMismatch green_ratio_statistics([1, 2], [1, 1]; sample_counts=[1])
    @test_throws DimensionMismatch green_ratio_statistics([1, 2], [1, 1]; transform=x -> x == 1.5 ? [x] : [x, x])
    @test_throws ArgumentError green_ratio_statistics(Float64[], Float64[])
end
