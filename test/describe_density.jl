import Measurements

@testset "Readable Density report" begin
    summary = Dict{String,Any}(
        "species" => [4, 2],
        "densityBySpecies" => Measurements.measurement.([0.5, 0.25], [0.01, 0.005]),
        "particleNumberBySpecies" => Measurements.measurement.([32.0, 16.0], [0.6, 0.3]),
        "totalParticleNumber" => Measurements.measurement(48.0, 0.8),
        "nblocks" => 100, "nsamples" => 1_000_000,
        "error_method" => "repeated blocking: maximum standard error",
    )
    original = deepcopy(summary)
    report = sprint(io -> describe_density(summary; io))
    @test occursin("Density (particles/site)", report)
    @test occursin("Particle number (particles)", report)
    @test occursin(r"(?m)^4\s+0\.5 ± 0\.01\s+32\.0 ± 0\.6$", report)
    @test occursin(r"(?m)^2\s+0\.25 ± 0\.005\s+16\.0 ± 0\.3$", report)
    @test occursin("Total particle number: 48.0 ± 0.8 particles", report)
    @test occursin(r"Saved blocks\s+100", report)
    @test occursin(r"Density samples\s+1,000,000", report)
    @test occursin("Uncertainty    Max repeated-blocking SE", report)
    for heading in ("OBSERVABLES", "MEASUREMENTS", "PRODUCTION STATISTICS", "ESTIMATOR CONVENTIONS")
        @test occursin(heading, report)
    end
    @test occursin(r"densityBySpecies\s+2\s+species\s+particles/site", report)
    @test occursin(r"totalParticleNumber\s+scalar", report)
    @test occursin("sampleCount-weighted block means", report)
    @test !occursin("error status", report)
    @test !occursin("converged", report)
    @test !occursin("Component", report)
    @test describe_density(summary; io=devnull) === nothing
    @test isequal(summary, original)
    @test sprint(io -> describe_density(missing; io)) == "No Density measurements available.\n"

    delete!(summary, "species")
    report = sprint(io -> describe_density(summary; io))
    @test occursin("Component 1", report) && occursin("Component 2", report)
    @test occursin("Species IDs unavailable", report)

    summary["densityBySpecies"] = Measurements.measurement.([0.5, 0.25], NaN)
    summary["particleNumberBySpecies"] = Measurements.measurement.([32.0, 16.0], NaN)
    summary["totalParticleNumber"] = Measurements.measurement(48.0, NaN)
    summary["nblocks"] = 1
    report = sprint(io -> describe_density(summary; io))
    @test occursin("0.5 (uncertainty unavailable)", report)
    @test occursin("48.0 (uncertainty unavailable)", report)
    @test !occursin("NaN", report)
    @test occursin(r"Saved blocks\s+1", report)

    delete!(summary, "error_method")
    delete!(summary, "nsamples")
    report = sprint(io -> describe_density(summary; io))
    @test occursin("Uncertainty    Unavailable", report)
    @test occursin(r"Density samples\s+Unavailable", report)

    summary["species"] = [1]
    @test_throws DimensionMismatch describe_density(summary; io=devnull)
    summary["species"] = [1, 1]
    @test_throws ArgumentError describe_density(summary; io=devnull)
    summary["species"] = [1, 2]
    summary["particleNumberBySpecies"] = summary["particleNumberBySpecies"][1:1]
    @test_throws DimensionMismatch describe_density(summary; io=devnull)
end
