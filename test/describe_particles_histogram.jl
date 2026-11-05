import Measurements

@testset "Readable ParticlesHistogram report" begin
    function make_summary(; nblocks=2, nbins=5)
        species_counts = zeros(Int64, 2, nbins)
        species_counts[:, 1:5] = [2 0 2 0 0; 1 0 0 0 3]
        histories = Dict(
            "particleNumberBySpeciesCounts" => repeat(species_counts; outer=(1, 1, nblocks)),
            "totalParticleNumberCounts" => repeat(Int64[1, 0, 0, 0, 1, 0, 2]; outer=(1, nblocks)),
            "zSectorSampleCount" => fill(Int64(4), nblocks),
        )
        metadata = Dict(
            "species" => [7, 2], "species_particle_number_bins" => collect(0:nbins-1),
            "total_particle_number_bins" => collect(0:6), "max_particle_number_by_species" => [2, 4],
        )
        return WormWeaver.particle_histogram_summary(histories, fill(4, nblocks), metadata; min_blocks=10)
    end
    summary = make_summary()
    original = deepcopy(summary)
    report = sprint(io -> describe_particles_histogram(summary; io))
    for heading in ("OBSERVABLES", "SPECIES AND BINS", "MOMENTS", "DISTRIBUTION OVERVIEW",
                    "PRODUCTION STATISTICS", "ESTIMATOR CONVENTIONS")
        @test occursin(heading, report)
    end
    @test occursin(r"particleNumberBySpeciesCounts\s+2×5\s+species × bin\s+counts", report)
    @test occursin(r"totalParticleNumberProbability\s+7\s+bin\s+dimensionless", report)
    @test occursin(r"Species 7\s+5\s+0–4\s+2", report)
    @test occursin(r"Species 2\s+5\s+0–4\s+4", report)
    @test first(findfirst("Species 7", report)) < first(findfirst("Species 2", report))
    moments = split(split(report, "\nMOMENTS\n")[2], "\nDISTRIBUTION OVERVIEW\n")[1]
    for s in 1:2
        row = only(filter(line -> startswith(line, "Species $(summary["species"][s])"), split(moments, '\n')))
        @test occursin(string(summary["meanParticleNumberBySpecies"][s]), row)
        @test occursin(string(summary["meanSquaredParticleNumberBySpecies"][s]), row)
    end
    total_row = only(filter(line -> startswith(line, "Total"), split(moments, '\n')))
    @test occursin(string(summary["meanTotalParticleNumber"]), total_row)
    @test occursin(string(summary["meanSquaredTotalParticleNumber"]), total_row)
    @test occursin(r"Species 7\s+2\s+0–2\s+0 \(tie: 2 bins\)\s+0\.5 ± 0\.0", report)
    @test occursin(r"Species 2\s+2\s+0–4\s+4\s+0\.75 ± 0\.0", report)
    @test occursin(r"Total\s+3\s+0–6\s+6\s+0\.5 ± 0\.0", report)
    @test occursin("Raw bin count / nsamples", report)
    @test occursin("Max repeated-blocking SE", report)
    @test occursin("no sample-count multiplication", report)
    @test occursin("Raw-count errors None stored", report)
    @test occursin("Second moments", report)
    @test !occursin("error status", report)
    @test occursin(r"Saved blocks\s+2", report)
    @test occursin(r"ParticlesHistogram samples\s+8", report)
    @test occursin(r"Z-sector configurations\s+8", report)
    @test describe_particles_histogram(summary; io=devnull) === nothing
    @test isequal(summary, original)
    @test sprint(io -> describe_particles_histogram(missing; io)) ==
        "No ParticlesHistogram measurements available.\n"

    # A wide stored axis must not expand into thousands of terminal rows.
    large = make_summary(; nbins=10_000)
    large_report = sprint(io -> describe_particles_histogram(large; io))
    @test occursin("10,000", large_report)
    @test occursin("0–9999", large_report)
    @test length(split(large_report, '\n')) == length(split(report, '\n'))

    single = make_summary(; nblocks=1)
    single_report = sprint(io -> describe_particles_histogram(single; io))
    @test occursin("0.5 (uncertainty unavailable)", single_report)
    @test !occursin("NaN", single_report)

    # Stored values are presented independently of the counts used for the overview.
    summary["meanTotalParticleNumber"] = Measurements.measurement(99.0, 0.1)
    summary["particleNumberBySpeciesProbability"][1, 1] = Measurements.measurement(0.123, 0.001)
    report = sprint(io -> describe_particles_histogram(summary; io))
    @test occursin("99.0 ± 0.1", report)
    @test occursin(r"0 \(tie: 2 bins\)\s+0\.123 ± 0\.001", report)

    partial = deepcopy(original)
    for key in ("species", "species_particle_number_bins", "total_particle_number_bins",
                "max_particle_number_by_species", "aggregation_method", "error_method")
        delete!(partial, key)
    end
    report = sprint(io -> describe_particles_histogram(partial; io))
    @test occursin("Species IDs unavailable", report)
    @test occursin(r"Component 1\s+2\s+Unavailable\s+Unavailable \(tie: 2 bins\)", report)
    @test occursin("Uncertainty      Unavailable", report)

    summary = deepcopy(original)
    summary["species"] = [7]
    @test_throws DimensionMismatch describe_particles_histogram(summary; io=devnull)
    summary["species"] = [7, 7]
    @test_throws ArgumentError describe_particles_histogram(summary; io=devnull)
    summary = deepcopy(original)
    summary["particleNumberBySpeciesProbability"] = summary["particleNumberBySpeciesProbability"][:, 1:2]
    @test_throws DimensionMismatch describe_particles_histogram(summary; io=devnull)
    summary = deepcopy(original)
    summary["species_particle_number_bins"] = [0, 2, 1, 3, 4]
    @test_throws ArgumentError describe_particles_histogram(summary; io=devnull)
    summary = deepcopy(original)
    summary["particleNumberBySpeciesCounts"][1, 5] = 1
    @test_throws ArgumentError describe_particles_histogram(summary; io=devnull)
    summary = deepcopy(original)
    summary["totalParticleNumberCounts"][1] = -1
    @test_throws ArgumentError describe_particles_histogram(summary; io=devnull)
end
