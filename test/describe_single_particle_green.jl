import Measurements

@testset "Readable SingleParticleGreenFunction report" begin
    summary = Dict{String,Any}(
        "gKTau" => fill(complex(Measurements.measurement(1.0, 0.1),
                                Measurements.measurement(-2.0, 0.2)), 2, 2, 2),
        "grEqualTime" => fill(Measurements.measurement(0.3, 0.01), 2, 3),
        "gREqualTimeHits" => reshape(collect(11:16), 2, 3),
        "species" => [7, 2], "momenta" => [1.0 2.0; 3.0 4.0],
        "displacements" => [-1 0; 0 0; 1 0],
        "displacement_vectors" => [-0.5 0.0; 0.0 0.0; 0.5 0.0],
        "tau_edges" => [0.25, 0.75, 3.0], "beta" => 64.0,
        "nblocks" => 6, "sampleCount" => 1000, "nsamples" => 1000,
        "zSectorVisitCount" => 111,
        "gKTauAggregationMethod" => "sample-count-weighted mean",
        "gKTauErrorMethod" => "repeated blocking: maximum standard error",
        "error_method" => "repeated blocking: maximum paired jackknife standard error",
        "error_status" => "ok", "normalization" => "EQUAL_TIME_ONLY_FORMULA",
        "equal_time_window" => 0.05, "equal_time_convention" => "one-sided inclusive window",
    )
    original = deepcopy(summary)
    report = sprint(io -> describe_single_particle_green(summary; io))
    for heading in ("OBSERVABLES", "COORDINATE GRIDS", "MEASUREMENT EXAMPLES",
                    "PRODUCTION STATISTICS", "ESTIMATOR CONVENTIONS")
        @test occursin(heading, report)
    end
    @test occursin(r"gKTau\s+2×2×2\s+species × k × τ\s+Raw weighted mean", report)
    @test occursin(r"grEqualTime\s+2×3\s+species × r\s+Normalized estimate", report)
    @test occursin(r"gREqualTimeHits\s+2×3\s+species × r\s+Raw hits", report)
    @test occursin(r"Species IDs\s+2\s+7\s+2", report)
    @test occursin(r"Momentum\s+2\s+\[1\.0000, 2\.0000\]\s+\[3\.0000, 4\.0000\]", report)
    @test occursin(r"Displacement\s+3\s+\[-1, 0\]\s+\[1, 0\]", report)
    @test occursin(r"τ bins\s+2\s+0\.2500\s+3\.0000", report)
    @test occursin("Coordinate grids do not imply measurement availability", report)
    for key in ("species", "momenta", "displacements", "displacement_vectors", "tau_edges")
        @test occursin("summary[\"$key\"]", report)
    end
    @test occursin("first stored entries, not typical values", report)
    @test occursin("gKTau[1,1,1]    species=7  k=[1.0000, 2.0000]  τ=[0.2500, 0.7500]", report)
    @test occursin("grEqualTime[1,1]    species=7  r=[-1, 0]", report)
    @test occursin("gREqualTimeHits[1,1]    species=7  r=[-1, 0]", report)
    @test occursin(r"Real\s+1\.0\s+± 0\.1", report)
    @test occursin(r"Imaginary\s+-2\.0\s+± 0\.2", report)
    @test occursin(r"Value\s+0\.3\s+± 0\.01\s+error status: ok", report)
    @test occursin(r"Hits\s+11\s+raw", report)
    @test occursin("Aggregation    sampleCount-weighted block means", report)
    @test occursin("Max repeated-blocking SE (Re/Im separately)", report)
    @test occursin("Max repeated-blocking paired jackknife SE", report)
    @test occursin("Window         0.05", report)
    @test occursin("Convention     one-sided inclusive window", report)
    @test occursin("Normalization  EQUAL_TIME_ONLY_FORMULA", report)
    @test occursin(r"Saved blocks\s+6", report)
    @test occursin(r"Green-function samples\s+1,000", report)
    @test occursin(r"Z-sector visits\s+111", report)
    separator_widths = [textwidth(line) for line in split(report, '\n')
                        if !isempty(line) && all(character -> character in ('═', '─'), line)]
    @test !isempty(separator_widths)
    @test all(==(first(separator_widths)), separator_widths)
    @test all(line -> textwidth(line) <= first(separator_widths), split(report, '\n'))
    @test describe_single_particle_green(summary; io=devnull) === nothing
    @test isequal(summary, original)
    @test sprint(io -> describe_single_particle_green(missing; io)) ==
        "No SingleParticleGreenFunction measurements available.\n"

    raw_only = deepcopy(summary)
    for key in ("grEqualTime", "gREqualTimeHits", "error_method", "normalization")
        delete!(raw_only, key)
    end
    report = sprint(io -> describe_single_particle_green(raw_only; io))
    @test occursin("gKTau[1,1,1]", report)
    @test count("Not present in this summary", report) == 2
    @test !occursin("paired jackknife", report)
    @test !occursin("EQUAL_TIME_ONLY_FORMULA", report)

    # Older aggregates can retain coordinate grids without any of these quantities.
    legacy = deepcopy(raw_only)
    delete!(legacy, "gKTau")
    report = sprint(io -> describe_single_particle_green(legacy; io))
    @test count("Not present in this summary", report) == 3
    @test occursin(r"Momentum\s+2\s+\[1\.0000, 2\.0000\]", report)
    @test !occursin("not measured", lowercase(report))

    for key in ("species", "momenta", "tau_edges", "displacements", "displacement_vectors")
        delete!(raw_only, key)
    end
    report = sprint(io -> describe_single_particle_green(raw_only; io))
    @test occursin(r"Species IDs\s+Unavailable\s+Unavailable\s+Unavailable", report)
    @test occursin(r"Momentum\s+Unavailable\s+Unavailable\s+Unavailable", report)
    @test occursin(r"τ bins\s+Unavailable\s+Unavailable\s+Unavailable", report)
    @test occursin("species=Unavailable", report)
    @test occursin("k=Unavailable", report)

    summary["gKTau"] .= complex(Measurements.measurement(1.0, NaN), Measurements.measurement(2.0, NaN))
    summary["grEqualTime"] .= Measurements.measurement(0.3, NaN)
    report = sprint(io -> describe_single_particle_green(summary; io))
    @test occursin(r"Real\s+1\.0 \(uncertainty unavailable\)", report)
    @test occursin(r"Imaginary\s+2\.0 \(uncertainty unavailable\)", report)
    @test occursin(r"Value\s+0\.3 \(uncertainty unavailable\)", report)
    @test !occursin("NaN", report)

    formatting = deepcopy(original)
    formatting["momenta"][1, 1] = 1.0e-8
    formatting["momenta"][end, end] = 1.0e7
    formatting["normalization"] =
        "nsites*beta/(C_r*equal_time_window*exp(production_log_eta_s)) * sum(hits_s_r)/sum(z_counts)"
    report = sprint(io -> describe_single_particle_green(formatting; io))
    @test occursin("1.0000e-08", report)
    @test occursin("1.0000e+07", report)
    @test occursin("production_log_eta_s", report)
    formatting["normalization"] = repeat("LONG_NORMALIZATION_TOKEN", 8)
    report = sprint(io -> describe_single_particle_green(formatting; io))
    @test occursin(formatting["normalization"], replace(report, r"\s" => ""))
    report_width = textwidth(split(report, '\n')[2])
    @test all(line -> textwidth(line) <= report_width, split(report, '\n'))

    summary["species"] = [7]
    @test_throws DimensionMismatch describe_single_particle_green(summary; io=devnull)
    summary["species"] = [7, 2]
    summary["tau_edges"] = [0.0, 1.0]
    @test_throws DimensionMismatch describe_single_particle_green(summary; io=devnull)
    summary["tau_edges"] = [0.25, 0.75, 3.0]
    summary["momenta"] = zeros(3, 2)
    @test_throws DimensionMismatch describe_single_particle_green(summary; io=devnull)
    summary["momenta"] = [1.0 2.0; 3.0 4.0]
    summary["gKTau"] = ones(2, 2)
    @test_throws DimensionMismatch describe_single_particle_green(summary; io=devnull)
end
