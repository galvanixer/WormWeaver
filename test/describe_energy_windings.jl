import Measurements

@testset "Readable Energy report" begin
    m = value -> Measurements.measurement(Float64(value), 0.01)
    summary = Dict{String,Any}(
        "diagonalEnergy" => m(10), "kineticEnergy" => m(-3), "totalEnergy" => m(99),
        "diagonalEnergyPerSite" => m(0.1), "kineticEnergyPerSite" => m(-0.03),
        "totalEnergyPerSite" => m(0.99), "nblocks" => 12, "nsamples" => 1200,
        "error_method" => "repeated blocking: maximum standard error",
    )
    original = deepcopy(summary)
    report = sprint(io -> describe_energy(summary; io))
    @test occursin("Energy (model energy)", report)
    @test occursin("Energy per site (model energy/site)", report)
    for (label, field) in (("Diagonal", "diagonalEnergy"), ("Kinetic", "kineticEnergy"), ("Total", "totalEnergy"))
        line = only(filter(line -> startswith(line, label), split(report, '\n')))
        @test occursin(string(summary[field]), line)
        @test occursin(string(summary[field * "PerSite"]), line)
    end
    @test occursin(r"Saved blocks\s+12", report)
    @test occursin(r"Energy samples\s+1,200", report)
    @test occursin("Uncertainty    Max repeated-blocking SE", report)
    for heading in ("OBSERVABLES", "MEASUREMENTS", "PRODUCTION STATISTICS", "ESTIMATOR CONVENTIONS")
        @test occursin(heading, report)
    end
    for field in ("diagonalEnergy", "kineticEnergy", "totalEnergy"), suffix in ("", "PerSite")
        @test occursin(Regex(field * suffix * raw"\s+scalar\s+—\s+model energy"), report)
    end
    @test occursin("sampleCount-weighted block means", report)
    @test !occursin("error status", report)
    @test describe_energy(summary; io=devnull) === nothing
    @test isequal(summary, original)
    @test sprint(io -> describe_energy(missing; io)) == "No Energy measurements available.\n"
    summary["error_method"] = "custom saved uncertainty method " * repeat("αβγ ", 40)
    report = sprint(io -> describe_energy(summary; io))
    @test occursin(replace(summary["error_method"], r"\s" => ""), replace(report, r"\s" => ""))
    width = textwidth(split(report, '\n')[2])
    @test all(line -> textwidth(line) <= width, split(report, '\n'))
    summary["totalEnergy"] = Measurements.measurement(99.0, NaN)
    report = sprint(io -> describe_energy(summary; io))
    @test occursin("99.0 (uncertainty unavailable)", report)
    @test !occursin("NaN", report)
    summary["totalEnergy"] = [m(99)]
    @test_throws ArgumentError describe_energy(summary; io=devnull)
end

@testset "Readable Windings report" begin
    m = value -> Measurements.measurement(Float64(value), 0.01)
    summary = Dict{String,Any}(
        "windingBySpecies" => m.([1 2; 3 4; 5 6]),
        "squaredWindingBySpecies" => m.([2 5; 10 17; 26 37]),
        "totalWinding" => m.([20, 30]), "squaredTotalWinding" => m.([400, 900]),
        "species_pairs" => [2 3; 1 2], "periodic_axes" => [true, false],
        "windingCorrelationBySpecies" => m.([15 24; 3 8]),
        "squaredCoflowBySpeciesPair" => m.([66 102; 18 38]),
        "squaredCounterflowBySpeciesPair" => m.([6 6; 6 6]),
        "nblocks" => 8, "nsamples" => 800,
        "error_method" => "repeated blocking: maximum standard error",
    )
    original = deepcopy(summary)
    report = sprint(io -> describe_windings(summary; io))
    @test occursin("Windings (dimensionless)", report)
    @test occursin("species IDs are not stored", report)
    @test occursin("Axis 1: periodic", report) && occursin("Axis 2: open", report)
    @test occursin(r"W₁\s+Axis 2\s+2\.0 ± 0\.01\s+5\.0 ± 0\.01", report)
    @test occursin(r"Axis 1\s+20\.0 ± 0\.01\s+400\.0 ± 0\.01", report)
    @test occursin(r"\(W₂, W₃\)\s+Axis 1\s+15\.0 ± 0\.01\s+66\.0 ± 0\.01\s+6\.0 ± 0\.01", report)
    @test first(findfirst("(W₂, W₃)", report)) < first(findfirst("(W₁, W₂)", report))
    @test occursin("⟨Wₐ⟩", report) && occursin("⟨Wₐ²⟩", report)
    @test occursin("Total over all components", report)
    @test occursin("⟨Σₐ Wₐ⟩", report) && occursin("⟨(Σₐ Wₐ)²⟩", report)
    @test occursin("⟨Wₐ Wᵦ⟩", report)
    @test occursin("Coflow: ⟨(Wₐ+Wᵦ)²⟩", report)
    @test occursin("Counterflow: ⟨(Wₐ−Wᵦ)²⟩", report)
    @test occursin("No stiffness normalization applied.", report)
    @test occursin(r"Saved blocks\s+8", report)
    @test occursin(r"Windings samples\s+800", report)
    @test occursin("Uncertainty    Max repeated-blocking SE", report)
    for heading in ("OBSERVABLES", "COMPONENTS AND AXES", "MEASUREMENTS",
                    "PRODUCTION STATISTICS", "ESTIMATOR CONVENTIONS")
        @test occursin(heading, report)
    end
    for field in ("windingBySpecies", "squaredWindingBySpecies", "totalWinding", "squaredTotalWinding",
                  "windingCorrelationBySpecies", "squaredCoflowBySpeciesPair", "squaredCounterflowBySpeciesPair")
        @test !occursin(field, report)
    end
    @test occursin(r"⟨\(Wₐ\+Wᵦ\)²⟩\s+2×2\s+pair × axis\s+dimensionless\s+Pair coflow", report)
    @test occursin("sampleCount-weighted block means", report)
    @test !occursin("error status", report)
    @test describe_windings(summary; io=devnull) === nothing
    @test isequal(summary, original)
    @test sprint(io -> describe_windings(missing; io)) == "No Windings measurements available.\n"

    many_components = deepcopy(summary)
    many_components["windingBySpecies"] = fill(m(1), 12, 2)
    many_components["squaredWindingBySpecies"] = fill(m(2), 12, 2)
    many_components["species_pairs"] = [2 12; 1 10]
    many_report = sprint(io -> describe_windings(many_components; io))
    @test occursin(r"W₁₂\s+Axis 2\s+1\.0 ± 0\.01", many_report)
    @test occursin("(W₂, W₁₂)", many_report) && occursin("(W₁, W₁₀)", many_report)

    delete!(summary, "periodic_axes")
    delete!(summary, "species_pairs")
    report = sprint(io -> describe_windings(summary; io))
    @test occursin("Axis 1: Boundary condition unknown", report)
    @test occursin("Pair row 1 — component mapping unavailable", report)
    @test occursin(string(m(15)), report)

    delete!(summary, "squaredCoflowBySpeciesPair")
    delete!(summary, "squaredCounterflowBySpeciesPair")
    report = sprint(io -> describe_windings(summary; io))
    @test occursin("Unavailable", report)
    @test occursin(r"⟨\(Wₐ\+Wᵦ\)²⟩\s+Not present in this summary", report)
    delete!(summary, "windingCorrelationBySpecies")
    @test occursin("Pair quantities: Unavailable", sprint(io -> describe_windings(summary; io)))

    # Two-component legacy summaries can contain coflow/counterflow without raw pair products.
    summary["windingBySpecies"] = summary["windingBySpecies"][1:2, :]
    summary["squaredWindingBySpecies"] = summary["squaredWindingBySpecies"][1:2, :]
    summary["species_pairs"] = [1 2]
    summary["squaredCoflowBySpeciesPair"] = reshape(m.([18, 38]), 1, 2)
    summary["squaredCounterflowBySpeciesPair"] = reshape(m.([6, 6]), 1, 2)
    report = sprint(io -> describe_windings(summary; io))
    @test occursin(r"\(W₁, W₂\)\s+Axis 1\s+Unavailable\s+18\.0 ± 0\.01", report)
    summary["windingBySpecies"][1, 1] = Measurements.measurement(1.0, NaN)
    report = sprint(io -> describe_windings(summary; io))
    @test occursin("1.0 (uncertainty unavailable)", report)
    @test !occursin("NaN", report)

    for key in ("squaredCoflowBySpeciesPair", "squaredCounterflowBySpeciesPair")
        summary[key] = Matrix{Measurements.Measurement{Float64}}(undef, 0, 2)
    end
    summary["species_pairs"] = zeros(Int, 0, 2)
    summary["windingBySpecies"] = summary["windingBySpecies"][1:1, :]
    summary["squaredWindingBySpecies"] = summary["squaredWindingBySpecies"][1:1, :]
    @test occursin("No component pairs", sprint(io -> describe_windings(summary; io)))

    summary["periodic_axes"] = [true]
    @test_throws DimensionMismatch describe_windings(summary; io=devnull)
    delete!(summary, "periodic_axes")
    summary["squaredTotalWinding"] = [m(1)]
    @test_throws DimensionMismatch describe_windings(summary; io=devnull)
    summary["squaredTotalWinding"] = m.([400, 900])
    summary["species_pairs"] = [1 2]
    @test_throws ArgumentError describe_windings(summary; io=devnull)
    summary["species_pairs"] = zeros(Int, 0, 2)
    summary["windingCorrelationBySpecies"] = fill(m(1), 1, 2)
    @test_throws DimensionMismatch describe_windings(summary; io=devnull)
end
