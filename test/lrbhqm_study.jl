import DataFrames
import JLD2
import YAML

function curated_crossing(t, mu; include=true, beta=8.0)
    (; model_id="alpha1p3", g2_mode="none", nsites=64, beta,
        hopping_alpha_1=1.3, t_over_U=t, mu_at_unit_density=mu,
        mu_lower=mu - 0.02, mu_upper=mu + 0.02, include)
end

@testset "LRBHQM unit-density estimates" begin
    data = DataFrames.DataFrame([
        curated_crossing(0.02, 0.2),
        curated_crossing(0.03, 0.3; include=false),
        curated_crossing(0.04, 0.4),
    ])
    model = (hopping_alpha_1=1.3,)

    sampled = WormWeaver.LRBHQM.estimate_mu_over_U(
        data; model, nsites=64, beta=8.0, t_over_U=0.02)
    @test sampled.method == :sampled
    @test sampled.mu_over_U == 0.2

    bridged = WormWeaver.LRBHQM.estimate_mu_over_U(
        data; model, nsites=64, beta=8.0, t_over_U=0.03)
    @test bridged.method == :interpolated
    @test bridged.mu_over_U ≈ 0.3
    @test bridged.supporting_t_over_U == (0.02, 0.04)

    beta16 = DataFrames.DataFrame([
        curated_crossing(0.02, 0.6; beta=16.0),
        curated_crossing(0.03, 0.7; beta=16.0),
    ])
    extrapolated = WormWeaver.LRBHQM.estimate_mu_over_U(
        vcat(data, beta16); model, nsites=64, beta=8.0, t_over_U=0.01)
    @test extrapolated.method == :extrapolated
    @test extrapolated.beta == 16.0
    @test extrapolated.mu_over_U ≈ 0.5
    fallback = WormWeaver.LRBHQM.estimate_mu_over_U(
        vcat(data, beta16); model, nsites=64, beta=12.0, t_over_U=0.025)
    @test fallback.beta == 16.0
    @test fallback.mu_over_U ≈ 0.65
    @test WormWeaver.LRBHQM._lookup_mode(Dict(), false) == "none"
    @test_throws ArgumentError WormWeaver.LRBHQM._lookup_mode(Dict(), true)
    @test_throws ArgumentError WormWeaver.LRBHQM.estimate_mu_over_U(
        data; model, nsites=64, beta=8.0, t_over_U=0.05)
end

@testset "LRBHQM study generation" begin
    mktempdir() do root
        template = Dict(
            "system" => Dict("beta" => 8.0, "ensemble" => "grand_canonical"),
            "lattice" => Dict("geometry" => "chain", "extent" => [64]),
            "model" => Dict(
                "family" => "bose_hubbard",
                "nspecies" => 2,
                "hopping" => Dict(
                    "kind" => "power_law",
                    "J" => [0.6, 0.6],
                    "alpha" => [1.3, 1.3],
                    "distance_cutoff" => [nothing, nothing],
                ),
                "onsite" => Dict(
                    "U" => [[20.0, 10.0], [10.0, 20.0]],
                    "mu" => [0.0, 0.0],
                ),
            ),
            "run" => Dict("seed" => 42),
            "sampler" => Dict(
                "multiworm_order" => 1,
                "moves" => [
                    Dict("source" => [0, 0], "requests" => [
                        Dict("kind" => "Open", "target" => [1, 0], "weight" => 0.5),
                        Dict("kind" => "Open", "target" => [0, 1], "weight" => 0.5),
                    ]),
                    Dict("source" => [1, 0], "requests" => [
                        Dict("kind" => "Close", "target" => [0, 0], "weight" => 1.0),
                    ]),
                    Dict("source" => [0, 1], "requests" => [
                        Dict("kind" => "Close", "target" => [0, 0], "weight" => 1.0),
                    ]),
                ],
            ),
        )
        template_path = joinpath(root, "template.yaml")
        YAML.write_file(template_path, template)
        definition = Dict{String,Any}(
            "append_date" => false,
            "project" => "LRBHQM",
            "template" => "./template.yaml",
            "output_dir" => "generated",
            "unit_density_data" => "alpha1p3_U12overU11_0p5.jld2",
            "sweep" => Dict("model.hopping.J" => [[0.6, 0.6]]),
        )
        source = joinpath(root, "study.yaml")
        YAML.write_file(source, definition)

        output = WormWeaver.LRBHQM.generate_study(source; match_sampling_mode=true)
        config = YAML.load_file(joinpath(output, "sim1", "config.yaml"))
        manifest = JLD2.load(joinpath(output, "manifest.jld2"), "manifest")
        expected_mu_over_U = 0.22261343788323623
        @test config["model"]["hopping"]["J"] == [0.6, 0.6]
        @test config["model"]["onsite"]["mu"] ≈ fill(20 * expected_mu_over_U, 2)
        @test manifest[1, "lrbhqm.unit_density_data"] ==
            "alpha1p3_U12overU11_0p5.jld2"
        @test manifest[1, "lrbhqm.model_id"] == "alpha1p3_U12overU11_0p5"
        @test manifest[1, "lrbhqm.lookup_mode"] == "none"
        @test manifest[1, "lrbhqm.match_sampling_mode"] === true
        @test manifest[1, "lrbhqm.t_over_U"] ≈ 0.03
        @test manifest[1, "lrbhqm.mu_over_U"] ≈ expected_mu_over_U
        @test manifest[1, "lrbhqm.estimation_method"] == "sampled"
        @test_throws ArgumentError generate_study(source)

        definition["output_dir"] = "invalid_mu"
        definition["overrides"] = Dict("model.onsite.mu" => [1.0, 1.0])
        YAML.write_file(source, definition)
        @test_throws ArgumentError WormWeaver.LRBHQM.generate_study(source)
        @test !ispath(joinpath(root, "invalid_mu"))
    end
end

@testset "Bundled LRBHQM unit-density study" begin
    mktempdir() do root
        source = joinpath(root, "study.yaml")
        definition = YAML.load_file(joinpath(
            @__DIR__, "..", "templates", "studies", "lrbhqm-unit-density.yaml"))
        definition["output_dir"] = "generated"
        definition["append_date"] = false
        definition["replicas"] = 1
        delete!(definition, "multiwormqmc_path")
        YAML.write_file(source, definition)

        output = WormWeaver.LRBHQM.generate_study(source)
        manifest = JLD2.load(joinpath(output, "manifest.jld2"), "manifest")
        @test DataFrames.nrow(manifest) == 4
        @test manifest[!, "model.hopping.J"] ==
            [[0.55, 0.55], [0.60, 0.60], [0.70, 0.70], [0.80, 0.80]]
        @test manifest[!, "lrbhqm.estimation_method"] == fill("sampled", 4)
        @test manifest[!, "lrbhqm.t_over_U"] ≈ [0.0275, 0.03, 0.035, 0.04]
    end
end
