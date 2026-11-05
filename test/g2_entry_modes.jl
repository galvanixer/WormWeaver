import YAML

@testset "Configured G2 entry modes" begin
    templates = joinpath(@__DIR__, "..", "templates", "multiwormqmc", "LRBHQM")
    for (variant, expected) in (("singleworm", ()), ("regular", (:regular,)),
                                ("aligned", (:aligned,)), ("anti_aligned", (:anti_aligned,)))
        path = joinpath(templates, "config_LRBHQM_$variant.yaml")
        bytes = read(path)
        config = YAML.load_file(path)
        original = deepcopy(config)
        @test g2_entry_modes(path) === expected
        @test g2_entry_modes(config) === expected
        @test config == original
        @test read(path) == bytes
    end
    config = YAML.load_file(joinpath(templates, "config_LRBHQM_regular.yaml"))
    request(kind; weight=0.2) = Dict("kind" => kind, "target" => [1, 1], "weight" => weight)

    # Both G1 source tables may independently offer multiple types of entry.
    mixed = deepcopy(config)
    mixed["sampler"]["moves"][2]["requests"] = [
        request("OpenColocatedAntiAligned"), request("Open"), request("OpenColocatedAligned")]
    @test g2_entry_modes(mixed) === (:regular, :aligned, :anti_aligned)
    generic = deepcopy(config)
    for table in generic["sampler"]["moves"][2:3]
        table["requests"] = [request("OpenColocated")]
    end
    @test g2_entry_modes(generic) === (:aligned, :anti_aligned)

    disabled = deepcopy(config)
    for table in disabled["sampler"]["moves"][2:3], item in table["requests"]
        sum(item["target"]) == 2 && (item["weight"] = 0.0)
    end
    @test g2_entry_modes(disabled) === () # Order two alone does not enable G2.
    disconnected = deepcopy(disabled)
    disconnected["sampler"]["moves"][1]["requests"][2]["weight"] = 0.0
    push!(disconnected["sampler"]["moves"][3]["requests"], request("OpenColocatedAligned"))
    @test g2_entry_modes(disconnected) === ()
    # Reachability is independent of the ordering of source tables.
    reverse!(mixed["sampler"]["moves"])
    @test g2_entry_modes(mixed) === (:regular, :aligned, :anti_aligned)

    invalid = Any[nothing, [], Dict(), Dict("sampler" => nothing)]
    for (key, value) in (("moves", nothing), ("moves", []), ("multiworm_order", 0),
                         ("multiworm_order", true), ("multiworm_order", 1))
        bad = deepcopy(config)
        bad["sampler"][key] = value
        push!(invalid, bad)
    end
    for (key, value) in (("weight", -1.0), ("weight", Inf), ("weight", NaN),
                         ("kind", "UnknownOpen"), ("target", [1]), ("target", [-1, 1]),
                         ("target", [true, false]), ("target", [1.0, 1.0]), ("target", [0, 0]))
        bad = deepcopy(config)
        bad["sampler"]["moves"][2]["requests"][end][key] = value
        push!(invalid, bad)
    end
    missing_table = deepcopy(config)
    pop!(missing_table["sampler"]["moves"])
    push!(invalid, missing_table)
    for bad in invalid
        result = @test_logs (:warn, r"Cannot interpret configuration") g2_entry_modes(bad)
        @test ismissing(result)
    end
    mktempdir() do root
        path = joinpath(root, "config.yaml")
        @test_logs (:warn, r"Cannot read configuration") @test ismissing(g2_entry_modes(path))
        write(path, "sampler: [unterminated")
        @test_logs (:warn, r"Cannot read configuration") @test ismissing(g2_entry_modes(path))
        # A scalar YAML root is not another config path to follow.
        YAML.write_file(path, path)
        @test_logs (:warn, r"YAML root must be a mapping") @test ismissing(g2_entry_modes(path))
        YAML.write_file(path, Dict("sampler" => Dict("moves" => [])))
        @test_logs (:warn, r"Cannot interpret configuration") @test ismissing(g2_entry_modes(path))
    end
end
