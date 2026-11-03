import DataFrames

@testset "Independent dashboard counts" begin
    df = DataFrames.DataFrame(lattice_size=Any[], beta=Float64[], status=String[],
                              failed=Bool[], log_source=String[])
    # Each tuple describes parameter group, stage, failure flag, source, count.
    rows = [([64], 64.0, "Finished", false, "local", 5),
            ([64], 64.0, "Finished", false, "scratch", 20),
            ([64], 64.0, "Checkpoint 3", false, "scratch", 5),
            ([128], 128.0, "Finished", false, "local", 5),
            ([128], 128.0, "Finished", false, "scratch", 15),
            ([128], 128.0, "Production Started", false, "scratch", 10),
            ([256], 128.0, "Finished", false, "local", 5),
            ([256], 128.0, "Finished", false, "scratch", 5),
            ([256], 128.0, "Tuning", true, "scratch", 5),
            ([256], 128.0, "Checkpoint 1", false, "scratch", 15),
            ([256], 256.0, "Finished", false, "local", 5),
            ([256], 256.0, "Thermalizing", true, "scratch", 10),
            ([256], 256.0, "No Log File", false, "none", 5),
            ([256], 256.0, "Checkpoint 2", false, "scratch", 10)]
    for (lattice, beta, stage, failed, source, count) in rows, _ in 1:count
        push!(df, (lattice, beta, stage, failed, source))
    end
    original = deepcopy(df)
    expected = DataFrames.combine(DataFrames.groupby(df,
        [:lattice_size, :beta, :status, :failed, :log_source]; sort=false),
        DataFrames.nrow => :n_simulations)
    data = status_summary_data(df)
    @test data.overall == (; total=120, finished=60, failed=15, remaining=45)
    @test isequal(data.summary, expected)
    @test data.has_source
    @test data.parameters.lattice_size == [[64], [128], [256], [256]]
    @test data.parameters.beta == [64.0, 128.0, 128.0, 256.0]
    @test data.parameters.total == fill(30, 4)
    @test data.parameters.finished == [25, 20, 10, 5]
    @test data.parameters.failed == [0, 0, 5, 10]
    @test data.parameters.remaining == [5, 10, 15, 15]
    @test sum(data.parameters.total) == sum(data.stages.n_simulations) == 120
    @test sum(data.stages.n_local) == 20
    @test sum(data.stages.n_scratch) == 95
    @test sum(data.stages.n_none) == 5
    @test sum(data.stages.n_unspecified) == 0
    @test data.stages.n_local .+ data.stages.n_scratch .+ data.stages.n_none .+
          data.stages.n_unspecified == data.stages.n_simulations
    @test isequal(df, original)

    io = IOBuffer()
    @test isequal(status_summary(df; io), expected)
    text = String(take!(io))
    @test occursin(r"Simulations\s+120", text)
    @test occursin(r"Finished\s+60\s+50\.0%", text)
    @test occursin(r"Failed\s+15\s+12\.5%", text)
    @test occursin(r"Remaining\s+45\s+37\.5%", text)
    @test occursin(r"TOTAL\s+120\s+60\s+15\s+45", text)
    @test occursin("PROGRESS BY PARAMETERS", text)
    @test occursin("PROGRESS BY STAGE", text)
    @test occursin(r"Local\s+Scratch\s+No log\s+Total", text)
    @test occursin("Failed (Tuning)", text)
    @test occursin("Failed (Thermalizing)", text)
    @test occursin("Log-derived progress; not Slurm job state.", text)
    @test isequal(status_summary(df; io, detailed=false), expected)
    compact = String(take!(io))
    @test occursin("PROGRESS BY PARAMETERS", compact)
    @test !occursin("PROGRESS BY STAGE", compact)
    @test !occursin('\e', text)

    # Color must never change the layout or numerical output.
    status_summary(df; io=IOContext(io, :color => true))
    colored = String(take!(io))
    @test occursin("\e[38;2;27;158;119m", colored)
    @test occursin("\e[38;2;217;95;2m", colored)
    @test replace(colored, r"\e\[[0-9;]*m" => "") == text
    status_summary(df; io=IOContext(io, :color => false))
    @test String(take!(io)) == text
    mktemp() do _, output
        redirect_stdout(output) do
            status_summary(df; io=devnull)
            status_summary_data(df)
        end
        flush(output)
        seekstart(output)
        @test read(output, String) == ""
    end

    legacy = DataFrames.select(df, DataFrames.Not(:log_source))
    legacy_data = status_summary_data(legacy)
    @test !legacy_data.has_source
    @test legacy_data.overall == data.overall
    @test sum(legacy_data.stages.n_unspecified) == 120
    @test !hasproperty(legacy_data.summary, :log_source)
    status_summary(legacy; io)
    legacy_text = String(take!(io))
    @test !occursin("Local", legacy_text)
    @test !occursin("Scratch", legacy_text)
    @test !occursin("Unspecified", legacy_text)
    @test occursin("PROGRESS BY STAGE", legacy_text)
    @test_throws ArgumentError status_summary_data(DataFrames.select(df, DataFrames.Not(:failed)))
end

@testset "Stage ordering and adaptive parameter display" begin
    known = ["No Log File", "Pre Tuning", "Tuning", "Thermalizing", "Production Started",
             "Checkpoint 1", "Checkpoint 2", "Checkpoint 3", "Checkpoint 4", "Finished"]
    # Reverse input order to ensure rendering order is independent of the
    # backwards-compatible first-appearance order of the returned DataFrame.
    labels = [reverse(known); "Tuning"; "Finished"; "Custom first"; "Custom second"]
    df = DataFrames.DataFrame(lattice_size=fill([16, 16], length(labels)),
        beta=fill(64.0, length(labels)), status=labels,
        failed=[falses(10); true; true; false; false],
        log_source=fill("scratch", length(labels)))
    original = deepcopy(df)
    data = status_summary_data(df)
    @test data.overall == (; total=14, finished=1, failed=2, remaining=11)
    @test data.summary.status[1] == "Finished"
    io = IOBuffer()
    summary = status_summary(df; io)
    text = String(take!(io))
    @test !occursin("PROGRESS BY PARAMETERS", text)
    @test length(findall("[16, 16]", text)) == 1
    @test occursin("Lattice [16, 16]  β 64.0", text)
    @test occursin(r"Stage\s+Local\s+Scratch\s+Total", text)
    @test !occursin("No log", text)
    @test !occursin("Unspecified", text)
    expected_order = [known[1:3]; "Failed (Tuning)"; known[4:end];
                      "Failed (Finished)"; "Custom first"; "Custom second"]
    lines = split(text, '\n')
    stage_lines = lines[findfirst(==("PROGRESS BY STAGE"), lines)+1:end]
    positions = [findfirst(line -> startswith(line, "$stage  "), stage_lines) for stage in expected_order]
    @test all(!isnothing, positions)
    @test issorted(positions)
    @test isequal(summary, data.summary)
    @test isequal(df, original)
    status_summary(df; io, detailed=false)
    @test !occursin("PROGRESS BY STAGE", String(take!(io)))

    # Missing parameter metadata and unknown source labels remain countable.
    mixed = DataFrames.DataFrame(lattice_size=Any[[8, 8, 8], missing, [8, 8, 8], [8, 8, 8]],
        beta=Any[32.0, missing, 32.0, 32.0],
        status=["Tuning", "No Log File", "Tuning", "Tuning"], failed=falses(4),
        log_source=Union{Missing,String}["local", "none", missing, "archive"])
    mixed_data = status_summary_data(mixed)
    @test mixed_data.parameters.total == [3, 1]
    @test ismissing(mixed_data.parameters.lattice_size[2])
    @test mixed_data.stages.n_unspecified == [2, 0]
    status_summary(mixed; io)
    mixed_text = String(take!(io))
    @test occursin("[8, 8, 8]", mixed_text)
    @test occursin("—", mixed_text)
    @test occursin(r"Local\s+Scratch\s+No log\s+Unspecified\s+Total", mixed_text)

    empty = mixed[1:0, :]
    empty_data = status_summary_data(empty)
    @test empty_data.overall == (; total=0, finished=0, failed=0, remaining=0)
    @test DataFrames.nrow(empty_data.parameters) == DataFrames.nrow(empty_data.stages) == 0
    result = status_summary(empty; io)
    empty_text = String(take!(io))
    @test startswith(empty_text, "STUDY STATUS\n")
    @test length(collect(eachmatch(r"0\.0%", empty_text))) == 3
    @test !occursin("NaN", empty_text)
    @test !occursin("PROGRESS BY", empty_text)
    @test isequal(result, empty_data.summary)
end
