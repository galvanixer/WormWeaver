import DataFrames

@testset "DataFrame column size report" begin
    table = DataFrames.DataFrame(
        small=fill(Int8(1), 100_000), large=collect(Int64, 1:100_000))
    output = IOBuffer()

    @test describe_column_sizes(table; io=output) === nothing
    lines = split(chomp(String(take!(output))), '\n')
    @test length(lines) == 2
    @test startswith(lines[1], "large")
    @test startswith(lines[2], "small")
    @test all(endswith(" KiB"), lines)

    mib_output = sprint(io -> describe_column_sizes(table; io, unit=:MiB))
    @test all(endswith(" MiB"), split(chomp(mib_output), '\n'))
    @test_throws ArgumentError describe_column_sizes(table; io=devnull, unit=:MB)

    @test sprint(io -> describe_column_sizes(DataFrames.DataFrame(); io)) == ""
end
