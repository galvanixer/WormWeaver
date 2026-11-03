import DataFrames
import CSV
import JLD2

@testset "Slurm array locations" begin
    mktempdir() do root
        fixture = joinpath(root, "queue.txt")
        arguments = joinpath(root, "arguments.txt")
        executable = joinpath(root, "squeue")
        write(executable, raw"""
        #!/bin/sh
        printf '%s\n' "$@" > "$SQUEUE_ARGUMENTS"
        if [ "$SQUEUE_FAIL" = "1" ]; then
            exit 1
        fi
        cat "$SQUEUE_FIXTURE"
        """)
        chmod(executable, 0o755)
        write(fixture, "18523211_3|18523211|RUNNING|hpc-n937\n" *
                       "18523211_2|18523213|RUNNING|hpc-n664\n" *
                       "18523211_10|18523214|PENDING|\n")

        withenv("PATH" => root * ":" * ENV["PATH"],
                "SQUEUE_FIXTURE" => fixture, "SQUEUE_ARGUMENTS" => arguments,
                "SQUEUE_FAIL" => "0") do
            table = slurm_array_tasks(18523211)
            @test table.index == [2, 3, 10]
            @test table.jobid == [18523213, 18523211, 18523214]
            @test table.task == ["18523211_2", "18523211_3", "18523211_10"]
            @test table.nodes[1:2] == ["hpc-n664", "hpc-n937"]
            @test ismissing(table.nodes[3])
            @test table.state == ["RUNNING", "RUNNING", "PENDING"]
            @test table.scratch_dir == ["/scratch/job.18523213",
                                       "/scratch/job.18523211",
                                       "/scratch/job.18523214"]

            sibling = slurm_array_tasks("18523211_2"; scratch_root="/custom/scratch")
            @test sibling.index == table.index
            @test "18523211" in readlines(arguments)
            @test !("18523211_2" in readlines(arguments))
            @test sibling.scratch_dir[1] == "/custom/scratch/job.18523213"
            @test_throws ArgumentError slurm_array_tasks("18523211; echo invalid")

            write(fixture, "")
            empty = slurm_array_tasks(18523211)
            @test DataFrames.nrow(empty) == 0
            @test names(empty) == names(table)

            withenv("SQUEUE_FAIL" => "1") do
                @test_throws ProcessFailedException slurm_array_tasks(18523211)
            end
        end
    end
end

@testset "Find manifest logs on scratch" begin
    mktempdir() do root
        bin, study = joinpath(root, "bin"), joinpath(root, "study")
        scratch = joinpath(root, "scratch 'quoted' \$(touch injected);")
        mkpath.([bin, study, scratch])
        queue, calls = joinpath(root, "queue"), joinpath(root, "calls")
        write(calls, "")
        write(queue, "18523211_3|18523211|RUNNING|hpc-n937\n" *
                     "18523211_2|18523213|RUNNING|hpc-n664\n" *
                     "18523211_1|18523214|PENDING|\n")
        scripts = Dict(
            "squeue" => raw"""
            #!/bin/bash
            printf 'squeue\n' >> "$HPC_CALLS"
            cat "$HPC_QUEUE"
            """,
            "scontrol" => raw"""
            #!/bin/bash
            printf 'scontrol\n' >> "$HPC_CALLS"
            printf 'hpc-n664\nhpc-n937\n'
            """,
            "srun" => raw"""
            #!/bin/bash
            jobid= node=
            while [[ "$1" == --* ]]; do
                case "$1" in
                    --jobid=*) jobid=${1#*=} ;;
                    --nodelist=*) node=${1#*=} ;;
                esac
                shift
            done
            printf 'srun|%s|%s\n' "$jobid" "$node" >> "$HPC_CALLS"
            [[ "$node" != "${HPC_FAIL_NODE:-}" ]] || exit 1
            exec "$@"
            """)
        for (name, script) in scripts
            path = joinpath(bin, name)
            write(path, script)
            chmod(path, 0o755)
        end
        quoted = "sim 'quoted' \$(touch injected);"
        manifest = DataFrames.DataFrame(simulation=["sim41", "sim1", "sim21", quoted, "sim41"],
                                        config=fill("unused/config.yaml", 5))
        native = joinpath(study, "manifest.jld2")
        JLD2.jldsave(native; manifest)
        # Only simulation labels determine scratch paths; config paths are not used.
        CSV.write(joinpath(study, "manifest.csv"), DataFrames.DataFrame(simulation=["wrong"]))
        paths = [joinpath(scratch, "job.$jobid", sim, "simulation.log")
                 for (jobid, sim) in ((18523211, "sim41"), (18523213, "sim21"),
                                      (18523213, quoted), (18523213, "sim41"))]
        for path in paths
            mkpath(dirname(path))
            write(path, "Do not parse this log.\n")
        end
        original = read.(paths)
        withenv("PATH" => bin * ":" * ENV["PATH"], "HPC_QUEUE" => queue,
                "HPC_CALLS" => calls, "HPC_FAIL_NODE" => "") do
            table = find_scratch_logs(study, 18523211; scratch_root=scratch)
            @test names(table) == ["simulation", "task", "jobid", "node", "log_path"]
            @test table.simulation == ["sim41", "sim41", "sim21", quoted]
            @test table.jobid == [18523213, 18523211, 18523213, 18523213]
            @test table.node == ["hpc-n664", "hpc-n937", "hpc-n664", "hpc-n664"]
            @test table.task == ["18523211_2", "18523211_3", "18523211_2", "18523211_2"]
            @test table.log_path == paths[[4, 1, 2, 3]]
            @test count(startswith("srun|"), readlines(calls)) == 2
            @test read.(paths) == original
            @test !ispath(joinpath(root, "injected"))
            @test isequal(find_scratch_logs(study, "18523211_2"; scratch_root=scratch), table)

            # Read complete text through the actual cat invocation. The quoted
            # path must remain one argument and text must be returned unchanged.
            content = "MultiWormQMC v0.2\nPRODUCTION\n  25/100 25% 1s ◆\nno final newline"
            write(paths[3], content)
            @test read_scratch_log(table[4, :]) == content
            @test study_log_status(IOBuffer(read_scratch_log(table[4, :]))) ==
                (; status="Checkpoint 1", completed_blocks=25, failed=false)
            @test read(paths[3], String) == content
            @test read_scratch_log(NamedTuple(table[4, :])) == content
            write(paths[3], "")
            @test read_scratch_log(table[4, :]) == ""
            @test !ispath(joinpath(root, "injected"))
            absent = (; jobid=18523213, node="hpc-n664", log_path=joinpath(root, "missing.log"))
            @test_throws ProcessFailedException read_scratch_log(absent)
            @test_throws ArgumentError read_scratch_log((; jobid=18523213))
            @test_throws ArgumentError read_scratch_log(merge(absent, (; jobid=0)))
            @test_throws ArgumentError read_scratch_log(merge(absent, (; node=missing)))
            @test_throws ArgumentError read_scratch_log(merge(absent, (; log_path="relative.log")))

            # Existence does not require opening/reading log contents.
            chmod(paths[1], 0o000)
            try
                @test isequal(find_scratch_logs(study, 18523211; scratch_root=scratch), table)
            finally
                chmod(paths[1], 0o600)
            end
            withenv("HPC_FAIL_NODE" => "hpc-n664") do
                @test_throws ProcessFailedException find_scratch_logs(study, 18523211; scratch_root=scratch)
                @test_throws ProcessFailedException read_scratch_log(table[1, :])
            end
            @test_throws ArgumentError find_scratch_logs(study, 18523211; scratch_root="relative")
            for bad in ("../sim1", "sim1/subdir", "", ".", "..", missing)
                JLD2.jldsave(native; manifest=DataFrames.DataFrame(simulation=[bad]))
                @test_throws ArgumentError find_scratch_logs(study, 18523211; scratch_root=scratch)
            end
            JLD2.jldsave(native; manifest=DataFrames.DataFrame(config=["sim1/config.yaml"]))
            @test_throws ArgumentError find_scratch_logs(study, 18523211; scratch_root=scratch)

            # CSV fallback and a node range both exercise the real Bash check.
            rm(native)
            CSV.write(joinpath(study, "manifest.csv"), manifest)
            write(queue, "18523211_3|18523211|RUNNING|hpc-n[664,937]\n")
            expanded = find_scratch_logs(study, 18523211; scratch_root=scratch)
            @test expanded.simulation == ["sim41", "sim41"]
            @test expanded.node == ["hpc-n664", "hpc-n937"]
            @test "scontrol" in readlines(calls)
            write(queue, "18523211_1|18523214|PENDING|\n")
            before = count(startswith("srun|"), readlines(calls))
            empty = find_scratch_logs(study, 18523211; scratch_root=scratch)
            @test DataFrames.nrow(empty) == 0
            @test names(empty) == names(table)
            @test count(startswith("srun|"), readlines(calls)) == before
        end
    end
end
