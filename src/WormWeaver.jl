module WormWeaver

include("study_generation/study.jl")
include("study_generation/init_study.jl")
include("study_generation/init_config.jl")
using .Study: generate_study

include("monitoring/study_status.jl")
include("monitoring/hpc.jl")
include("monitoring/status_summary.jl")
include("aggregation/statistics.jl")
include("aggregation/ratio_statistics.jl")
include("aggregation/particle_histogram.jl")
include("aggregation/windings.jl")
include("aggregation/single_particle_green.jl")
include("aggregation/observable_io.jl")
include("aggregation/aggregate.jl")
include("aggregation/observable_table.jl")

export generate_study, init_study, init_config, aggregate, save_aggregate, observable_table
export green_ratio_statistics, study_status, study_log_status, status_summary, status_manifest
export status_summary_data
export slurm_array_tasks, find_scratch_logs, read_scratch_log

end
