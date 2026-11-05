module WormWeaver

import Printf

include("g2_entry_modes.jl")

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
include("aggregation/static_structure_factor.jl")
include("aggregation/counterflow_green.jl")
include("aggregation/observable_io.jl")
include("aggregation/aggregate.jl")
include("aggregation/observable_table.jl")
include("aggregation/describe_observables.jl")

export generate_study, init_study, init_config, aggregate, save_aggregate, observable_table
export g2_entry_modes
export describe_density, describe_energy, describe_windings, describe_single_particle_green
export describe_particles_histogram
export green_ratio_statistics, study_status, study_log_status, status_summary, status_manifest
export status_summary_data
export slurm_array_tasks, find_scratch_logs, read_scratch_log

end
