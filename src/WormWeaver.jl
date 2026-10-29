module WormWeaver

include("study_generation/study.jl")
include("study_generation/init_study.jl")
include("study_generation/init_config.jl")
using .Study: generate_study

include("aggregation/statistics.jl")
include("aggregation/ratio_statistics.jl")
include("aggregation/particle_histogram.jl")
include("aggregation/windings.jl")
include("aggregation/single_particle_green.jl")
include("aggregation/observable_io.jl")
include("aggregation/aggregate.jl")
include("aggregation/observable_table.jl")

export generate_study, init_study, init_config, aggregate, observable_table
export green_ratio_statistics

end
