module WormWeaver

include("study_generation/study.jl")
include("study_generation/init_study.jl")
include("study_generation/init_config.jl")
using .Study: generate_study

export generate_study, init_study, init_config

end
