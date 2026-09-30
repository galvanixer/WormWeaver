module WormWeaver

include("study.jl")
include("init_study.jl")
using .Study: generate_study

export generate_study, init_study

end
