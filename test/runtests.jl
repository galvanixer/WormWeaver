using Test
using WormWeaver

include("g2_entry_modes.jl")

include("study.jl")

include("init_study.jl")
include("init_config.jl")
include("aggregation.jl")
include("density_correlation.jl")
include("static_structure_factor.jl")
include("counterflow_green.jl")
include("describe_density.jl")
include("describe_energy_windings.jl")
include("describe_single_particle_green.jl")
include("describe_particles_histogram.jl")
include("green_statistics.jl")
include("single_particle_green.jl")
include("partial_aggregation.jl")
include("study_status.jl")
include("status_summary.jl")
include("hpc.jl")
