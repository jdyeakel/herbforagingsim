module herbforagingsim

using DataFrames
using Distributions
using LinearAlgebra
using Base.Threads

include("trait_and_rate_functions.jl")
include("dailyforage.jl")
include("withindaysim_v2.jl")
include("acrossdaysim.jl")
include("calculate_histogram.jl")
include("find_bin_index.jl")
include("gc_sample.jl")
include("decompose_to_integer_components.jl")
include("herbsim_individuals.jl")
include("findrhomin.jl")
include("smartpath.jl")

export 

plant_digestibility,
apparent_digestibility,
generate_trait_array,
find_velocity,
generate_rate_array,
bite_size_allo,
bite_rate_allo,
alpha_allo,
number_of_chews,
chew_rate_allo,
chew_allo,
mean_retention_time,
max_retention_time,
min_retention_time,
gut_capacity_g,
mean_particle_mass,
outflow_rate,
maxfatstorage,
metabolic_cost,
indperarea,
fatrecoveryrate,
expectedlifetime,
foragingtime,
gut_turnover_time,
dailyfoodintake,

dailyforage,
withindaysim,
acrossdaysim,
herbsim_individuals,

calculate_histogram,
find_bin_index,
gc_sample,
decompose_to_integer_components,
findrhomin,
smartpath



end # module herbforagingsim
