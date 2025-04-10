module herbforagingsim

using DataFrames
using Distributions
using LinearAlgebra
using Base.Threads

include("trait_and_rate_functions.jl")
include("dailyforage.jl")
include("withindaysim_v2.jl")
include("withindaysim_v3.jl")
include("acrossdaysim.jl")
include("calculate_histogram.jl")
include("find_bin_index.jl")
include("gc_sample.jl")
include("decompose_to_integer_components.jl")
include("herbsim_individuals.jl")
include("findrhomin.jl")
include("smartpath.jl")
include("dailyforage_expectedvalues.jl")
include("breakpoint_find.jl")
include("breakpoint_compare_models.jl")
include("breakpoint_fit.jl")
include("dailyforage_experiment.jl")
include("dailyforage_trim.jl")

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
reactionwidth,
reactionheight,

dailyforage,
withindaysim,
withindaysim_split,
acrossdaysim,
herbsim_individuals,

calculate_histogram,
find_bin_index,
gc_sample,
decompose_to_integer_components,
findrhomin,
smartpath,
dailyforage_expectedvalues,
breakpoint_find,
breakpoint_compare_models,
breakpoint_fit,
dailyforage_experiment,
dailyforage_trim



end # module herbforagingsim
