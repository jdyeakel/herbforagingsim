#// Description

using Revise

using herbforagingsim

using DataFrames
using Distributions
using LinearAlgebra
using Base.Threads

using UnicodePlots
using JLD2
using StatsPlots
using ProgressMeter

# Load variables from saved file
filename = smartpath("data/simdata/massrhozeta.jld2")
@load filename vars_to_save
for (name, value) in vars_to_save
    @eval Main $(Symbol(name)) = $value
end




#Extract rhomin
rhomin = Array{Float64}(undef,l_zetavec)

for j=1:l_zetavec
    rhomin[j] = findrhomin(rhoexpvec,survival[:,j]);
end
