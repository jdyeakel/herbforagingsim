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


#HERBIVORE
#Define mass of herbivore
mass = 100;
#Define tooth and gut type of herbivore
teeth = "bunodont"; # bunodont, acute/obtuse lophs, lophs and non-flat, lophs and flat
gut_type = "caecum"; # caecum, colon, non-rumen foregut, rumen foregut

#RESOURCE
#Set richness
# rho = 1*10^-9; #This is very small, because we have set mu = 1
# rhoexp = -6.19;
# rhoexp = -7.07;

#Define resource traits
mu = 1;
alpha = 3
edensity = 18.2;
# zeta = 1.5;
p_bad = 0.05;
configurations = 20000; #works fine - 100000 for more perfect distributions
runs = 100;

rhoexpvec = collect(-7.2:0.05:-6.8);
l_rhoexpvec = length(rhoexpvec);
zetavec = collect(1:0.25:2);
l_zetavec = length(zetavec);

survival = Array{Float64}(undef,l_rhoexpvec,l_zetavec);
fatmean = Array{Float64}(undef,l_rhoexpvec,l_zetavec);
fatCV = Array{Float64}(undef,l_rhoexpvec,l_zetavec);

@showprogress 1 "Computing..." for i=1:l_rhoexpvec
    for j=1:l_zetavec

        rhoexp = rhoexpvec[i];
        zeta = zetavec[j]

        #Run within-day and across day sims
        gains_inds, 
        costs_inds, 
        gut_inds, 
        fat_inds, 
        fatsynth_inds = herbsim_individuals(mass,teeth,gut_type,rhoexp,mu,alpha,edensity,zeta,configurations,p_bad,runs);
        relfat_inds = fat_inds/maxfatstorage(mass,37.7)[1];

        #Calculate metrics
        survival[i,j] = sum(gut_inds[:,end] .> 0) / runs;
        fatmean[i,j] = mean(mean(gut_inds[:,end-100:end],dims=2));
        fatCV[i,j] = mean(std(gut_inds[:,end-100:end],dims=2) ./ mean(gut_inds[:,end-100:end],dims=2));

    end
    
end

# path = joinpath(homedir(), "Dropbox/PostDoc/2024_Anna_bse/herbforaging/herbforagingsim/data/simdata/rhozeta.jld2")
# @save path survival fatmean fatCV

filename = smartpath("data/simdata/rhozeta.jld2")
vars_to_save = Dict(n => getfield(Main, n) for n in setdiff(names(Main), [:Base, :Core, :Main]))
@save filename vars_to_save


UnicodePlots.heatmap(survival)

UnicodePlots.heatmap((fatCV))

UnicodePlots.heatmap((fatmean))


#Extract rhomin
rhomin = Array{Float64}(undef,l_zetavec)
for j=1:l_zetavec
    rhomin[j] = findrhomin(rhoexpvec,survival[:,j]);
end
UnicodePlots.lineplot(zetavec,rhomin)