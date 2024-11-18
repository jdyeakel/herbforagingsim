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
massexpvec = collect(1.5:0.05:4.4);
massvec = 10 .^ massexpvec;
l_massvec = length(massvec);
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
runs = 200;

rhoexpvec = collect(-7.5:0.01:-6.7);
l_rhoexpvec = length(rhoexpvec);
zetavec = [1.0,1.75,2.0]; #collect(1:0.5:2);
l_zetavec = length(zetavec);

survival = Array{Float64}(undef,l_massvec,l_rhoexpvec,l_zetavec);
fatmean = Array{Float64}(undef,l_massvec,l_rhoexpvec,l_zetavec);
fatCV = Array{Float64}(undef,l_massvec,l_rhoexpvec,l_zetavec);

@showprogress 1 "Computing..." for m=1:l_massvec
    for i=1:l_rhoexpvec
        for j=1:l_zetavec

            mass = massvec[m];
            rhoexp = rhoexpvec[i];
            zeta = zetavec[j]

            #Run within-day and across day sims
            gains_inds, 
            costs_inds, 
            gut_inds, 
            fat_inds, 
            fatsynth_inds = herbsim_individuals(mass,teeth,gut_type,rhoexp,mu,alpha,edensity,zeta,configurations,p_bad,runs)
            # relfat_inds = fat_inds/maxfatstorage(mass,37.7)[1];

            #Calculate metrics
            survival[m,i,j] = sum(gut_inds[:,end] .> 0) / runs;
            fatmean[m,i,j] = mean(mean(gut_inds[:,end-100:end],dims=2));
            fatCV[m,i,j] = mean(std(gut_inds[:,end-100:end],dims=2) ./ mean(gut_inds[:,end-100:end],dims=2));

        end
    end
end

filename = smartpath("data/simdata/massrhozeta_maxgut3.jld2")
# @save filename survival fatmean fatCV
vars_to_save = Dict(n => getfield(Main, n) for n in setdiff(names(Main), [:Base, :Core, :Main]))
@save filename vars_to_save
