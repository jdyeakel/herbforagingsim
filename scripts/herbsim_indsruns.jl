using Revise

using herbforagingsim

using DataFrames
using Distributions
using LinearAlgebra
using Base.Threads

using UnicodePlots
using JLD2
using StatsPlots



#HERBIVORE
#Define mass of herbivore
mass = 30;
#Define tooth and gut type of herbivore
teeth = "bunodont"; # bunodont, acute/obtuse lophs, lophs and non-flat, lophs and flat
gut_type = "caecum"; # caecum, colon, non-rumen foregut, rumen foregut

#RESOURCE
#Set richness
# rho = 1*10^-9; #This is very small, because we have set mu = 1
# rhoexp = -6.19;
rhoexp = -6.7;

#Define resource traits
mu = 1;
alpha = 3
edensity = 18.2;
zeta = 1.5;
p_bad = 0.05;
configurations = 20000; #works fine - 100000 for more perfect distributions
runs = 50;

#Run within-day and across day sims
gains_inds, 
costs_inds, 
gut_inds, 
fat_inds, 
fatsynth_inds = herbsim_individuals(mass,teeth,gut_type,rhoexp,mu,alpha,edensity,zeta,configurations,p_bad,runs);
relfat_inds = fat_inds/maxfatstorage(mass,37.7)[1];

# PLOT fat reserves for individuals across time
alpha = 0.6;
color = :lightblue;
p = Plots.plot(collect(1:length(relfat_inds[1,:]))./365,relfat_inds[1,:],
            xlabel="Time (years)", 
            ylabel="Relative Fat Reserves",
            legend=false,
            color = color, 
            alpha = alpha);
for i=2:runs
    Plots.plot!(p, collect(1:length(relfat_inds[i,:]))./365,relfat_inds[i,:],
    color = color, 
    alpha = alpha);
end
p

