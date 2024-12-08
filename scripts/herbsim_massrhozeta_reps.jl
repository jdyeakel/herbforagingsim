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

#Saving as individual files, so we could do more
reps = collect(1:10);

#HERBIVORE
#Define mass of herbivore
massexpvec = collect(1.5:0.05:4.4);
massvec = 10 .^ massexpvec;
l_massvec = length(massvec);
#Define tooth and gut type of herbivore
teeth = "all"; # all, bunodont, acute/obtuse lophs, lophs and non-flat, lophs and flat
gut_type = "rumen foregut"; # caecum, colon, non-rumen foregut, rumen foregut

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
zetavec = [1.0,2.0]; #collect(1:0.5:2);
l_zetavec = length(zetavec);

@showprogress 1 "Computing..." for r=reps

    survival = Array{Float64}(undef, l_massvec, l_rhoexpvec, l_zetavec)
    fatmean = Array{Float64}(undef, l_massvec, l_rhoexpvec, l_zetavec)
    fatCV = Array{Float64}(undef, l_massvec, l_rhoexpvec, l_zetavec)

    # Compute total iterations
    total_iterations = l_massvec * l_rhoexpvec * l_zetavec

    # Parallel loop over combined indices
    @threads for index in 1:total_iterations
        # Compute m, i, j from the flattened index
        m = div(index - 1, l_rhoexpvec * l_zetavec) + 1
        i = div(mod(index - 1, l_rhoexpvec * l_zetavec), l_zetavec) + 1
        j = mod(index - 1, l_zetavec) + 1

        mass = massvec[m]
        rhoexp = rhoexpvec[i]
        zeta = zetavec[j]

        # Run within-day and across-day sims
        gains_inds, 
        costs_inds, 
        gut_inds, 
        fat_inds, 
        fatsynth_inds = herbsim_individuals(mass, teeth, gut_type, rhoexp, mu, alpha, edensity, zeta, configurations, p_bad, runs)

        # Calculate metrics
        survival[m, i, j] = sum(gut_inds[:, end] .> 0) / runs
        fatmean[m, i, j] = mean(mean(gut_inds[:, end-100:end], dims=2))
        fatCV[m, i, j] = mean(std(gut_inds[:, end-100:end], dims=2) ./ mean(gut_inds[:, end-100:end], dims=2))
    end

    filename = smartpath("data/simdata/massrhozetareps_foregut/massrhozeta_rep.jld2", [r])
    # @save filename survival fatmean fatCV
    vars_to_save = Dict(n => getfield(Main, n) for n in setdiff(names(Main), [:Base, :Core, :Main]))
    @save filename vars_to_save

    # GC.gc()  # Force garbage collection after each rep
end

# #Validate index scheme NOTE: appears to work!

# # Generate all possible combinations of (m, i, j)
# all_combinations = Set{Tuple{Int, Int, Int}}([(m, i, j) for m in 1:l_massvec, i in 1:l_rhoexpvec, j in 1:l_zetavec])

# # Flattened index to (m, i, j) mapping and save combinations
# visited_combinations = Set{Tuple{Int, Int, Int}}()
# for index in 1:total_iterations
#     m = div(index - 1, l_rhoexpvec * l_zetavec) + 1
#     i = div(mod(index - 1, l_rhoexpvec * l_zetavec), l_zetavec) + 1
#     j = mod(index - 1, l_zetavec) + 1

#     # Save the combination
#     push!(visited_combinations, (m, i, j))
# end

# # Check if all combinations are visited
# missing_combinations = setdiff(all_combinations, visited_combinations)
# if !isempty(missing_combinations)
#     error("Missing combinations: $missing_combinations")
# else
#     println("All combinations visited successfully!")
# end

# @showprogress 1 "Computing..." for r=reps

#     survival = Array{Float64}(undef,l_massvec,l_rhoexpvec,l_zetavec);
#     fatmean = Array{Float64}(undef,l_massvec,l_rhoexpvec,l_zetavec);
#     fatCV = Array{Float64}(undef,l_massvec,l_rhoexpvec,l_zetavec);

#     @threads for m=1:l_massvec
#         for i=1:l_rhoexpvec
#             for j=1:l_zetavec

#                 mass = massvec[m];
#                 rhoexp = rhoexpvec[i];
#                 zeta = zetavec[j]

#                 #Run within-day and across day sims
#                 gains_inds, 
#                 costs_inds, 
#                 gut_inds, 
#                 fat_inds, 
#                 fatsynth_inds = herbsim_individuals(mass,teeth,gut_type,rhoexp,mu,alpha,edensity,zeta,configurations,p_bad,runs)
#                 # relfat_inds = fat_inds/maxfatstorage(mass,37.7)[1];

#                 #Calculate metrics
#                 survival[m,i,j] = sum(gut_inds[:,end] .> 0) / runs;
#                 fatmean[m,i,j] = mean(mean(gut_inds[:,end-100:end],dims=2));
#                 fatCV[m,i,j] = mean(std(gut_inds[:,end-100:end],dims=2) ./ mean(gut_inds[:,end-100:end],dims=2));

#             end
#         end
#     end


#     filename = smartpath("data/simdata/massrhozetareps_foregut/massrhozeta_rep.jld2",[r])
#     # @save filename survival fatmean fatCV
#     vars_to_save = Dict(n => getfield(Main, n) for n in setdiff(names(Main), [:Base, :Core, :Main]))
#     @save filename vars_to_save

#     # GC.gc()  # Force garbage collection after each rep

    
# end


