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
massexpvec = collect(1.5:0.05:2.5);
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
runs = 100;

rhoexpvec = collect(-7.25:0.025:-6.7);
l_rhoexpvec = length(rhoexpvec);
zetavec = [1.0,2.0]; #collect(1:0.5:2);
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





# filename = smartpath("data/simdata/massrhozeta_maxgut4_split.jld2")
# # @save filename survival fatmean fatCV
# vars_to_save = Dict(n => getfield(Main, n) for n in setdiff(names(Main), [:Base, :Core, :Main]))
# @save filename vars_to_save



#Extract rhomin
rhomin = Array{Float64}(undef,l_massvec,l_zetavec);
for m=1:l_massvec
    mass_survival = survival[m,:,:];
    for j=1:l_zetavec
        rhomin[m,j] = findrhomin(rhoexpvec,mass_survival[:,j],0.5);
    end
end


#What are the breakpoints (slope and mass) as a function of survival threshold?

surv_thresholdvec = collect(0:0.1:1);
l_surv_thresholdvec = length(surv_thresholdvec);
#Extract rhomin
rhomin_surv = Array{Float64}(undef,l_surv_thresholdvec,l_massvec,l_zetavec);

for i=1:l_surv_thresholdvec
    surv_threshold = surv_thresholdvec[i];
    for m=1:l_massvec
        mass_survival = survival[m,:,:];
        for j=1:l_zetavec
            rhomin_surv[i,m,j] = findrhomin(rhoexpvec,mass_survival[:,j],surv_threshold);
        end
    end
end

break_mass_surv = Array{Float64}(undef,l_surv_thresholdvec,l_zetavec)
break_slope_surv = Array{Float64}(undef,l_surv_thresholdvec,l_zetavec)
for j=1:l_zetavec
    for i=1:l_surv_thresholdvec
        zetaindex=j;
        break_mass, index_break = breakpoint_find(massvec,rhomin_surv[i,:,zetaindex])
        a1, b1, a2, b2, mass_range1, rhomin_fit1, mass_range2, rhomin_fit2 = breakpoint_fit(massvec,rhomin_surv[i,:,:],break_mass,index_break,zetaindex)
        error_all, error_two_piece, F_value, p_value = breakpoint_compare_models(massvec, rhomin_surv[i,:,zetaindex], index_break)

        break_mass_surv[i,j] = break_mass;
        break_slope_surv[i,j] = b1;
    end
end

mean(break_mass_surv,dims=1)

mean(break_slope_surv,dims=1)


# Find the function corresponding to the fit.
zetaindex = 1

break_mass, index_break = breakpoint_find(massvec,rhomin[:,zetaindex])
a1, b1, a2, b2, mass_range1, rhomin_fit1, mass_range2, rhomin_fit2 = breakpoint_fit(massvec,rhomin,break_mass,index_break,zetaindex)

error_all, error_two_piece, F_value, p_value = breakpoint_compare_models(massvec, rhomin[:,zetaindex], index_break)

# Original data
pb = Plots.scatter(massvec, rhomin[:,zetaindex], 
    label="Data", 
    xlabel="Mass", 
    ylabel="Rhomin", 
    xscale=:log10, 
    yscale=:log10,
    framestyle=:box)  
# Fitted function
plot!(pb,mass_range1, rhomin_fit1, label="Fit1", lw=2, xscale=:log10, yscale=:log10)
plot!(pb,mass_range2, rhomin_fit2, label="Fit1", lw=2, xscale=:log10, yscale=:log10)

for i=l_zetavec
    zetaindex = i

    break_mass, index_break = breakpoint_find(massvec,rhomin[:,zetaindex])
    a1, b1, a2, b2, mass_range1, rhomin_fit1, mass_range2, rhomin_fit2 = breakpoint_fit(massvec,rhomin,break_mass,index_break,zetaindex)

    error_all, error_two_piece, F_value, p_value = breakpoint_compare_models(massvec, rhomin[:,zetaindex], index_break)

    # Original data
    Plots.scatter!(pb,massvec, rhomin[:,zetaindex], label="Data", xlabel="Mass", ylabel="Rhomin", xscale=:log10, yscale=:log10)
    # Fitted function
    plot!(pb,mass_range1, rhomin_fit1, label="Fit1", lw=2, xscale=:log10, yscale=:log10)
    plot!(pb,mass_range2, rhomin_fit2, label="Fit1", lw=2, xscale=:log10, yscale=:log10)
end
pb

