#// Description

using Revise

using herbforagingsim

using DataFrames
using Distributions
using LinearAlgebra
using Base.Threads

using Plots
using UnicodePlots
using StatsPlots
using ColorSchemes
using JLD2
using CSV
using ProgressMeter

# Load variables from saved file
# NOTE: Both of these seem to work the same
filename = smartpath("data/simdata/massrhozeta_maxgut3.jld2")
filename = smartpath("data/simdata/massrhozeta_maxgut4_split.jld2")
# filename = smartpath("data/simdata/massrhozetareps/massrhozeta_rep_9.jld2")

@load filename vars_to_save
for (name, value) in vars_to_save
    @eval Main $(Symbol(name)) = $value
end
println("Loaded variables: ", collect(keys(vars_to_save)))


#Extract rhomin
rhomin = Array{Float64}(undef,l_massvec,l_zetavec);
for m=1:l_massvec
    mass_survival = survival[m,:,:];
    for j=1:l_zetavec
        rhomin[m,j] = findrhomin(rhoexpvec,mass_survival[:,j],0.8);
    end
end

filename=smartpath("data/simdata/rhovalues_maxgut.csv")
maxgutdata = CSV.read(filename, DataFrame)
mass_maxgut = Array(maxgutdata[!,:mass])
rhomin_maxgut = Array(maxgutdata[!,[:zeta1,:zeta175,:zeta2]])



# PLOT fat reserves for individuals across time
line_width=2;
# Sample colors at evenly spaced intervals for more distinct differences
num_colors = min(10, l_zetavec)  # Adjust the number of colors for distinctiveness
color_indices = range(0, stop=1, length=num_colors)
color_scheme = [ColorSchemes.roma[c] for c in color_indices]  # Choose distinct colors

p = Plots.scatter(massvec,rhomin[:,1],
            xlabel="Mass (kg)", 
            ylabel="Minimum rho",
            xlims=(10^1.5,10^4.5),
            xscale=:log10,
            yscale=:log10,
            legend=false,
            color = color_scheme[1], 
            linewidth=line_width,
            framestyle=:box);
for i=2:l_zetavec
    Plots.scatter!(p, massvec,rhomin[:,i],
    xscale=:log10,
    yscale=:log10,
    color = color_scheme[i], 
    linewidth=line_width);
end; p
Plots.savefig(p, string(homedir(),"/Dropbox/PostDoc/2024_herbforaging/figures/fig_massrhozeta_maxgut2.pdf"))


# Strength of selection = Delta rhomin / Delta mass
Delta_mass = diff(massvec)
Delta_rhomin = diff(rhomin[:,1])
Delta_rhomin[Delta_rhomin .== 0] .= NaN
selstrength = log.(abs.(Delta_rhomin)) ./ log.(Delta_mass)
pos_noNaN = findall(x -> !isnan(x), selstrength)

ps = Plots.scatter(massvec[2:end][pos_noNaN],abs.(Delta_rhomin[pos_noNaN]),
    xlabel="Mass (kg)", 
    ylabel="Strength of Selection",
    xlims=(10^1.5,10^4.5),
    xscale=:log10,
    yscale=:log10,
    legend=false,
    color = color_scheme[1], 
    linewidth=line_width,
    framestyle=:box  # Ensures all four sides are framed
    );
Plots.plot!(ps,massvec[2:end][pos_noNaN],abs.(Delta_rhomin[pos_noNaN]),
    linewidth=line_width,
    color = color_scheme[1])

for j=2:l_zetavec
    Delta_rhomin = diff(rhomin[:,j])
    Delta_rhomin[Delta_rhomin .== 0] .= NaN
    selstrength = log.(abs.(Delta_rhomin)) ./ log.(Delta_mass)
    pos_noNaN = findall(x -> !isnan(x), selstrength)
    
    Plots.scatter!(ps,massvec[2:end][pos_noNaN],abs.(Delta_rhomin[pos_noNaN]),
    xscale=:log10,
    yscale=:log10,
    color = color_scheme[j], 
    linewidth=line_width);
    
    Plots.plot!(ps,massvec[2:end][pos_noNaN],abs.(Delta_rhomin[pos_noNaN]),
    linewidth=line_width,
    color = color_scheme[j])
end
ps

Plots.savefig(ps, string(homedir(),"/Dropbox/PostDoc/2024_herbforaging/figures/fig_massrhozeta_strength_maxgut2.pdf"))



# Find the function corresponding to the fit.
zetaindex = 1

break_mass, index_break = breakpoint_find(massvec[1:50],rhomin[1:50,zetaindex])
a1, b1, a2, b2, mass_range1, rhomin_fit1, mass_range2, rhomin_fit2 = breakpoint_fit(massvec[1:50],rhomin[1:50,:],break_mass,index_break,zetaindex)

error_all, error_two_piece, F_value, p_value = breakpoint_compare_models(massvec[1:50], rhomin[1:50,zetaindex], index_break)

# Original data
pb = Plots.scatter(massvec, rhomin[:,zetaindex], 
    label="Data", 
    xlabel="Mass", 
    ylabel="Rhomin", 
    xscale=:log10, 
    yscale=:log10,
    framestyle=:box,
    legend=false)  
# Fitted function
plot!(pb,mass_range1, rhomin_fit1, label="Fit1", lw=2, xscale=:log10, yscale=:log10)
plot!(pb,mass_range2, rhomin_fit2, label="Fit1", lw=2, xscale=:log10, yscale=:log10)
plot!(pb,mass_maxgut,rhomin_maxgut[:,1],color=:black)

for i=3
    zetaindex = i

    break_mass, index_break = breakpoint_find(massvec[1:50],rhomin[1:50,zetaindex])
    a1, b1, a2, b2, mass_range1, rhomin_fit1, mass_range2, rhomin_fit2 = breakpoint_fit(massvec[1:50],rhomin[1:50,:],break_mass,index_break,zetaindex)

    error_all, error_two_piece, F_value, p_value = breakpoint_compare_models(massvec[1:50], rhomin[1:50,zetaindex], index_break)

    # Original data
    Plots.scatter!(pb,massvec, rhomin[:,zetaindex], label="Data", xlabel="Mass", ylabel="Rhomin", xscale=:log10, yscale=:log10)
    # Fitted function
    plot!(pb,mass_range1, rhomin_fit1, label="Fit1", lw=2, xscale=:log10, yscale=:log10)
    plot!(pb,mass_range2, rhomin_fit2, label="Fit1", lw=2, xscale=:log10, yscale=:log10)
    plot!(pb,mass_maxgut,rhomin_maxgut[:,i],color=:black)
end
pb

Plots.savefig(pb, string(homedir(),"/Dropbox/PostDoc/2024_herbforaging/figures/fig_massrhozeta_breakpoints_maxgut4_split.pdf"))


#What are the breakpoints (slope and mass) as a function of survival threshold?

surv_thresholdvec = collect(0.1:0.1:0.9);
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

break_mass_surv = Array{Float64}(undef,l_surv_thresholdvec,l_zetavec);
break_slope_surv = Array{Float64}(undef,l_surv_thresholdvec,l_zetavec);
for j=1:3
    for i=1:l_surv_thresholdvec
        zetaindex=j;
        break_mass, index_break = breakpoint_find(massvec[1:50],rhomin_surv[i,1:50,zetaindex])
        a1, b1, a2, b2, mass_range1, rhomin_fit1, mass_range2, rhomin_fit2 = breakpoint_fit(massvec[1:50],rhomin_surv[i,1:50,:],break_mass,index_break,zetaindex)
        error_all, error_two_piece, F_value, p_value = breakpoint_compare_models(massvec[1:50], rhomin_surv[i,1:50,zetaindex], index_break)

        break_mass_surv[i,j] = break_mass;
        break_slope_surv[i,j] = b1;
    end
end

mean(break_mass_surv,dims=1)
std(break_mass_surv,dims=1)




#Examination of mean mass values pre and post Eocene-Oligocene Transition

filename = smartpath("data/Alroy_bodysize1998.csv");
massdata = CSV.read(filename, DataFrame) 
# Convert all entries to lowercase
massdata[!, :"Limb morphology"] .= lowercase.(massdata[!, :"Limb morphology"])
UnicodePlots.scatterplot(-massdata[!,:"Last appearance (mya)"],log.(massdata[!,:"Mass (ln g)"]))
# Subset the DataFrame
subsetdata = filter(row -> row[:"Limb morphology"] in ["graviportal", "protounguligrade","unguligrade"], massdata)

fa_mya = subsetdata[!,:"First appearance (mya)"];
la_mya = subsetdata[!,:"Last appearance (mya)"];
sp_mass = (exp.(subsetdata[!,:"Mass (ln g)"]))/1000;

pre_eot_pos = findall(x -> x < 56  && x > 34 ,fa_mya);
post_eot_pos = findall(x -> x < 34 && x > 10 ,fa_mya);

exp(mean(log.(sp_mass[pre_eot_pos])))
exp(mean(log.(sp_mass[post_eot_pos])))

"graviportal", "protounguligrade", "unguligrade","digitigrade","protodigitigrade"