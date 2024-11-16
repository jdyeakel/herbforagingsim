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
using ProgressMeter

# Load variables from saved file
filename = smartpath("data/simdata/massrhozeta_maxgut2.jld2")
@load filename vars_to_save
for (name, value) in vars_to_save
    @eval Main $(Symbol(name)) = $value
end




#Extract rhomin
rhomin = Array{Float64}(undef,l_massvec,l_zetavec);
for m=1:l_massvec
    mass_survival = survival[m,:,:];
    for j=1:l_zetavec
        rhomin[m,j] = findrhomin(rhoexpvec,mass_survival[:,j]);
    end
end


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
Plots.savefig(p, "/Users/justinyeakel/Dropbox/PostDoc/2024_herbforaging/figures/fig_massrhozeta_maxgut2.pdf")


# Strength of selection = Delta rhomin / Delta mass
Delta_mass = diff(massvec)
Delta_rhomin = diff(rhomin[:,1])
Delta_rhomin[Delta_rhomin .== 0] .= NaN
selstrength = Delta_rhomin ./ Delta_mass
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
    selstrength = Delta_rhomin ./ Delta_mass
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

Plots.savefig(ps, "/Users/justinyeakel/Dropbox/PostDoc/2024_herbforaging/figures/fig_massrhozeta_strength_maxgut2.pdf")

