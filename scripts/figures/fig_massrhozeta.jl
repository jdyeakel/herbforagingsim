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
filename = smartpath("data/simdata/massrhozeta.jld2")
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

p = Plots.plot(massvec,rhomin[:,1],
            xlabel="Mass (kg)", 
            ylabel="Minimum rho",
            xscale=:log10,
            yscale=:log10,
            legend=false,
            color = color_scheme[1], 
            linewidth=line_width);
for i=2:l_zetavec
    Plots.plot!(p, massvec,rhomin[:,i],
    xscale=:log10,
    yscale=:log10,
    color = color_scheme[i], 
    linewidth=line_width);
end; p
Plots.savefig(p, "/Users/jdyeakel/Dropbox/PostDoc/2024_herbforaging/figures/fig_massrhozeta.pdf")
