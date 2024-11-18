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
zetaindex = 3
break_mass, index_break = find_breakpoint(massvec[1:35],rhomin[1:35,zetaindex])
massvec1 = massvec[1:index_break]
rhomin1 = rhomin[1:index_break,zetaindex]
massvec2 = massvec[index_break+1:end]
rhomin2 = rhomin[index_break+1:end,zetaindex]

log_massvec1 = log.(massvec1)
log_rhomin1 = log.(rhomin1)
# Fit a line
X1 = hcat(ones(length(log_massvec1)), log_massvec1)
coeffs1 = X1 \ log_rhomin1  # Linear regression
log_a1, b1 = coeffs1
a1 = exp(log_a1)  # Convert log(a1) to a1

log_massvec2 = log.(massvec2)
log_rhomin2 = log.(rhomin2)
# Fit a line
X2 = hcat(ones(length(log_massvec2)), log_massvec2)
coeffs2 = X2 \ log_rhomin2  # Linear regression
log_a2, b2 = coeffs2
a2 = exp(log_a2)  # Convert log(a2) to a2

mass_range1 = range(minimum(massvec), stop=break_mass, length=100)
rhomin_fit1 = a1 .* mass_range1 .^ b1;
mass_range2 = range(break_mass, stop= maximum(massvec))
rhomin_fit2 = a2 .* mass_range2 .^ b2;

# Original data
Plots.scatter(massvec, rhomin[:,zetaindex], label="Data", xlabel="Mass", ylabel="Rhomin", xscale=:log10, yscale=:log10)

# Fitted function
plot!(mass_range1, rhomin_fit1, label="Fit1", lw=2, xscale=:log10, yscale=:log10)
plot!(mass_range2, rhomin_fit2, label="Fit1", lw=2, xscale=:log10, yscale=:log10)


error_all, error_two_piece, F_value, p_value = compare_breakpoint_models(massvec[1:35], rhomin[1:35,zetaindex], index_break)