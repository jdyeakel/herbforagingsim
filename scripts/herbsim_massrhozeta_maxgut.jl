using Revise

using herbforagingsim

using DataFrames
using Distributions
using LinearAlgebra
using Base.Threads

using Plots
using StatsPlots
using UnicodePlots
using ColorSchemes
using JLD2
using ProgressMeter
using GLM


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

rhoexpvec = collect(-8.0:0.01:-5.7);
l_rhoexpvec = length(rhoexpvec);
zetavec = [1.0,1.75,2.0]; #collect(1:0.5:2);
l_zetavec = length(zetavec);

reps = 500;

gains = Array{Float64}(undef,reps,l_massvec,l_rhoexpvec,l_zetavec);
gaindiff = Array{Float64}(undef,reps,l_massvec,l_rhoexpvec,l_zetavec);
gainremainder = Array{Float64}(undef,reps,l_massvec,l_rhoexpvec,l_zetavec);
costs = Array{Float64}(undef,reps,l_massvec,l_rhoexpvec,l_zetavec);
encounters = Array{Float64}(undef,reps,l_massvec,l_rhoexpvec,l_zetavec);
cons_maxgut = Array{Float64}(undef,l_massvec);

@showprogress 1 "Computing..." for r=1:reps
    @threads for m=1:l_massvec
        for i=1:l_rhoexpvec
            for j=1:l_zetavec

                mass = massvec[m];
                rhoexp = rhoexpvec[i];
                zeta = zetavec[j];

                rho = 10^rhoexp;


                beta = bite_size_allo(mass); # mass in kg, bite size in g/bite
                chewrate = chew_allo(mass, teeth); #g/s
                t_chewgram = 1/chewrate; #s/g
                tchew = t_chewgram * beta; #s/g * g/bite = s/bite
                maxgut = gut_capacity_g(mass, gut_type) #grams
                bcost_kJps, fcost_kJps = metabolic_cost(mass);
                velocity = find_velocity(mass); # m/s
                tmax_bout, _ = foragingtime(mass) .* (60*60) #hours * 60 min/hour * 60 sec/min = secs
                #Consumer population density: individuals/m^2
                n = indperarea(mass); #individuals/m^2


                # rho * mu * (1/beta) = bite encounters/m^2 ~ bite encounter rate
                m_res = rho*mu*(1/beta);

                #Adjusted for competitive landscape
                mprime = m_res/n;
                alphaprime = alpha*n^(zeta-2);

                configurations = 200000;
                gammadist = Gamma(alphaprime, mprime / alphaprime);

                gains_daily, costs_daily, encounters_daily = dailyforage(gammadist,tchew, beta, maxgut, velocity, bcost_kJps, fcost_kJps, edensity, tmax_bout)


                gains[r,m,i,j] = gains_daily
                costs[r,m,i,j] = costs_daily
                encounters[r,m,i,j] = encounters_daily
                cons_maxgut[m] = maxgut*edensity

                gaindiff[r,m,i,j] = gains[r,m,i,j] - cons_maxgut[m]

                #Calculate remainder gain - cost 
                gainremainder[r,m,i,j] = gains[r,m,i,j] - costs[r,m,i,j]

                # #Run within-day and across day sims
                # gains_inds, 
                # costs_inds, 
                # gut_inds, 
                # fat_inds, 
                # fatsynth_inds = herbsim_individuals(mass,teeth,gut_type,rhoexp,mu,alpha,edensity,zeta,configurations,p_bad,runs)
                # # relfat_inds = fat_inds/maxfatstorage(mass,37.7)[1];

                # #Calculate metrics
                # survival[m,i,j] = sum(gut_inds[:,end] .> 0) / runs;
                # fatmean[m,i,j] = mean(mean(gut_inds[:,end-100:end],dims=2));
                # fatCV[m,i,j] = mean(std(gut_inds[:,end-100:end],dims=2) ./ mean(gut_inds[:,end-100:end],dims=2));

            end
        end
    end
end

# UnicodePlots.lineplot(rhoexpvec,gaindiff[4,:,1])


#Extract rhomin
rhomin = Array{Float64}(undef,reps,l_massvec,l_zetavec);
for r=1:reps
    for m=1:l_massvec
        mass_gain = gaindiff[r,m,:,:];
        for j=1:l_zetavec
            rhomin[r,m,j] = findrhomin(rhoexpvec,mass_gain[:,j],0.0);
        end
    end
end

rhomin_mean = mean(rhomin, dims=1)

rhomin_mean = dropdims(rhomin_mean, dims=1)

# PLOT fat reserves for individuals across time
line_width=2;
# Sample colors at evenly spaced intervals for more distinct differences
num_colors = min(10, l_zetavec)  # Adjust the number of colors for distinctiveness
color_indices = range(0, stop=1, length=num_colors)
color_scheme = [ColorSchemes.roma[c] for c in color_indices]  # Choose distinct colors

p = Plots.plot(massvec,rhomin_mean[:,1],
            xlabel="Mass (kg)", 
            ylabel="Gutmax rho",
            xlims=(10^1.5,10^4.5),
            xscale=:log10,
            yscale=:log10,
            legend=false,
            color = color_scheme[1], 
            linewidth=line_width,
            framestyle=:box);
for i=2:l_zetavec
    Plots.plot!(p, massvec,rhomin_mean[:,i],
    xscale=:log10,
    yscale=:log10,
    color = color_scheme[i], 
    linewidth=line_width);
end; p


mean_gainremainder = vec(mean(gainremainder,dims=1)[:,:,40,1])
positive_gainremainder = findall(x->x>0,mean_gainremainder)
# plot(log.(massvec[positive_gainremainder])[20:37],log.(mean_gainremainder[positive_gainremainder])[20:37])


# Extract the log-transformed values
x = log.(massvec[positive_gainremainder])[end-17:end]
y = log.(mean_gainremainder[positive_gainremainder])[end-17:end]
# Create a DataFrame for the linear model
df = DataFrame(x = x, y = y)
# Fit a linear model: y ~ x
model = lm(@formula(y ~ x), df)
# Print the model summary to see the slope
println(coef(model))

intercept, slope = coef(model)
# Generate the fitted line
fitted_y = intercept .+ slope .* x;
plot(log.(massvec[positive_gainremainder]),log.(mean_gainremainder[positive_gainremainder]))
plot!(x, fitted_y, label="Fitted Line", linewidth=2, color=:red)




# Evaluate the SLOPE across Rho
remainder_rho_slope = Array{Float64}(undef,l_rhoexpvec)
minmass = Array{Float64}(undef,l_rhoexpvec)
for i=1:l_rhoexpvec
    mean_gainremainder = vec(mean(gainremainder,dims=1)[:,:,i,1])
    positive_gainremainder = findall(x->x>0,mean_gainremainder)
    if length(positive_gainremainder) > 0
        minmass[i] = massvec[positive_gainremainder][1]
    else
        minmass[i] = NaN
    end
    # plot(log.(massvec[positive_gainremainder])[20:37],log.(mean_gainremainder[positive_gainremainder])[20:37])
    # Extract the log-transformed values
    if length(positive_gainremainder)>17
        x = log.(massvec[positive_gainremainder])[end-17:end]
        y = log.(mean_gainremainder[positive_gainremainder])[end-17:end]
        # Create a DataFrame for the linear model
        df = DataFrame(x = x, y = y)
        # Fit a linear model: y ~ x
        model = lm(@formula(y ~ x), df)
        # Print the model summary to see the slope
        # println(coef(model))
        intercept, slope = coef(model)

        remainder_rho_slope[i] = slope
    else
        remainder_rho_slope[i] = NaN
    end
end



pslope = plot(rhoexpvec,remainder_rho_slope,
    xlabel="Richness",
    ylabel="Energetic remainder exponent",
    ylims=(0.8,1.3),
    width=2,
    frame=:box,
    legend=false)
#Define a Too-Low-Rho
#Where remainder < 0 for all masses 
minrho = rhoexpvec[findall(x->x<100,minmass)][1]
minrho_yvalues = collect(0.8:0.1:1.3)
minrho_xvalues = repeat([minrho], outer=length(minrho_yvalues))
# Add the line to the plot
Plots.plot!(pslope, minrho_xvalues, minrho_yvalues,
    color=:orange,
    width=2,
    linestyle = :dash)
minrho = rhoexpvec[findall(x->x<500,minmass)][1]
minrho_yvalues = collect(0.8:0.1:1.3)
minrho_xvalues = repeat([minrho], outer=length(minrho_yvalues))
# Add the line to the plot
Plots.plot!(pslope, minrho_xvalues, minrho_yvalues,
    color=:red,
    width=2,
    linestyle = :dash)


#SAME THING BUT ACROSS GUT TYPES

# Define your gut types
gut_types = ["caecum", "colon", "non-rumen foregut", "rumen foregut"]
n_gut = length(gut_types)

# Your pre-existing parameter definitions
massexpvec = collect(1.5:0.05:4.4);
massvec = 10 .^ massexpvec;
l_massvec = length(massvec);

teeth = "all";  # remains constant here

# RESOURCE parameters
mu = 1;
alpha = 3;
edensity = 18.2;
p_bad = 0.05;
configurations = 20000;
runs = 200;

rhoexpvec = collect(-8.0:0.01:-5.7);
l_rhoexpvec = length(rhoexpvec);
zetavec = [1.0, 1.75, 2.0];
l_zetavec = length(zetavec);

reps = 500

# Pre-allocate an array to store slopes for each gut type and each rhoexp value:
remainder_rho_slope_all = Array{Float64}(undef, n_gut, l_rhoexpvec);

gains_array = Array{Array{Float64}}(undef,n_gut);
costs_array = Array{Array{Float64}}(undef,n_gut);
encounters_array = Array{Array{Float64}}(undef,n_gut);
cons_maxgut_array = Array{Array{Float64}}(undef,n_gut);
gaindiff_array = Array{Array{Float64}}(undef,n_gut);
gainremainder_array = Array{Array{Float64}}(undef,n_gut);

# Loop over each gut type
for g in 1:n_gut
    current_gut_type = gut_types[g]
    println("Processing gut type: $current_gut_type")
    
    # Pre-allocate arrays for this gut type
    gains = Array{Float64}(undef, reps, l_massvec, l_rhoexpvec, l_zetavec)
    costs = Array{Float64}(undef, reps, l_massvec, l_rhoexpvec, l_zetavec)
    encounters = Array{Float64}(undef, reps, l_massvec, l_rhoexpvec, l_zetavec)
    cons_maxgut = Array{Float64}(undef, l_massvec)
    gaindiff = Array{Float64}(undef, reps, l_massvec, l_rhoexpvec, l_zetavec)
    gainremainder = Array{Float64}(undef, reps, l_massvec, l_rhoexpvec, l_zetavec)
    
    # Run simulation over reps, mass, rhoexp, and zeta
    @showprogress 1 "Computing for gut type $current_gut_type ..." for r in 1:reps
        @threads for m in 1:l_massvec
            for i in 1:l_rhoexpvec
                for j in 1:l_zetavec
                    mass = massvec[m]
                    rhoexp = rhoexpvec[i]
                    zeta = zetavec[j]
                    rho = 10^rhoexp
                    
                    # Compute consumer and resource properties
                    beta = bite_size_allo(mass)
                    chewrate = chew_allo(mass, teeth)
                    t_chewgram = 1 / chewrate
                    tchew = t_chewgram * beta
                    maxgut = gut_capacity_g(mass, current_gut_type)
                    bcost_kJps, fcost_kJps = metabolic_cost(mass)
                    velocity = find_velocity(mass)
                    tmax_bout, _ = foragingtime(mass) .* (60 * 60)
                    n = indperarea(mass)
                    
                    m_res = rho * mu * (1 / beta)
                    mprime = m_res / n
                    alphaprime = alpha * n^(zeta - 2)
                    
                    configurations = 200000  # note: adjust if needed
                    gammadist = Gamma(alphaprime, mprime / alphaprime)
                    
                    # Compute daily foraging outcomes
                    gains_daily, costs_daily, encounters_daily = dailyforage(gammadist, tchew, beta, maxgut, velocity, bcost_kJps, fcost_kJps, edensity, tmax_bout)
                    
                    gains[r,m,i,j] = gains_daily
                    costs[r,m,i,j] = costs_daily
                    encounters[r,m,i,j] = encounters_daily
                    cons_maxgut[m] = maxgut*edensity

                    gaindiff[r,m,i,j] = gains[r,m,i,j] - cons_maxgut[m]

                    #Calculate remainder gain - cost 
                    gainremainder[r,m,i,j] = gains[r,m,i,j] - costs[r,m,i,j]
                end
            end
        end
    end

    gains_array[g] = gains
    costs_array[g] = costs
    encounters_array[g] = encounters
    cons_maxgut_array[g] = cons_maxgut
    gaindiff_array[g] = gaindiff
    gainremainder_array[g] = gainremainder
    
    # # Calculate slope for each rhoexp value
    # remainder_rho_slope = Array{Float64}(undef, l_rhoexpvec)
    # minmass = Array{Float64}(undef, l_rhoexpvec)
    # for i in 1:l_rhoexpvec
    #     # Take mean across reps (and using the first zeta value)
    #     mean_gainremainder = vec(mean(gainremainder, dims=1)[:, :, i, 1])
    #     positive_gainremainder = findall(x -> x > 0, mean_gainremainder)
    #     if !isempty(positive_gainremainder)
    #         minmass[i] = massvec[positive_gainremainder][1]
    #     else
    #         minmass[i] = NaN
    #     end
        
    #     # Only perform the linear fit if there are enough positive points
    #     if length(positive_gainremainder) > 17
    #         x = log.(massvec[positive_gainremainder])[end-17:end]
    #         y = log.(mean_gainremainder[positive_gainremainder])[end-17:end]
    #         df = DataFrame(x = x, y = y)
    #         model = lm(@formula(y ~ x), df)
    #         # Extract coefficients: intercept and slope
    #         intercept, slope = coef(model)
    #         remainder_rho_slope[i] = slope
    #     else
    #         remainder_rho_slope[i] = NaN
    #     end
    # end
    
    # # Save the slope curve for the current gut type
    # remainder_rho_slope_all[g, :] = remainder_rho_slope
end


#EXAMINE SLOPE EXTRACTION
remainder_rho_slope_all = Array{Float64}(undef, n_gut, l_rhoexpvec);

for g in 1:n_gut
    gains = gains_array[g]
    costs = costs_array[g]
    encounters = encounters_array[g]
    cons_maxgut = cons_maxgut_array[g]
    gaindiff = gaindiff_array[g]
    gainremainder = gainremainder_array[g]

    
    # Calculate slope for each rhoexp value
    remainder_rho_slope = Array{Float64}(undef, l_rhoexpvec)
    minmass = Array{Float64}(undef, l_rhoexpvec)
    for i in 1:l_rhoexpvec
        # Take mean across reps (and using the first zeta value)
        mean_gainremainder = vec(mean(gainremainder, dims=1)[:, :, i, 1])
        positive_gainremainder = findall(x -> x > 0, mean_gainremainder)
        if !isempty(positive_gainremainder)
            minmass[i] = massvec[positive_gainremainder][1]
        else
            minmass[i] = NaN
        end
        
        # Only perform the linear fit if there are enough positive points
        num_est = 20
        if length(positive_gainremainder) > num_est
            x = log.(massvec[positive_gainremainder])[end-num_est:end]
            y = log.(mean_gainremainder[positive_gainremainder])[end-num_est:end]
            df = DataFrame(x = x, y = y)
            model = lm(@formula(y ~ x), df)
            # Extract coefficients: intercept and slope
            intercept, slope = coef(model)
            remainder_rho_slope[i] = slope
        else
            remainder_rho_slope[i] = NaN
        end
    end
    # Save the slope curve for the current gut type
    remainder_rho_slope_all[g, :] = remainder_rho_slope
end

filename=smartpath("data/simdata/gainsremainder_acrossgut.jld2")
@save filename gut_types n_gut massexpvec massvec l_massvec teeth mu alpha edensity p_bad configurations runs rhoexpvec l_rhoexpvec zetavec l_zetavec reps gains_array costs_array encounters_array cons_maxgut_array gaindiff_array gainremainder_array remainder_rho_slope_all minmass 

@load filename gut_types n_gut massexpvec massvec l_massvec teeth mu alpha edensity p_bad configurations runs rhoexpvec l_rhoexpvec zetavec l_zetavec reps gains_array costs_array encounters_array cons_maxgut_array gaindiff_array gainremainder_array remainder_rho_slope_all minmass 


gut_types_cap = ["Caecum","Colon","Non-rumen foregut","Rumen foregut"]
gut_type_maxgutslope = [0.860,0.919,0.881,0.897];
using LaTeXStrings
pslope = plot(
    xlabel = L"Richness, $\rho$",
    ylabel = L"Forage excess exponent, $g_R \propto M^\gamma$",
    ylims = (0.8, 1.3),
    legend = :topright,
    foreground_color_legend = nothing
              )
colors = palette(:auto, n_gut)  # Get exactly n_gut colors
for g in 1:n_gut
    col = colors[g]
    gutslope = gut_type_maxgutslope[g]
    x = log.(massvec)
    y = log.(vec(mean(costs_array[g],dims=1)[:,:,end,1]))
    df = DataFrame(x = x, y = y)
    model = lm(@formula(y ~ x), df)
    intc, _ = coef(model)
    x = log.(massvec)
    y = log.(vec(mean(gains_array[g],dims=1)[:,:,end,1]))
    df = DataFrame(x = x, y = y)
    model = lm(@formula(y ~ x), df)
    intg, _ = coef(model)
    int = exp(intg)/exp(intc)
    exp_slope = 0.75 + ((int*(gutslope-0.75))*massvec[end]^(gutslope-0.75))/(int*massvec[end]^(gutslope-0.75) - 1)
    plot!(rhoexpvec, remainder_rho_slope_all[g, :], label=gut_types_cap[g], lw=2,color=col)
    scatter!(pslope, [rhoexpvec[end]], [exp_slope],markersize=5,legend=false,color=col)
end
minrho = rhoexpvec[findall(x->x<100,minmass)][1]
minrho_yvalues = collect(0.8:0.1:1.3)
minrho_xvalues = repeat([minrho], outer=length(minrho_yvalues))
# Add the line to the plot
Plots.plot!(pslope, minrho_xvalues, minrho_yvalues,
    color=:orange,
    width=2,
    linestyle = :dash)
minrho = rhoexpvec[findall(x->x<500,minmass)][1]
minrho_yvalues = collect(0.8:0.1:1.3)
minrho_xvalues = repeat([minrho], outer=length(minrho_yvalues))
# Add the line to the plot
Plots.plot!(pslope, minrho_xvalues, minrho_yvalues,
    color=:red,
    width=2,
    linestyle = :dash)
display(pslope)
Plots.savefig(pslope, string(homedir(),"/Dropbox/PostDoc/2024_herbforaging/figures/fig_gainsremainder_slope.pdf"))




guti = 4
gainsi = vec(mean(gains_array[guti],dims=1)[:,:,end,1])
costsi = vec(mean(costs_array[guti],dims=1)[:,:,end,1])
plot(log.(massvec),log.(gainsi))
plot!(log.(massvec),log.(edensity.*gut_capacity_g.(massvec, gut_types[guti])))

plot(log.(massvec),log.(gainsi-costsi))
plot!(log.(massvec),log.(edensity.*gut_capacity_g.(massvec, gut_types[guti])))

plot(log.(massvec),log.(gainsi-costsi))
x = log.(massvec)
y = log.(costsi)
df = DataFrame(x = x, y = y)
model = lm(@formula(y ~ x), df)
println(coef(model))


gut_capacity_g.(massvec, gut_types[g])














gut_types_cap = ["Caecum","Colon","Non-rumen foregut","Rumen foregut"]
gut_type_maxgutslope = [0.860,0.919,0.881,0.897];
using LaTeXStrings
pslope = plot(
    xlabel = L"Richness, $\rho$",
    ylabel = L"Forage excess exponent, $g_R \propto M^\gamma$",
    ylims = (0.8, 1.3),
    legend = :topright,
    foreground_color_legend = nothing
              )
colors = palette(:auto, n_gut)  # Get exactly n_gut colors
for g in 1:n_gut
    col = colors[g]
    plot!(rhoexpvec, remainder_rho_slope_all[g, :], label=gut_types_cap[g], lw=2,color=col)
    scatter!(pslope, [rhoexpvec[end]], [gut_type_maxgutslope[g]],markersize=5,legend=false,color=col)
end
minrho = rhoexpvec[findall(x->x<100,minmass)][1]
minrho_yvalues = collect(0.8:0.1:1.3)
minrho_xvalues = repeat([minrho], outer=length(minrho_yvalues))
# Add the line to the plot
Plots.plot!(pslope, minrho_xvalues, minrho_yvalues,
    color=:orange,
    width=2,
    linestyle = :dash)
minrho = rhoexpvec[findall(x->x<300,minmass)][1]
minrho_yvalues = collect(0.8:0.1:1.3)
minrho_xvalues = repeat([minrho], outer=length(minrho_yvalues))
# Add the line to the plot
Plots.plot!(pslope, minrho_xvalues, minrho_yvalues,
    color=:red,
    width=2,
    linestyle = :dash)
display(pslope)












# scatter!(pslope, [rhoexpvec[1]], [remainder_rho_slope[1]],
#     markersize=5, markercolor=:black)
# Plots.annotate!(pslope, rhoexpvec[1] + 0.08, remainder_rho_slope[1] + 0.01,
#     text(string(round(remainder_rho_slope[1];digits=3)), :black, 10))





Plots.savefig(p, string(homedir(),"/Dropbox/PostDoc/2024_herbforaging/figures/fig_massrhozeta_maxgutgains.pdf"))

# dataexport = [massvec rhomin];
# Create a DataFrame
df = DataFrame(
    mass = massvec,
    zeta1 = rhomin_mean[:, 1],
    zeta175 = rhomin_mean[:, 2],
    zeta2 = rhomin_mean[:, 3]
)
filename=smartpath("data/simdata/rhovalues_foregut_maxgut.csv")
CSV.write(filename, df)

# Find the function corresponding to the fit.
zetaindex = 1
break_mass_surv = Array{Float64}(undef,2);
break_slope_surv = Array{Float64}(undef,2);

break_mass, index_break = breakpoint_find(massvec[1:50],rhomin[1:50,zetaindex])
a1, b1, a2, b2, mass_range1, rhomin_fit1, mass_range2, rhomin_fit2 = breakpoint_fit(massvec[1:50],rhomin[1:50,:],break_mass,index_break,zetaindex)

error_all, error_two_piece, F_value, p_value = breakpoint_compare_models(massvec[1:50], rhomin[1:50,zetaindex], index_break)

break_mass_surv[1] = break_mass;
break_slope_surv[1] = b1;

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



for i=3
    zetaindex = i

    break_mass, index_break = breakpoint_find(massvec[1:50],rhomin[1:50,zetaindex])
    a1, b1, a2, b2, mass_range1, rhomin_fit1, mass_range2, rhomin_fit2 = breakpoint_fit(massvec[1:50],rhomin[1:50,:],break_mass,index_break,zetaindex)

    error_all, error_two_piece, F_value, p_value = breakpoint_compare_models(massvec[1:50], rhomin[1:50,zetaindex], index_break)

    break_mass_surv[i-1] = break_mass;
    break_slope_surv[i-1] = b1;

    # Original data
    Plots.scatter!(pb,massvec, rhomin[:,zetaindex], label="Data", xlabel="Mass", ylabel="Rhomin", xscale=:log10, yscale=:log10)
    # Fitted function
    plot!(pb,mass_range1, rhomin_fit1, label="Fit1", lw=2, xscale=:log10, yscale=:log10)
    plot!(pb,mass_range2, rhomin_fit2, label="Fit1", lw=2, xscale=:log10, yscale=:log10)
end
pb

Plots.savefig(pb, string(homedir(),"/Dropbox/PostDoc/2024_herbforaging/figures/fig_massrhozeta_breakpoints_maxgutgains.pdf"))











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




# filename = smartpath("data/simdata/massrhozeta_maxgut5_split.jld2")
# # @save filename survival fatmean fatCV
# vars_to_save = Dict(n => getfield(Main, n) for n in setdiff(names(Main), [:Base, :Core, :Main]))
# @save filename vars_to_save
