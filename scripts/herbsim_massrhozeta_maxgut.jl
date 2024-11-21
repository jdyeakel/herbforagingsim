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

reps = 100;

gains = Array{Float64}(undef,reps,l_massvec,l_rhoexpvec,l_zetavec);
gaindiff = Array{Float64}(undef,reps,l_massvec,l_rhoexpvec,l_zetavec);
costs = Array{Float64}(undef,reps,l_massvec,l_rhoexpvec,l_zetavec);
encounters = Array{Float64}(undef,reps,l_massvec,l_rhoexpvec,l_zetavec);
cons_maxgut = Array{Float64}(undef,l_massvec);

@showprogress 1 "Computing..." for r=1:reps
    for m=1:l_massvec
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

UnicodePlots.lineplot(rhoexpvec,gaindiff[4,:,1])

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

p = Plots.scatter(massvec,rhomin[:,1],
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
    Plots.scatter!(p, massvec,rhomin[:,i],
    xscale=:log10,
    yscale=:log10,
    color = color_scheme[i], 
    linewidth=line_width);
end; p
Plots.savefig(p, string(homedir(),"/Dropbox/PostDoc/2024_herbforaging/figures/fig_massrhozeta_maxgutgains.pdf"))

# dataexport = [massvec rhomin];
# Create a DataFrame
df = DataFrame(
    mass = massvec,
    zeta1 = rhomin_mean[:, 1],
    zeta175 = rhomin_mean[:, 2],
    zeta2 = rhomin_mean[:, 3]
)
filename=smartpath("data/simdata/rhovalues_maxgut.csv")
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
