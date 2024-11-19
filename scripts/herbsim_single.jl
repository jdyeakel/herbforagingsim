using Revise

using herbforagingsim

using DataFrames
using Distributions
using LinearAlgebra

using Plots
using UnicodePlots
using StatsPlots
using JLD2
using CSV



#HERBIVORE
#Define mass of herbivore
mass = 10^1.5;
#Define tooth and gut type of herbivore
teeth = "bunodont"; # bunodont, acute/obtuse lophs, lophs and non-flat, lophs and flat
gut_type = "caecum"; # caecum, colon, non-rumen foregut, rumen foregut

#NOTE: Need to think more about this! See Owen-Smith 1988
# tmax_bout = 2*60*60; # Set at 1/2 day (6) hours (43200 seconds)
tmax_bout, tmax_bout_upper95 = foragingtime(mass) .* (60*60) #hours * 60 min/hour * 60 sec/min = secs

#RESOURCE
#Set richness
# rho = 1*10^-9; #This is very small, because we have set mu = 1
# rhoexp = -6.19;
rhoexp = -7.07#-7.07; #-6.7 #
rho = 1*10^rhoexp

#Define resource traits
res_traits = (mu = 1, alpha = 3, edensity = 18.2);
zeta = 2.;
mu = res_traits[:mu];
alpha = res_traits[:alpha];
edensity = res_traits[:edensity]; 

configurations = 20000; #works fine - 100000 for more perfect distributions

# This provides a return distribution
@time gains, costs, prob = withindaysim(rho,alpha,mu,zeta,edensity,mass,teeth,gut_type,tmax_bout,configurations);

# @time gains, costs, prob = withindaysim_split(rho,alpha,mu,zeta,edensity,mass,teeth,gut_type,tmax_bout,configurations);


#Confirmed - these are the gains
UnicodePlots.lineplot(gains,sum(prob,dims=1)')
#Confirmed - these are the costs
UnicodePlots.lineplot(costs,sum(prob,dims=2))

UnicodePlots.heatmap(prob)



cyears = expectedlifetime(mass); #Mean expected lifespan (yrs); Calder 1984; 10;
p_bad = 0.05; #Probability of a bad day (no autorcorrelation)
gdaily, cdaily, cgut, cfat, cfat_synthesizing = acrossdaysim(gains,costs,prob,edensity,mass,gut_type,cyears,p_bad);
relfatres = cfat/maxfatstorage(mass,37.7)[1];
Plots.plot(collect(1:length(relfatres))./365,relfatres,
            xlabel="Time (years)", 
            ylabel="Relative Fat Reserves",
            legend=false)





            
# Compare simulated vs. expected for daily simulations

# SIMULATED
# This provides a return distribution
@time gains, costs, prob = withindaysim(rho,alpha,mu,zeta,edensity,mass,teeth,gut_type,tmax_bout,configurations);


#Confirmed - these are the gains
UnicodePlots.lineplot(gains,sum(prob,dims=1)')
#Confirmed - these are the costs
UnicodePlots.lineplot(costs,sum(prob,dims=2))

UnicodePlots.heatmap(prob)

# EXPECTED NOTE: I don't think this works!
E_gains, Var_gains, E_costs, Var_costs = dailyforage_expectedvalues(mass, gut_type, teeth, rho, mu, alpha, zeta, edensity, tmax_bout)

#How does the mean gain compare to the simulated gain?
[E_gains, dot(gains,sum(prob,dims=1)')]

#How does the mean cost compare to the simulated cost?
[E_costs, dot(costs,sum(prob,dims=2))]




#Comparing stochasticity with within day splits
# This provides a return distribution
nosplit = Array{Float64}(undef,50,2);
split = Array{Float64}(undef,50,2);
for i=1:50
    @time gains, costs, prob = withindaysim(rho,alpha,mu,zeta,edensity,mass,teeth,gut_type,tmax_bout,configurations);

    nosplit[i,:] = [mean(gains),mean(costs)]

    @time gains, costs, prob = withindaysim_split(rho,alpha,mu,zeta,edensity,mass,teeth,gut_type,tmax_bout,configurations);

    split[i,:] = [mean(gains),mean(costs)]

end

mean(nosplit,dims=1)
mean(split,dims=1)


std(nosplit,dims=1)
std(split,dims=1)