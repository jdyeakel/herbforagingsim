using Plots
using UnicodePlots
using Distributions

mass = 100*1000 # ALL IN GRAMS
B0 = 0.047 #Wg-3/4
Emprime = 7000 #J/g
aprime = B0/Emprime
epsilon_lam = 0.95
eta = 3/4

fatmass =  (0.0202*(mass)^1.19)

start_perc = (mass - fatmass)/mass

epsilonvec = collect(start_perc:0.001:1)
t_recover = Array{Float64}(undef,length(epsilonvec));

for i=1:length(epsilonvec)
    epsilon_sig = start_perc;
    epsilon_lam = epsilonvec[i];
    t_recover[i] = log(((1-(epsilon_sig * epsilon_lam)^(1-eta))/(1-epsilon_lam^(1-eta))))*(mass^(1-eta)/(aprime*(1-eta)))
end
#Take geometric mean
exp(mean(log.(t_recover[1:end-5])))

UnicodePlots.scatterplot(epsilonvec,(t_recover))



start_perc = 0.9 #(mass - fatmass)/mass
epsilon_sig = start_perc;
epsilon_lam = start_perc+0.01;
t_recover = log(((1-(epsilon_sig * epsilon_lam)^(1-eta))/(1-epsilon_lam^(1-eta))))*(mass^(1-eta)/(aprime*(1-eta)))