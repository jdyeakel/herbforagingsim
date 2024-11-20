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


#Consumer attributes
mass = 100;
#Define tooth and gut type of herbivore
teeth = "bunodont"; # bunodont, acute/obtuse lophs, lophs and non-flat, lophs and flat
gut_type = "caecum"; # caecum, colon, non-rumen foregut, rumen foregut
beta = bite_size_allo(mass); # mass in kg, bite size in g/bite
chewrate = chew_allo(mass, teeth); #g/s
t_chewgram = 1/chewrate; #s/g
t_chew = t_chewgram * beta; #s/g * g/bite = s/bite
maxgut = gut_capacity_g(mass, gut_type) #grams
bcost_kJps, fcost_kJps = metabolic_cost(mass);
velocity = find_velocity(mass); # m/s
tmax_bout, _ = foragingtime(mass) .* (60*60) #hours * 60 min/hour * 60 sec/min = secs
#Consumer population density: individuals/m^2
n = indperarea(mass); #individuals/m^2

#Resource attributes
rhoexp = -6.5;
rho = 10^rhoexp;
#Define resource traits
mu = 1;
alpha = 3
edensity = 18.2;
zeta = 1;

# rho * mu * (1/beta) = bite encounters/m^2 ~ bite encounter rate
m = rho*mu*(1/beta);

#Adjusted for competitive landscape
mprime = m/n;
alphaprime = alpha*n^(zeta-2);


# ########################################

# #Examine without competition
# #Specify gamma distribution characterizing the distribution of Lmabda
# #We don't need to broadcast operations because we just have a single resource type
# gammadist = Gamma(alpha, m / alpha);

# # Sample from gammadist once and compute the inverse
# gamma_samples = rand(gammadist, 100000);
# scales = 1.0 ./ gamma_samples;

# # Draw from Exponential directly using the scale parameter
# distance_draws = rand.(Exponential.(scales));

# #Estimate mean and sd from the drawn distribution
# est_mean = mean(distance_draws);
# est_sd = std(distance_draws);

# #Calculated mean and sd from our equations
# calc_mean = alpha/(m*(alpha-1));
# calc_sd = sqrt(alpha^3/(((alpha-1)^2)*(alpha-2)*m^2));

# # Compare the estimated and calculated means and variances
# comparison_df = DataFrame(
#     Metric = ["Mean", "Std"],
#     Estimated = [est_mean, est_sd],
#     Calculated = [calc_mean, calc_sd]
# )

# # Display the table
# println(comparison_df)


########################################

#Now examine WITH effects of competition
#Specify gamma distribution characterizing the distribution of Lmabda
#We don't need to broadcast operations because we just have a single resource type
gammadist = Gamma(alphaprime, mprime / alphaprime);

# Sample from gammadist once and compute the inverse
gamma_samples = rand(gammadist, 100000);
scales = 1.0 ./ gamma_samples;

# Draw from Exponential directly using the scale parameter
distance_draws = rand.(Exponential.(scales));

#Estimate mean and sd from the drawn distribution
est_mean = mean(distance_draws);
est_sd = std(distance_draws);


#Calculated mean and sd from our equations
calc_mean = (alpha*n^(zeta-1)) / (m*(alpha*n^(zeta-2)-1));
# calc_sd = sqrt((((n^2)*(alpha*n^(zeta-2))^3) / ((alpha*n^(zeta-2)-1)^2)*(alpha*n^(zeta-2)-2)*m^2));
calc_sd = sqrt(((alpha*n^(zeta-2))^3) / (((alpha*n^(zeta-2)-1)^2)*(alpha*n^(zeta-2)-2)*(m/n)^2))

# Compare the estimated and calculated means and variances
comparison_df = DataFrame(
    Metric = ["Mean", "Std"],
    Estimated = [est_mean, est_sd],
    Calculated = [calc_mean, calc_sd]
)

# Display the table
println(comparison_df)

#Expectation and Variance of travel time
ExpT = (1/velocity)*((alpha*n^(zeta-1)) / (m*(alpha*n^(zeta-2)-1)))
VarT = (1/velocity^2)*(((alpha*n^(zeta-2))^3) / (((alpha*n^(zeta-2)-1)^2)*(alpha*n^(zeta-2)-2)*(m/n)^2))

exp_encounters = (tmax_bout/(ExpT + t_chew)) + (1/2)*((2*tmax_bout)/(ExpT+t_chew)^3)*VarT

# From my mathematica notebook, I get this expression for the Expected number of encounters in tbout:
# exp_encounters = (tmax_bout * (t_chew^2-(2 * n^(1+zeta) * t_chew * alpha)/(m * n^2 * velocity - m * n^zeta * velocity * alpha) + (2 * n^(2 + 2 * zeta) * alpha^2) / (m^2 * velocity^2 * (2 * n^4 - 3 * n^(2+zeta) * alpha + n^(2 * zeta) * alpha^2)))) / (t_chew - (n^(1 + zeta) * alpha) / (m * n^2 * velocity - m * n^zeta * velocity * alpha))^3

# Simulation of the daily forage
configurations = 20000;
gains = Vector{Float64}(undef, configurations);
costs = Vector{Float64}(undef, configurations);
encounters = Vector{Float64}(undef, configurations);
travel = Vector{Float64}(undef, configurations);
chew = Vector{Float64}(undef, configurations);
rest = Vector{Float64}(undef, configurations);
mean_travel_perbite = Vector{Float64}(undef, configurations);
var_travel_perbite = Vector{Float64}(undef, configurations);
cum_t = Vector{Float64}(undef, configurations);
@threads for i in 1:configurations
    gains_daily, costs_daily, encounters_daily, t_travel_daily, t_chew_daily, t_rest_daily, travel_perbite, cum_t_daily = dailyforage_experiment(mass, gut_type, teeth, rho, mu, alpha, zeta, edensity, tmax_bout)
    
    # Store results directly in preallocated arrays
    gains[i] = gains_daily
    costs[i] = costs_daily
    encounters[i] = encounters_daily
    travel[i] = t_travel_daily
    chew[i] = t_chew_daily
    rest[i] = t_rest_daily
    cum_t[i] = cum_t_daily

    mean_travel_perbite[i] = mean(travel_perbite)
    var_travel_perbite[i] = var(travel_perbite)
end

comparison_ne = DataFrame(
    Metric = ["MeanT", "VarT","Mean Enc"],
    Simulated = [mean(mean_travel_perbite), mean(var_travel_perbite),mean(encounters)],
    Calculated = [ExpT, VarT, exp_encounters]
)



# Now I've written a new 'trimmed' version of the dailyforage simulation (dailyforage_v2.jl)

configurations = 200000;
gammadist = Gamma(alphaprime, mprime / alphaprime);

# Preallocate arrays with the size of configurations
gains = Vector{Float64}(undef, configurations)
costs = Vector{Float64}(undef, configurations)
encounters = Vector{Float64}(undef, configurations)
@time @threads for i in 1:configurations
    gains_daily, costs_daily, encounters_daily = dailyforage(gammadist,tchew, beta, maxgut, velocity, bcost_kJps, fcost_kJps, edensity, tmax_bout)
    
    # Store results directly in preallocated arrays
    gains[i] = gains_daily
    costs[i] = costs_daily
    encounters[i] = encounters_daily
end




# Preallocate arrays with the size of configurations
gains_trim = Vector{Float64}(undef, configurations)
costs_trim = Vector{Float64}(undef, configurations)
encounters_trim = Vector{Float64}(undef, configurations)
@time @threads for i in 1:configurations
    gains_daily, costs_daily, encounters_daily = dailyforage_trim(gammadist,tchew, beta, maxgut, velocity, bcost_kJps, fcost_kJps, edensity, tmax_bout)
    
    # Store results directly in preallocated arrays
    gains_trim[i] = gains_daily
    costs_trim[i] = costs_daily
    encounters_trim[i] = encounters_daily
end

