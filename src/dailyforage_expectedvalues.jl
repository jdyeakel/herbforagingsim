function dailyforage_expectedvalues(mass, gut_type, teeth, rho, mu, alpha, zeta, edensity, tmax_bout)
    # Metabolic costs (kJ per second)
    bcost_kJps, fcost_kJps = metabolic_cost(mass)
    
    # Foraging velocity (m/s)
    velocity = find_velocity(mass)
    
    # Bite size (grams per bite)
    beta = bite_size_allo(mass)
    
    # Chewing rate (grams per second)
    chewrate = chew_allo(mass, teeth)
    tchewgram = 1 / chewrate  # Time to chew one gram (s/g)
    tchew = tchewgram * beta  # Time to chew one bite (s)
    
    # Maximum gut capacity (grams)
    maxgut = gut_capacity_g(mass, gut_type)
    
    # Consumer population density (individuals per m²)
    n = indperarea(mass)
    
    # Resource encounter rate per area (bites per m²)
    m = rho * mu * (1 / beta)
    
    # Adjusted for competitive landscape
    mprime = m / n
    alphaprime = alpha * n^(zeta - 2)
    
    # Gamma distribution parameters (resource encounter rate)
    theta = mprime / alphaprime  # Scale parameter θ of Gamma distribution
    # Ensure alphaprime > 2 for finite variance
    if alphaprime <= 2
        error("alphaprime must be greater than 2 for variance to be finite.")
    end
    
    # Lomax (Pareto Type II) distribution parameters
    # Scale parameter λ = θ, shape parameter k = alphaprime
    # Expected distance to resource
    E_distance = theta / (alphaprime - 1)
    Var_distance = (theta^2 * alphaprime) / ((alphaprime - 1)^2 * (alphaprime - 2))
    
    # Expected travel time per foraging cycle
    E_ttravel = E_distance / velocity
    Var_ttravel = Var_distance / velocity^2
    
    # Expected time per foraging cycle (travel + chew)
    E_ti = E_ttravel + tchew
    Var_ti = Var_ttravel  # Variance due to travel time only (chewing time is constant)
    
    # Expected number of foraging cycles (bites)
    N_time_limit = tmax_bout / E_ti
    N_gut_limit = maxgut / beta
    if N_time_limit < N_gut_limit
        # Time limit is the limiting factor
        E_N = N_time_limit
        Var_N = (Var_ti / E_ti^3) * tmax_bout
    else
        # Gut capacity is the limiting factor
        E_N = N_gut_limit
        Var_N = 0  # N is deterministic when gut capacity is limiting
    end
    
    # Expected total foraging time
    E_T = E_N * E_ti
    Var_T = E_N * Var_ti + Var_N * E_ti^2
    
    # Expected total travel and chewing times
    E_t_travel_total = E_N * E_ttravel
    Var_t_travel_total = E_N * Var_ttravel + Var_N * E_ttravel^2
    E_t_chew_total = E_N * tchew  # Chewing time is deterministic
    Var_t_chew_total = 0  # No variance in chewing time
    
    # Expected rest time
    E_t_rest = (24 * 60 * 60) - E_T
    Var_t_rest = Var_T  # Variance in rest time equals variance in foraging time
    
    # Expected gains (kJ)
    E_gains = E_N * beta * edensity
    Var_gains = Var_N * (beta * edensity)^2
    
    # Expected costs (kJ)
    E_costs = fcost_kJps * E_t_travel_total + bcost_kJps * (E_t_chew_total + E_t_rest)
    Var_costs = (fcost_kJps^2) * Var_t_travel_total + (bcost_kJps^2) * (Var_t_chew_total + Var_t_rest)
    
    return E_gains, Var_gains, E_costs, Var_costs
end
