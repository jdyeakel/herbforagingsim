function withindaysim(
    rho,
    alpha, # Resource dispersion
    mu,  # Resource mean
    zeta, # Resource variability scaling
    edensity, # Resource energy density kJ/gram
    mass, # KILOGRAMS
    teeth, # bunodont, acute/obtuse lophs, lophs and non-flat, lophs and flat
    tmax_bout,
    configurations
    )

    # probability = Array{Float64}(undef,kmax+1);
    # kinfo = Array{Float64}(undef,kmax+1);
    # # probability = SharedArray{Float64}(kmax+1);
    # # kinfo = SharedArray{Float64}(kmax+1);
    # data = zeros(Float64,configurations);

    # gains = Float64[]
    # costs = Float64[]
 
    # for _ in 1:configurations
    #     gains_daily, costs_daily = dailyforage(mass,teeth,rho,mu,alpha,zeta,edensity,tmax_bout);
    #     push!(gains,gains_daily);
    #     push!(costs,costs_daily);
    # end

    # Preallocate arrays with the size of configurations
    gains = Vector{Float64}(undef, configurations)
    costs = Vector{Float64}(undef, configurations)
    @threads for i in 1:configurations
        gains_daily, costs_daily = dailyforage(mass, teeth, rho, mu, alpha, zeta, edensity, tmax_bout)
        
        # Store results directly in preallocated arrays
        gains[i] = gains_daily
        costs[i] = costs_daily
    end

    nbins = 50;

    # Calculate histograms
    gains_bin_edges, gains_bin_counts = calculate_histogram(gains, nbins)
    costs_bin_edges, costs_bin_counts = calculate_histogram(costs, nbins)

    # Midpoints for better representation
    gains_bin_midpoints = (gains_bin_edges[1:end-1] + gains_bin_edges[2:end]) / 2
    costs_bin_midpoints = (costs_bin_edges[1:end-1] + costs_bin_edges[2:end]) / 2

    # Joint probability matrix setup
    gain_bins = range(minimum(gains), maximum(gains), length=nbins+1)
    cost_bins = range(minimum(costs), maximum(costs), length=nbins+1)
    joint_histogram = zeros(Int, nbins, nbins)

    # Fill in joint histogram
    for i in 1:length(gains)
        gain_bin_index = find_bin_index(gains[i], gain_bins)
        cost_bin_index = find_bin_index(costs[i], cost_bins)
        joint_histogram[gain_bin_index, cost_bin_index] += 1
    end

    # Normalize to create joint probabilities
    joint_probabilities = joint_histogram / sum(joint_histogram)

    return gains_bin_midpoints, costs_bin_midpoints, joint_probabilities

end
