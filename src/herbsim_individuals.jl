function herbsim_individuals(mass,teeth,gut_type,rhoexp,mu,alpha,edensity,zeta,configurations,p_bad,runs)

    #Best fit
    tmax_bout, _ = foragingtime(mass) .* (60*60) #hours * 60 min/hour * 60 sec/min = secs

    # #Upper 95CI
    # _, tmax_bout = foragingtime(mass) .* (60*60) #hours * 60 min/hour * 60 sec/min = secs

    rho = 1*10^rhoexp
    
    # res_traits = (mu, alpha, edensity);
    
    # Use a single withinday distribution
    gains, costs, probs = withindaysim_split(rho,alpha,mu,zeta,edensity,mass,teeth,gut_type,tmax_bout,configurations);

    cyears = expectedlifetime(mass); #Mean expected lifespan (yrs); Calder 1984; 10;
    daysincyears = Int64(floor(365*cyears));

    gains_inds = Array{Float64}(undef,runs,daysincyears);
    costs_inds = Array{Float64}(undef,runs,daysincyears);
    gut_inds = Array{Float64}(undef,runs,daysincyears);
    fat_inds = Array{Float64}(undef,runs,daysincyears);
    fatsynth_inds = Array{Float64}(undef,runs,daysincyears);

    # Many across-day simulations
    @threads for i=1:runs
    # for i=1:runs
        gdaily, cdaily, cgut, cfat, cfat_synthesizing = acrossdaysim(gains,costs,probs,edensity,mass,gut_type,cyears,p_bad);
        # relfatres = cfat/maxfatstorage(mass,37.7)[1];

        gains_inds[i,:] = gdaily;
        costs_inds[i,:] = cdaily;
        gut_inds[i,:] = cgut;
        fat_inds[i,:] = cfat;
        fatsynth_inds[i,:] = cfat_synthesizing;

    end

    return gains_inds, costs_inds, gut_inds, fat_inds, fatsynth_inds

end
