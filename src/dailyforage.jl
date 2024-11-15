function dailyforage(mass,gut_type,teeth,rho,mu,alpha,zeta,edensity,tmax_bout)

    bcost_kJps, fcost_kJps = metabolic_cost(mass);

    velocity = find_velocity(mass); # m/s

    # CALCULATE tchew
    # NOTE: BITE SIZE SEEMS SMALL
    # bite_rate = bite_rate_allo(mass); # mass in kg, bite/s
    beta = bite_size_allo(mass); # mass in kg, bite size in g/bite

    chewrate = chew_allo(mass, teeth); #g/s
    tchewgram = 1/chewrate; #s/g
    tchew = tchewgram * beta; #s/g * g/bite = s/bite

    maxgut = gut_capacity_g(mass, gut_type) #grams

    #Consumer population density: individuals/m^2
    n = indperarea(mass); #individuals/m^2

    # rho * mu * (1/beta) = bite encounters/m^2 ~ bite encounter rate
    m = rho*mu*(1/beta);

    #Adjusted for competitive landscape
    mprime = m/n;
    alphaprime = alpha*n^(zeta-2);

    #Define Gamma Distribution for resource availability
    gammadist = Gamma(alphaprime,mprime/alphaprime); #mean = alpha * m/alpha

    #Initialize
    t_travel = 0.0;
    t_chew = 0.0;
    bites = 0;
    gut = 0.0
    t=0.0;

    #Simulate daily returns and costs
    while t < tmax_bout && gut < maxgut
        
        # tic += 1

        # Draw distance to next resource
        # First draw encounter rate
        rg = rand(gammadist);
        distance_to_resource = rand(Exponential(1.0/rg));
            
        #The forager will move towards the closest resource
        ttravel = distance_to_resource/velocity; # distance / velcity = time
        
        t += ttravel; # time
        t_travel += ttravel; # time
        
        # If the forager successfully reaches the resource
        # Each resource is equal to 1 mouthful
        if (t + ttravel) < tmax_bout
            # number_of_successes += 1;
            
            # Time passes while chewing
            t += tchew; #time
            t_chew += tchew; #time

            # Pass Mouth-Unit to Gut Total (boundary conditions in across-day sim)
            # resgain is kJ per volume (mouth volume)
            gut += beta; #grams/bite
            bites += 1;
            
        end

    end


    #kJ gains
    gains = gut * edensity; #grams * kJ/gram = kJ returns

    #Non-forage time
    #NOTE: changed tmax_bout to t, because you can stop if your gut fills, such that t < tmax_bout, in which case, the remainder is rest.
    #11/15/2024
    t_rest = (24*60*60) - t; #tmax_bout;

    #kJ costs
    costs = fcost_kJps*t_travel + bcost_kJps*(t_chew + t_rest);

    return gains, costs
    
end