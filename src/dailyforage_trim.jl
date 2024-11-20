function dailyforage_trim(gammadist, tchew, beta, maxgut, velocity, bcost_kJps, fcost_kJps, edensity, tmax_bout)

    max_enc = 100000;

    # Sample from gammadist once and compute the inverse
    gamma_samples = rand(gammadist, max_enc);
    scales = 1.0 ./ gamma_samples;

    # Draw from Exponential directly using the scale parameter
    distance_draws = rand.(Exponential.(scales));
    
    encounter_time = (distance_draws ./ velocity) .+ tchew;

    cum_encounter_time = cumsum(encounter_time);
    final_encounter = findfirst(x -> x > tmax_bout, cum_encounter_time);

    # total_encounter_time = cum_encounter_time[final_encounter - 1]

    per_encounter_gains = repeat([beta],final_encounter - 1);

    cum_gains = cumsum(per_encounter_gains)

    final_gain = something(findfirst(x -> x > maxgut, cum_gains), final_encounter - 1);
    
    #stop at reaching maxgut or stop at reaching tmax_bout?
    last_encounter = minimum([final_encounter - 1,final_gain - 1])

    rest_time = (24*60*60) - cum_encounter_time[last_encounter];

    total_cost = (cum_encounter_time[last_encounter] - (last_encounter * tchew))*fcost_kJps + ((last_encounter * tchew) + rest_time)*bcost_kJps;

    total_gain = cum_gains[last_encounter]*edensity
    
    return total_gain, total_cost, last_encounter
    
    
    # travel_perbite = Array{Float64}(undef,0)

    

    # #Initialize
    # t_travel = 0.0;
    # tchew = 0.0;
    # bites = 0;
    # gut = 0.0
    # t=0.0;

    # #Simulate daily returns and costs
    # while t < tmax_bout && gut < maxgut
        
    #     # tic += 1

    #     # Draw distance to next resource
    #     # First draw encounter rate
    #     rg = rand(gammadist);
    #     distance_to_resource = rand(Exponential(1.0/rg));
            
    #     #The forager will move towards the closest resource
    #     ttravel = distance_to_resource/velocity; # distance / velcity = time
        
    #     t += ttravel; # time
    #     t_travel += ttravel; # time

    #     push!(travel_perbite,ttravel);
        
    #     # If the forager successfully reaches the resource
    #     # Each resource is equal to 1 mouthful
    #     if (t + ttravel) < tmax_bout
    #         # number_of_successes += 1;
            
    #         # Time passes while chewing
    #         t += tchew; #time
    #         tchew += tchew; #time

    #         # Pass Mouth-Unit to Gut Total (boundary conditions in across-day sim)
    #         # resgain is kJ per volume (mouth volume)
    #         gut += beta; #grams/bite
    #         bites += 1;
            
    #     end

    # end


    # #kJ gains
    # gains = gut * edensity; #grams * kJ/gram = kJ returns

    # #Non-forage time
    # #NOTE: changed tmax_bout to t, because you can stop if your gut fills, such that t < tmax_bout, in which case, the remainder is rest.
    # #11/15/2024
    # t_rest = (24*60*60) - t; #tmax_bout;

    # #kJ costs
    # costs = fcost_kJps*t_travel + bcost_kJps*(tchew + t_rest);

    # return gains, costs, bites, t_travel, tchew, t_rest, travel_perbite, t
    
end