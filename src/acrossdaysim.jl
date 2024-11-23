function acrossdaysim(
    gains, 
    costs,
    probs,
    edensity, # Resource energy density kJ/gram
    mass, # KILOGRAMS
    gut_type,
    cyears,
    probbad) #consumer years

    #Second in a day
    secday = 60*60*24;

    # Maximum gut capacity
    # Think about this...
    maxgut = gut_capacity_g(mass, gut_type) * edensity; # grams * kJ/gram = kJ

    #Currently 4 days for a 500 Kg mammal
    #Correct according to Muller
    mrt = mean_retention_time(mass, gut_type); # seconds / PARTICLE

    # Partical mass = gram / particle (cube)
    # Paricle length = mm
    #NOTE: Not using this right now
    particle_mass, particle_length = mean_particle_mass(mass, gut_type); 
    
    # Calculate surface area per gram of ingesta
    #NOTE: Not using this right now
    ingesta_SA_per_gram = (6*particle_length^2)/particle_mass; #SA/particle * (particle/gram) mm^2/gram
    ingesta_SA_per_kj = ingesta_SA_per_gram * edensity; #mm^2 per kJ

    # Passage rate of food (rate of flow from gut to body)
    # passrate = (1/mrt) * particle_mass * edensity; #particle/s * gram / particle * kJ/gram = kJ/s

    # Single day gut passage
    # passage rate of a kJ within a single particle (needs to be multiplied by kJ in stomach to get total kJ passing in a day)
    # passrate_day = passrate * secday; # grams/s * s/day * kJ/gram = kJ/day

    # Metabolic loss per day
    # activehrs = 12;
    # maxfatstorage, cost_basalrate, cost_fieldrate = find_metabolism(mass,activehrs);
    #maxfatstorage = kJ
    #costs = kJ/day

    #Metabolic efficiency of fat synthesis
    Eprime = 7; #kJ cost of synthesizing 1 gram of fat (Yeakel et al. 2018)
    edensity_fat = 37.7; #kJ per gram of fat (Prange et al . AmNat 1979)

    maxfat_kj, maxfat_kg = maxfatstorage(mass,edensity_fat); # kJ

    #Time to rebuild entire fat storage:
    t_max_recovery = 1/fatrecoveryrate(mass,(mass-maxfat_kg)/mass,0.95) #seconds
    t_max_recovery_days = t_max_recovery/secday; # sec / (sec/day) = days

    #Efficiencies
    #Digestive Efficiencies
    eta_d = 0.6;

    
    eta_m = edensity_fat/(edensity_fat + Eprime); #(%)kJ cost of synthesizing 1 kJ of fat storage


    #Active time
    # field_cost = cost_fieldrate;    #kJ/day
    # rest_cost = cost_basalrate; #KJ/day

    # How long does it take to empty gut from maxgut to 1/2 maxgut?
    # twait = (maxgut/2)*(1/passrate); #grams * s/grams = seconds

     #### NEEDS WORK ####
    # Digestibility
    # NOTE What is the relationship between retention time and digestibility????
    # Closer to 1 for ruminants
    # Less than 1 for simple digestive systems

    # Because MRT is positively linked to the digestive efficiency of a herbivore (Foose, 1982; Udén and Van Soest, 1982; Clauss et al., 2007b),
    # epsilon = (max_retention_time(mass) - mean_retention_time(mass, gut_type))/(max_retention_time(mass) - (min_retention_time(mass)));
    # epsilon = 0.1;

    daysincyears = Int64(floor(365*cyears));  

    cgut = zeros(daysincyears);
    cfat = zeros(daysincyears);
    gain_day = zeros(daysincyears);
    cost_day = zeros(daysincyears);

    cfat_synthesizing = zeros(daysincyears);

    #Starting conditions
    cgut[1] = maxgut;
    # cabsorpta.
    cfat[1] = maxfat_kj;
    
    for t=2:daysincyears

        if cfat[t-1] > 0.0

             # We could build in environmental stress here
            # (Prevalence of Good <-> Bad days w/ autocorrelation)
            gooddaydraw = rand();

            if gooddaydraw < (1-probbad)
                
                #good day
                #Draw daily return (allowed to be > stomach size)
                gain_day[t], cost_day[t] = gc_sample(gains,costs,probs); #kJ, kJ

            else

                #bad day
                #Draw daily return and cost, but reset return to zero
                gain_day[t], cost_day[t] = gc_sample(gains,costs,probs); #kJ, kJ
                gain_day[t] = 0.0; #or some proportion?

            end

            
            #kj digested per gut cycle
            # (1/SA*time) * SA * time/cycle * kJ
            #Biophysical option
            #Gut surface area (from previous time interval - not just-updated interval)
            # gut_SA = particle_SA_kj*cgut[t-1];
            # absorpta = absorb_rate_perSA * gut_SA * mrt * cgut[t-1];
            
            #Shortcut
            absorpta = eta_d * cgut[t-1];

            #Updated: 10/31/2024 - happy Halloween!
            #Absorpta: absorbed for use from prior day in kJ
            # absorpta = eta_d*passrate_day*cgut[t-1]; #NOTE: This is very low
            
            poo = (1 - eta_d)*cgut[t-1];

            #Change in yesterday's gut
            deltagut_pregain =  -absorpta - poo;

            cgut[t] = clamp(cgut[t-1] + deltagut_pregain, 0.0, maxgut);

            #Change in today's gut
            #At the end of the day, the gut should have what we gained during the day
            deltagut_postgain = gain_day[t];

            cgut[t] = clamp(cgut[t] + deltagut_postgain, 0.0, maxgut);


            #Change in fat 
            deltafat = absorpta - cost_day[t];

            #Fat can be gained or lost, but at inverse rates
            #NOTE: need to think about this inverse rate idea
            if deltafat > 0

                prop_fat_gained = deltafat/maxfat_kj;
                t_synthesis_day = (prop_fat_gained*t_max_recovery_days); #days
                synthesis_rate_perday = 1/t_synthesis_day;

                #Days over which we deposit delta fat
                #First account for time retained in the gut
                t_retention = Int64(round(mrt/secday)); #days
                

                #### Adding gains to fat post-synthesis. Two options:

                #### OPTION 1: distribute gains across synthesis days ####
                #NOTE: this seems more realistic to me.
                #Second account for the synthesis time
                propday_future = decompose_to_integer_components(t_synthesis_day)
                t_future = (t + t_retention) .+ cumsum(propday_future .> 0)
                #Trim periods beyond max time step
                idx = searchsortedfirst(t_future, daysincyears + 1)  # Find first index where element > daysincyears
                t_future = copy(view(t_future, 1:idx-1))  # Create a view to avoid copying
                propday_future = copy(view(propday_future, 1:idx-1))
                #Add to synthesizing fat vector
                if length(t_future) > 0
                    cfat_synthesizing[t_future] = cfat_synthesizing[t_future] .+ 
                                                ((synthesis_rate_perday * deltafat) .* propday_future)
                end

                # #### OPTION 2: add gains at the end of synthesis days ####
                # t_future = t + t_retention + Int64(round(t_synthesis_day));
                # if t_future < daysincyears
                #     cfat_synthesizing[t_future] = cfat_synthesizing[t_future] + deltafat;
                # end
                

                

                #No use of current fat storage
                cfat[t] = cfat[t-1];

                #  #If there is net gain, add delta fat at rate eta_m
                # cfat[t] = clamp(cfat[t-1] + synthesis_rate_perday * deltafat, 0.0, maxfat_kj);
            
            else

                #If there is net cost, use current fat with 1/eta_m usage penalty
                #If storing fat is inefficient (low eta_m), the usage fee will be higher (high 1/eta_m)
                cfat[t] = clamp(cfat[t-1] + (1/eta_m) * deltafat, 0.0, maxfat_kj);

            end

            #Store synthesized fat, regardless of deltafat
            cfat[t] = clamp(cfat[t] + eta_m * cfat_synthesizing[t], 0.0, maxfat_kj);


            #Older code
            #Change in stomach contents
            # deltagut = (gain_day[t] - cgut[t-1]*passrate_day);
            # cgut[t] = clamp(cgut[t-1] + deltagut, 0.0, maxgut);
            # cgut[t] = minimum([maximum([cgut[t-1] + deltagut,0.0]),maxgut]);
            # deltafat = (epsilon*cgut[t-1]*passrate_day - dailycost);
            # deltafat = (cgut[t-1]*passrate_day - dailycost);
            #We might want an incorporation compartment... where gut contents are passed to a compartment, and added to fat at a constant rate...
            # cfat[t] = clamp(cfat[t-1] + deltafat, 0.0, maxfat);
            # cfat[t] = minimum([maximum([cfat[t-1] + deltafat,0.0]),maxfat]);
        end
    end
        
    # ctraits = tuple(maxfat)

    return gain_day, cost_day, cgut, cfat, cfat_synthesizing #, ctraits


end