# most but not all of these are currently in use.  

# plant digestibility stuff
# from nutiritional content of savanna plants
# codron et al 2007
# tree foliage: %N 2.1, %NDF 38.1, %ADF 26.2, %ADL 12.1
# grasses:      %N 1.0, %NDF 68.2, %ADF 36.3, %ADL 4.9
# cell wall digestibility CWD = 100 - (ADL/NDF*100)

function plant_digestibility(adl, ndf)
    # takes a plant's acid detergent lignin and neutral detergent fiber percentages
    # returns a % digestibility rating (cell wall digestibility)
    cwd = 100 - (ADL/NDF*100)
    return cwd
end

function apparent_digestibility(cwd, mrt)
    # temporary solution
    # how much digestion actually occurs is product of cell wall properties and retention time
    # it looks like 96hrs is used as max in digestion papers, so we'll say thats max possible
    # missing effect of mean particle size/reduction
    digestion = mrt/(96*60*60) # mrt [s] / [hrs->min->s]
    return digestion
end


function generate_trait_array(species_list)
    # takes the species list and returns an array of alpha and beta rates
    trait_array = []
    for i = 1:number_species
        a,b,c,d = find_metabolism(species_list[i,1])
        push!(trait_array, [
            bite_size_allo(species_list[i,1], "graze") # don't need full mouth, just full normal bite
            gut_volume_g(species_list[i,1], species_list[i,3])
            mean_retention_time(species_list[i,1], species_list[i,3])
            mean_particle_mass(species_list[i,1], species_list[i,3])
            find_velocity(species_list[i,1])
            a
            b
            c
            d
            ])
    end
    # mouth size, gut size, mrt, mps, velocity, initial_energy_state, basal, field, fat_max
    return trait_array 
end


function find_velocity(mass)
    # mass in kg
    # from kramer 2010 "allometric scaling of resource acquisition"
    #Consumer Velocity (meters/second)
    velocity = (0.5 * mass^0.13);
    return velocity
end



function generate_rate_array(species_list)
    # takes the species list and returns an array of alpha and beta rates
    rate_array = []
    for i = 1:number_species
        push!(rate_array, [alpha_allo(species_list[i,1], "graze")
            alpha_allo(species_list[i,1], "browse")
            beta_allo(species_list[i,1], species_list[i,2])])
    end
    return rate_array # alpha graze, alpha browse, beta
end


function bite_size_allo(mass)
    # from shipley 94 "the scaling of intake rate"
    # mass in kg
    # bite size in g
    #if plant == "browse"
    #    bite_size = (0.057* (mass)^0.63); # [g]
    #elseif plant == "graze"
    #    bite_size = (0.026 * (mass)^0.59); # [g]
    #end

    # Why elephants have trunks Pretorius 2015
    bite_size = (0.002 * (mass)^0.969); # [g]
    return bite_size

end



function bite_rate_allo(mass)
    # Why elephants have trunks Pretorius 2015
    bite_rate = 0.37 * mass^(-0.024) #(bites/s)  
    return bite_rate
end



function alpha_allo(mass)
    bite_rate = bite_rate_allo(mass); # mass in kg, bite/s
    bite_size = bite_size_allo(mass) # mass in kg, bite size in g
    alpha = bite_rate * bite_size # bite/s * g/bite
    return alpha
end



function number_of_chews(mass)
    # shipley 94
   #chew/g (processed to mean particle size (allo) and swallowed)
   # mass in kg
    
   #uses the correction factor reported in Shipley
   chews_per_gram = 343.71* mass^(-0.83); 
    

   return chews_per_gram

end



function chew_rate_allo(mass, teeth)
   # from "dental functional morphology predicts scaling"
   # mass in kg
   # duration in ms
   if teeth == "bunodont"
       chewing_cycle_duration = (228.0* (mass)^0.246) / 1000; # 2.358 [ms -> s]

   elseif teeth == "acute/obtuse lophs"
       chewing_cycle_duration = (299.2 * (mass)^0.173) / 1000; # 2.476 [ms -> s]

   elseif teeth == "lophs and non-flat"
       chewing_cycle_duration = (320.6 * (mass)^0.154) / 1000; # 2.506 [ms -> s]

   elseif teeth == "lophs and flat"
       chewing_cycle_duration = (262.4* (mass)^0.207) / 1000; # 2.419 [ms -> s]0

   end

   return 1 / (chewing_cycle_duration )  #[s/chew -> chews/s]

end


function chew_allo(mass, teeth)
    # allometric function for mouth->gut rate
   chew_rate = chew_rate_allo(mass, teeth); # chews/s
   chews_per_gram = number_of_chews(mass)   # chew/g
   chew = chew_rate / chews_per_gram        # chew/s / chew/g -> g/s
   return chew
end



function mean_retention_time(mass, gut_type)
    # from Muller 2013 - Assessing the Jarmain-Bell Principle
    # mrt in hrs
    # bm in kg
    # mean retention of a particle in the gut [s]

    if gut_type == "caecum"
        mean_retention_time = (21.7 * (mass)^0.271)

    elseif gut_type == "colon"
        mean_retention_time = (47.1 * (mass)^0.0) #not significant, double check using zero is right

    elseif gut_type == "non-rumen foregut"
        mean_retention_time = (30.3 * (mass)^0.109)

    elseif gut_type == "rumen foregut"
        mean_retention_time = (27.6 * (mass)^0.133)
    end
    

    return mean_retention_time * 60 * 60 # [hr -> s]

end

function max_retention_time(mass)
    # from Muller 2013 - Assessing the Jarmain-Bell Principle
    # mrt in hrs
    # bm in kg
    # mean retention of a particle in the gut [s]

    mean_retention_time_v = Array{Float64}(undef,4);

    # "caecum"
    mean_retention_time_v[1] = (21.7 * (mass)^0.271)

    # "colon"
    mean_retention_time_v[2] = (47.1 * (mass)^0.0) #not significant, double check using zero is right

    # "non-rumen foregut"
    mean_retention_time_v[3] = (30.3 * (mass)^0.109)

    # "rumen foregut"
    mean_retention_time_v[4] = (27.6 * (mass)^0.133)

    return maximum(mean_retention_time_v) * 60 * 60 # [hr -> s]

end

function min_retention_time(mass)
    # from Muller 2013 - Assessing the Jarmain-Bell Principle
    # mrt in hrs
    # bm in kg
    # mean retention of a particle in the gut [s]

    mean_retention_time_v = Array{Float64}(undef,4);

    # "caecum"
    mean_retention_time_v[1] = (21.7 * (mass)^0.271)

    # "colon"
    mean_retention_time_v[2] = (47.1 * (mass)^0.0) #not significant, double check using zero is right

    # "non-rumen foregut"
    mean_retention_time_v[3] = (30.3 * (mass)^0.109)

    # "rumen foregut"
    mean_retention_time_v[4] = (27.6 * (mass)^0.133)

    return minimum(mean_retention_time_v) * 60 * 60 # [hr -> s]

end




function gut_capacity_g(mass, gut_type)
    # from Muller 2013 - Assessing the Jarmain-Bell Principle
    # dry mass of guts in kg
    # bm in kg

    if gut_type == "caecum"
        capacity = (0.025*(mass)^0.860);

    elseif gut_type == "colon"
        capacity = (0.029 * (mass)^0.919);

    elseif gut_type == "non-rumen foregut"
        capacity = (0.030 * (mass)^0.881);

    elseif gut_type == "rumen foregut"
        capacity = (0.041 * (mass)^0.897);

    end

    return capacity * 1000 #[kg -> g]

end



function mean_particle_mass(mass, gut_type)
    # from "comparative chewing efficiency"
    # mean particle size in mm
    # ingesta particle size
    # mass in g (!!!!!) This might need to be KG
    # "We investigated faecal particle size, which reflects the size of ingesta particles after both mastication and specialized processes such as rumination."

    if gut_type == "rumen foregut"
        
        #Ruminant
        particle_length = (0.177)*mass^(0.15)
    
    elseif gut_type == "non-rumen foregut"

        #Nonruminant foregut fermenter
        particle_length = (0.176)*mass^(0.57)

    else
        
        #Hindgut fermenter = presume to be caecum or colon
        particle_length = (0.478)*mass^(0.146)

    end

    # OLD - Do NOT think these are correct!
    # mass_g = mass * 1000; # convert mass kg -> g
    # mean_particle_size = 0.0;
    # if gut_type == "rumen foregut" # ruminant
    #     mean_particle_size = (7.74 * (mass_g)^0.22)

    # elseif gut_type == "colon" # hindgut
    #     mean_particle_size = (6.61 * (mass_g)^0.26)

    #     else mean_particle_size = (6.61 * (mass_g)^0.26) # print("what are these guts?")

    # end
    #Sphere
    # volume = (4/3) * pi * (1/2 * mean_particle_size)^3; # [mm^3]
    #Cube
    volume = particle_length^3; # [mm^3]

    #Water density = 0.001 g/mm^3
    # density = 0.001;
    #Median particle density from Clauss et al. 2011 (basically waterish!)
    density = 0.0012  # g/mm^3

    particle_mass = density * volume; # [g/mm^3 * mm^3 = g]

    return particle_mass, particle_length

end


function outflow_rate(gut_fill_max, mrt)
    # function to capture processing rate of gut
    # gut_fill [g], mrt [s]
    gamma =  gut_fill_max / mrt #gut_fill_max  #[g/s]
    return gamma
end





# revised metabolism
# originally from Prins and Van Langevelde (2008)
# I'm getting from kramer & prins 2010
# cb = 3.4 * mass^-0.75 # j/s
# cf = 1.54 * cb # multip factor travel, cropping, or handling
# cd = 1.24 # multip factor for digestion




function metabolic_cost(mass)
    # takes a terrestrial mammal body mass (kg), returns storage masses and metabolic rates
    #Convert to grams
    mass_g = mass*1000;

    b0_bmr = 0.018; #watts g^-0.75
    b0_fmr = 0.047; #watts g^-0.75
    bcost_watts = (b0_bmr*(mass_g^0.75)); #Cost in Watts
    fcost_watts = (b0_fmr*(mass_g^0.75)); #Cost in Watts
    
    #Watt = 0.001 kJ/s
    kJps = 0.001;

    bcost_kJps = bcost_watts*kJps; #kJ per s
    fcost_kJps = fcost_watts*kJps; #kJ per s

    # mass_g = mass * 1000; #[kj]->[g]
    
    # function for setting initial energy and metabolic costs from body mass
    
    #Joules per gram
    # joules_per_gram = 20000; #varies between 7000 to 36000 [int] [J/g]
    # kjoules_per_gram = joules_per_gram / 1000;   # [int/int=int] [kJ/g]

    # initial_energy_state = mass_g * kjoules_per_gram;


    # from Yeakel et al 2018 Dynamics of starvation and recovery
    #mass at which you die, no fat or muscle
    # fat_mass =  0.02*mass_g^1.19           #[g]
    #muscle_mass = 0.383*mass_g^1.0         #[g]
    # you starve when muscle and fat stores have been depleted
    #mass_starve = round(mass_g - (fat_mass + muscle_mass));  #[g]

    #how many kj units does this organism have?
    #convert grams to kJ [int-float*int=float] [g-g*kJ/g=kJ]
    # storage_kj = (fat_mass) * kjoules_per_gram;

    # Metabolic constants for the basal and field metabolic rate
   # b0_basal_met_rate = 0.018; #[watts] g^-0.75,
   # b0_field_met_rate = 0.047; #[watts] g^-0.75,

    #costs: f/df + sleeping over active hours
    #cost_wh_basal = (b0_basal_met_rate * (mass_g^0.75)); #watt*hour, cost of basal metabolic rate in watt hours
    #cost_wh_field = (b0_field_met_rate * (mass_g^0.75)); #watt*hour, cost of field metabolic rate in watt hours

    #Convert to kiloJoules
    #watt_hour_to_kJ = 3.6;  #  [float], [kJ/watthour]

    #metabolic costs per hour
    #cost_basal_hr = cost_wh_basal * watt_hour_to_kJ; # [float], [wh*kJ/wh=kJ/hr]
    #cost_field_hr = cost_wh_field * watt_hour_to_kJ; # [float], [wh*kJ/wh=kJ/hr]

    # metabolic costs per second
    #cost_basal = cost_basal_hr / 60 / 60;   # [float], [kJ/s]
    #cost_field = cost_field_hr / 60 / 60;   # [float], [kJ/s]
    
    # cost_basal =  3.4/1000 * mass^-0.75 # j/s -> kj/s
    # cost_field = 1.54 * cost_basal # kj/s

    #Metabolic costs from Sevilleta papers
    #These euations are for mass in grams, hence (mass*1000)
     #Metabolic constants for the basal and field metabolic rate
    # b0_bmr = 0.018; #watts g^-0.75
    # b0_fmr = 0.047; #watts g^-0.75
    # nonactivehours = 24-activehours;
    # cwh_df = (b0_bmr*((mass*1000)^0.75))*activehours + (b0_bmr*((mass*1000)^0.75))*nonactivehours; #watt*hour
    # cwh_f = (b0_fmr*((mass*1000)^0.75))*activehours + (b0_bmr*((mass*1000)^0.75))*nonactivehours; #watt*hour
    # #Convert to kiloJoules
    # whrkJ = 3.6;

    # cost_basal = cwh_df*whrkJ; #kJ per day
    # cost_field = cwh_f*whrkJ; #kJ per day

    return bcost_kJps, fcost_kJps

end

function indperarea(mass)
    #Enter mass in kg
    #Convert to grams
    massg = mass*1000;
    popdensity = (0.0116)*massg^-0.776; #inds/area (from Damuth)
    return popdensity
end


function maxfatstorage(mass,edensity_fat)
    mass_g = mass * 1000; #[kg]->[g]

    # from Yeakel et al 2018 Dynamics of starvation and recovery
    #mass at which you die, no fat or muscle
    fat_mass =  0.02*mass_g^1.19;           #[g]

    storage_kg = fat_mass/1000;

    #Joules per gram
    # joules_per_gram = 20000; #varies between 7000 to 36000 [int] [J/g]
    # kjoules_per_gram = joules_per_gram / 1000;   # [int/int=int] [kJ/g]
    kjoules_per_gram = edensity_fat; #kJ/g

    storage_kj = (fat_mass) * kjoules_per_gram;

    return storage_kj, storage_kg
end


function fatrecoveryrate(mass,start_perc,end_perc)
    #From Yeakel et al. 2018
    #mass should be in grams

    mass_g = mass*1000;

    B0 = 0.047; #Wg-3/4
    Emprime = 7000; #J/g
    aprime = B0/Emprime;
    eta = 3/4;

    epsilon_sig = start_perc;
    epsilon_lam = end_perc;
    t_recover = log(((1-(epsilon_sig * epsilon_lam)^(1-eta))/(1-epsilon_lam^(1-eta))))*(mass_g^(1-eta)/(aprime*(1-eta))) #seconds

    fatsynthrate = 1/t_recover; #1/sec

    return fatsynthrate
end

function expectedlifetime(mass)
    #Calder 1984;
    #mass in kg
    #expected lifetime in years
    explifetime = 5.08*mass^0.35;
    return explifetime
end

function foragingtime(mass)
    #mass in kg
    #foraging time in % of 24 hours

    forageperc = 21.0885*mass^0.124
    forageperc_upper95 = 27*mass^0.17

    foragehrs = (forageperc/100)*24
    foragehrs_upper95 = (forageperc_upper95/100)*24

    return foragehrs, foragehrs_upper95 # hours
end

function gut_turnover_time(mass,gut_type,edensity)
    # mass in kg
    # edensity of food in kJ/g

    # maxgut
    maxgut = gut_capacity_g(mass, gut_type) * edensity; #kJ

    #metabolic demand (rate) kJ per day
    #From Calder book pg 126
    dailymetrate = 504*mass^0.76;

    turnovertime = maxgut / dailymetrate; # days
    turnovertime_hrs = turnovertime * 24; # days * 24 hours/day = hrs

    return turnovertime_hrs

end

function dailyfoodintake(mass)
    intake_kJday = 971*mass^0.73
    return intake_kJday
end



# # Bite chew time allometry
# massvec = [10^i for i=collect(0:0.1:5)];
# teeth = "bunodont";
# betavec = bite_size_allo.(massvec);
# chewratevec = chew_allo.(massvec,teeth);
# tchewgramvec = 1 ./ chewratevec;
# tchewvec = tchewgramvec .* betavec
# namespace = smartpath("figures/tchew_allometry.pdf")
# R"""
# pdf($namespace,width=4,height=4)
# plot($massvec,$tchewvec,type='l',lwd=2,xlab='Mass (kg)',ylab='Bite chew time (s)')
# dev.off()
# """
