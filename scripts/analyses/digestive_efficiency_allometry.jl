
using herbforagingsim

#HERBIVORE
#Define mass of herbivore
mass = 100;
#Define tooth and gut type of herbivore
teeth = "bunodont"; # bunodont, acute/obtuse lophs, lophs and non-flat, lophs and flat
gut_type = "caecum"; # caecum, colon, non-rumen foregut, rumen foregut

edensity = 18.2

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
particle_mass, particle_length = mean_particle_mass(mass, gut_type); 

#particle density
d_p = 0.0012  # g/mm^3

gut_SA = ((maxgut/edensity)/d_p)/particle_length #mm^2

# Calculate surface area per gram of ingesta
ingesta_SA_per_gram = (6*particle_length^2)/particle_mass; #SA/particle * (particle/gram) mm^2/gram
ingesta_SA_per_kj = ingesta_SA_per_gram * edensity; #mm^2 per kJ

#Digestive efficiency
eta_d = 0.6;



gut_t = 0.5*maxgut;

r_a = (1/6)*(eta_d*d_p*edensity)*(particle_length/mrt)*(1/gut_t)

#Or rearrange to solve for eta_d
#Set r_a
r_a = 6.294*10^-13; #arithmetic mean
eta_d = (6*r_a*mrt)/(d_p*particle_length*edensity)*gut_t


b0_fmr = 0.047; #watts g^-0.75
fcost_watts = (b0_fmr*((mass*1000)^0.75)); #Cost in Watts
kJps = 0.001;
fcost_kJps = fcost_watts*kJps; #kJ per s

#Minimum r_a to make up for the metabolic cost
min_r_a = fcost_kJps/gut_SA

eta_d = (6*min_r_a*mrt)/(d_p*particle_length*edensity)*0.035
