
# coding: utf-8
import numpy as np

posteriors_file = "/Users/cv/exonailer/results/228735255_full_K2_white_quadratic/posterior_parameters.dat"

posteriors = {}
with open(posteriors_file,"r") as f:
    for line in f:
        if line[0] == "#" or line[0] == " ":
            pass
        else:
            ls = line.split()
            #print type(ls[0])
            posteriors[ls[0]] = np.float(ls[1])
            posteriors[ls[0] + "_up"] = np.float(ls[2])
            posteriors[ls[0] + "_down"] = np.float(ls[3])

K = posteriors["K"]*1000. #Semi-amplitude m/s
e = posteriors["ecc"]
P = posteriors["P"]*24.*60.*60. #Period in seconds
i = posteriors["inc"]*2*np.pi/360. #inclination in radians
G = 6.67408e-11 # Gravitational constant in m^3 kg^-1 s^-2
M_sun = 1.98855e30 # Sun Mass in kg
M_star = 1.005*M_sun #Star mass in kg

M_jup = 1.89813e27 #Jupiter mass in kg

M_star_up = 0.021*M_sun
M_star_down = 0.020*M_sun

M_p = (K*np.sqrt(1.-e**2.)/np.sin(i))*(P/(2.*np.pi*G))**(1./3.)*(M_star**(2./3.))/M_jup # Planet mass in jupiter mass

K_up = posteriors["K_up"]*1000.
e_up = posteriors["ecc_up"]
P_up = posteriors["P_up"]*24.*60.*60.
i_up = posteriors["inc_up"]*2*np.pi/360.

#using dK as (delta K / K )
dK_up = K_up/K
de_up = (1./2.)*(1/(1-e**2))*e_up
di_up = (np.cos(i)/np.sin(i))*i_up
dP_up = (1./3.)*(P_up/P)
dM_star_up = (2./3.)*(M_star_up/M_star)

M_p_up = M_p * np.sqrt(dK_up**2+de_up**2+di_up**2+dP_up**2+dM_star_up**2)


K_down = posteriors["K_down"]*1000.
e_down = posteriors["ecc_down"]
P_down = posteriors["P_down"]*24.*60.*60.
i_down = posteriors["inc_down"]*2*np.pi/360.

dK_down = K_down/K
de_down = (1./2.)*(1/(1-e**2))*e_down
di_down = (np.cos(i)/np.sin(i))*i_down
dP_down = (1./3.)*(P_down/P)
dM_star_down = (2./3.)*(M_star_down/M_star)

M_p_down = M_p * np.sqrt(dK_down**2+de_down**2+di_down**2+dP_down**2+dM_star_down**2)


p = posteriors["p"]
p_up = posteriors["p_up"]
p_down = posteriors["p_down"]


R_sun = 6.95700e8 # Sun radius in m
R_jup =  7.1492e7 # Jupiter radius in m

R_star = 0.987 * R_sun 
R_star_err = 0.011 * R_sun

R_p = R_star * p / R_jup # Planet radius in Jupiter radius


dp_up = p_up/p
dR_star = R_star_err / R_star

R_p_up = R_p * np.sqrt(dp_up**2+dR_star**2)

dp_down = p_down/p

R_p_down = R_p * np.sqrt(dp_down**2+dR_star**2)

a = np.float(posteriors["a"])
a_up = np.float(posteriors["a_up"])
a_down = np.float(posteriors["a_down"])

rho_jup = 1.33e3 #Jupiter density in kg m^-3


rho_p = M_p / ( R_p**3)


rho_sun = 1.408e3 #Sun density in kg m^-3

#check small term
#(p**3/(rho_p*rho_jup))/rho_sun

rho_star = ((((3*np.pi)/(G*P**2))*a**3))/1000.#-(p**3/(rho_p*rho_jup))/rho_sun

dP_up = 2*P_up/P
da_up = 3*a_up/a

rho_star_up = rho_star * np.sqrt(dP_up**2+da_up**2)

dP_down = 2*P_down/P
da_down = 3*a_down/a

rho_star_down = rho_star * np.sqrt(dP_down**2+da_down**2)