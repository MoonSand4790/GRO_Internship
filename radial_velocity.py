import numpy as np

# Radial velocity due to Doppler shift
def rel_vel(f_obs):
    del_f = f_int - f_obs
    vr = c*float(del_f/f_int)
    return vr

# Observed frequencies for Neutral Hydrogen Cloud
f_obs = float(input('Observed frequency of radio signal (MHz):'))

# Intrinsic 21cm Frequency
f_int = 1420.405751 #(MHz)

c = 3e8 # m/s
    
res = np.round(rel_vel(f_obs)/1000,2)

print('Radial velocity of the gas cloud is:',res,'km/s')
