# Short program to calculate t_crunch

import numpy as np

# Units
ly = 9.46e15 # m
M_sun = 2e30 # kg

# Numerical setup
M_tot = 1000*M_sun
R = 20*ly
rho = M_tot/(4*np.pi*R**3/3.0)
G = 6.678e-11

t_crunch = np.sqrt(3*np.pi/(32*G*rho))
sec_to_years = 1./(365*24*60*60)

print "t_crunch = %e years" % (t_crunch*sec_to_years)
print "t_crunch = %.2f million years" % (t_crunch*sec_to_years*1e-6)