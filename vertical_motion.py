from numpy import *
from math import *
from matplotlib.pyplot import *
from Solver_functions import *
from My_user_interface import *
import nose.tools as nt 

d =  0.5                # diameter of body m
#rho_air = 1.2          # density of air kg/m3
#rho_water = 1000       # density of water kg/m3
rho_body = 1003         # Copper ball density 8910
                        # human body density 1003    
#V = (4./3)*pi*(d/2)**3  # Volume of sphere
V = 0.08                # Volume of human body
g = 9.81                # gravitational constan
mu = 1.8*10**(-5)       # dynamic viscosity of fluid
                        # air = 1.8*10**(-5) Pa
                        # water = 8.9*10**(-4) Pa
#A = pi*(d/2)**2         # area    
A = 0.9                 # Area of human body
#Cd = 0.45               # drag coefficient
    # Sphere = 0.45
    # semi-sphere = 0.42
    # Cube = 1.05
    # long Cylinder (long alnog vertical direction) = 0.82
    # Rocket = 0.75
    # man in Upright position = 1.0 - 1.3
    # flat plate perpendicular to flow = 1.3
    # streamlined body = 0.04
Cd = 1.2
m = 80    

# Read variables from command line        
rho_air, v0, T, dt = read_command_line()

# Calculate coefficients of differential equation
a_Sd = 3*pi*d*mu/rho_body*V
a_Qd = 0.5*Cd*(rho_air*A/rho_body*V)
b = g*(rho_air/rho_body - 1)
Re_v = rho_air*d/mu  

# numeric solution output
v, t ,Re ,N = solver_vertical_motion(Re_v, v0, a_Sd, a_Qd, b, T, dt)

# exact solution output
w = exact_solution_vertical_motion(Re, v0, a_Sd, a_Qd, b, T, dt)

# Nose-test 
diff = abs(w - v).max()
nt.assert_almost_equal(diff, 0, delta=1E-1)

# Maximum error between exact and numeric solution
print 'diff = %0.10f' % diff

#Plot velocity by numeric_solution exact_solution agains time-shots number
figure(1) 
plot(t,v,'r-',t, w, 'bo')
show()


