from numpy import *
from math import *
from matplotlib.pyplot import *
from Solver_functions import *
from My_user_interface import *




d =  0.3                # diameter of body m
#rho_air = 1.2          # density of air kg/m3
#rho_water = 1000       # density of water kg/m3
rho_body = 8910        # Copper ball density 8910
                        # human body density 1003    
V = (4./3)*pi*(d/2)**3  # Volume of sphere
g = 9.81                # gravitational constan
mu = 1.8*10**(-5)       # dynamic viscosity of fluid
                        # air = 1.8*10**(-5) Pa
                        # water = 8.9*10**(-4) Pa
A = pi*(d/2)**2         # area    

Cd = 0.45               # drag coefficient
    # Sphere = 0.45
    # semi-sphere = 0.42
    # Cube = 1.05
    # long Cylinder (long alnog vertical direction) = 0.82
    # Rocket = 0.75
    # man in Upright position = 1.0 - 1.3
    # flat plate perpendicular to flow = 1.3
    # streamlined body = 0.04
    

    
rho_air, v0, T, dt = read_command_line()
a_Sd = 3*pi*d*mu/rho_body*V
a_Qd = 0.5*Cd*(rho_air*A/rho_body*V)
b = g*(rho_air/rho_body - 1)
Re_v = rho_air*d/mu  
v, t ,R ,N = solver_vertical_motion(Re_v, v0, a_Sd, a_Qd, b, T, dt)

Fd = zeros(N+1)
for n in range(0,N):
    if R[n]<1 :
        Fd[n] = 3*pi*d*mu*v[n] #Fd_s
    else:
        Fd[n] = -0.5*Cd*A*abs(v[n])*v[n]  #Fd_q
        
Fg = -1.0*V*rho_body
Fb = rho_air*g*V
print Fd

figure(2)
plot(R,Fd,'--')
show()
