# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 16:10:01 2012

@author: puya
"""
from numpy import *
from math import *
def solver_vertical_motion(Re_v, v0, a_Sd, a_Qd, b, T, dt):
    """ Solve v' = -a*|v|*v + b, v(0)=v0, for t in (0,T] with steps of dt. 
        a = 0.5*Cd*(rho*a/rho_b*V)
        b = g*(rho/rho_b - 1)
    """
    N = int(T/dt)               # no of time intervals
    T = N*dt                    # adjust T to fit time step dt
    Re = zeros(N+1)              # values of Re varies with v0
    v = zeros(N+1)              # array of u[n] values
    t = linspace(0, T, N+1)     # time mesh
    v[0] = v0    
    Re[0] = Re_v*abs(v0)               # assign initial conditiol	
    for n in range(0, N):
        #Re_v = rho_air*d/mu 
        if Re[n] < 1 :
            v[n+1] = ((1 - dt*0.5*a_Sd)*v[n] + dt*b)/(1 + dt*0.5*a_Sd)
             # Stokes drag solver
        else:
            v[n+1] = (v[n]+dt*b)/(1+dt*a_Qd*abs(v[n]))
            # Quadratic drag solver
        Re[n+1] = Re_v*abs(v[n+1])
    return v, t, Re, N  




def exact_solution_vertical_motion(Re, v0, a_Sd, a_Qd, b, T, dt):

	N = int(T/dt)
        T = N*dt
	w = []
	g= 9.81
	for n in range(0,N+1):
		
		if Re[n] < 1 :
			w.append((b/a_Qd)*(1-exp(-1*a_Qd*n)))
		else:
			w.append(-1*sqrt(abs(b/a_Qd))*tanh(sqrt(abs(a_Qd*b))*n))
				
	return w

def Leapfrog_secondOrder(a, b, I, T, dt):
        """Solves u'=-a(t)*u+b(t), u(0)=I"""

	dt = float(dt)
	# avoid integer division
	N = int(round(T/dt))
	# no of time intervals
	T = N*dt
	# adjust T to fit time step dt
	u = zeros(N+1)
	# array of u[n] values
	t = linspace(0, T, N+1) # time mesh
	u[0] = I
	# assign initial condition
	for n in range(0, N):
		# n=0,1,...,N-1
 		u[n+1] = u[n-1] + 2*dt*(-a(t[n])*u[n]+b(t[n]))
	return u, t

def Leapfrog_firstOrder(a, b, I, T, dt):
        """Solves u'=-a(t)*u+b(t), u(0)=I"""

	dt = float(dt)
	# avoid integer division
	N = int(round(T/dt))
	# no of time intervals
	T = N*dt
	# adjust T to fit time step dt
	u = zeros(N+1)
	# array of u[n] values
	t = linspace(0, T, N+1) # time mesh
	u[0] = I
	# assign initial condition
	for n in range(0, N):
		# n=0,1,...,N-1
 		u[n+1] = u[n] + dt*(-a(t[n])*u[n]+b(t[n]))
	return u, t




