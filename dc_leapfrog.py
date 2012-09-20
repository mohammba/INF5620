from numpy import *
from matplotlib.pyplot import *
from Solver_functions import *
from My_user_interface import *


def a(t):
        return 1

def b(t):
        return 1

def exact_solution(t):
        return 1-exp(-t)

# solution of numerical of fist order and second order in dt
numerical_2th, t = Leapfrog_secondOrder( a, b, I=0, T=10, dt=0.5)
numerical_1th, t = Leapfrog_firstOrder( a, b, I=0, T=10, dt=0.5)
# plot exact solution and numeric for compare the result
"""
figure(1)
plot(t,numerical_2th,'b')
plot(t, exact_solution(t), 'r')
legend(['numerical_2th','exact'])
xlabel('time')
show()

figure(2)
plot(t,numerical_1th,'b')
plot(t, exact_solution(t), 'r')
legend(['numerical_1th','exact'])
xlabel('time')
show()
"""
# copy from page.32 main_decay.pdf
dt_values = linspace(0.00001, 0.01, 10)
print len(dt_values)
r = {}
# estimated convergence rates
E_values = []
for dt in dt_values:
	numerical,t =  Leapfrog_secondOrder(a, b, I=0, T=10, dt=dt)
        E = numerical - exact_solution(t)
        E = sqrt(dt*sum(e**2))
        E_values.append(E)
	# Compute convergence rates
	m = len(dt_values)
	r = [log(E_values[i-1]/E_values[i])/
	log(dt_values[i-1]/dt_values[i]) for i in range(1, m, 1)]

print r
