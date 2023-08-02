import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, diff
#normalization function
def Norm(_x1, _x2):
    return (_x1**2+_x2**2)**(1/2)

#initial condition
x_initial = np.array([1, 0])
v_initial = np.array([0, 1])
X0 = np.array([x_initial[0], x_initial[1], v_initial[0], v_initial[1]])
M = 1   ; m = 1 ; G = 9.8 

#variables
x1, x2, v1, v2 = symbols('x1 x2 v1 v2')

#energys
E = (1/2)*m*Norm(v_initial[0], v_initial[1])**2 \
    -  G*M*m/Norm(x_initial[0], x_initial[1])
T = E + G*M*m/Norm(x1, x2)

#geodesic equation
F1 = -(diff(T, x1)*v1 + diff(T, x2)*v2)*v1 \
    + 1/(2*T)*diff(T, x1)*Norm(v1, v2)**2
F2 = -(diff(T, x1)*v1 + diff(T, x2)*v2)*v2 \
    + 1/(2*T)*diff(T, x2)*Norm(v1, v2)**2
def geod(X, t):
    return np.array([X[2], X[3], F1.subs({x1:X[0], x2:X[1], v1:X[2], v2:X[3]}),\
        F2.subs({x1:X[0], x2:X[1], v1:X[2], v2:X[3]})])

#garph setting
t_final = 5
divide_num = t_final*(3*10**2)
t = np.linspace(0, t_final, divide_num)

#fourth order Runge-Kutta method
def RK4(func, X0, t):
  dt = t[1] - t[0]
  nt = len(t)
  X  = np.zeros([nt, len(X0)])
  X[0] = X0
  for i in range(nt-1):
    k1 = func(( X[i]             ), ( t[i]        ))
    k2 = func(( X[i] + dt/2. * k1), ( t[i] + dt/2.))
    k3 = func(( X[i] + dt/2. * k2), ( t[i] + dt/2.))
    k4 = func(( X[i] + dt    * k3), ( t[i] + dt   ))
    X[i+1] = X[i] + (dt / 6.) * ((k1) + (2. * k2) + (2. * k3) + (k4))
  return X

RK4_result = RK4(geod, X0, t)

#show 2-Dimension graph
plt.plot(RK4_result[:,0], RK4_result[:,1], linewidth = 2)
plt.title("1-body gravitational field")
plt.show()