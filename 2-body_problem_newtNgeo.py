import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, diff
from matplotlib.animation import FuncAnimation
#normalization function
def Norm(_X):
    sum = 0
    for i in range(len(_X)):
        sum += _X[i] ** 2
    return sum**(1/2)

#initial condition
x_initial = np.array([5, 0])
v_initial = np.array([0, 1.5])
X0 = np.array([x_initial[0], x_initial[1], v_initial[0], v_initial[1]])
M = 10   ; m = 1 ; G = 10
_type = "" ; bar = 0
foo = input("Choose geodesic or newtonian. Write g or n : ")
if (foo == "g" or foo == "G"):
    bar = 1
    _type = "geodesic"
elif(foo == "n" or foo == "N"):
    bar = 0
    _type = "newtonian"
else:
    print("please write the correct function type")
    exit()


#variables
x1, x2, v1, v2 = symbols('x1 x2 v1 v2')

#energys
E = (1/2)*m*Norm([v_initial[0], v_initial[1]])**2 \
    -  G*M*m/Norm([x_initial[0], x_initial[1]])
T = E + G*M*m/Norm([x1, x2])

#geodesic equation
F1 = -(diff(T, x1)*v1 + diff(T, x2)*v2)*v1 \
    + 1/(2*T)*diff(T, x1)*Norm([v1, v2])**2
F2 = -(diff(T, x1)*v1 + diff(T, x2)*v2)*v2 \
    + 1/(2*T)*diff(T, x2)*Norm([v1, v2])**2
def geod(X, t):
    return np.array([X[2], X[3], F1.subs({x1:X[0], x2:X[1], v1:X[2], v2:X[3]}),\
        F2.subs({x1:X[0], x2:X[1], v1:X[2], v2:X[3]})])
    
#newtonian equation
def F1n (X):
    return -G*M/Norm(X[0:2])**3*X[0]
def F2n (X):
    return -G*M/Norm(X[0:2])**3*X[1]
def newtonian(X, t):
    return np.array([X[2], X[3], F1n(X), F2n(X)])

#garph setting
t_final = 4
divide_num = (10**3)
t = np.linspace(0, t_final, divide_num)

#fourth order Runge-Kutta method
dt = t[1] - t[0]
nt = len(t)
X  = np.zeros([nt, len(X0)])
X[0] = X0
plt.title(f"2B gravi field <{_type}> (final : {t_final}, div : {divide_num}, \n [{X0[0]}, {X0[1]}, {X0[2]}, {X0[3]}], M={M}, m={m}, G={G})")
# plt.xlim(-1.5, 1.5 )
# plt.ylim(-1.5, 1.5)
for i in range(nt-1):
    if (bar): 
        k1 = newtonian(( X[i]             ), ( t[i]        ))
        k2 = newtonian(( X[i] + dt/2. * k1), ( t[i] + dt/2.))
        k3 = newtonian(( X[i] + dt/2. * k2), ( t[i] + dt/2.))
        k4 = newtonian(( X[i] + dt    * k3), ( t[i] + dt   ))
        X[i+1] = X[i] + (dt / 6.) * ((k1) + (2. * k2) + (2. * k3) + (k4))
    else:
        k1 = geod(( X[i]             ), ( t[i]        ))
        k2 = geod(( X[i] + dt/2. * k1), ( t[i] + dt/2.))
        k3 = geod(( X[i] + dt/2. * k2), ( t[i] + dt/2.))
        k4 = geod(( X[i] + dt    * k3), ( t[i] + dt   ))
        X[i+1] = X[i] + (dt / 6.) * ((k1) + (2. * k2) + (2. * k3) + (k4))
    plt.plot(X[0:i+1,0]*(m/(M+m)), X[0:i+1,1]/2, linewidth=2, color='b')
    plt.plot(-X[0:i+1,0]*(M/(M+m)), -X[0:i+1,1]/2, linewidth=2, color='r')
    plt.pause(0.1)
plt.show()

print("End")