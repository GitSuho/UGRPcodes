import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, diff
from matplotlib.animation import FuncAnimation
#normalization function
def Norm(_x1, _x2):
    return (_x1**2+_x2**2)**(1/2)

#initial condition
x_initial = np.array([1, 0])
v_initial = np.array([0, 3])
X0 = np.array([x_initial[0], x_initial[1], v_initial[0], v_initial[1]])
M = 1   ; m = 1 ; G = 9.8 

#variables
x1, x2, v1, v2 = symbols('x1 x2 v1 v2')

#energys
E = (1/2)*m*Norm(v_initial[0], v_initial[1])**2 \
    -  G*M*m/Norm(x_initial[0], x_initial[1])
T = E + G*M*m/Norm(x1, x2)

# geodesic equation
F1 = (-G*M*m*(x1*v1+x2*v2)/(Norm(x1, x2)**3))/T*v1 \
    + 1/T*diff(T, x2)*v2*v1 + 1/(2*T)*diff(T, x1)*(v2)**2
F2 = (-G*M*m*(x1*v1+x2*v2)/(Norm(x1, x2)**3))/T*v2 \
    + 1/T*diff(T, x1)*v2*v1 + 1/(2*T)*diff(T, x2)*(v1)**2    



# F1 = -(diff(T, x1)*v1 + diff(T, x2)*v2)*v1 \
#     + 1/(2*T)*diff(T, x1)*Norm(v1, v2)**2
# F2 = -(diff(T, x1)*v1 + diff(T, x2)*v2)*v2 \
#     + 1/(2*T)*diff(T, x2)*Norm(v1, v2)**2
def geod(X, t):
    return np.array([X[2], X[3], F1.subs({x1:X[0], x2:X[1], v1:X[2], v2:X[3]}),\
        F2.subs({x1:X[0], x2:X[1], v1:X[2], v2:X[3]})])

#garph setting
t_final = 3
divide_num = t_final*(10**2)
t = np.linspace(0, t_final, divide_num)

#fourth order Runge-Kutta method
dt = t[1] - t[0]
nt = len(t)
X  = np.zeros([nt, len(X0)])
X[0] = X0
plt.title(f"1B gravi field (f_final : {t_final}, divide num : {divide_num})")
plt.xlim(-1.5, 1.5 )
plt.ylim(-1.5, 1.5)
j = 0
for i in range(nt-1):
    # if ((X[i][2] < 0.1 and X[i][2] > -0.1 ) and (X[i][3] < 0.1 and X[i][3] > -0.1 )):
    #     X[i][2] *= 10
    #     X[i][3] *= 10
    j += 1
    if (j > 100):
        print("working")
        j = 0
    
    k1 = geod(( X[i]             ), ( t[i]        ))
    k2 = geod(( X[i] + dt/2. * k1), ( t[i] + dt/2.))
    k3 = geod(( X[i] + dt/2. * k2), ( t[i] + dt/2.))
    k4 = geod(( X[i] + dt    * k3), ( t[i] + dt   ))
    X[i+1] = X[i] + (dt / 6.) * ((k1) + (2. * k2) + (2. * k3) + (k4))
    # if ((X[i+1][2] < 0.1 and X[i+1][2] > -0.1 ) or (X[i+1][3] < 0.1 and X[i+1][3] > -0.1 )):
    #     print("initial : " , X[i+1][2],X[i+1][3] )
    #     X[i+1][2] = X[i+1][2]/Norm(X[i+1][2],X[i+1][2] ) /10
    #     X[i+1][3] = X[i+1][3]/Norm(X[i+1][3], X[i+1][3]) /10
    #     print(X[i+1][2],X[i+1][3] )
    plt.plot(X[0:i+1,0], X[0:i+1,1], linewidth=2, color='b')
    plt.pause(0.1)
plt.show()

print("End")
