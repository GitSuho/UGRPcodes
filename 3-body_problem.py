import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, diff
#normalization geodtion
def Norm(_list):
    result = 0
    for i in _list:
        result += i**2
    return result**(1/2)

#inner product geodtion
def DotPro(_lis1, _lis2):
    _len = len(_lis1)
    if not(_len == len(_lis2)):
        print("ERROR:2889638hfiu2979378egf");exit()
    result = 0
    for i in range(_len):
        result += _lis1[i]*_lis2[i]
    return result

#initial condition
x_initial = np.array([[-1/2, -1/3], [1/2, -1/3], [0, 2/3]])
v_initial = np.array([[1/3, -2/3], [-1/3, 2/3], [-1, 0]])
m = np.array([1, 1, 1]) ; G = 9.8
lis = [[],[],[]]
for i in range(3):
    for j in range(2):
        lis[i].append(x_initial[i][j])
for i in range(3):
    for j in range(2):  
        lis[i].append(v_initial[i][j])
initial_val = np.array(lis)

#variables
xv_symbol = symbols('xv_symbol:12')
xv_symbol = np.array(xv_symbol).reshape((3, 4))


#energys
E = 0
for i in range(3):
    j = (i+1)%3
    E += (1/2)*m[i]*Norm(v_initial[i])**2 \
        - G*m[i]*m[j]/Norm(x_initial[i]-x_initial[j])
T = E
for i in range(3):
    j = (i+1)%3
    T += G*m[i]*m[j]\
        / Norm([xv_symbol[i][0] - xv_symbol[j][0] ,
                xv_symbol[i][1] - xv_symbol[j][1]])
T_diff = 0
for i in range(3):
    j = (i+1)%3
    T_diff += G*m[i]*m[j]*DotPro(xv_symbol[i][0:2]-xv_symbol[j][0:2], xv_symbol[i][2:4]-xv_symbol[j][2:4])\
                             / Norm(xv_symbol[i][0:2]-xv_symbol[j][0:2])**2

#geodesic equation
F = [[0, 0] for _ in range(3)]
for i in range(3):
    for j in range(2):
        F[i][j] = T_diff/T*xv_symbol[i,j+2]
        for k in range(3):
            for l in range(2):
                F[i][j] = F[i][j] - 1/T*diff(T, xv_symbol[i][j])*xv_symbol[k][l+2] \
                                  + 1/(2*T)*diff(T, xv_symbol[i][j])*xv_symbol[k][l+2]**2

print(F)
exit()   

def geod(_X, t):
    _var = {xv_symbol[0, 0]: _X[0, 0], xv_symbol[0, 1]: _X[0, 1], xv_symbol[0, 2]: _X[0, 2], xv_symbol[0, 3]: _X[0, 3],
            xv_symbol[1, 0]: _X[1, 0], xv_symbol[1, 1]: _X[1, 1], xv_symbol[1, 2]: _X[1, 2], xv_symbol[1, 3]: _X[1, 3],
            xv_symbol[2, 0]: _X[2, 0], xv_symbol[2, 1]: _X[2, 1], xv_symbol[2, 2]: _X[2, 2], xv_symbol[2, 3]: _X[2, 3]}
    return np.array([[_X[0, 2], _X[0, 3], F[0][0].subs(_var), F[0][1].subs(_var)],
                     [_X[1, 2], _X[1, 3], F[1][0].subs(_var), F[1][1].subs(_var)],
                     [_X[2, 2], _X[2, 3], F[2][0].subs(_var), F[2][1].subs(_var)]])

#graph setting
t_final = 5
divide_num = t_final*(100)
t = np.linspace(0, t_final, divide_num)
plt.title("3-body gravitational field")

#Runge-Kutta method
dt = t[1] - t[0]


nt = len(t)
X  = np.array([initial_val for _ in range(nt) ])
for i in range(nt-1):
    k1 = geod(( X[i]             ), ( t[i]        ))
    k2 = geod(( X[i] + dt/2. * k1), ( t[i] + dt/2.))
    k3 = geod(( X[i] + dt/2. * k2), ( t[i] + dt/2.))
    k4 = geod(( X[i] + dt    * k3), ( t[i] + dt   ))
    X[i+1] = X[i] + (dt / 6.) * ((k1) + (2. * k2) + (2. * k3) + (k4))

#plot the graph in live
    for j in range(3):
        plt.plot(X[0:i+1,j, 0], X[0:i+1,j ,1], linewidth = 2)
    plt.pause(0.1)
plt.show()
print("end")