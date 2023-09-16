import matplotlib.pyplot as plt
import numpy as np

def Norm(lis):
    sum = 0
    for i in lis:
        sum += i**2
    return sum**0.5

E = 1/2
ecc = 0.75

M = 9
m = 1
G = 1
v_init = [0, (E/(-m/2+m/(ecc+1)))**0.5]
x_init = [(ecc+1)*G*M/(v_init[1]**2), 0]

rp = x_init[0]
ra = rp/((1-ecc)/(1+ecc))

a =  (ra + rp)/2
b = (1-ecc**2)**0.5*a



#analytic solution
def R_anly(theta):
    l = x_init[0]*v_init[1]-x_init[1]*v_init[0]
    k = G*M
    A = 1/Norm(x_init) - k/(m*l**2)
    return (m*l**2/k) / (1 + (m*l**2*A/k)*np.cos(theta))




#Error function
def Error(theta1, plot_interval):
    x1, y1 = R_anly(theta1)*np.cos(theta1), R_anly(theta1)*(np.sin(theta1))
    
    delta_x = plot_interval/((1 + (x1**2*b**4)/(y1**2*a**4))**0.5)
    delta_y = delta_x*Norm([x1*b**2/(y1*a**2)])
    
    # print(Norm([delta_x, delta_y]))
    
    
    if y1 > 0:
        x2 = x1 - delta_x
    else:
        x2 = x1 + delta_x 
        
    if x1 > 0:
        y2 = y1 + delta_x*Norm([x1*b**2/(y1*a**2)])        
    else:
        y2 = y1 - delta_x*Norm([x1*b**2/(y1*a**2)])
    # y2 = b**2/y1 - x1*b**2/(y1*a**2)*x2

    if y2 > 0:
        theta2 = np.arccos(x2/Norm([x2, y2]))
    else:
        theta2 = 2*np.pi - np.arccos(x2/Norm([x2, y2]))
        
    plt.plot([x1, x2], [y1, y2], 'r-')
    print(Norm([x1-x2, y1-y2]))

    return Norm([x2, y2]) - R_anly(theta2)    


#plot an analytic Solution's trajectory
count = 10000
orbit_x = []
orbit_y = []
theta_arr = [ i*2*np.pi/count for i in range(count)  ]
for i in theta_arr:
    rad = R_anly(i)
    orbit_x.append(rad*np.cos(i))
    orbit_y.append(rad*np.sin(i))
    
plt.plot(orbit_x, orbit_y, 'g-')
plt.plot(0, 0, 'g')



plot_interval = 0.01
theta_lis = [ i*(2*np.pi-0.01)/1000 for i in range(1000)]
Err_lis = []
for theta in theta_lis:
    Err_lis.append(Error(theta, plot_interval))

plt.show()

plt.plot(theta_lis, Err_lis)
plt.show()    