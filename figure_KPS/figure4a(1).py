import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import math
from scipy.optimize import fsolve


def Cal_DiscreteCurvature(d, theta):
    return 2 * math.cos(theta) / d

def Dis_2points(p1, p2):
    return ((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2) ** 0.5

def Cal_ContinuousCurvature(_t, _a, _b):
    return _a * _b / (((_a * math.sin(_t)) ** 2 + ((_b * math.cos(_t)) ** 2)) ** 1.5)

a = 5
b = 5

dis_curvature_points = []
con_curvature_points = []
plot_interval = 0.1

div_num = 1
for i in range(div_num):
    ### given angle
    t = 2 * math.pi * i / div_num

    ### given point
    pi1 = [a * math.cos(t), b * math.sin(t)]


    
    ### find points
    x, y = sp.symbols('x y')
    eq1 = sp.Eq((x - pi1[0]) ** 2 + (y - pi1[1]) ** 2, plot_interval ** 2)
    eq2 = sp.Eq(x ** 2 / a ** 2 + y ** 2 / b ** 2, 1)
    solutions = sp.solve((eq1, eq2), (x, y))
    

    ### back and front points
    pi0 = solutions[0]
    pi2 = solutions[1]
    
    ### find theta by using cosine law
    p_middle = np.array([(pi0[0] + pi2[0]) / 2, (pi0[1] + pi2[1]) / 2])
    dis_A = Dis_2points(pi1, p_middle)
    dis_B = Dis_2points(pi1, pi0)
    dis_C = Dis_2points(pi0, p_middle)
    theta = math.acos((dis_A ** 2 + dis_B ** 2 - dis_C ** 2) / (2 * dis_A * dis_B))

    ### store curvatures
    dis_curvature = Cal_DiscreteCurvature(plot_interval, theta)
    con_curvature = Cal_ContinuousCurvature(t, a, b)
    dis_curvature_points.append(dis_curvature)
    con_curvature_points.append(con_curvature)

plt.scatter(dis_curvature_points, con_curvature_points)
print(dis_curvature_points)
print(con_curvature_points)
plt.plot([0, 5], [0, 5])
plt.show()

