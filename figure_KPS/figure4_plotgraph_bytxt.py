import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import math
import matplotlib.colors as mcl
from matplotlib.colors import LinearSegmentedColormap

"""
def Cal_DiscreteCurvature(d, theta):
    return 2 * math.cos(theta) / d


def Dis_2points(p1, p2):
    return ((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2) ** 0.5


def Cal_ContinuousCurvature(_t, _a, _b):
    return _a * _b / (((_a * math.sin(_t)) ** 2 + ((_b * math.cos(_t)) ** 2)) ** 1.5)


def Cal_ContinuousCurvature_Rate(t, a, b):
    return - (3 * a * b * (a ** 2 - b ** 2) * math.cos(t) * math.sin(t)) / (b ** 2 * math.cos(t) ** 2 + a ** 2 * math.sin(t) ** 2) ** 2.5


ab_list = [ [4 , 2], [3.9 , 8/3.9], [3.8 , 8/3.8 ]]


con_curvature_points = [[],[],[]]
a_errors = [[],[],[]]
thetas = [ math.pi / 2 * i / 100 for i in range(100)] 


file_list = ["CurPerAbsErr_e0.866_a4b2_PlotInt0.001.txt", "CurPerAbsErr_e0.850_a3.9b2.0512820512820515_PlotInt0.001.txt", "CurPerAbsErr_e0.832_a3.8b2.1052631578947367_PlotInt0.001.txt"]
for i in range(3):
    with open(file_list[i], 'r') as rf:
        foo = rf.readlines()
        for line in foo:
            var = line.strip().split()
            con_curvature_points[i].append(float(var[0]))
            a_errors[i].append(float(var[1]))



plt.scatter(con_curvature_points[0], a_errors[0], c=thetas, cmap='Reds_r', label = "e = 0.866" )
plt.scatter(con_curvature_points[1], a_errors[1], c=thetas, cmap='Greens_r', label = "e = 0.850")
plt.scatter(con_curvature_points[2], a_errors[2], c=thetas, cmap='Blues_r', label = "e = 0.832" )


plt.xlabel("Curvature (m^-1)", fontsize=11)
plt.ylabel("Absolute error of curvature and discrete curvature (m^-1)", fontsize=11)
plt.xlim(0.0, 1.2)
plt.legend(loc = "lower left")

# plt.colorbar().set_label("Angle (rad)", fontsize=11)

plt.grid()
current_values = plt.gca().get_yticks()
plt.gca().set_yticklabels(['{:.6f}'.format(x) for x in current_values])


plt.show()
"""

ab_list = [ [4 , 2], [3.9 , 8/3.9], [3.8 , 8/3.8 ]] #"e = 0.866""e = 0.850""e = 0.832"


