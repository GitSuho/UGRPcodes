import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import math
import matplotlib.colors as mcl
from matplotlib.colors import LinearSegmentedColormap


def Cal_DiscreteCurvature(d, theta):
    return 2 * math.cos(theta) / d


def Dis_2points(p1, p2):
    return ((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2) ** 0.5


def Cal_ContinuousCurvature(_t, _a, _b):
    return _a * _b / (((_a * math.sin(_t)) ** 2 + ((_b * math.cos(_t)) ** 2)) ** 1.5)


def Cal_ContinuousCurvature_Rate(t, a, b):
    return - (3 * a * b * (a ** 2 - b ** 2) * math.cos(t) * math.sin(t)) / (b ** 2 * math.cos(t) ** 2 + a ** 2 * math.sin(t) ** 2) ** 2.5


ab_list = [ [4 , 2], [3.9 , 8/3.9], [3.8 , 8/3.8 ]] #"e = 0.866""e = 0.850""e = 0.832"
ecc_list = ["0.866", "0.850", "0.832"]


for ab in range(3):
    a ,b = ab_list[ab][0] , ab_list[ab][1]
    
    file_name = f"CurPerAbsErr_e{ecc_list[ab]}_a{a}b{b}_PlotInt0.001.txt"
    with open(file_name,"w") as wf:

        
        
        rate_curvature_points = []
        thetas = []

        plot_interval = 0.01
        f = plot_interval

        div_num = 100
        for i in range(div_num):
            print(f'{ab}, {i}')
            # given angle
            t = math.pi / 2 * i / div_num
            thetas.append(t)

            # given point
            pi1 = [a * math.cos(t), b * math.sin(t)]
            c = pi1[0]
            d = pi1[1]

            # find points
            x, y = sp.symbols('x y')

            eq1 = sp.Eq((x - pi1[0]) ** 2 + (y - pi1[1]) ** 2, plot_interval ** 2)
            eq2 = sp.Eq(x ** 2 / a ** 2 + y ** 2 / b ** 2, 1)
            solutions = sp.solve((eq1, eq2), (x, y))

            # back and front points
            pi0 = solutions[0]
            pi2 = solutions[1]

            # find theta by using cosine law
            p_middle = np.array([(pi0[0] + pi2[0]) / 2, (pi0[1] + pi2[1]) / 2])
            dis_A = Dis_2points(pi1, p_middle)
            dis_B = Dis_2points(pi1, pi0)
            dis_C = Dis_2points(pi0, p_middle)

            # theta = float(np.arccos((dis_A**2 + dis_B**2 - dis_C**2)/(2*dis_A*dis_B)))
            theta = math.acos((dis_A ** 2 + dis_B ** 2 - dis_C ** 2) / (2 * dis_A * dis_B))

            # store curvatures
            dis_curvature = Cal_DiscreteCurvature(plot_interval, theta)
            con_curvature = Cal_ContinuousCurvature(t, a, b)
            rate_curvature = Cal_ContinuousCurvature_Rate(t, a, b)

            a_error = dis_curvature - con_curvature
            
            
            wf.write(str(con_curvature) + " " + str(a_error) + '\n')
    