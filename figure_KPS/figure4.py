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


ab_list = [ [4 , 2], [3.9 , 8/3.9], [3.8 , 8/3.8 ]]


dis_curvature_points = [[],[],[]]
con_curvature_points = [[],[],[]]
"""
r_errors = [[],[],[]]
"""
a_errors = [[],[],[]]


for ab in range(3):
    a ,b = ab_list[ab][0] , ab_list[ab][1]
    
    
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

        '''
        sol1x = 1 / 2 * (c + d - (-c ** 2 + 2 * c * d - d ** 2 + 2 * f ** 2) ** (1 / 2))
        sol1y = (b * (2 * a ** 2 - 2 * c * d - f ** 2 + c * (-c ** 2 + 2 * c * d - d ** 2 + 2 * f ** 2) ** (1 / 2) + d * (-c ** 2 + 2 * c * d - d ** 2 + 2 * f ** 2) ** (1 / 2)) ** (1 / 2))/(2 ** (1 / 2) * a)
        sol1 = [sol1x, sol1y]
        sol2x = 1 / 2 * (c + d + (-c ** 2 + 2 * c * d - d ** 2 + 2 * f ** 2) ** (1 / 2))
        sol2y = (b * (2 * a ** 2 - 2 * c * d - f ** 2 - c * (-c ** 2 + 2 * c * d - d ** 2 + 2 * f ** 2) ** (1 / 2) - d * (-c ** 2 + 2 * c * d - d ** 2 + 2 * f ** 2) ** (1 / 2)) ** (1 / 2))/(2 ** (1 / 2) * a)
        sol2 = [sol2x, sol2y]
        solutions = [sol1, sol2]
        print(solutions)
        '''
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
        dis_curvature_points[ab].append(dis_curvature)
        con_curvature_points[ab].append(con_curvature)
        rate_curvature_points.append(rate_curvature)

        """
        r_error = (dis_curvature - con_curvature) / con_curvature * 100
        r_errors[ab].append(r_error)
        """
        a_error = dis_curvature - con_curvature
        a_errors[ab].append(a_error)
    '''
    plt.scatter(dis_curvature_points, con_curvature_points)
    plt.plot([0, 5], [0, 5])
    plt.show()
    plt.plot(thetas, dis_curvature_points)
    plt.show()
    plt.plot(thetas, con_curvature_points)
    plt.show()
    '''


# print(thetas)
# print(dis_curvature_points)
# print(con_curvature_points)
# print(r_errors)
# print(a_errors)

###########################################################################

"""
plt.plot([0, 1.2], [0, 1.2], color = 'black')


plt.scatter(con_curvature_points[0], dis_curvature_points[0] , color = 'red')
plt.scatter(con_curvature_points[1], dis_curvature_points[1] , color = 'green')
plt.scatter(con_curvature_points[2], dis_curvature_points[2] , color = 'blue')


plt.xlabel("Curvature (rad/m)", fontsize=11)
plt.ylabel("Discrete curvature (rad/m)", fontsize=11)
plt.xlim(0.0, 1.2)
plt.ylim(0.0, 1.2)
plt.grid()
plt.show()
"""
###########################################################################

"""
plt.scatter(con_curvature_points, r_errors, c=thetas, cmap='Blues_r')
plt.xlabel("Curvature (rad/m)", fontsize=11)
plt.ylabel("Relative error of curvature and discrete curvature (%)", fontsize=11)
plt.xlim(0.0, 1.2)
plt.ylim(-0.0020, 0.0020)
plt.colorbar().set_label("Angle (rad)", fontsize=11)
plt.grid()
plt.show()
"""

###########################################################################
plt.scatter(con_curvature_points[0], a_errors[0], c=thetas, cmap='Reds_r' , label = "e = 0.866")
plt.scatter(con_curvature_points[1], a_errors[1], c=thetas, cmap='Greens_r', label = "e = 0.850")
plt.scatter(con_curvature_points[2], a_errors[2], c=thetas, cmap='Blues_r' , label = "e = 0.832")


plt.xlabel("Curvature (rad/m)", fontsize=11)
plt.ylabel("Absolute error of curvature and discrete curvature (rad/m)", fontsize=11)
plt.xlim(0.0, 1.2)
#plt.ylim(-0.00002, 0.00001)
plt.legend(loc = "lower left")

plt.colorbar().set_label("Angle (rad)", fontsize=11)

plt.grid()
current_values = plt.gca().get_yticks()
plt.gca().set_yticklabels(['{:.6f}'.format(x) for x in current_values])


plt.show()