import numpy as np
import matplotlib.pyplot as plt

# 초기 조건 설정
x0 = 0.0  # 초기 x 값
y0 = 0.0  # 초기 y 값
h = 2   # 스텝 크기
end_x = 2.0  # 종료 x 값

# 미분 방정식을 정의합니다. 이 경우, y' = x^2 입니다.
def f(x, y):
    return 2*x

# Runge-Kutta 4차수치해법을 구현합니다.
def runge_kutta_4th(x0, y0, h, end_x):
    x_values = [x0]
    y_values = [y0]

    while x_values[-1] < end_x:
        print("while")
        x = x_values[-1]
        y = y_values[-1]

        k1 = h * f(x, y)
        k2 = h * f(x + h/2, y + k1/2)
        k3 = h * f(x + h/2, y + k2/2)
        k4 = h * f(x + h, y + k3)
        
        
        

        x_new = x + h
        y_new = y + (k1 + 2*k2 + 2*k3 + k4) / 6

        x_values.append(x_new)
        y_values.append(y_new)
        
        plt.plot(x_values, [y, k1], label='k1')
        plt.plot(x_values, [y, k2], label='k2')
        plt.plot(x_values, [y, k3], label='k3', linestyle=':')
        plt.plot(x_values, [y, k4], label='k4')
        # plt.scatter(x_new, k2)
        # plt.scatter(x_new, k3)
        plt.scatter(x_new, k4)
        plt.scatter(x_new, k1)

    return x_values, y_values

# Runge-Kutta 메소드를 사용하여 미분 방정식을 풉니다.
x_values, y_values = runge_kutta_4th(x0, y0, h, end_x)

# 결과를 시각화합니다.
plt.plot(x_values, y_values, label='Numerical Solution')
plt.plot(np.linspace(0, end_x, 1000), [(x**2) for x in np.linspace(0, end_x, 1000)], label='Analytical Solution', linestyle='--')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.grid(True)
plt.show()
