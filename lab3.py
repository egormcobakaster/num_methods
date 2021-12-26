import numpy as np
import matplotlib.pyplot as plt
import math


a = 0
b = 1
I = (2 + math.pi) / 4


def func(x):
    return math.sqrt(2 - x**2)


def rect(n, func=func):
    h = (b - a) / (n - 1)
    X = np.linspace(a, b, n)
    res = 0
    for i in range(1, n):
        res += func(X[i])
    return res * h


def trapez(n, func=func):
    h = (b - a) / (n - 1)
    X = np.linspace(a, b, n)
    res = (func(X[0]) + func(X[n-1])) / 2
    for i in range(1, n - 1):
        res += func(X[i])
    return res * h


def Simpson(n, func=func):
    h = (b - a) / (2 * n)
    X = np.linspace(a, b, (2 * n + 1))
    res = func(X[0]) + func(X[2 * n])
    for i in range(1, 2 * n + 1, 2):
        res += 4 * func(X[i])
    for j in range(2, 2 * n, 2):
        res += 2 * func(X[j])
    return res * h / 3


h_list = []
rect_list = []
trapez_list = []
simpson_list = []
for i in range(50, 300):
    h_list.append(math.log((b - a) / (i-1)))
    rect_list.append(math.log(abs(rect(i) - I)))
    trapez_list.append(math.log(abs(trapez(i) - I)))
    simpson_list.append(math.log(abs(Simpson(i) - I)))
    print(Simpson(i), rect(i), trapez(i))
X = np.array(h_list)
Y1 = np.array(rect_list)
Y2 = np.array(trapez_list)
Y3 = np.array(simpson_list)
A1 = np.concatenate([np.ones(len(h_list)).reshape(-1,1), X.reshape(-1,1)], axis=1)

w1 = np.linalg.inv(A1.T @ A1) @ A1.T @ Y1
w2 = np.linalg.inv(A1.T @ A1) @ A1.T @ Y2
w3 = np.linalg.inv(A1.T @ A1) @ A1.T @ Y3
print(w1)
print(w2)
print(w3)
plt.grid(True)
plt.plot(h_list, rect_list, color="r", linewidth=2, label="Rectangles method")
plt.plot(h_list, trapez_list, color="b", linewidth=2, label="Trapezes method")
plt.plot(h_list, simpson_list, color="g", linewidth=2, label="Simpson method")
plt.legend(loc='lower right', fontsize=8)
plt.show()