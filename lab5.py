import numpy as np
import math
from matplotlib import pyplot as plt

a = 0
b = 1
u_0 = 2
u1_0 = -2


def u2(x, u, u1):
    return -x / (1 + x*x) * u1 + 1 / (1 + x*x) * u + (3 - 2*x + 4*x**2) / (1 + x*x) * math.exp(-2*x)


def u_real(x):
    return np.power((1 + x*x), 1/2) + math.exp(-2*x)


def Runge_Kutta(X, F, u0, u1, h):
    u = np.zeros(len(X))
    w = np.zeros(len(X))
    u[0] = u0
    w[0] = u1

    for i in range(len(X) - 1):
        coef_u = np.zeros(4)
        coef_w = np.zeros(4)
        coef_w[0] = F(X[i], u[i], w[i])
        coef_u[0] = w[i]
        coef_w[1] = F(X[i] + h / 2, u[i] + h / 2 * coef_u[0], w[i] + h / 2 * coef_w[0])
        coef_u[1] = w[i] + h / 2 * coef_w[0]
        coef_w[2] = F(X[i] + h / 2, u[i] + h / 2 * coef_u[1], w[i] + h / 2 * coef_w[1])
        coef_u[2] = w[i] + h / 2 * coef_w[1]
        coef_w[3] = F(X[i] + h, u[i] + h * coef_u[2], w[i] + h * coef_w[2])
        coef_u[3] = w[i] + h * coef_w[2]
        w[i+1] = w[i] + h / 6 * (coef_w[0]+ 2 * (coef_w[1] + coef_w[2]) + coef_w[3])
        u[i+1] = u[i] + h / 6 * (coef_u[0]+ 2 * (coef_u[1] + coef_u[2]) + coef_u[3])
    return u


def Adams(X, F, u0, u1, h):
    u = np.zeros(len(X))
    w = np.zeros(len(X))
    u[0] = u0
    w[0] = u1

    for i in range(2):
        coef_u = np.zeros(4)
        coef_w = np.zeros(4)
        coef_w[0] = F(X[i], u[i], w[i])
        coef_u[0] = w[i]
        coef_w[1] = F(X[i] + h / 2, u[i] + h / 2 * coef_u[0], w[i] + h / 2 * coef_w[0])
        coef_u[1] = w[i] + h / 2 * coef_w[0]
        coef_w[2] = F(X[i] + h / 2, u[i] + h / 2 * coef_u[1], w[i] + h / 2 * coef_w[1])
        coef_u[2] = w[i] + h / 2 * coef_w[1]
        coef_w[3] = F(X[i] + h, u[i] + h * coef_u[2], w[i] + h * coef_w[2])
        coef_u[3] = w[i] + h * coef_w[2]
        w[i + 1] = w[i] + h / 6 * (coef_w[0] + 2 * (coef_w[1] + coef_w[2]) + coef_w[3])
        u[i + 1] = u[i] + h / 6 * (coef_u[0] + 2 * (coef_u[1] + coef_u[2]) + coef_u[3])

    for i in range(len(X) - 3):
        F_0 = F(X[i], u[i], w[i])
        F_1 = F(X[i + 1], u[i + 1], w[i + 1])
        F_2 = F(X[i + 2], u[i + 2], w[i + 2])
        w[i+3] = (w[i + 2] + h * (23 / 12 * F_2 - 4 / 3 * F_1 + 5 / 12 * F_0))
        u[i+3] = u[i + 2] + h * (23 / 12 * w[i + 2] - 4 / 3 * w[i + 1] + 5 / 12 * w[i])

    return u


def Euler(X, F, u0, u1, h):
    u = np.zeros(len(X))
    w = np.zeros(len(X))
    u[0] = u0
    w[0] = u1

    for i in range(len(X) - 1):
        w[i+1] = w[i] + h * F(X[i], u[i], w[i])
        u[i+1] = u[i] + h * w[i]

    return u


def Runge_correction(u_h, u_h_2, p):
    runge_corr = []
    for i in range(len(u_h)):
        n = 2 * i
        runge_elem = (u_h_2[n] - u_h[i]) / (2**p - 1)
        runge_corr.append(runge_elem)

    return runge_corr


h1 = 0.01
x1 = np.linspace(a, b, round((b - a) / h1) + 1)
u1 = [u_real(x) for x in x1]
RK = Runge_Kutta(x1, u2, u_0, u1_0, h1)
AD = Adams(x1, u2, u_0, u1_0, h1)
Eul = Euler(x1, u2, u_0, u1_0, h1)
fig, ax = plt.subplots(figsize=(14, 8))
ax.plot(x1, u1, color="r", label='u')
ax.plot(x1, RK, color="g", label='RK')
ax.plot(x1, AD, color="b", label='Ad')
ax.plot(x1, Eul, color="brown", label='Eul')
ax.grid(which="major", linewidth=1.2)
ax.grid(which="minor", linestyle="--", color="gray", linewidth=0.5)

ax.set_xlabel("x", fontsize=14)
ax.set_ylabel("u", fontsize=14)
plt.legend()
plt.show()
h2 = 0.1
n1=5
n2=100
errors1 = []
errors2 = []
errors3 = []
steps = []
for n in range(n1,n2):
    x, h = np.linspace(a, b, n, retstep=True)
    RK = Runge_Kutta(x, u2, u_0, u1_0, h)
    AD = Adams(x, u2, u_0, u1_0, h)
    Eul = Euler(x, u2, u_0, u1_0, h)
    u_h = [u_real(x) for x in x]
    steps.append(np.log(h))
    max_errorRK = max([abs(rk - u) for rk, u in zip(RK, u_h)])
    max_errorAD = max([abs(rk - u) for rk, u in zip(AD, u_h)])
    max_errorEul = max([abs(rk - u) for rk, u in zip(Eul, u_h)])

    errors1.append(np.log(max_errorRK))
    errors2.append(np.log(max_errorAD))
    errors3.append(np.log(max_errorEul))
fig, ax = plt.subplots(figsize=(14, 8))
ax.plot(steps, errors1, color="g", label='RK')
ax.plot(steps, errors2, color="b", label='Ad')
ax.plot(steps, errors3, color="brown", label='Eul')
ax.grid(which="major", linewidth=1.2)
ax.grid(which="minor", linestyle="--", color="gray", linewidth=0.5)

ax.set_xlabel("x", fontsize=14)
ax.set_ylabel("u", fontsize=14)
plt.legend()
plt.show()
X = np.array(steps)
Y1 = np.array(errors1)
Y2 = np.array(errors2)
Y3 = np.array(errors3)
A1 = np.concatenate([np.ones(len(steps)).reshape(-1,1), X.reshape(-1,1)], axis=1)
w1 = np.linalg.inv(A1.T @ A1) @ A1.T @ Y1
w2 = np.linalg.inv(A1.T @ A1) @ A1.T @ Y2
w3 = np.linalg.inv(A1.T @ A1) @ A1.T @ Y3
print(w1)
print(w2)
print(w3)

errors1 = []
errors2 = []
steps = []
steps = []


for n in np.arange(n1, n2 + 1, step=10):

    x0, h = np.linspace(a, b, n, retstep=True)

    x1 = []
    for i in range(len(x0) - 1):
        x1.append(x0[i])
        x1.append((x0[i] + x0[i + 1]) / 2)
    x1.append(x0[-1])

    RK0 = Runge_Kutta(x0, u2, u_0, u1_0, h)
    RK1 = Runge_Kutta(x1, u2, u_0, u1_0, h/2)
    u = [u_real(x) for x in x0]

    runge_range = Runge_correction(RK0, RK1, 4)

    RK_inc = []
    for i in range(len(runge_range)):

        RK_inc.append(RK1[2*i] + runge_range[i])

    max_error0 = max([abs(rk - u) for rk, u in zip(RK0, u)])
    max_error1 = max([abs(rk - u) for rk, u in zip(RK_inc, u)])

    steps.append(np.log(h))
    errors1.append(np.log(max_error0))
    errors2.append(np.log(max_error1))
X = np.array(steps)
Y1 = np.array(errors1)
Y2 = np.array(errors2)

A1 = np.concatenate([np.ones(len(steps)).reshape(-1,1), X.reshape(-1,1)], axis=1)
w1 = np.linalg.inv(A1.T @ A1) @ A1.T @ Y1
w2 = np.linalg.inv(A1.T @ A1) @ A1.T @ Y2

print(w1)
print(w2)

fig, ax = plt.subplots(figsize=(14, 8))
ax.plot(steps, errors1, color="g", label='RK')
ax.plot(steps, errors2, color="b", label='RK_cor')

ax.grid(which="major", linewidth=1.2)
ax.grid(which="minor", linestyle="--", color="gray", linewidth=0.5)

ax.set_xlabel("x", fontsize=14)
ax.set_ylabel("u", fontsize=14)
plt.legend()
plt.show()








































