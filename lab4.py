import numpy as np
import matplotlib.pyplot as plt
import math
a = 0
b = 10

def func(x):
    return x**5 - 7 * x**3 - 3 * x - 2


def dif_func(x):
    return 5 * x**4 - 21 * x**2 - 3

def dichotomy(func, a, b, eps, vec_x):
    if (func(a) * func(b) < 0):
        c = (b + a) / 2
        vec_x.append(c)
        if abs(b - a) > eps:
            if func(c) == 0:
                return c, 0
            elif func(a) * func(c) < 0:
                x, n = dichotomy(func, a, c, eps, vec_x)
                return x, n +1
            else:
                x, n = dichotomy(func, c, b, eps, vec_x)
                return x, n +1
        else:
            return c, 0
    elif (func(a) * func(b) == 0):
        if func(a) == 0:
            return a, 0
        else:
            return b, 0

    else:
        return None, 0


def Newton(func, x, eps, vec_x):
    n = 0
    while (dif_func(x) == 0):
        x += 1e-5
    while (abs(func(x)) > eps):
        x = x - func(x) / dif_func(x)
        vec_x.append(x)
        n = n + 1
    return x, n


def sec(x0, x1, func, eps, vec_x):
    x2 = x0 - func(x0) * (x0 - x1) / (func(x0) - func(x1))
    vec_x.append(x2)
    if abs(x2 - x1) > eps:
        x, n = sec(x1, x2, func, eps, vec_x)
        return x, n + 1

    else:
        return x2, 1

x = np.linspace(a, b, 100)
y = [func(x) for x in x]
y1 = [dif_func(x) for x in x]

fig, ax  = plt.subplots(figsize=(14, 5))

ax.plot(x, y, color="red")
ax.plot(x, y1, color="g")
ax.grid(which="major", linewidth=1.2)
ax.grid(which="minor", linestyle="--", color="gray", linewidth=0.5)


plt.show()

roots = []
steps = []

epss = [1e-3, 1e-6, 1e-9]

vec_x = [[],[],[]]
i=0
for eps in epss:
    root,  step_val = dichotomy(func, a, b, eps, vec_x[i])
    roots.append(root)
    steps.append(step_val)
    i = i + 1
print('dich')
print(roots)
print(steps)
print(vec_x[2])
X = vec_x[2]
n = len(X)
print(math.log(abs((X[n-1] - X[n-2]) / (X[n-2] - X[n-3]))) / math.log(abs((X[n-2] - X[n-3]) / (X[n-3] - X[n-4]))))


roots = []
steps = []
vec_x1 = [[],[],[]]
i=0
for eps in epss:
    x, n = Newton(func, (a + b)/2, eps, vec_x1[i])
    roots.append(x)
    steps.append(n)
    i = i + 1
print('Newton')
print(roots)
print(steps)
print(vec_x1[0])
X = vec_x1[2]
n = len(X)
print(math.log(abs((X[n-1] - X[n-2]) / (X[n-2] - X[n-3]))) / math.log(abs((X[n-2] - X[n-3]) / (X[n-3] - X[n-4]))))


roots = []
steps = []
vec_x2 = [[],[],[]]
i=0
for eps in epss:
    x,  n = sec(b, b - eps, func, eps, vec_x2[i])
    roots.append(x)
    steps.append(n)
    i = i + 1
print('sec')
print(roots)
print(steps)
print(vec_x2[0])
X = vec_x2[2]
n = len(X)
print(math.log(abs((X[n-1] - X[n-2]) / (X[n-2] - X[n-3]))) / math.log(abs((X[n-2] - X[n-3]) / (X[n-3] - X[n-4]))))