import matplotlib.pyplot as plt
import math
import numpy as np
from numpy import linalg
def func(x):
    return np.exp(-x) * np.sin(x)
def dif(x):
    return np.exp(-x) * (np.cos(x) - np.sin(x))


def dudif(x):
    return -np.exp(-x) * 2 * np.cos(x)

def numdif_1(Y, n, h):
    ND = []
    for i in range(0, n):
        if i < n-1:
            ND.append((Y[i+1] - Y[i]) / h)
        else:
            ND.append((Y[i] - Y[i-1]) / h)
    return np.array(ND)

def numdif(Y, n, h):
    ND = []
    for i in range(0, n):
        if i == 0:
            ND.append((-3 / 2 * Y[0] + 2 * Y[1] - 1 / 2 * Y[2]) / h)
        elif i == n-1:
            ND.append((3 / 2 * Y[n-2] - 2 * Y[n-2] + 1 / 2 * Y[n-1]) / h)
        else:
            ND.append((Y[i+1] - Y[i-1]) / (2 * h))
    return np.array(ND)
def numdudif(Y_1, n, h):
    ND = []
    X = np.linspace(a - 2 * h, b + 2 * h, n + 4)
    Y = func(X)
    for i in range(2,n+2):
        ND.append(( -Y[i - 2] + 16 * Y[i - 1] - 30 * Y[i] + 16 * Y[i + 1] - Y[i +2]) / (12 * h*h))
    return np.array(ND)


a = 0.
b = 10.
sruv = []
sruv1 = []
sruv2 = []
h = []
for n in range(0, 100):

    k = 100 + n * 5
    h.append((b - a) / (k-1))
    X = np.linspace(a, b, k)
    Y = func(X)
    NDD = numdudif(Y, k, h[n])
    ND = numdif(Y, k, h[n])
    DD = dudif(X)
    D = dif(X)
    ND1 = numdif_1(Y, k, h[n])
    delt = []
    delt1=[]
    delt2 = []
    for i in range(0, k):
        delt.append(np.abs(DD[i] - NDD[i]))
    sruv.append(np.array(delt).max())
    for i in range(0, k):
        delt1.append(np.abs(D[i] - ND[i]))
    sruv1.append((np.array(delt1).max()))
    for i in range(0, k):
        delt2.append(np.abs(D[i] - ND1[i]))
    sruv2.append((np.array(delt2).max()))
print(NDD)
print(ND - D)
plt.plot(X, NDD, color = 'red')
plt.plot(X, ND, color = 'g')
plt.plot(X, Y)

plt.show()
plt.plot(np.log(h), np.log(sruv1), color = 'red')
plt.plot(np.log(h), np.log(sruv), color = 'green')
plt.plot(np.log(h), np.log(sruv2))
X = np.log(h)
Y1 = np.log(sruv1)
Y2 = np.log(sruv)
Y3 = np.log(sruv2)
A1 = np.concatenate([np.ones(len(h)).reshape(-1,1), X.reshape(-1,1)], axis=1)

w1 = np.linalg.inv(A1.T @ A1) @ A1.T @ Y1
w2 = np.linalg.inv(A1.T @ A1) @ A1.T @ Y2
w3 = np.linalg.inv(A1.T @ A1) @ A1.T @ Y3
print(w1)
print(w2)
print(w3)
plt.grid(True)
plt.show()

