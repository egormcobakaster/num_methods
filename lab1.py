import matplotlib.pyplot as plt
import math
import numpy as np


def func(x):
    return np.sin(np.exp(x/3)/10)

def cheb(a, b, n):
    X=[]
    for i in range(1,n+1):
        X.append((a+b)/2 + (a-b) / 2 * math.cos((2*i -1)/(2 * n) * math.pi))
    return np.array(X)
def interpol(Xnew, X, Y, n):
    Ynew = np.zeros((n - 1))
    l = len(Xnew)
    for k in range(n - 1):

        S = 0.

        for i in range(n):
            P = 1.
            for j in range(n):
                if j != i:
                    P = P * (Xnew[k] - X[j]) / (X[i] - X[j])
            S = S + P * Y[i]
        Ynew[k] = S

    return Ynew


a=np.array([])
np.append(a,0)
a=0.
b=10.
d=20
Nmas=[]
print(Nmas)
for i in range(0,d):
    Nmas.append(10 + 10*i)
h=[]
deltmas = []
for n in Nmas:
    h.append(10/n)
    Xmas = cheb(a, b, n)
    Fmas = func(Xmas)

    Xnew = np.linspace(a, b, num=n -1)
    Fnew = interpol(Xnew, Xmas, Fmas, n)
    Fmas = func(Xnew)

    Fdelt = np.abs(Fnew - Fmas)

    delt = Fdelt.mean()
    deltmas.append(delt)
plt.plot(h, deltmas, color = 'red')
plt.show()
print(deltmas)
h=[]
deltmas = []
for n in Nmas:
    h.append(10/n)
    Xmas = np.linspace(a, b, num=n)
    Fmas = func(Xmas)

    Xnew = np.linspace(a, b, num=n -1)
    Fnew = interpol(Xnew, Xmas, Fmas, n)
    Fmas = func(Xnew)

    Fdelt = np.abs(Fnew - Fmas)

    delt = Fdelt.mean()
    deltmas.append(delt)
plt.plot(h, np.log(deltmas), color = 'red')
plt.show()
print(deltmas)
deltmas1 = []


n=100
Xmas = cheb(a, b, n)

Fmas = func(Xmas)
Xnew=np.zeros(n-1)
for i in range(0,n-1):
    Xnew[i] = Xmas[i]/2+Xmas[i+1]/2
print(Xmas)
print(Xnew)
Fnew = interpol(Xnew, Xmas, Fmas, n)
Fmas = func(Xnew)
print(Fmas)
print(Fnew)
plt.plot(Xnew, Fnew)
plt.plot(Xnew, Fmas, color = 'red')
plt.show()
print((np.abs(Fmas - Fnew).mean()))

