from scipy import optimize
from numpy import log
import matplotlib.pyplot as plt

a = 2.
b = 1.
c = 1.
d = 1.

N = 10**2
dt = 0.7

X = [1.]*(N+1)
Y = [1.]*(N+1)


def discrete_variational(delta=10**-3):
    for i in range(N):
        def f(x): return (x[0]-X[i]-(x[0]+X[i])*(x[1]+Y[i])/4*(a*(log(x[1])-log(Y[i]))/(x[1]-Y[i])-b)*dt, x[1]-Y[i]-(x[0]+X[i])*(x[1]+Y[i])/4*(c-d*(log(x[0])-log(X[i]))/(x[0]-X[i]))*dt)
        res = optimize.root(f, [X[i]+delta, Y[i]+delta])  # 0徐算を防ぐため微少量だけずらす
        X[i+1], Y[i+1] = res.x[0], res.x[1]


discrete_variational()

plt.plot(X, Y)
plt.title(f'discrete variational(dt={dt}, N={N})')
plt.grid()
plt.show()
