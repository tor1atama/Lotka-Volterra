from scipy import optimize
from numpy import log
import matplotlib.pyplot as plt

a = 2.
b = 1.
c = 1.
d = 1.

M = 10**5
X = [1.]*M
Y = [1.]*M


def discrete_variational(N=10**4, dt=0.7, delta=10**-3):
    for i in range(N):
        def f(x): return (x[0]-X[i]-(x[0]+X[i])*(x[1]+Y[i])/4*(a*(log(x[1])-log(Y[i]))/(x[1]-Y[i])-b)* dt, x[1]-Y[i]-(x[0]+X[i])*(x[1]+Y[i])/4*(c-d*(log(x[0])-log(X[i]))/(x[0]-X[i]))*dt)
        res = optimize.root(f, [X[i]+delta, Y[i]+delta])  # 0徐算を防ぐため微少量だけずらす
        X[i+1], Y[i+1] = res.x[0], res.x[1]
    plt.plot(X[:N+1], Y[:N+1], label=f'discrete variational(dt={dt}, N={N})')


def RK(N=10**3, dt=10**-2):
    def f(x, y): return (a*x-b*x*y, c*x*y-d*y)
    for i in range(N):
        k1 = f(X[i], Y[i])[0]*dt, f(X[i], Y[i])[1]*dt
        k2 = f(X[i]+k1[0]/2, Y[i]+k1[1]/2)[0] * dt, f(X[i]+k1[0]/2, Y[i]+k1[1]/2)[1]*dt
        k3 = f(X[i]+k2[0]/2, Y[i]+k2[1]/2)[0] * dt, f(X[i]+k2[0]/2, Y[i]+k2[1]/2)[1]*dt
        k4 = f(X[i]+k3[0], Y[i]+k3[1])[0]*dt, f(X[i]+k3[0], Y[i]+k3[1])[1]*dt

        X[i+1] = X[i]+(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6
        Y[i+1] = Y[i]+(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6
    plt.plot(X[:N+1], Y[:N+1], label=f'RK4(dt={dt}, N={N})')


RK()
discrete_variational()

plt.grid()
plt.legend()
plt.show()