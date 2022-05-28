import matplotlib.pyplot as plt

a = 2.
b = 1.
c = 1.
d = 1.


def f(x, y): return a*x-b*x*y, c*x*y-d*y


N = 10**3
dt = 0.01

X = [1.]*(N+1)
Y = [1.]*(N+1)


def RK():
    for i in range(N):
        k1 = f(X[i], Y[i])[0]*dt, f(X[i], Y[i])[1]*dt
        k2 = f(X[i]+k1[0]/2, Y[i]+k1[1]/2)[0]*dt, f(X[i]+k1[0]/2, Y[i]+k1[1]/2)[1]*dt
        k3 = f(X[i]+k2[0]/2, Y[i]+k2[1]/2)[0]*dt, f(X[i]+k2[0]/2, Y[i]+k2[1]/2)[1]*dt
        k4 = f(X[i]+k3[0], Y[i]+k3[1])[0]*dt, f(X[i]+k3[0], Y[i]+k3[1])[1]*dt

        X[i+1] = X[i]+(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6
        Y[i+1] = Y[i]+(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6


RK()
plt.plot(X, Y)
plt.title(f'RK4(dt={dt}, N={N})')
plt.grid()
plt.show()
