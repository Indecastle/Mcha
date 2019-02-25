import numpy as np
import matplotlib.pyplot as plt
import math

f = lambda x: math.exp(x) + math.log(x, math.e) - x
fp = lambda x: math.exp(x) + 1/x
f1 = lambda x: math.asin(0.1 + 1.2*x) - x
f21 = lambda x: math.sqrt(0.5 - 0.5 * x**2)
f22 = lambda x: -math.sqrt(0.5 - 0.5 * x**2)

def plots(y=f, y1=f1, y21=f21, y22=f22):
    fig, (ax1, ax2) = plt.subplots(
        nrows=1, ncols=2,
        figsize=(8, 4)
    )

    X = np.arange(0.01, 3, 0.01)
    Y = [y(np.real(x)) for x in X]
    ax1.plot(X, Y, 'g')

    X = np.arange(-0.91, 0.75, 0.01)
    Y = [y1(np.real(x)) for x in X]
    ax2.plot(X, Y, 'g')

    X = np.append( np.arange(-1, 1, 0.01), 1)
    Y = [y21(np.real(x)) for x in X]
    ax2.plot(X, Y, 'b')

    X = np.append( np.arange(-1, 1, 0.01), 1)
    Y = [y22(np.real(x)) for x in X]
    ax2.plot(X, Y, 'b')

    ax1.axis([0, 3, -2,16])
    #ax1.set_aspect('equal')
    ax1.grid(True, which='both')
    ax1.axhline(y=0, color='k')
    ax1.axvline(x=0, color='k')

    ax2.axis([-1.1,1.1, -1.1,1.1])
    #ax2.set_aspect('equal')
    ax2.grid(True, which='both')
    ax2.axhline(y=0, color='k')
    ax2.axvline(x=0, color='k')

    plt.show()


plots()


def hord(a, b, E: float, f):
    x1, x0, fb = a, 0, f(b)
    counter = 0
    while True:
        counter += 1
        x0 = x1
        fx0 = f(x1)
        x1 = x0 - fx0 / (fb - fx0) * (b - x0)
        if abs(x1 - x0) < E:
            break
    print(f"x[{counter}] = {x1}")


def kasat(a, E, f, fp):
    xk1, xk = a, a
    counter = 0
    while True:
        counter += 1
        xk = xk1
        xk1 = xk - f(xk)/fp(xk)
        if abs(xk1 - xk) < E:
            break
    print(f"x[{counter}] = {xk1}")


#def force():



E = 0.00001
hord(0.1, 3.0, E, f)
kasat(0.1, E, f, fp)