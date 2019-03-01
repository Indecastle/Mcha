import numpy as np
import sympy as sm
import matplotlib.pyplot as plt
import math
import timeit

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

    #X = np.append( np.arange(-1, 1, 0.01), 1)
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
    return counter, x1


def kasat(a, E, f, fp):
    xk1, xk = a, a
    counter = 0
    while True:
        counter += 1
        xk = xk1
        xk1 = xk - f(xk)/fp(xk)
        if abs(xk1 - xk) < E:
            break
    return counter, xk1


def newton(FF, Symb, X0, E=0.001):
    n = len(Symb)
    Xk0 = np.array(X0)
    m = {}
    for i in range(n):
        m[Symb[i]] = Xk0[i]

    J = [[0]*n for i in range(n)]
    for i in range(n):  # Матрица Якоби
        for j in range(n):
            J[i][j] = sm.diff(FF[i], Symb[j])
    #print(*J, sep='\n')

    # Jk = [[0]*n for i in range(n)]
    Jk = np.eye(n)
    counter = 1
    while True:
        counter += 1
        for i in range(n):
            m[Symb[i]] = np.real(Xk0[i])
        Fxk = np.array([float(fi.subs(m).n()) for fi in FF], float)

        for i in range(n): # Матрица Якоби сj значениями
            for j in range(n):
                Jk[i][j] = float(J[i][j].subs(m))

        Xk1 = Xk0 - np.linalg.inv(Jk) @ Fxk
        if all(abs(Xk1[i] - Xk0[i]) < E for i in range(n)):
            break
        Xk0 = Xk1
    return counter , Xk1


def newton_mod(FF, Symb, X0, E=0.001):
    n = len(Symb)
    Xk0 = np.array(X0)
    m = {}
    for i in range(n):
        m[Symb[i]] = Xk0[i]

    J = [[0]*n for i in range(n)]
    for i in range(n):  # Матрица Якоби
        for j in range(n):
            J[i][j] = sm.diff(FF[i], Symb[j])
    #print(*J, sep='\n')

    Jk = np.eye(n)
    for i in range(n):  # Матрица Якоби с начальными значениями
        for j in range(n):
            Jk[i][j] = float(J[i][j].subs(m))
    invJk = np.linalg.inv(Jk)

    counter = 1
    while True:
        counter += 1
        for i in range(n):
            m[Symb[i]] = np.real(Xk0[i])
        Fxk = np.array([float(fi.subs(m).n()) for fi in FF], float)

        Xk1 = Xk0 - invJk @ Fxk
        if all(abs(Xk1[i] - Xk0[i]) < E for i in range(n)):
            break
        Xk0 = Xk1
    return counter , Xk1


def MPI_s(FF: list, F0: list, E):
    n = len(F0)
    F1 = []
    F0 = F0.copy()
    counter = 0
    while True:
        counter += 1
        F1 = [FF[i](*F0) for i in range(n)]
        if all(abs(F1[i] - F0[i]) < E for i in range(n)):
            break
        F0 = F1
    return counter, F1 #[round(dec, 3) for dec in F1]


#x, y, z = sm.symbols('x y z')
#f1 = x**2 + y**2 + z**2 - 1
#f2 = 2*x**2 + y**2 - 4*z
#f3 = 3*x**2 - 4*y + z**2
#X = newton((f1,f2,f3), (x,y,z), (1,0.5,0.5), 0.1)
#print(f"[{X[0]}] - {X[1]}")


f2x = lambda x,y: 5/6*math.sin(x+y) - 1/12
f2y = lambda x,y: (1 - x**2) / (2*y)
#f2y = lambda x,y: math.sqrt(0.5 - 0.5 * x**2)
#f2y2 = lambda x,y: -math.sqrt(0.5 - 0.5 * x**2)
E = 0.001


def main():
    #plots()

    counter, x1 = hord(0.1, 3.0, E, f)
    print(f"Метод хорд: x[{counter}] = {x1}")
    counter, x1 = kasat(0.1, E, f, fp)
    print(f"Метод касательных: x[{counter}] = {x1}")

    x, y = sm.symbols('x y')
    f1 = sm.sin(x+y) - 1.2*x - 0.1
    f2 = x**2 + 2*y**2 - 1
    X = newton((f1,f2), (x,y), (0.5,0.5), E)
    print(f"Метод Ньютона: [{X[0]}] - {X[1]}")

    X = newton_mod((f1, f2), (x, y), (0.5, 0.5), E)
    print(f"Метод Ньютона(мод.): [{X[0]}] - {X[1]}")

    counter, F1 = MPI_s([f2x, f2y], [0.5, 0.5], E)
    print(f"МПИ: {[counter]} - {F1}")


if __name__ == '__main__':
    main()
    input()
    #print(min(timeit.repeat("f()", "from __main__ import f", number=50)))
