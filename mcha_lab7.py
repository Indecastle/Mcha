import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import math
import timeit
from tabulate import tabulate

x, y = sp.symbols('x y')


def error_estimation(method, func, a, b, eps):
    n = 1
    while True:
        if abs(method(func, a, b, 2 * n) - method(func, a, b, n))/3 < eps:
            break
        n *= 2
    #print(abs(method(func, a, b, n * 2) - method(func, a, b, n)))
    return n * 2


def trapez(func, a, b, n):
    h = (b - a) / n
    return h * sum([(func(a + h * i) + func(a + h * (i + 1))) / 2 for i in range(n)])


def simpson(func, a, b, n):
    h = (b - a) / n
    return h / 3 * sum(
        [func(a + h * i) + 4 * func(a + h * (i + 1)) + func(a + h * (i + 2)) for i in range(0, n - 1, 2)])


def newton_leibniz(A, a, b):
    h = sp.integrate(A, x)
    return sp.simplify(h.subs(x, b) - h.subs(x, a))



def de_error_estimation(ode_method, func, a, b, x0, y0, eps):
    n = 2
    while True:
        k, y1 = ode_method(func, a, b, x0, y0, n)
        k, y2 = ode_method(func, a, b, x0, y0, n//2)
        if abs(y2[-1] - y1[-1]) < eps:
            break
        n *= 2
    return n


def runge(func, a, b, x0, y0, n):
    h = (b - a) / n
    x = np.empty(n + 1, float)
    y = np.empty(n + 1, float)
    x[0] = x0
    y[0] = y0
    for i in range(n):
        F1 = h * func(x[i], y[i])
        F2 = h * func(x[i] + h / 2, y[i] + F1 / 2)
        F3 = h * func(x[i] + h / 2, y[i] + F2 / 2)
        F4 = h * func(x[i] + h, y[i] + F3)
        y[i + 1] = y[i] + 1 / 6 * (F1 + F4 + 2 * (F2 + F3))
        x[i + 1] = x[i] + h
    return x, y


def adams(f, a, b, x0, y0, n):
    h = (b - a) / n
    x = np.empty(n + 1)
    y = np.empty(n + 1)
    x[0] = x0
    y[0] = y0
    x[1] = x[0] + h
    y[1] = y[0] + h * f(x[0], y[0])
    for i in range(1, n):
        p = y[i] + h / 2 * (3 * f(x[i], y[i]) - f(x[i - 1], y[i - 1]))
        x[i + 1] = x[i] + h
        y[i + 1] = y[i] + h / 2 * (f(x[i], y[i]) + f(x[i + 1], p))
    return x, y


def eiler(func, a, b, x0, y0, n):
    h = (b - a) / n
    x = np.empty(n + 1, float)
    y = np.empty(n + 1, float)
    x[0] = x0
    y[0] = y0
    for i in range(n):
        y[i + 1] = y[i] + h * func(x[i], y[i])
        x[i + 1] = x[i] + h
    return x, y



def main():
    plt.grid()

    A = x * 3**-x
    func = sp.lambdify(x, A)
    a_S = 0
    b_S = 1.5
    eps = 0.0001

    n = error_estimation(trapez, func, a_S, b_S, eps)
    print('шаг:', (b_S - a_S) / n)
    print('метод трапеций:')
    t = trapez(func, a_S, b_S, n)
    print('n:  ', t)
    t2 = trapez(func, a_S, b_S, n // 2)
    print('n/2:', t2)
    print('метод Симпcон:')
    s = simpson(func, a_S, b_S, n)
    print('n:  ', s)
    s2 = simpson(func, a_S, b_S, n // 2)
    print('n/2:', s2)
    f = newton_leibniz(A, a_S, b_S)
    print('метод Ньютона:', f, '=', f.n(16))
    sym = sp.integrate(A, (x, a_S, b_S))
    print('точное решение:', sp.simplify(sym), '=', sym.n(16))

    print()
    
    diff_A = x*y**2 - y
    #diff_A = x * y ** 2 + y
    diff_func = sp.lambdify((x, y), diff_A)
    diff_aS = 0
    diff_bS = 2
    x0 = 0
    y0 = 1

    nd = de_error_estimation(runge, diff_func, diff_aS, diff_bS, x0, y0, 0.0001)
    runge_x, runge_y = runge(diff_func, diff_aS, diff_bS, x0, y0, nd)
    runge_x_2, runge_y_2 = runge(diff_func, diff_aS, diff_bS, x0, y0, nd // 2)
    plt.plot(runge_x, runge_y, label='Рунге')
    table1 = [['xᵢ', 'yᵢ', 'ỹᵢ', '∆ᵢ = |yᵢ - ỹᵢ|']]
    for i in range(len(runge_x)):
        if i % 2 == 0:
            table1.append([runge_x[i], runge_y[i], runge_y_2[i // 2], abs(runge_y_2[i // 2] - runge_y[i])])
        #else:
        #    table.append([runge_x[i], runge_y[i], None, None])
    print("Метод Рунге:")
    print(tabulate(table1, tablefmt='fancy_grid'))

    adams_x, adams_y = adams(diff_func, diff_aS, diff_bS, x0, y0, nd)
    adams_x_2, adams_y_2 = adams(diff_func, diff_aS, diff_bS, x0, y0, nd // 2)
    plt.plot(adams_x, adams_y, label='Адамс')
    table2 = [['xᵢ', 'yᵢ', 'ỹᵢ', '∆ᵢ = |yᵢ - ỹᵢ|']]
    for i in range(len(adams_x)):
        if i % 2 == 0:
            table2.append([adams_x[i], adams_y[i], adams_y_2[i // 2], abs(adams_y_2[i // 2] - adams_y[i])])
    print("Метод Адамса:")
    print(tabulate(table2, tablefmt='fancy_grid'))

    eiler_x, eiler_y = eiler(diff_func, diff_aS, diff_bS, x0, y0, nd)
    eiler_x_2, eiler_y_2 = eiler(diff_func, diff_aS, diff_bS, x0, y0, nd // 2)
    plt.plot(eiler_x, eiler_y, label='Эйлер')
    table3 = [['xᵢ', 'yᵢ', 'ỹᵢ', '∆ᵢ = |yᵢ - ỹᵢ|']]
    for i in range(len(eiler_x)):
        if i % 2 == 0:
            table3.append([eiler_x[i], eiler_y[i], eiler_y_2[i // 2], abs(eiler_y_2[i // 2] - eiler_y[i])])
    print("Метод Эйлера:")
    print(tabulate(table3, tablefmt='fancy_grid'))

    f = sp.Function('f')
    C1 = sp.Symbol('C1')
    solve = sp.simplify( sp.dsolve(sp.diff(f(x), x) + f(x) - x * f(x) ** 2, f(x)) )
    print('решение:', solve)
    CC = sp.solve(solve.subs({x:x0, f(x):y0}), C1)[0]
    #print('C1 =', CC)
    print('задача Коши:', solve.rhs.subs(C1, CC))
    solution = sp.lambdify(x, solve.rhs.subs(C1, CC))
    xL = np.linspace(diff_aS, diff_bS)
    plt.plot(xL, solution(xL), 'r--', label='точное решение')
    
    legend = plt.legend(loc='upper center', shadow=True, fontsize='x-large', frameon=True)
    #legend.get_frame().set_facecolor('C0')
    plt.show()
    

if __name__ == '__main__':
    main()


    #x = sm.symbols('x')
    #a = sm.plot((sm.sqrt(x),(x,-5,5)), (x**2, (x,-2,2)), show=False, legend=True)
    #a[0].line_color = 'green'
    #a[0].label = '$\\frac{dy}{dx}$'
    #a[1].label = '$y$ y'
    #a.show()

    
    #print(min(timeit.repeat("f()", "from __main__ import f", number=50)))
