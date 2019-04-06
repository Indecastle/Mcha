import numpy as np
import sympy as sm
import matplotlib.pyplot as plt
import math
import timeit
from tabulate import tabulate

X = (0.234, 0.649, 1.382, 2.672, 2.849)
Y = (0.511, 0.982, 2.411, 3.115, 4.184)
x = sm.symbols('x')
#plt.style.use('ggplot')
#print(plt.style.available)

def inter(X, Y):
    n = len(X)
    def li(i):
        result = 1
        for j in [j for j in range(n) if j != i]:
            result *= (x-X[j])/(X[i]-X[j])
        return result
    L = 0
    for i in range(n):
        L += Y[i]*li(i)
    return sm.expand(L)
        

def konechn(X, Y):
    n = len(X)+1
    table = [['']*n for _ in range(n)]
    table[0][0] = 'x'
    table[0][1] = 'y'
    table[0][2] = '∆y'
    for i in range(3, n):
        table[0][i] = f'∆{i-1}yk'
    for i in range(1, n):
        table[i][0] = X[i-1]
        table[i][1] = Y[i-1]
    for j in range(2, n):
        for i in range(1, n-j+1):
            table[i][j] = table[i+1][j-1] - table[i][j-1]
    return table


def razdel(X, Y):
    n = len(X)
    m = (n+n-1) + 1
    n+=1
    table = [['']*(n) for _ in range(m)]
    table[0][0] = 'x'
    table[0][1] = 'y'
    table[0][2] = '1-го порядка'
    for j in range(3, n):
        table[0][j] = f'{j-1}-го порядка'
    for i in range(1, n):
        table[i*2-1][0] = X[i-1]
        table[i*2-1][1] = Y[i-1]

    pol = table[1][1]
    for j in range(2, n):
        for i in range(1+j-1, m-j+1, 2):
            h = j - 1
            table[i][j] = (table[i+1][j-1] - table[i-1][j-1]) / (table[i+h][0] - table[i-h][0])

        mul=1
        for i in range(j-1):
            mul *= (x-X[i])
        pol += table[j][j] * mul

    
    return table, pol
    

def ksn(X, Y):
    n = len(X)
    a, b = 0, 0
    fx = []
    for i in range(n-1):
        A = np.array([[X[i], 1], [X[i+1], 1]], float)
        B = np.array([Y[i], Y[i+1]], float)
        solve = np.linalg.solve(A, B)
        fx.append([solve[0]*x + solve[1]])

    args = []
    for i, f in enumerate(fx):
        args.append((*f, (x, X[i], X[i+1])))
    #print(args)
    plot = sm.plot(*args, show=False)
    return fx, plot


def ksk(X, Y):
    n = len(X)
    a, b, c = 0, 0, 0
    fx = []
    for i in range(1, n//2 + 1):
        A = np.array([[X[2*i-2]**2, X[2*i-2], 1], 
                      [X[2*i-1]**2, X[2*i-1], 1], 
                      [X[2*i]**2, X[2*i], 1]], float)
        B = np.array([Y[2*i-2], Y[2*i-1], Y[2*i]], float)
        solve = np.linalg.solve(A, B)
        fx.append([solve[0]*x**2 + solve[1]*x + solve[2]])

    args = []
    for i, f in enumerate(fx):
        args.append((*f, (x, X[2*i], X[2*i+2])))
    #print(args)
    plot = sm.plot(*args, show=False)
    return fx, plot


def splain(X, Y):
    n = len(X)
    h = [0.0] * n
    I = [0.0] * n
    b = [0.0] * n
    lam = [0.0] * n
    C = [0.0] * n

    for i in range(1,n):
        h[i] = X[i] - X[i-1]
        I[i] = (Y[i] - Y[i-1]) / h[i]
    b[1] = -0.5 * (h[2] / (h[1] + h[2]))
    lam[1] = 1.5 * ((I[2] - I[1]) / (h[1] + h[2]))
    for i in range(3, n):
        b[i-1] = -(h[i] / (2 * h[i-1] + 2 * h[i] + h[i-1] * b[i-2]))
        lam[i-1] = (2 * I[i] - 3 * I[i-1] - h[i-1] * lam[i-2]) / \
            (2 * h[i-1] + 2*h[i] + h[i-1]*b[i-2])
    for i in range(n-1, 1, -1):
        C[i-1] = b[i-1] * C[i] + lam[i-1]
    A = [0.0] * n
    A[0] = Y[0]
    B = [0.0] * n
    D = [0.0] * n
    for i in range(1, n):
        A[i] = Y[i]
        B[i] = I[i] + (2/3) * (C[i] * h[i]) + (1/3) * (h[i] * C[i-1])
        D[i] = (C[i] - C[i-1]) / (3 * h[i])

    G = [[A[k] + B[k]*(x - X[k]) + C[k]*(x - X[k])**2 + D[k]*(x-X[k])**3] for k in range(1, n)]

    args = []
    for i, f in enumerate(G):
        args.append((*f, (x, X[i], X[i+1])))
    #print(args)
    plot = sm.plot(*args, show=False)
    return G, plot
    


def main():
    #X = [-1.5, -0.75, 0, 0.75, 1.5]
    #y = [-14,1014, -0.931596, 0, 0.931596, 14.1014]
    #Y = [sm.tan(-1.5), sm.tan(-0.75), sm.tan(0), sm.tan(0.75), sm.tan(1.5)]

    #X = [0, 1, 2, 3]
    #Y = [-2,-5,0,-4]

    #X = [-1, 0, 1, 2]
    #Y = [4, 2, 0, 1]

    #X = [0,1,2,3]
    #Y = [math.sin(30*i/180*math.pi) for i in X]

    #X = [0.1, 0.3, 0.5, 0.6, 0.9]
    #Y = [0.4, 0.6, 0.75, 0.75, 0.83]

    polinom = inter(X, Y)
    print("Полином Лежандра:", polinom)
    print("L4(x1+x2) =", polinom.subs(x, X[1]+X[2]))
    #sm.plot(polinom, (x,X[0], X[-1]))

    polinom2 = konechn(X,Y)
    print(tabulate(polinom2, tablefmt = 'fancy_grid'))

    polinom3, pol = razdel(X,Y)
    print(tabulate(polinom3, tablefmt = 'fancy_grid'))
    print("Полином Ньютона:", pol.n(2))
    print("N4(x1+x2) =", pol.subs(x, X[1]+X[2]))
    #sm.plot(pol, (x,X[0], X[-1]))

    fx, plot = ksn(X, Y)
    print(tabulate(fx, tablefmt='fancy_grid'))
    #plot.show()
    fx2, plot = ksk(X, Y)
    #plot.show()
    print(tabulate(fx2, tablefmt='fancy_grid'))

    G, plot = splain(X, Y)
    #plot.show()
    print(tabulate(G, tablefmt='fancy_grid'))



    args = []
    args.append((polinom, (x,X[0], X[-1])))
    args.append((pol, (x,X[0], X[-1])))
    for i, f in enumerate(fx):
        args.append((*f, (x, X[i], X[i+1])))
    for i, f in enumerate(fx2):
        args.append((*f, (x, X[2*i], X[2*i+2])))
    for i, f in enumerate(G):
        args.append((*f, (x, X[i], X[i+1])))
    #print(args)
    plot = sm.plot(*args, show=False, legend=True)
    plot[0].line_color = 'b'
    plot[1].line_color = 'b'
    plot[0].label = 'Полином Лежандра'
    plot[1].label = 'Полином Ньютона'
    counter = 2
    for i in range(len(fx)):
        plot[counter].line_color = 'g'
        plot[counter].label = ''
        counter += 1
    plot[counter-1].label = 'Кусочно-непрерывная'
    for i in range(len(fx2)):
        plot[counter].line_color = 'y'
        plot[counter].label = ''
        counter += 1
    plot[counter-1].label = 'Кусочно-квадратичная'
    for i in range(len(G)):
        plot[counter].line_color = 'r'
        plot[counter].label = ''
        counter += 1
    plot[counter-1].label = 'Кубический сплайн'
    
    plot.show()



if __name__ == '__main__':
    main()


    #x = sm.symbols('x')
    #a = sm.plot((sm.sqrt(x),(x,-5,5)), (x**2, (x,-2,2)), show=False, legend=True)
    #a[0].line_color = 'green'
    #a[0].label = '$\\frac{dy}{dx}$'
    #a[1].label = '$y$ y'
    #a.show()

    
    #print(min(timeit.repeat("f()", "from __main__ import f", number=50)))
