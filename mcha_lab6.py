import numpy as np
import sympy as sm
import matplotlib.pyplot as plt
import math
import timeit
from tabulate import tabulate

#x = (0.234, 0.649, 1.382, 2.672, 2.849)
#y = (0,511, 0.982, 2.411, 3.115, 4.184)
x = sm.symbols('x')
plt.style.use('ggplot')
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
    




def main():
    X = [-1.5, -0.75, 0, 0.75, 1.5]
    #y = [-14,1014, -0.931596, 0, 0.931596, 14.1014]
    Y = [sm.tan(-1.5), sm.tan(-0.75), sm.tan(0), sm.tan(0.75), sm.tan(1.5)]

    #X = [0, 1, 2, 3]
    #Y = [-2,-5,0,-4]

    #X = [-1, 0, 1, 2]
    #Y = [4, 2, 0, 1]

    #X = [0,1,2,3]
    #Y = [math.sin(30*i/180*math.pi) for i in X]

    polinom = inter(X, Y)
    print(polinom)
    print("L4(x1+x2) =", polinom.subs(x, X[1]+X[2]))
    #sm.plot(polinom, (x,-1.5, 1.5), ylim=(-15,15))

    table = konechn(X,Y)
    print(tabulate(table, tablefmt = 'fancy_grid'))

    table, pol = razdel(X,Y)
    print(tabulate(table, tablefmt = 'fancy_grid'))
    print("Полином Ньютона:", pol.n(2))
    print("N4(x1+x2) =", pol.subs(x, X[1]+X[2]))
    #sm.plot(pol, (x,-1.5, 1.5), ylim=(-15,15))



if __name__ == '__main__':
    main()


    #x = sm.symbols('x')
    #a = sm.plot((sm.sqrt(x),(x,-5,5)), (x**2, (x,-2,2)), show=False, legend=True)
    #a[0].line_color = 'green'
    #a[0].label = '$\\frac{dy}{dx}$'
    #a[1].label = '$y$ y'
    #a.show()

    
    #print(min(timeit.repeat("f()", "from __main__ import f", number=50)))
