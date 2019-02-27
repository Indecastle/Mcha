import numpy as np
import math
import random
"""
a = np.array([[3, 1, 1],
              [1, 11, 1],
              [1, 0, 5]], float)
b = np.array([21, 21, 48], float)
"""
"""
a = np.array([[3.910, 0.129, 0.283, 0.107],
              [0.217, 4.691, 0.279, 0.237],
              [0.201, 0.372, 2.987, 0.421],
              [0.531, 0.196, 0.236, 5.032]], float)
b = np.array([0.395, 0.432, 0.127, 0.458], float)
"""
#n = len(b)


def getU(A):
    n = np.shape(A)[0]
    U = np.eye(n)
    U[(0,0)] = math.sqrt(A[(0,0)])
    U[0,1:] = [A[0][j] / U[0][0] for j in range(1, n)]

    for i in range(n):
        ss = sum([U[k][i]**2 for k in range(i)])
        U[i][i] = math.sqrt(A[i][i] - ss)
        for j in range(i+1, n):
            ss = sum([U[k][i] * U[k][j] for k in range(i)])
            U[i][j] = (A[i][j] - ss) / U[i][i]
    return U


def Eigenvalue(A: np.array, B: np.array):
    n = len(B)
    A = A.copy()
    B = B.copy()
    V = np.eye(n)
    rad = 0
    Uv = np.eye(n)
    
    OK = True
    while OK:
        i, j = 0, 1
        OK = False
        for ii in range(n):
            for jj in range(ii+1, n):
                if abs(A[i][j]) <= abs(A[ii][jj]) and abs(np.real(A[ii][jj])) > 0.01e-3:
                    i, j = ii, jj
                    OK = True

        U = np.eye(n)
        if A[j][j]-A[i][i] == 0:
            rad = math.pi/4
        else:
            rad = 0.5 * math.atan(2 * A[i][j] / (A[j][j]-A[i][i]))
        s = math.sin(rad)
        c = math.cos(rad)
        U[i][i] = U[j][j] = c
        U[i][j] = s
        U[j][i] = -s
        K = np.dot(np.transpose(U), A)
        A = np.dot(K, U)
        Uv = np.dot(Uv,U)
        #print(np.around(A, decimals=2))
        #input()
    print(np.around(A, decimals=2))
    eiges = [round(np.real(A[i][i]), 8) for i in range(n)]
    return (eiges, Uv)


def main():
    a = np.array([[3.910, 0.129, 0.283, 0.107],
                  [0.217, 4.691, 0.279, 0.237],
                  [0.201, 0.372, 2.987, 0.421],
                  [0.531, 0.196, 0.236, 5.032]], float)
    b = np.array([0.395, 0.432, 0.127, 0.458], float)
    U = getU(a)
    a_ = U.T @ U

    eiges, V = Eigenvalue(a_, b)

    print("собственные значения", eiges)
    print("собственные вектора:\n", V)

    print("\nпроверка:", *np.linalg.eig(a), sep='\n')



if __name__ == "__main__":
    main()
