import numpy as np
import math
import random

a = np.array([[3, 1, 1, -3],
              [1, 11, 1, 2],
              [1, 0, 5, 7],
              [-3, 2, 7, 11]], float)
b = np.array([21, 21, 48, 61], float)
"""
a = np.array([[3.910, 0.129, 0.283, 0.107],
              [0.217, 4.691, 0.279, 0.237],
              [0.201, 0.372, 2.987, 0.421],
              [0.531, 0.196, 0.236, 5.032]], float)
b = np.array([0.395, 0.432, 0.127, 0.458], float)
"""
n = len(b)


def radfromgrad(grad):
    return math.pi*grad/180


def Eigenvalue(A: np.array, B: np.array):
    n = len(B)
    A = A.copy()
    B = B.copy()
    V = np.eye(n)
    rad = 0
    Uv = np.eye(n)
    for XXX in range(20): # для точности
        for i in range(n):
            for j in range(i+1, n):
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
    eiges = [round(np.real(A[i][i]), 8) for i in range(n)]
    return (eiges, Uv)


def main():
    a = np.array([[3.910, 0.129, 0.283, 0.107],
                  [0.217, 4.691, 0.279, 0.237],
                  [0.201, 0.372, 2.987, 0.421],
                  [0.531, 0.196, 0.236, 5.032]], float)
    b = np.array([0.395, 0.432, 0.127, 0.458], float)

    eiges, V = Eigenvalue(a, b)

    print("собственные значения", eiges)
    print("собственные вектора:\n", V)

    print("\nпроверка:", *np.linalg.eig(a), sep='\n')



if __name__ == "__main__":
    main()
