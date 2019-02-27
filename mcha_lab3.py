import numpy as np
import math
"""
a = np.array([[3, 0, 0, -3],
              [0, 11, 1, 2],
              [1, 0, 5, 0],
              [-1, 2, 7, 11]], float)
b = np.array([21, 21, 48, 61], float)
"""
a = np.array([[3.910, 0.129, 0.283, 0.107],
              [0.217, 4.691, 0.279, 0.237],
              [0.201, 0.372, 2.987, 0.421],
              [0.531, 0.196, 0.236, 5.032]], float)
b = np.array([0.395, 0.432, 0.127, 0.458], float)

n = len(b)

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

bg = False
for i in range(n):
    for j in range(n):
        if a[i][j] != a[j][i]:
            bg = True
at = a.transpose()
a_ = a
b_ = b
if bg == True and 0:
    #at = a.transpose()
    #a_ = np.dot(at,a)
    u = getU(a)
    a_ = u.T @ u
    #b_ = np.dot(at,b)
print("a*:", a_, sep='\n')
print("b*:", b_, sep='\n', end='\n\n')


u = getU(a)
ut = u.T
y = np.linalg.solve(ut, b_)
x = np.linalg.solve(u, y)

e = np.eye(n)
inv = np.eye(n)

det = 1
for i in range(n):
    det *= u[i][i] ** 2

#for i in range(n):
#    inv[:,i] =  np.linalg.solve(a, e[:,i]) # тоже работает

for i in range(n):
    y = np.linalg.solve(ut, e[:,i])
    inv[:, i] = np.linalg.solve(u, y)

#print("U =", ut@u, sep='\n', end='\n\n')
print("X =", x)
print("Det =", det)
#print(np.linalg.det(a))
print("inv:", inv, sep='\n')
#print(np.linalg.inv(a))