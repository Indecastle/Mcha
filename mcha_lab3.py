import numpy as np
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

at = a.transpose()
a_ = np.dot(at,a)
b_ = np.dot(at,b)
print("a*:", a_, sep='\n')
print("b*:", b_, sep='\n', end='\n\n')

u = np.eye(n)

u[(0,0)] = a_[(0,0)]**(1/2)
for j in range(1, n):
    u[0][j] = a_[0][j] / u[0][0]

for i in range(n):
    sum = 0
    for k in range(i):
        sum += u[k][i]**2
    u[i][i] = (a_[i][i] - sum) ** 0.5
    for j in range(i+1, n):
        sum = 0
        for k in range(i):
            sum += u[k][i] * u[k][j]
        u[i][j] = (a_[i][j] - sum) / u[i][i]

ut = np.transpose(u)
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

print("X =", x)
print("Det =", det)
#print(np.linalg.det(a_))
print("inv:", inv, sep='\n')
#print(np.linalg.inv(a_))