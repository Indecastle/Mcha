import numpy as np
import sys
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
a = np.array([[3, 1, 1],
              [1, 4, 1],
              [1, 1, 5]], float)
b = np.array([1, -1, 1], float)

n = len(b)
E = 0.001
if len(sys.argv) > 1:
    E = float(sys.argv[1])


def MPI(a: np.array, b: np.array, E: float):
    a = a.copy()
    b = b.copy()
    n = len(b)
    x1 = list(range(n))
    x01 = [b[i] for i in range(n)]
    h = 1  # counter iteration

    #for i in range(n):
    #    sum1 = 0
    #    for num in a[i]:
    #        sum1 += num
    #    sum1 -= a[i][i]
    #    if a[i][i] < sum1:
    #        print("не выполняется условие сходимости!")
    #        return tuple(), 0
    #print("выполняется условие сходимости!")

    for i in range(n):
        ii = a[i][i]
        for j in range(n):
            a[(i, j)] /= ii
        b[i] /= ii
        a[(i, i)] = 0

    #eig = [(np.real(f) < 1.0) for f in np.linalg.eig(a)[0]]
    #if all(eig):
    #    print("выполняется условие сходимости!")
    #else:
    #    print("не выполняется условие сходимости!")
    #    return tuple(), 0

    if not (np.linalg.norm(a, np.inf) < 1 or np.linalg.norm(a, 1) < 1 or np.linalg.norm(a) < 1):
        print("не выполняется условие сходимости!")
        print("Нормы:")
        print(np.linalg.norm(a, np.inf), np.linalg.norm(a, 1), np.linalg.norm(a))
        return tuple(), 0
    else:
        print("выполняется условие сходимости!")

    while True:
        for i in range(n):
            x1[i] = b[i]
            for k in range(n):
                x1[i] -= a[i][k] * x01[k]
        # print(x1)
        if all(abs(x1[i] - x01[i]) < E for i in range(n)):
            break
        x01 = list(x1)
        h += 1
    # print(h)
    return (x1, h)


def Zeid(a: np.array, b: np.array, E: float):
    a = a.copy()
    b = b.copy()
    n = len(b)
    x1 = list(range(n))
    x01 = [b[i] for i in range(n)]
    h = 1  # counter iteration

    for i in range(n):
        sum1 = 0
        for num in a[i]:
            sum1 += abs(num)
        sum1 -= abs(a[i][i])
        if a[i][i] < sum1:
            print("не выполняется условие сходимости!")
            return tuple(), 0
    else:
        print("выполняется условие сходимости!")

    for i in range(n):
        ii = a[i][i]
        for j in range(n):
            a[(i, j)] /= ii
        b[i] /= ii
        a[(i, i)] = 0

    while True:
        for i in range(n):
            x1[i] = b[i]
            for k in range(n):
                if k > i:
                    x1[i] -= a[i][k] * x01[k]
                else:
                    x1[i] -= a[i][k] * x1[k]
        # print(x1)
        if all(abs(x1[i] - x01[i]) < E for i in range(n)):
            break
        x01 = list(x1)
        h += 1
    #print(h)
    return (x1, h)

print("\nМПИ: ")
solve = MPI(a, b, E)
print("Итераций: " + str(solve[1]))
for t in enumerate(solve[0]):
    print(f"x{t[0]} = {t[1]}")

print("\nМетод Зейделя: ")
solve = Zeid(a, b, E)
print("Итераций: " + str(solve[1]))
for t in enumerate(solve[0]):
    print(f"x{t[0]} = {t[1]}")