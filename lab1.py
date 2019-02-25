import matplotlib.pyplot as plt
import numpy
import math

xr = numpy.arange(0.25,0.46,0.01)
plt.plot(xr,[0.85,0.825,0.8,0.775,0.75,0.725,0.675,0.65,0.6,0.575,0.55,0.5,0.46,0.44,0.42,0.41,0.4,0.39,0.38,0.37,0.36],label='15V')
plt.plot(xr,[1.275,1.25,1.25,1.225,1.2,1.15,1.15,1.1,1.05,1,0.9,0.875,0.825,0.775,0.7,0.67,0.65,0.625,0.6,0.6,0.575],label='20V')
plt.plot(xr,[1.75,1.725,1.7,1.7,1.7,1.65,1.65,1.6,1.575,1.5,1.45,1.375,1.3,1.225,1.15,1.1,1.025,0.95,0.9,0.85,0.825],label='25V')
plt.xlabel('Ток соленоида', size=14)
plt.ylabel('Ток анода', size=14)


u = 1
u0 = 4*math.pi * 10**(-7)
l = 0.07
N = 2000
Ic = 0.34
ra = 0.31 * 10**(-3)
Bk = (Ic * N) / (2*l) * u * u0 * 1.6
m = 9.1*10**(-31)

q = lambda u: m*8*u/(ra**2*Bk**2)

print('Элементарный заряд при 15, 20 и 25 вольт соответственно\n')
print(q(15))
print(q(20))
print(q(25))




def diff(b):
    a = []
    for i in range(20):
        a.append(abs(b[i+1] - b[i]))
    return a

a = [0.85,0.825,0.8,0.775,0.75,0.725,0.675,0.65,0.6,0.575,0.55,0.5,0.46,0.44,0.42,0.41,0.4,0.39,0.38,0.37,0.36]
b = [1.275,1.25,1.25,1.225,1.2,1.15,1.15,1.1,1.05,1,0.9,0.875,0.825,0.775,0.7,0.67,0.65,0.625,0.6,0.6,0.575]
c = [1.75,1.725,1.7,1.7,1.7,1.65,1.65,1.6,1.575,1.5,1.45,1.375,1.3,1.225,1.15,1.1,1.025,0.95,0.9,0.85,0.825]
xr = numpy.arange(0.25,0.45,0.01)

plt.plot(xr,diff(a), label = 'ΔIa')
plt.plot(xr,diff(b))
plt.plot(xr,diff(c))
plt.legend()
plt.savefig('graph.png',format='png')
plt.show()