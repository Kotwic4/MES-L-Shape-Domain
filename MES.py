import math
import operator
import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import scipy.interpolate as interp
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

def g(x, y):
    r = math.sqrt(x ** 2 + y ** 2)
    cos = y/r  # sin(fi + pi/2) = cos(fi)
    if y <= 0 or x <= 0:
        cos = -cos
    return (r ** 2) ** (1/3) * (cos ** 2) ** (1/3)

p = int(input("Enter p: "))
if p < 1:
    p = 1
    print("p must be >= 1")
    print("Assuming p = 1")

N = 2 * p + 1
h = (1/p)

MatrixA = [[0 for x in range(N*N)] for y in range(N*N)]

Matrix = [[2/3, -1/6, -1/3, -1/6],
          [-1/6, 2/3, -1/6, -1/3],
          [-1/3, -1/6, 2/3, -1/6],
          [-1/6, -1/3, -1/6, 2/3]]
MatrixB = [0 for x in range(N*N)]

for i in range(N-1):
    for j in range(N-1):
        _4 = i*N + j
        _3 = _4 + 1
        _1 = _4 + N
        _2 = _1 + 1
        temp = [[_4, 3], [_3, 2], [_2, 1], [_1, 0]]
        for x in temp:
            for y in temp:
                val = Matrix[x[1]][y[1]]*p*p
                MatrixA[x[0]][y[0]] += val

for i in range(p, N):
    for j in range(0, p+1):
        x = i * N + j
        for k in range(0, N*N):
            MatrixA[x][k] = 0
        MatrixA[x][x] = 1

for i in range(1, N-1):
    MatrixB[i] = (1/2) * h * (g(-1+(i*h), 1) + g(-1+((i+1)*h), 1))
for i in range(1, N-1):
    MatrixB[i*N+N-1] = (1/2) * h * (g(1, 1-(i*h)) + g(1, 1-((i+1)*h)))
for i in range(1, p):
    MatrixB[i*N] = (1/2) * h * (g(-1, 1-(i*h)) + g(-1, 1-((i+1)*h)))
for i in range(p+1, N-1):
    MatrixB[(N-1)*N+i] = (1/2) * h * (g(-1+(i*h), -1) + g(-1+((i+1)*h), -1))
MatrixB[0] = (1/2) * h * (g(-1+h, 1) + g(-1, 1-h))
MatrixB[N-1] = (1/2) * h * (g(1-h, 1) + g(1, 1-h))
MatrixB[N*N-1] = (1/2) * h * (g(-1+h, -1) + g(-1, -1+h))
MatrixZ = np.linalg.solve(MatrixA, MatrixB)
X = [(x % N)*h-1 for x in range(N*N)]
Y = [1-(y // N)*h for y in range(N*N)]
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

plotx,ploty = np.meshgrid(np.linspace(np.min(X),np.max(X),50),
                           np.linspace(np.min(Y),np.max(Y),50))

plotz = interp.griddata((X,Y),
                        MatrixZ.tolist(),
                        (plotx,ploty),
                        method='linear')

surf = ax.plot_surface(
    plotx, ploty, plotz,
    cstride=1, rstride=1,
    cmap=cm.coolwarm,
    antialiased=True,
    linewidth=0)

plt.show()



