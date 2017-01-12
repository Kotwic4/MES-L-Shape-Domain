import math
import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import scipy.interpolate as interp
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

test = False
e = 0.00001


def g(x, y):
    if test:
        if (abs(y + 1.0) < e) and (x >= 0):
            return -x
        if (abs(x + 1.0) < e) and (y >= 0):
            return -y
        if abs(y - 1.0) < e:
            return x
        if abs(x - 1.0) < e:
            return y
    else:
        r = math.sqrt(x ** 2 + y ** 2)
        cos = x/r  # sin(fi + pi/2) = cos(fi)
        if y <= 0 or x <= 0:
            cos = -cos
        return (r ** 2) ** (1/3) * (cos ** 2) ** (1/3)


def getx(i):
    return (i % n_a)*h_a-1


def gety(i):
    return 1-(i // n_a)*h_b


def read_parameter(label):
    x = int(input("Enter " + label + ": "))
    if x < 1:
        x = 1
        print(label + "must be >= 1")
        print("Assuming " + label + " = 1")
    n = 2 * x + 1
    h = 1 / x
    d = h / 2
    return x, n, h, d


a, n_a, h_a, d_a = read_parameter("a")
b, n_b, h_b, d_b = read_parameter("b")
nab = n_a * n_b

MatrixA = [[0 for x in range(nab)] for y in range(nab)]
a1 = 2/3 * (a * b)
a2 = -1/6 * (a * b)
a3 = -1/3 * (a * b)
Matrix = [[a1, a2, a3, a2],
          [a2, a1, a2, a3],
          [a3, a2, a1, a2],
          [a2, a3, a2, a1]]
MatrixB = [0 for x in range(nab)]

for i in range(n_b-1):
    for j in range(n_a-1):
        _4 = i*n_a + j
        _3 = _4 + 1
        _1 = _4 + n_a
        _2 = _1 + 1
        temp = [[_4, 3], [_3, 2], [_2, 1], [_1, 0]]
        for x in temp:
            for y in temp:
                val = Matrix[x[1]][y[1]]
                MatrixA[x[0]][y[0]] += val

for i in range(b, n_b):
    for j in range(0, a+1):
        x = i * n_a + j
        for k in range(0, nab):
            MatrixA[x][k] = 0
            MatrixA[k][x] = 0
        MatrixA[x][x] = 1

for i in range(1, n_a-1):
    index = i
    x = getx(index)
    y = gety(index)
    MatrixB[index] = (1/2) * h_b * (g(x-d_a, y) + g(x+d_a, y))
for i in range(1, n_b-1):
    index = (i + 1) * n_a - 1
    x = getx(index)
    y = gety(index)
    MatrixB[index] = (1/2) * h_a * (g(x, y+d_b) + g(x, y-d_b))
for i in range(1, b):
    index = i*n_a
    x = getx(index)
    y = gety(index)
    MatrixB[index] = (1/2) * h_a * (g(x, y+d_b) + g(x, y-d_b))
for i in range(a+1, n_a-1):
    index = nab-n_a+i
    x = getx(index)
    y = gety(index)
    MatrixB[index] = (1/2) * h_b * (g(x-d_a, -1) + g(x+d_a, -1))
MatrixB[0] = (1/2) * (g(-1+d_a, 1) * h_b + g(-1, 1-d_b) * h_a)
MatrixB[n_a-1] = (1/2) * (g(1-d_a, 1) * h_b + g(1, 1-d_b) * h_a)
MatrixB[nab-1] = (1/2) * (g(1-d_a, -1) * h_b + g(1, -1+d_b) * h_a)
MatrixZ = np.linalg.solve(MatrixA, MatrixB)
for i in range(0,nab):
    MatrixZ[i] *= a * b
    if test:
        x = getx(i)
        y = gety(i)
        if x >= 0 or y >= 0:
            MatrixZ[i] -= x * y
X = [getx(i) for i in range(nab)]
Y = [gety(i) for i in range(nab)]
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')


plotx, ploty = np.meshgrid(np.linspace(np.min(X), np.max(X), n_a),
                           np.linspace(np.min(Y), np.max(Y), n_b))


plotz = interp.griddata((X, Y),
                        MatrixZ.tolist(),
                        (plotx, ploty),
                        method='linear')

surf = ax.plot_surface(
    plotx, ploty, plotz,
    cstride=1, rstride=1,
    cmap=cm.coolwarm,
    antialiased=True,
    linewidth=0)

plt.show()
