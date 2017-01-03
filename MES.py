import math


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


