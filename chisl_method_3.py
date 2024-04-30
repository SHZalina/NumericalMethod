import numpy as np
alpha , beta, gamma = 3, 8, 2
def fx(x):
    return -12 * x**3 * (1 - x)**8 + 64 * x**4 * (1 - x)**7 - 6 * x * (1 - x)**8 + 48 * x**2 * (1 - x)**7 - 56 * x**3 * (1 - x)**6 - 56 * x**5 * (1 - x)**6 + (x**3 + x**4) * (1 - x)**8


def qx(x):
    return 1 + x

def ux(x):
  return x**3 * (1 - x)**8

def px(x):
    return 1 + x**2 

def Zeydel(n):
    preduu = np.zeros(n + 1)
    usl = 0
    kol = 0

    for i in range(n + 1):
        preduu[i] = 0

    while True:
        for i in range(n + 1):
            preduu[i] = y_iter[i]
        usl = -1

        for i in range(1, n):
            y_iter[i] = (fi[i] + ai[i + 1] * preduu[i + 1] + ai[i] * y_iter[i - 1]) / (ai[i + 1] + ai[i] + gi[i] * h * h)

        kol += 1

        for i in range(1, n):
            current_usl = abs(-ai[i + 1] * y_iter[i + 1] + ai[i + 1] * y_iter[i] + ai[i] * y_iter[i] - ai[i] * y_iter[i - 1] + gi[i] * h * h * y_iter[i] - fi[i])
            if current_usl > usl:
                usl = current_usl

        if usl <= 0.000001:
            break

    print("Метод Зейделя")
    print("{:<20} {:<25} {:<25} {:<25}".format("i*h", "yi(zeydel)", "yi^k", "|yi^k-y(zeydel)|"))

    for i in range(n + 1):
        print("{0:20} {1:25} {2:25} {3:25}".format(i * h, y_iter[i], arrY[i], arrY[i] - y_iter[i] ))

    print("Количество итерраций =", kol)
    print('\n')


def SOR(n):
    print("Метод верхней релаксации")
    predu = np.zeros(n + 1)
    usl = -1
    kol_iter = 100000
    w_temp = 0

    for i in range(n + 1):
        predu[i] = 0

    for om in np.arange(1.05, 2, 0.05):
        y_iter = np.zeros(n + 1)
        kol = 0

        while True:
            for i in range(n + 1):
                predu[i] = y_iter[i]
            usl = -1

            for i in range(1, n):
                y_iter[i] = (fi[i] + ai[i + 1] * predu[i + 1] + ai[i] * y_iter[i - 1]) / (ai[i + 1] + ai[i] + gi[i] * h * h) * om + (1 - om) * predu[i]

            for i in range(1, n):
                if abs(-ai[i + 1] * y_iter[i + 1] + ai[i + 1] * y_iter[i] + ai[i] * y_iter[i] - ai[i] * y_iter[i - 1] + gi[i] * h * h * y_iter[i] - fi[i]) > usl:
                    usl = abs(-ai[i + 1] * y_iter[i + 1] + ai[i + 1] * y_iter[i] + ai[i] * y_iter[i] - ai[i] * y_iter[i - 1] + gi[i] * h * h * y_iter[i] - fi[i])

            kol += 1

            if usl <= 0.000001:
                break

        if kol < kol_iter:
            w_temp = om
            kol_iter = kol

        print("{0:3} {1:3}".format("w =", om), end=", ")
        print("{0:12} {1:2}".format("iterations =", kol))

    y_iter = np.zeros(n + 1)
    kol = 0

    while True:
        for i in range(n + 1):
            predu[i] = y_iter[i]
        usl = -1

        for i in range(1, n):
            y_iter[i] = (fi[i] + ai[i + 1] * predu[i + 1] + ai[i] * y_iter[i - 1]) / (ai[i + 1] + ai[i] + gi[i] * h * h) * w_temp + (1 - w_temp) * predu[i]

        kol += 1

        for i in range(1, n):
            if abs(-ai[i + 1] * y_iter[i + 1] + ai[i + 1] * y_iter[i] + ai[i] * y_iter[i] - ai[i] * y_iter[i - 1] + gi[i] * h * h * y_iter[i] - fi[i]) > usl:
                usl = abs(-ai[i + 1] * y_iter[i + 1] + ai[i + 1] * y_iter[i] + ai[i] * y_iter[i] - ai[i] * y_iter[i - 1] + gi[i] * h * h * y_iter[i] - fi[i])

        if usl <= 0.000001:
            break

    print("Наименьшее количество итерраций при ω =", w_temp)
    print("Количество итераций =", kol)
    print("{0:20} {1:25} {2:25} {3:25}".format("i*h", "yi(relax)", "yi^k", "|yi^k-yi(relax)|"))
    for i in range(n + 1):
        print("{0:20} {1:25} {2:25} {3:25}".format(i * h, y_iter[i], arrY[i], arrY[i] - y_iter[i]))
    print('\n')

print("Введите n:")
n = int(input())
h = 1 / n
eps = h ** 6
x = np.zeros(n + 1)
ui = np.zeros(n + 1)
ai = np.zeros(n + 1)
gi = np.zeros(n + 1)
fi = np.zeros(n + 1)
alpha = np.zeros(n + 1)
betta = np.zeros(n + 1)
A = np.zeros(n + 1)
B = np.zeros(n + 1)
C = np.zeros(n + 1)
A1 = np.zeros((n, n))
arrY = np.zeros(n + 1)
y_iter = np.zeros(n + 1)
alpha[1] = 0
betta[1] = 0

arrY[0] = 0
arrY[n] = 0

for i in range(n + 1):
    ai[i] = px(i * h)
    gi[i] = qx(i * h)
    fi[i] = fx(i * h) * h ** 2

for i in range(n):
    x[i] = i * h
    ui[i] = ux(x[i])

for i in range(1, n):
    A[i] = -ai[i]
    B[i] = -(ai[i] + ai[i + 1] + h ** 2 * gi[i])
    C[i] = -ai[i + 1]

for i in range(n - 1):
    A1[i + 1][i] = A[i]

for i in range(n):
    A1[i][i] = B[i]

for i in range(n - 1):
    A1[i][i + 1] = C[i]

for i in range(1, n):
    alpha[i + 1] = C[i] / (B[i] - A[i] * alpha[i])
    betta[i + 1] = (A[i] * betta[i] - fi[i]) / (B[i] - A[i] * alpha[i])

print('alpha')
print (alpha)
for i in range(n - 1, 0, -1):
    arrY[i] = alpha[i + 1] * arrY[i + 1] + betta[i + 1]

print("Метод прогонки")
print("{:<20} {:<25} {:<25} {:<25}".format("i*h", "yi", "u(ih)", "|yi-u(ih)|"))

for i in range(n + 1):
    print("{0:20} {1:25} {2:25} {3:25}".format(i * h, arrY[i], ui[i], abs(arrY[i] - ui[i] )))
print('\n')


Zeydel(n)
SOR(n)
print(B)
print('\n')
print(C)