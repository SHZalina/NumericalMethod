import numpy as np
import math
import matplotlib.pyplot as plt

def Taylor(x, e):
    n = 0
    a = -x**2 / 4
    sum = a
    while abs(a) >= e:
        n += 1
        a *= -x**2 * n / ((n+1)*(2*n+2)*(2*n+1))
        sum += a
    return sum

def f(t):
    if t == 0: return 0
    else: return (math.cos(t) - 1) / t

#в эту функцию будут передаваться только n = 1 и n = 1024, так как в данной составной квадратурной формуле мы останавливаем вычисления при n = 1024 (можно не рассматривать остальные n, как я поняла)
def left_riemann_sum(a, b, n):
    h = (b - a) / n
    x_i = [a + i * h for i in range(0, n+1)]
    y_n = 0
    if n == 1:
        y_n = (b - a) * f(a)
    else:
        for i in range(n):
             y_n += f(x_i[i]) * h
    return y_n

def center_riemann_sum(a, b, n):
    if n == 1:
        y_n = (b - a) * f((a + b) / 2)
        return y_n
    
    h = (b - a) / n
    y_n = 0
    e = 10 ** (-6)
    h = (b - a) / n 
    x_i = [a + i * h for i in range(0, n+1)]
    for i in range(0, n):
         y_n += f((x_i[i] + x_i[i+1]) / 2) * h

    y_2n = 0
    h_dop = (b - a) / (2 * n)
    x_i_dop = [a + i * h_dop for i in range(0, 2*n+1)]
    for i in range(0, 2*n):
         y_2n += f((x_i_dop[i] + x_i_dop[i+1]) / 2) * h_dop

    else:
        while (abs(y_n - y_2n) > e) and (n != 1024):
            y_n = y_2n
            n *= 2
            h = (b - a) / (2 * n)
            x_i = [a + i * h for i in range(0, 2 * n + 1)]
            y_2n = 0
            for i in range(0, 2 * n):
                y_2n += f((x_i[i] + x_i[i+1]) / 2) * h

        return y_n, n

# квадратурная ф-ла симпсона по трем узлам + составная кв.ф.симпсона
def simpsons_sum(a, b, n):
    if n == 1:
        y_n = (b - a)/6 * (f(a) + 4 * f((a + b) / 2) + f(b))
        return y_n
    
    h = (b - a) / n
    y_n = 0
    e = 10 ** (-6)
    h = (b - a) / n
    x_i = [a + i * h for i in range(0, n+1)]
    for i in range(0, n):
         y_n += (f(x_i[i]) + 4 * f((x_i[i] + x_i[i+1]) / 2) + f(x_i[i+1])) * h/6

    y_2n = 0
    h_dop = (b - a) / (2 * n)
    x_i_dop = [a + i * h_dop for i in range(0, 2*n+1)]
    for i in range(0, 2*n):
         y_2n += (f(x_i_dop[i]) + 4 * f((x_i_dop[i] + x_i_dop[i+1]) / 2) + f(x_i_dop[i+1])) * h_dop/6

    else:
        while (abs(y_n - y_2n) > e) and (n != 1024):
            y_n = y_2n
            n *= 2
            h = (b - a) / (2 * n)
            x_i = [a + i * h for i in range(0, 2 * n + 1)]
            y_2n = 0
            for i in range(0, 2 * n):
                y_2n += (f(x_i[i]) + 4 * f((x_i[i] + x_i[i+1]) / 2) + f(x_i[i+1])) * h/6

        return y_n, n

# квадратураня ф-ла гаусса + составная кв.ф.гаусса
def gauss_quadrature(a, b, n):
    if n == 1:
        x1 = a + (b - a) / 2 * (1 - 1 / math.sqrt(3))
        x2 = a + (b - a) / 2 * (1 + 1 / math.sqrt(3))
        y_n = (b - a) / 2 * (f(x1) + f(x2))
        return y_n
    
    z = 1 / math.sqrt(3)
    h = (b - a) / n
    y_n = 0
    e = 10 ** (-6)
    h = (b - a) / n
    x_i = [a + i * h for i in range(0, n+1)]
    for i in range(0, n):
         y_n += (f(x_i[i] + h / 2 * (1 - z)) + f(x_i[i] + h / 2 * (1 + z))) * h/2

    y_2n = 0
    h_dop = (b - a) / (2 * n)
    x_i_dop = [a + i * h_dop for i in range(0, 2*n+1)]
    for i in range(0, 2*n):
         y_2n +=  (f(x_i_dop[i] + h / 2 * (1 - z)) + f(x_i_dop[i] + h / 2 * (1 + z))) * h_dop/2

    else:
        while (abs(y_n - y_2n) > e):
            y_n = y_2n
            n *= 2
            h = (b - a) / (2 * n)
            x_i = [a + i * h for i in range(0, 2 * n + 1)]
            y_2n = 0
            for i in range(0, 2 * n):
                y_2n += (f(x_i[i] + h / 2 * (1 - z)) + f(x_i[i] + h / 2 * (1 + z))) * h/2
        return y_n, n

a = 0.4
b = 4
n = 11
h = (b - a) / (n - 1)
e = 10 ** (-6)
x_values = np.empty(n, dtype = float)

for i in range(n):
    x_values[i] = a + i * h
print('Иксы:', x_values)

f_values = np.zeros(shape=n) #массив значений вычисления ряда Тейлора
for i in range(n):
    f_values[i] = round(Taylor(x_values[i], e), 6)
print('Значения ряда Тейлора в иксах:', f_values)

#Таблица для формулы левых прямоугольников
lrs_1 = []
lrs_N = []
for i in range(n):
    lrs = left_riemann_sum(0, x_values[i], 1)
    lrs_1.append(round(lrs, 15))        
    lrs_S = left_riemann_sum(0, x_values[i], 1024)
    lrs_N.append(round(lrs_S, 15))

print ('Таблица для формулы левых прямоугольников')
print ('x', '\t', '   S(x) ', '\t', ' Y1  ', '\t', '|S - Y1| ', '\t\t', 'Yn', '\t\t', '    |S - Yn|', '\t\t', ' n')
for i in range(n):
    print(round(x_values[i], 2), '\t', f_values[i], '\t', lrs_1[i], '\t', abs(f_values[i] - lrs_1[i]), '\t', lrs_N[i], '\t', round(abs(f_values[i] - lrs_N[i]), 15), '\t', '1024')
print('\n')

#Таблица для формулы центральных прямоугольников
сrs_1 = []
сrs_N = []
crs_counts = []

for i in range(0, n):
    сrs = center_riemann_sum(0, x_values[i], 1)
    сrs_1.append(round(сrs, 15))      
    
    сrs_S, crs_n = center_riemann_sum(0, x_values[i], 2)
    сrs_N.append(round(сrs_S, 15))
    crs_counts.append(crs_n)
   
print ('Таблица для формулы центральных прямоугольников')
print ('x', '\t', '   S(x) ', '\t\t', '  Y1  ', '\t', '   |S - Y1| ', '\t\t\t', 'Yn', '\t\t', '  |S - Yn|', '\t\t', 'n')
for i in range(0, n):
    print(round(x_values[i], 2), '\t', f_values[i], '\t', сrs_1[i], '\t', round(abs(f_values[i] - сrs_1[i]), 15), '\t', сrs_N[i], '\t', round(abs(f_values[i] - сrs_N[i]), 15), '\t', crs_counts[i])
print('\n')

#Таблица для формулы Симпсона
sim_1 = []
sim_N = []
sim_counts = []

for i in range(0, n):
    sim = simpsons_sum(0, x_values[i], 1)
    sim_1.append(round(sim, 15))      
    
    sim_S, sim_n = simpsons_sum(0, x_values[i], 2)
    sim_N.append(round(sim_S, 15))
    sim_counts.append(sim_n)
   
print ('Таблица для формулы Симпсона')
print ('x', '\t', '   S(x) ', '\t\t', '  Y1  ', '\t', '   |S - Y1| ', '\t\t\t', 'Yn', '\t\t', '  |S - Yn|', '\t\t', 'n')
for i in range(0, n):
    print(round(x_values[i], 2), '\t', f_values[i], '\t', sim_1[i], '\t', round(abs(f_values[i] - sim_1[i]), 15), '\t', sim_N[i], '\t', round(abs(f_values[i] - sim_N[i]), 15), '\t', sim_counts[i])
print('\n')

#Таблица для формулы Гаусса
gau_1 = []
gau_N = []
gau_counts = []

for i in range(0, n):
    gau = gauss_quadrature(0, x_values[i], 1)
    gau_1.append(round(gau, 15))      
    
    gau_S, gau_n = gauss_quadrature(0, x_values[i], 2)
    gau_N.append(round(gau_S, 15))
    gau_counts.append(gau_n)
   
print ('Таблица для формулы Гаусса')
print ('x', '\t', '   S(x) ', '\t\t', '  Y1  ', '\t', '   |S - Y1| ', '\t\t\t', 'Yn', '\t\t', '  |S - Yn|', '\t\t', 'n')
for i in range(0, n):
    print(round(x_values[i], 2), '\t', f_values[i], '\t', gau_1[i], '\t', round(abs(f_values[i] - gau_1[i]), 15), '\t', gau_N[i], '\t', round(abs(f_values[i] - gau_N[i]), 15), '\t', gau_counts[i])
print('\n')