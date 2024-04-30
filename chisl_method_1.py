import numpy as np
import math
import matplotlib.pyplot as plt

def basic_polynom(x, i, x_values):
    divider = 1
    res = 1
    for j in range(len(x_values)):
        if j != i:
            res *= (x-x_values[j])
            divider *= (x_values[i]-x_values[j])
    return res/divider

def Lagrange_polynom(x, x_values, f_values):
    res = 0
    for i in range(len(x_values)):
        res += f_values[i]*basic_polynom(x, i, x_values)
    return res

def Taylor(x, e):
    n = 0
    a = -x**2 / 4
    sum = a
    while abs(a) >= e:
        n += 1
        a *= -x**2 * n / ((n+1)*(2*n+2)*(2*n+1))
        sum += a
    return sum

a = 0.4
b = 4
n = 6
h0 = (b - a)/(n - 1)
e = 10**(-6)
x_values = np.empty(n, dtype = float)

for i in range(n):
    x_values[i] = a+i*h0

print('Равнораспределенные узлы:', x_values)

f_values = np.zeros(shape=n) #массив значений вычисления ряда Тейлора
for i in range(n):
    f_values[i] = round(Taylor(x_values[i], e), 6)
print('Значения ряда Тейлора в равнораспределленых узлах:', f_values)


x1 = 0.8
x2 = 2.0
x3 = 3.2
print ('x1:', x1,', f1:', round(Taylor(x1, e), 6))
print ('x2:', x2,', f2:', round(Taylor(x2, e), 6))
print ('x3:', x3,', f1:', round(Taylor(x3, e), 6))

n_Lagrange = 11
h = (b - a)/(n_Lagrange - 1)
x_Lagrange = np.empty(n_Lagrange, dtype = float)

for i in range(n_Lagrange):
    x_Lagrange[i] = a+i*h
print('x_Lagrange:', x_Lagrange)

f_Lagrange = np.zeros(shape=n_Lagrange) #массив значений вычисления ряда Тейлора для новых узлов при n=11
for i in range(n_Lagrange):
    f_Lagrange[i] = round(Taylor(x_Lagrange[i], e), 6)
print('f_Lagrange:', f_Lagrange)

plt.plot(x_Lagrange, f_Lagrange)
plt.xlabel(r'$xi$')
plt.ylabel(r'$f(xi)$')
plt.title('Табулирование Ci(x) с использованием ряда Тейлора')
plt.grid(True)
plt.show()

lag_pol1 = []
for x in x_Lagrange:
    lag_pol1.append(Lagrange_polynom (x, x_values, f_values))
print('Ln:', lag_pol1)

for i in range(n_Lagrange):
    print('x_Lagrange:', round(x_Lagrange[i], 2), ', Ln:', lag_pol1[i])

pogr = []
for i in range(n_Lagrange):
    pogr.append(math.fabs(f_Lagrange[i] - lag_pol1[i]))

print('Погрешности интерполирования при равнораспределленых узлах:', pogr)
print('Максимальная погрешность интерполирования? при равнораспределленых узлах:', max(pogr))
print('\n')

plt.plot(x_Lagrange, pogr)
plt.xlabel(r'$xi$')
plt.ylabel(r'$pogr$')
plt.title('Погрешность между функциями Ci(x) и L_n(x)')
plt.grid(True)
plt.show()


x_Cheb = np.empty(n, dtype = float)
for i in range(n):
    x_Cheb[i] = (b - a)/2 * math.cos((2*i + 1)/(2*n)*math.pi) + (a + b)/2
print('Узлы Чебышёва:', x_Cheb)

f_Cheb = np.zeros(shape=n) #массив значений вычисления ряда Тейлора
for i in range(n):
    f_Cheb[i] = round(Taylor(x_Cheb[i], e), 6)
print('Значения ряда Тейлора в узлах Чебышёва:', f_Cheb)

lag_pol2 = []
for x in x_Lagrange:
    lag_pol2.append(Lagrange_polynom (x, x_Cheb, f_Cheb))
print('Ln для узлов Чебышёва:', lag_pol2)

for i in range(n_Lagrange):
    print('x_Lagrange для узлов Чебышёва:', round(x_Lagrange[i], 2), ', Ln:', lag_pol2[i])

pogr2 = []
for i in range(n_Lagrange):
    pogr2.append(math.fabs(f_Lagrange[i] - lag_pol2[i]))

print('Погрешности интерполирования при узлах Чебышёва:', pogr2)
print('Максимальная погрешность интерполирования? при узлах Чебышёва:', max(pogr2))
print('\n')


plt.plot(x_Lagrange, pogr2)
plt.xlabel(r'$xi$')
plt.ylabel(r'$pogr2$')
plt.title('Погрешность между функциями Ci(x) и L_n(x)')
plt.grid(True)
plt.show()

pogr_ravn = []
pogr_cheb = []
nodes = [6, 12, 15, 18, 25, 36, 42, 50, 53]
for n0 in nodes:
    h0 = (b - a)/(n0 - 1)
    x_nodes = np.empty(n0, dtype = float) #хранятся новые значения n узлов
    for i in range(n0):
        x_nodes[i] = a+i*h0

    #print ('n0:', n0, '  x:', x_nodes)

    f_nodes = np.zeros(shape=n0) #массив значений вычисления ряда Тейлора для значений n узлов
    for i in range(n0):
        f_nodes[i] = round(Taylor(x_nodes[i], e), 6)
    #print ('n0:', n0, '  f:', f_nodes)

    lag_pol_nodes = []
    for x in x_Lagrange:
        lag_pol_nodes.append(Lagrange_polynom (x, x_nodes, f_nodes))
    #print('n0:', n0, 'Ln:', lag_pol_nodes)

    pogr_nodes = []
    for i in range(n_Lagrange):
        pogr_nodes.append(math.fabs(f_Lagrange[i] - lag_pol_nodes[i]))
   
    #print("Погрешности интерполирования при равнораспределленых узлах: ", pogr_nodes)
    #print("Максимальная погрешность интерполирования при равнораспределленых узлах: ", max(pogr_nodes))

    pogr_ravn.append(max(pogr_nodes))

    x_cheb = np.empty(n0, dtype = float)
    for i in range(n0):
        x_cheb[i] = (b - a)/2 * math.cos((2*i + 1)/(2*n0)*math.pi) + (a + b)/2
    #print('n0:', n0, '   Узлы Чебышёва:', x_cheb)

    f_cheb = np.zeros(shape=n0) #массив значений вычисления ряда Тейлора для узлов Чебышева
    for i in range(n0):
        f_cheb[i] = round(Taylor(x_cheb[i], e), 6)
    #print('n0:', n0, 'Значения ряда Тейлора в узлах Чебышёва:', f_cheb)

    lag_pol_cheb = []
    for x in x_Lagrange:
        lag_pol_cheb.append(Lagrange_polynom (x, x_cheb, f_cheb))
    #print('Ln для узлов Чебышёва:', lag_pol_cheb)

    pogr_ch = []
    for i in range(n_Lagrange):
        pogr_ch.append(math.fabs(f_Lagrange[i] - lag_pol_cheb[i]))

    #print('Погрешности интерполирования при узлах Чебышёва:', pogr_ch)
    #print('Максимальная погрешность интерполирования? при узлах Чебышёва:', max(pogr_ch))
    pogr_cheb.append(max(pogr_ch))
    #print('\n')


print('Погрешности интерполирования при n равнораспределленых узлах: ')
for i in range (len(nodes)):
    print("n = ", nodes[i], " e(x)= ", pogr_ravn[i])

plt.plot(nodes, pogr_ravn)
plt.xlabel(r'$nodes$')
plt.ylabel(r'$pogr_ravn$')
plt.grid(True)
plt.show()

print('\n')

print('Погрешности интерполирования при узлах Чебышёва: ')
for i in range (len(nodes)):
    print("n = ", nodes[i], " e(x)= ", pogr_cheb[i])

plt.plot(nodes, pogr_cheb)
plt.xlabel(r'$nodes$')
plt.ylabel(r'$pogr_cheb$')
plt.grid(True)
plt.show()