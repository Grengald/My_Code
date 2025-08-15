import numpy as np
import matplotlib.pyplot as plt

def function(x):
    return 1 / (2 + x**3)

def dif_function(x):
    return -(3*x**2)/(x**6+4*x**3+4)

def dif_dif_function(x):
    return (12*(x**4-x))/(x**9+6*x**6+12*x**3+8)

a = -1
b = 1
n_values = [10 * 5**i for i in range(8)]
err_sqrt = []
err_trap = []
err_trap_corrected = []  # ошибка для уточнённого метода трапеций
err_Simps = []
h_values = []

theor_int = 1.04154051729347627

for n in n_values:
    h = (b - a) / n
    h_values.append(h)
    x_source = np.linspace(a, b, n+1)
    y_source = function(x_source)
    
    sqrt_int = 0
    for z in range(len(y_source)-1):
        sqrt_int += (y_source[z+1]) * h
    err_sqrt.append(abs(sqrt_int - theor_int))
    
    trap_int = 0
    for x in range(len(y_source)-1):
        trap_int += h * (y_source[x] + y_source[x+1]) / 2
    err_trap.append(abs(trap_int - theor_int))
    
    # уточнённый метод трапеций
    # вычисляем среднее значение второй производной
    dif_dif_values = dif_dif_function(x_source)
    mean_dif_dif = np.mean(dif_dif_values)
    # вычитаем поправку
    correction = (b - a) / 12 * h**2 * mean_dif_dif
    trap_corrected_int = trap_int - correction
    err_trap_corrected.append(abs(trap_corrected_int - theor_int))
    
    Simps_int = 0
    for c in range(1, len(y_source)-1, 2):
        Simps_int += h * (y_source[c-1] + 4*y_source[c] + y_source[c+1]) / 3
    err_Simps.append(abs(Simps_int - theor_int))

plt.figure(figsize=(10, 6))
plt.loglog(h_values, err_sqrt, 'x-', label='Метод прямоугольников')
plt.loglog(h_values, err_trap, 'x-', label='Метод трапеций')
plt.loglog(h_values, err_trap_corrected, 'x-', label='Уточнённый метод трапеций')
plt.loglog(h_values, err_Simps, 'x-', label='Метод Симпсона')
plt.xlabel('Шаг сетки h')
plt.ylabel('Максимальная погрешность')
plt.legend()
plt.title('Зависимость погрешности от шага сетки')
plt.show()
