import numpy as np
import matplotlib.pyplot as plt

def plot3(name, axisX, axisY, x1, y1, label1, x2, y2, label2, x3, y3, label3):
    plt.plot(x1, y1, color='blue', label=label1)
    plt.plot(x2, y2, color='red', label=label2)
    plt.plot(x3, y3, color='green', label=label3)
    plt.xlabel(axisX)
    plt.ylabel(axisY)
    plt.title(name)
    plt.grid(True)
    plt.legend()
    plt.show()

def function(x):
    return x / (1 + np.tan(x)**2)

def errordif(x, y, teordif):
    error = np.max(np.abs(y - teordif))
    return error

def errordif2(x, y, teordif2):
    error = np.max(np.abs(y - teordif2))
    return error

a = -1.5
b = 1.5
n_values = [10 * 5**i for i in range(10)]
err_right = []
err_cent = []
err_2 = []
err_4 = []
h_values = []

for n in n_values:
    h = (b - a) / n
    h_values.append(h)
    x = np.linspace(a - 2*h, b + 2*h, n + 5)
    y = function(x)

    def true_derivative(x):
        return (1 + np.tan(x)**2 - 2 * x * np.tan(x) * (1/np.cos(x))**2) / (1 + np.tan(x)**2)**2

    def true_second_derivative(x):
        cos_x = np.cos(x)
        sin_x = np.sin(x)
        tan_x = np.tan(x)
        sec_x2 = 1 / cos_x**2  # sec^2(x)
        tan_x2_1 = tan_x**2 + 1
    
        numerator = -2 * (cos_x**2 * (tan_x + x * sec_x2) * tan_x2_1**2
                        - x * tan_x * (4 * tan_x * tan_x2_1
                        - 2 * cos_x * sin_x * tan_x2_1**2))
        denominator = cos_x**4 * tan_x2_1**4
    
        term1 = (2 * tan_x) / (cos_x**2 * tan_x2_1**2)
        term2 = -(4 * tan_x**3) / (cos_x**2 * tan_x2_1**3)
        term3 = -(4 * tan_x) / (cos_x**2 * tan_x2_1**3)
    
        return numerator / denominator + term1 + term2 + term3

    teordif = true_derivative(x)
    teordif2 = true_second_derivative(x)

    right_dif = (y[3:n+4] - y[2:n+3]) / h
    err_right.append(errordif(x[2:n+3], right_dif, teordif[2:n+3]))
  
    centr_dif = (y[3:n+4] - y[1:n+2]) / (2*h)
    err_cent.append(errordif(x[2:n+3], centr_dif, teordif[2:n+3]))
  
    second_dif2 = (y[3:n+4] - 2*y[2:n+3] + y[1:n+2]) / h**2
    err_2.append(errordif2(x[2:n+3], second_dif2, teordif2[2:n+3]))
  
    second_dif4 = (16*y[3:n+4] + 16*y[1:n+2] - y[0:n+1] - y[4:n+5] - 30*y[2:n+3]) / (12*h**2)
    err_4.append(errordif2(x[2:n+3], second_dif4, teordif2[2:n+3]))
  
    #plot3(f'График первой производной при N={n} узлах', 'X', "Y'", x[2:n+3], teordif[2:n+3], 'Производная', x[2:n+3], right_dif, 'Правая разность', x[2:n+3], centr_dif, 'Центральная разность')
    #plot3(f'График второй производной при N={n} узлах', 'X', "Y''", x[2:n+3], teordif2[2:n+3], '2-ая производная', x[2:n+3], second_dif2, '2-ой порядок точности', x[2:n+3], second_dif4, '4-ый порядок точности')

plt.figure(figsize=(10, 6))
plt.loglog(h_values, err_right, 'x-', label='Правая разность', linestyle='dashed')
plt.loglog(h_values, err_cent, 'x-', label='Центральная разность', linestyle='dotted')
plt.loglog(h_values, err_2, 'x-', label='Вторая производная (2-й порядок)', linestyle='dashdot')
plt.loglog(h_values, err_4, 'x-', label='Вторая производная (4-й порядок)', linestyle='solid')
plt.xlabel('Шаг сетки h')
plt.ylabel('Максимальная погрешность')
plt.legend()
plt.title('Зависимость погрешности от шага сетки')
plt.axvline(x=h_values[np.argmin(err_4)], color='red', linestyle='dashed', label='Оптимальный шаг')
plt.legend()
plt.show()


