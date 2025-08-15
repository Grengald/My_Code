import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return np.cos(np.exp(x / 2) / 25)

def Lagrange(x, y, x0):
    L = 0
    n = len(x)
    for i in range(n):
        l_i = 1
        for j in range(n):
            if i != j:
                l_i *= (x0 - x[j]) / (x[i] - x[j])
        L += y[i] * l_i
    return L

def error(x1, y1):
    errors = [abs(y1[i] - f(x1[i])) for i in range(len(x1))]
    return max(errors)

a, b = 0, 10
N_values = [5, 10, 50, 100, 150]
uniform_errors = []
chebyshev_errors = []

for N in N_values:
    x = np.linspace(a, b, N)
    y = f(x)
    x1 = (x[:-1] + x[1:]) / 2
    y1 = np.array([Lagrange(x, y, xi) for xi in x1])
    max_error_uniform = error(x1, y1)
    uniform_errors.append(max_error_uniform)
    
    k = np.arange(1, N + 1)
    x_cheb = 0.5 * (a + b) + 0.5 * (b - a) * np.cos((2 * k - 1) * np.pi / (2 * N))
    y_cheb = f(x_cheb)
    x1_cheb = (x_cheb[:-1] + x_cheb[1:]) / 2
    y1_cheb = np.array([Lagrange(x_cheb, y_cheb, xi) for xi in x1_cheb])
    max_error_cheb = error(x1_cheb, y1_cheb)
    chebyshev_errors.append(max_error_cheb)
    

    x_t = np.linspace(a, b, 1000)
    y_t = f(x_t)
    #plt.ylim(-1, 1)
    plt.figure(figsize=(8, 5))
    #plt.ylim(-1, 1)
    plt.plot(x_t, y_t, label=r'$f(x) = \cos\left(\frac{e^{x/2}}{25}\right)$', color='b')
    plt.scatter(x1, y1, color='g', label='Узлы интерполяции', zorder=3)
    plt.plot(x1, y1, 'r-', label='Интерполированная функция', zorder=2)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(f'Интерполяция Лагранжа при N={N} (равномерная сетка)')
    plt.legend()
    plt.grid(True)
    plt.show()
    
    plt.figure(figsize=(8, 5))
    plt.plot(x_t, y_t, label=r'$f(x) = \cos\left(\frac{e^{x/2}}{25}\right)$', color='b')
    plt.scatter(x1_cheb, y1_cheb, color='g', label='Узлы интерполяции (Чебышев)', zorder=3)
    plt.plot(x1_cheb, y1_cheb, 'm-', label='Интерполированная функция (Чебышев)', zorder=2)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(f'Интерполяция Лагранжа при N={N} (узлы Чебышева)')
    plt.legend()
    plt.grid(True)
    plt.show()


plt.figure(figsize=(10, 6))
plt.plot(N_values, uniform_errors, label='Равномерная сетка')
plt.plot(N_values, chebyshev_errors, label='Узлы Чебышева')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Количество узлов N')
plt.ylabel('Максимальная погрешность')
plt.title('Зависимость максимальной погрешности интерполяции от количества узлов')
plt.legend()
plt.grid(True)
plt.show()




        
