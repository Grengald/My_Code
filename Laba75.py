import numpy as np
import matplotlib.pyplot as plt

# Исходные параметры
h = 0.05
n = int(1 / h)
tao = 0.05
T = 1
M = int(T / tao)
aa = 0.5

x = np.linspace(0, 1, n + 2)
t = np.linspace(0, T, M + 1)
X, T_grid = np.meshgrid(x, t)

# Точные и вспомогательные функции
uo = lambda x, t: x**2 + np.sinh(x * t)
f = lambda x, t: -2 + (2 * x**2 - t**2) * np.sinh(x * t)
u1 = lambda x: x**2
u2 = lambda x: x
u3 = lambda t: t
u4 = lambda t: 3 + np.sinh(t) + t * np.cosh(t)

# Проверка условия устойчивости
stability = (tao**2) / (h**2)
is_stable = stability <= 1

# Решение первого порядка
def solve_first_order():
    a = np.array([u1(i * h) for i in range(n + 2)])
    b = np.array([a[i] + tao * u2(i * h) for i in range(n + 2)])
    result = np.zeros((M + 1, n + 2))
    result[0] = a
    result[1] = b

    for k in range(2, M + 1):
        c = np.zeros(n + 2)
        for e in range(1, n + 1):
            c[e] = (
                2 * b[e]
                - a[e]
                + aa * (tao / h) ** 2 * (b[e + 1] - 2 * b[e] + b[e - 1])
                + tao**2 * f(e * h, (k - 1) * tao)
            )
        # Левое граничное условие 1 порядка
        c[0] = (u3(k * tao) + c[1] / h) / (1 + 1 / h)
        # Правое граничное условие 1 порядка
        c[-1] = (u4(k * tao) - c[-2] / h) / (1 - 1 / h)
        a, b = b, c
        result[k] = c
    return result

# Точное решение
def analytical_solution():
    return uo(X, T_grid)

# Норма Чебышёва
def chebyshev_norm(numerical, analytical):
    return np.max(np.abs(numerical - analytical))

# Решения
sol_first = solve_first_order()
u_exact = analytical_solution()
error = sol_first - u_exact
max_error = chebyshev_norm(sol_first, u_exact)

# Построение 3D-графиков
fig = plt.figure(figsize=(15, 6))

# Численное решение
ax1 = fig.add_subplot(121, projection='3d')
ax1.plot_surface(X, T_grid, sol_first, cmap='viridis')
ax1.set_title('Численное решение (1-й порядок)')
ax1.set_xlabel('x')
ax1.set_ylabel('t')
ax1.set_zlabel('u')

# Ошибка
ax2 = fig.add_subplot(122, projection='3d')
ax2.plot_surface(X, T_grid, error, cmap='inferno')
ax2.set_title('Ошибка (1-й порядок)')
ax2.set_xlabel('x')
ax2.set_ylabel('t')
ax2.set_zlabel('Ошибка')

plt.tight_layout()
plt.show()

is_stable, stability, max_error

