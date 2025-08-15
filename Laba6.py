import numpy as np
import matplotlib.pyplot as plt

def u0(x):
    return x * np.sin(x) + np.cos(x)
def du0(x):
    return x * np.cos(x)

def p(x):
    return np.tan(x)

def q(x):
    return x

def f(x):
    return (1 + x) * np.cos(x) + x**2 * np.sin(x)

def f2(a_arr, b_arr, c_arr, f_arr):
    length = len(a_arr)
    D = np.zeros(length)
    B = np.zeros(length)

    D[0] = f_arr[0] / b_arr[0]
    B[0] = -c_arr[0] / b_arr[0]

    for i in range(1, length):
        a_k = a_arr[i]; b_k = b_arr[i]; c_k = c_arr[i]; f_k = f_arr[i]
        denom = b_k + a_k * B[i-1]
        D[i] = (f_k - a_k * D[i-1]) / denom
        B[i] = -c_k / denom

    B[-1] = 0
    res = np.zeros(length)
    res[-1] = D[-1]
    for i in range(length - 2, -1, -1):
        res[i] = D[i] + B[i] * res[i+1]
    return res

# Граничные условия
# Левое условие: u(0) + u'(0) = 1 (смешанное)
alfa_1 = 1; betta_1 = 1; gamma1 = 1
# Правое условие: u(1) = 1.3818 (Дирихле)
alfa_2 = 1; betta_2 = 0; gamma2 = u0(1)

N = [10, 20, 40, 80, 160]
steps = []
x_arrays = []
solutions1 = []
solutions2 = []
errors_1 = []
errors_2 = []

def get_subplots(n_rows=1, n_cols=1, title=''):
    _fig, _ax = plt.subplots(nrows=n_rows, ncols=n_cols, tight_layout=True)
    _fig.set_figheight(9)
    _fig.set_figwidth(16)
    _ax.grid()
    if title != '':
        _fig.suptitle(title, weight='bold', fontsize=14)
    return _fig, _ax

for n in N:
    h_step = (1 - 0) / n
    x_array = np.linspace(0, 1, n + 1)
    p_array = p(x_array)
    q_array = q(x_array)
    f_array = f(x_array)

    a_array = np.zeros_like(x_array)
    b_array = np.zeros_like(x_array)
    c_array = np.zeros_like(x_array)

    for i in range(1, n):
        a_array[i] = 1 / h_step**2 - p_array[i] / (2 * h_step)
        b_array[i] = -2 / h_step**2 + q_array[i]
        c_array[i] = 1 / h_step**2 + p_array[i] / (2 * h_step)

    # Первый порядок точности (оставляем без изменений)
    # Левое условие (смешанное)
    b_array[0] = alfa_1 - betta_1 / h_step
    c_array[0] = betta_1 / h_step
    f_array[0] = gamma1

    # Правое условие (Дирихле)
    a_array[-1] = 0
    b_array[-1] = 1
    c_array[-1] = 0
    f_array[-1] = gamma2

    steps.append(h_step)
    x_arrays.append(x_array)
    solutions1.append(f2(a_array.copy(), b_array.copy(), c_array.copy(), f_array.copy()))
    errors_1.append(np.max(np.abs(solutions1[-1] - u0(x_array))))

    # Второй порядок точности (новая аппроксимация согласно изображению)
    # Сначала заполняем внутренние точки (как было)
    a2 = np.zeros_like(x_array)
    b2 = np.zeros_like(x_array)
    c2 = np.zeros_like(x_array)
    f2 = np.zeros_like(x_array)

    for i in range(1, n):
        a2[i] = 1 / h_step**2 - p_array[i] / (2 * h_step)
        b2[i] = -2 / h_step**2 + q_array[i]
        c2[i] = 1 / h_step**2 + p_array[i] / (2 * h_step)
        f2[i] = f_array[i]

    A0 = alfa_1 - 3 * betta_1 / (2 * h_step)
    B0 = 4 * betta_1 / (2 * h_step)
    C0 = -1 * betta_1 / (2 * h_step)
    F0 = gamma1

    A0_new = A0 -  a2[1] * C0 / c2[1]
    B0_new = B0 - b2[1] * C0 / c2[1]
    C0_new = 0
    F0_new = F0 - f2[1] * C0 / c2[1]

    a2[0] = 0
    b2[0] = A0_new
    c2[0] = B0_new
    f2[0] = F0_new
    a2[-1] = 0
    b2[-1] = 1
    c2[-1] = 0
    f2[-1] = gamma2

    solutions2.append(f7(a2.copy(), b2.copy(), c2.copy(), f2.copy()))
    errors_2.append(np.max(np.abs(solutions2[-1] - u0(x_array))))

plt.figure(figsize=(14, 6))
plt.loglog(steps, errors_1, 'g-', label='1 порядок', markersize=8)
plt.loglog(steps, errors_2, 'r-', label='2 порядок', markersize=8)
plt.xlabel('Шаг сетки h', fontsize=14)
plt.ylabel('Погрешность', fontsize=14)
plt.title('График зависимости погрешности от шага сетки', fontsize=16)
plt.grid(True, which='both', linestyle='-', alpha=0.7)
plt.legend(fontsize=12)

plt.show()


