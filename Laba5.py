import numpy as np
import matplotlib.pyplot as plt

def f(y, x):
    u, v = y
    du = v
    dv = -np.cos(x)*v - np.sin(x)*u + x*np.sin(x)
    return np.array([du, dv])

def solution(x):
    return x + np.cos(x)

def el_method(f, y0, x_start, x_end, N):
    step = (x_end - x_start) / N
    x_arr = np.linspace(x_start, x_end, N+1)
    y_arr0 = np.zeros((N+1, len(y0)))
    y_arr0[0,:] = y0
    y_arr=np.zeros((N+1,1))
    y_arr[0,0]=1
    for i in range(N):
        y_arr0[i+1] = y_arr0[i] + step * f(y_arr0[i], x_arr[i])
        y_arr[i+1]=y_arr0[i+1,0]
    return (x_arr, y_arr)

def rk4_method(f, y0, x_start, x_end, N):
    step = (x_end - x_start) / N
    x_arr = np.linspace(x_start, x_end, N+1)
    y_arr0 = np.zeros((N+1, len(y0)))
    y_arr0[0,:] = y0
    y_arr=np.zeros((N+1,1))
    y_arr[0,0]=1
    for i in range(N):
        k1 = f(y_arr0[i], x_arr[i])
        k2 = f(y_arr0[i] + step/2*k1, x_arr[i] + step/2)
        k3 = f(y_arr0[i] + step/2*k2, x_arr[i] + step/2)
        k4 = f(y_arr0[i] + step*k3, x_arr[i] + step)
        y_arr0[i+1] = y_arr0[i] + step/6 * (k1 + 2*k2 + 2*k3 + k4)
        y_arr[i+1]=y_arr0[i+1,0]
    return (x_arr, y_arr)

def ad3_method(f, y0, x_start, x_end, N):
    step = (x_end - x_start) / N
    x_arr = np.linspace(x_start, x_end, N+1)
    y_arr0 = np.zeros((N+1, len(y0)))
    y_arr0[0,:] = y0
    y_arr=np.zeros((N+1,1))
    y_arr[0,0]=1
    for i in range(2):
        k1 = f(y_arr0[i], x_arr[i])
        k2 = f(y_arr0[i] + step/2*k1, x_arr[i] + step/2)
        k3 = f(y_arr0[i] + step/2*k2, x_arr[i] + step/2)
        k4 = f(y_arr0[i] + step*k3, x_arr[i] + step)
        y_arr0[i+1] = y_arr0[i] + step/6 * (k1 + 2*k2 + 2*k3 + k4)
        y_arr[i+1]=y_arr0[i+1,0]
    for i in range(2, N):
        y_arr0[i+1] = y_arr0[i] + step/12 * (23*f(y_arr0[i], x_arr[i]) - 16*f(y_arr0[i-1], x_arr[i-1]) + 5*f(y_arr0[i-2], x_arr[i-2]))
        y_arr[i+1]=y_arr0[i+1,0]
    return (x_arr, y_arr)

def runge_rule(y_h, y_h2, p):
    return np.abs((y_h2 - y_h) / (2**p - 1))

def rk_general(f, y0, x_start, x_end, N, order=4):
    step = (x_end - x_start) / N
    x_arr = np.linspace(x_start, x_end, N + 1)
    y_arr0 = np.zeros((N + 1, len(y0)))
    y_arr0[0, :] = y0
    y_arr = np.zeros((N + 1, 1))
    y_arr[0] = y0[0]
    
    if order == 1:
        a = []
        b = [1]
        c = [0]
    elif order == 2:
        a = [[1]]
        b = [0.5, 0.5]
        c = [0, 1]
    elif order == 3:
        a = [[0.5], [-1, 2]]
        b = [1/6, 4/6, 1/6]
        c = [0, 0.5, 1]
    elif order == 4:
        a = [[0.5], [0, 0.5], [0, 0, 1]]
        b = [1/6, 1/3, 1/3, 1/6]
        c = [0, 0.5, 0.5, 1]
    else:
        raise ValueError("Error")
    
    for i in range(N):
        x = x_arr[i]
        y = y_arr0[i]
        k = []
        k.append(f(y, x))
        for j in range(1, order):
            y_j = y.copy()
            for l in range(j):
                y_j += step * a[j-1][l] * k[l]
            k.append(f(y_j, x + c[j] * step))        
        y_next = y.copy()
        for j in range(order):
            y_next += step * b[j] * k[j]
        
        y_arr0[i + 1] = y_next
        y_arr[i + 1] = y_next[0]
    
    return (x_arr, y_arr)


    
x_start = 0
x_end = 1
N = 20
y0 = np.array([1.0, 1.0])

el_solution = el_method(f, y0, x_start, x_end, N)
rk4_solution = rk4_method(f, y0, x_start, x_end, N)
ad3_solution = ad3_method(f, y0, x_start, x_end, N)
y_exact = solution(el_solution[0])

plt.figure(figsize=(12, 6))
plt.plot(el_solution[0] , el_solution[1], 'b-', label='Метод Эйлера')
plt.plot(rk4_solution[0] , rk4_solution[1], 'r--', label='Метод Рунге-Кутта 4')
plt.plot(ad3_solution[0] , ad3_solution[1], 'g-.', label='Метод Адамса 3')
plt.plot(el_solution[0], y_exact, 'k:', label='Точное решение')
plt.title('Сравнение численных методов (h=0.05)')
plt.xlabel('x')
plt.ylabel('u')
plt.legend()
plt.grid()
plt.show()

steps=[0.01]
for i in range(1,10):
  steps.append(0.05/(2**i))
max_errors_el = []
max_errors_rk4 = []
max_errors_ad3 = []
max_errors_rk2 = []
max_errors_rk3 = []
for step in steps:
    N = int((x_end - x_start) / step)

    x, y = el_method(f, y0, x_start, x_end, N)
    errors_el = []
    for i in range(len(x)):
      error = np.abs(y[i] - solution(x[i]))
      errors_el.append(error)
    max_errors_el.append(max(errors_el))

    x, y = rk4_method(f, y0, x_start, x_end, N)
    errors_rk4 = []
    for i in range(len(x)):
      error = np.abs(y[i] - solution(x[i]))
      errors_rk4.append(error)
    max_errors_rk4.append(max(errors_rk4))

    x, y = ad3_method(f, y0, x_start, x_end, N)
    errors_ad3 = []
    for i in range(len(x)):
      error = np.abs(y[i] - solution(x[i]))
      errors_ad3.append(error)
    max_errors_ad3.append(max(errors_ad3))
    
    x, y = rk_general(f, y0, x_start, x_end, N, order=2)
    errors_rk2 = []
    for i in range(len(x)):
      error = np.abs(y[i] - solution(x[i]))
      errors_rk2.append(error)
    max_errors_rk2.append(max(errors_rk2))
    
    x, y = rk_general(f, y0, x_start, x_end, N, order=3)
    errors_rk3 = []
    for i in range(len(x)):
      error = np.abs(y[i] - solution(x[i]))
      errors_rk3.append(error)
    max_errors_rk3.append(max(errors_rk3))

    
plt.figure(figsize=(12, 6))
plt.loglog(steps, max_errors_el, 'bo-', label='Метод Эйлера (RK1)')
plt.loglog(steps, max_errors_rk2, 'go-', label='Метод Рунге-Кутта 2')
plt.loglog(steps, max_errors_rk3, 'co-', label='Метод Рунге-Кутта 3')
plt.loglog(steps, max_errors_rk4, 'ro-', label='Метод Рунге-Кутта 4')
plt.loglog(steps, max_errors_ad3, 'mo-', label='Метод Адамса 3')
plt.title('Зависимость максимальной ошибки от шага интегрирования')
plt.xlabel('Шаг')
plt.ylabel('Максимальная ошибка')
plt.legend()
plt.grid()
plt.show()


h = 0.1
N_h = int((x_end - x_start) / h)
N_h2 = int((x_end - x_start) / (h/2))
error_el_runge=[]
error_el_real=[]
error_rk4_runge=[]
error_rk4_real=[]
error_ad3_runge=[]
error_ad3_real=[]
x_h, y_h = el_method(f, y0, x_start, x_end, N_h)
x_h2, y_h2 = el_method(f, y0, x_start, x_end, N_h2)
for i in range(len(x_h)):
  error_el_runge.append(runge_rule(y_h[i], y_h2[2*i], 1))
  error_el_real.append(np.abs(y_h2[2*i] - solution(x_h2[2*i])))
print("Метод Эйлера:")
print(f"Оценка погрешности по Рунге: {np.max(error_el_runge)}")
print(f"Фактическая погрешность: {np.max(error_el_real)}\n")

x_h, y_h = rk4_method(f, y0, x_start, x_end, N_h)
x_h2, y_h2 = rk4_method(f, y0, x_start, x_end, N_h2)
for i in range(len(x_h)):
  error_rk4_runge.append(runge_rule(y_h[i], y_h2[2*i], 4))
  error_rk4_real.append(np.abs(y_h2[2*i] - solution(x_h2[2*i])))
print("Метод Рунге-Кутта 4:")
print(f"Оценка погрешности по Рунге: {np.max(error_rk4_runge)}")
print(f"Фактическая погрешность: {np.max(error_rk4_real)}\n")

x_h, y_h = ad3_method(f, y0, x_start, x_end, N_h)
x_h2, y_h2 = ad3_method(f, y0, x_start, x_end, N_h2)
for i in range(len(x_h)):
  error_ad3_runge.append(runge_rule(y_h[i], y_h2[2*i], 3))
  error_ad3_real.append(np.abs(y_h2[2*i] - solution(x_h2[2*i])))
print("Метод Адамса 3:")
print(f"Оценка погрешности по Рунге: {np.max(error_ad3_runge)}")
print(f"Фактическая погрешность: {np.max(error_ad3_real)}")

