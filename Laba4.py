import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return x**3-2*x+2

def df(x):
    return 3*x**2-2

def bisection(f, x0, x1, delta, max_iters=1000):
    if f(x0) * f(x1) >= 0:
        raise ValueError("Функция должна менять знак на концах интервала")

    iterations = 0
    while abs(x1 - x0) > delta:
        x2 = (x0 + x1) / 2
        if f(x2) == 0 or abs(x1 - x0) < delta:
            break
        elif f(x0) * f(x2) < 0:
            x1 = x2
        else:
            x0 = x2
        
        iterations += 1
        if iterations > max_iters:
            raise RuntimeError("Превышено максимальное количество итераций")

    return x2, iterations, abs(f(x2))

def newton(function, diff_function, x0, max_iterations=100, min=1/(10**6)):
    x = x0
    iterations = []
    for i in range(max_iterations):
        fx = function(x)
        dfx = diff_function(x)
        if abs(dfx) < min:
            print("Производная близка к нулю. Метод не может продолжаться.")
            return None, iterations
        x_new = x - fx / dfx
        iterations.append((x, x_new))
        if abs(x_new - x) < min:
            print(f"Найден корень: {x_new} за {i+1} итераций.")
            return x_new, iterations
        x = x_new
    print(f"Метод не сошелся за {max_iterations} итераций. Последнее приближение: {x}")
    return None, iterations

a, b = 0, 10
x0_newton = 0
delta_values = [10**(-3), 10**(-6), 10**(-9)]
c=4.79277761106522851037625755855006058199844786980422
for delta in delta_values:
    try:
        point_d, iterations_d, error_d = bisection(f, a, b, delta)
        print(f"Метод дихотомии (Δ={delta}): x={point_d:.10f}, невязка={error_d:.2e}, итераций={iterations_d}")
    except Exception as e:
        print(f"Метод дихотомии (Δ={delta}): Ошибка - {e}")

    try:
        point_n, iterations_n, error_n = newton(f, df, x0_newton, c, delta)
        print(f"Метод Ньютона (Δ={delta}): x={point_n:.10f}, невязка={error_n:.2e}, итераций={iterations_n}")
    except Exception as e:
        print(f"Метод Ньютона (Δ={delta}): Ошибка - {e}")

    print()


x_vals = np.linspace(0, 10, 1000)
y_vals = f(x_vals)

plt.plot(x_vals, y_vals, label="4*sin(x/2) + cos(x)*th(x) - x + 2")
plt.axhline(0, color='black', linewidth=0.8)
plt.xlabel("x")
plt.ylabel("f(x)")
plt.legend()
plt.title("График функции")
plt.grid()
plt.show()
