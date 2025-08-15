import numpy as np
import multiprocessing as mp
import time
start_time = time.time()
a = 10
b = 10
R = 1
x0 = 5
y0 = 5
N = 10000000

def estimate_area(seed):
    np.random.seed(seed)  
    x = np.random.uniform(0, a, N)
    y = np.random.uniform(0, b, N)
    L = np.sum(np.sqrt((x - x0)**2 + (y - y0)**2) < R)
    S0 = a * b
    S = S0 * L / N
    return S

if __name__ == '__main__':
    with mp.Pool(processes=1) as pool:
        seeds = list(range(10))
        Ss = pool.map(estimate_area, seeds)

    for i, s in enumerate(Ss, 1):
        print(f"Попытка {i}: {s}")

    print("Мат. ожидание:", np.mean(Ss))
    print("Дисперсия:", np.var(Ss))
    end_time = time.time()
    num_cores = mp.cpu_count()
    print(f"Число доступных ядер: {num_cores}")
    elapsed_time = end_time - start_time
    print(f"The task took {elapsed_time:.2f} seconds to complete.")
