import numpy as np
import time
start_time = time.time()
a=10
b=10
R=1
x0=5
y0=5
N=100000000
Ss=[]

for i in range(10):
    L=0
    x=np.random.uniform(0, a, N)
    y=np.random.uniform(0, b, N)
    
    L=np.sum(np.sqrt((x - x0)** 2+(y - y0)**2)< R)
    
    S0=a * b
    S=S0 * L / N
    Ss.append(S)
    print(S)

print("Мат. ожидание:", np.mean(Ss))
print("Дисперсия:", np.var(Ss))
end_time = time.time()
elapsed_time = end_time - start_time
print(f"The task took {elapsed_time:.2f} seconds to complete.")
