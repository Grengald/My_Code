import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

N = 10000  
x_vals = []
y_vals = []
z_vals = []
for i in range(N):
    while True:
        l = np.random.uniform(-1, 1)
        m = np.random.uniform(-1, 1)
        n = np.random.uniform(-1, 1)
        length = np.sqrt(l**2 + m**2 + n**2)

        if length > 1:
            continue
        l = l / length
        m = m / length
        n = n / length
        x_vals.append(l)
        y_vals.append(m)
        z_vals.append(n)
        break

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x_vals, y_vals, z_vals, s=5)
ax.set_box_aspect([1, 1, 1])  
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
ax.set_title("Изотропный источник частиц")

plt.show()
