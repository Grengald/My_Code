import random
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def generate_direction_spherical():
    phi = 2 * math.pi * random.random()
    cos_theta = 2 * random.random() - 1
    theta = math.acos(cos_theta)
    x = math.sin(theta) * math.cos(phi)
    y = math.sin(theta) * math.sin(phi)
    z = math.cos(theta)
    return theta, phi, x, y, z

N = 1000

results = [generate_direction_spherical() for i in range(N)]

theta_vals, phi_vals, x_vals, y_vals, z_vals = zip(*results)

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')

ax.scatter(x_vals, y_vals, z_vals, color='blue', marker='o', s=20, alpha=1)
ax.set_box_aspect([1, 1, 1])
ax.set_title("Распределение направлений частиц (сферические координаты)")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")

plt.show()
