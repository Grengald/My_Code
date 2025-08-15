import numpy as np

def is_inside_triangle_ray_casting(px, py, T1, T2, T3):
    
    def intersects(p1, p2, x, y):
        x1, y1 = p1
        x2, y2 = p2
        if y1 > y2:
            x1, y1, x2, y2 = x2, y2, x1, y1
        if y < y1 or y >= y2:
            return False
        
        x_intersect = x1 + (y - y1) * (x2 - x1) / (y2 - y1)
        return x < x_intersect

    edges = [(T1, T2), (T2, T3), (T3, T1)]
    count = sum(intersects(p1, p2, px, py) for p1, p2 in edges)
    return count % 2 == 1

print('Введите ширину прямоугольника:')
A = int(input())
print('Введите длину прямоугольника:')
B = int(input())
T1 = (0, 0)
T2 = (1, 1)
T3 = (1, 0)
print('Введите число точек:')
N = int(input())
S0 = A * B
Ss = []
for i in range(10):
    L = 0

    for n in range(N):
        x = np.random.uniform(0, A)
        y = np.random.uniform(0, B)
        if is_inside_triangle_ray_casting(x, y, T1, T2, T3):
            L += 1

    S = S0 * L / N
    Ss.append(S)
    print(S)

print("Мат. ожидание:", np.mean(Ss))
print("Дисперсия:", np.var(Ss))
