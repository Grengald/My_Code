import numpy as np
import matplotlib.pyplot as plt

def napravlenie():
    tetta = np.arccos(2 * np.random.random() - 1)
    fi = 2 * np.pi * np.random.random()
    x_dir = np.sin(tetta) * np.cos(fi)
    y_dir = np.sin(tetta) * np.sin(fi)
    z_dir = np.cos(tetta)

    return np.array([x_dir, y_dir, z_dir])

def peresechenie(particle_origin, ray_direction, cube_center, cube_side):
    half_side = cube_side / 2
    cube_min = cube_center - half_side
    cube_max = cube_center + half_side
    t_enter, t_exit = -np.inf, np.inf
    for axis in range(3):
        if np.abs(ray_direction[axis]) <= 1e-7:
            if not (cube_min[axis] <= particle_origin[axis] <= cube_max[axis]):
                return False
        else:
            t1 = (cube_min[axis] - particle_origin[axis]) / ray_direction[axis]
            t2 = (cube_max[axis] - particle_origin[axis]) / ray_direction[axis]
            t_enter = max(t_enter, min(t1, t2))
            t_exit = min(t_exit, max(t1, t2))

    return t_exit > max(t_enter, 0)


def tests():
    cube_side = 1.0
    cube_center = np.array([0.0, 0.0, 0.0])
    test_sources = {
        "центр грани": np.array([0.0, 0.0, 0.5]),
        "середина ребра": np.array([0.5, 0.5, 0.0]),
        "вершина": np.array([0.5, 0.5, 0.5])
    }

    for source_name, source_pos in test_sources.items():
        hit_count = 0
        total_tests = 100000

        for _ in range(total_tests):
            direction = napravlenie()
            if peresechenie(source_pos, direction, cube_center, cube_side):
                hit_count += 1
        actual_percent = 100 * hit_count / total_tests
        print(f"Тест '{source_name}':")
        print(f"  Процентр попаданий: {actual_percent:.2f}%")
        print("----------------------------------")
def popadanie(source_start_pos=np.array([0.0, 0.0, 0.0])):
    cube_side = 1.0
    cube_center = np.array([0.0, 0.0, 0.0])
    particles_per_step = 10000
    distance_steps = 100
    max_distance = 5
    distance_range = np.linspace(0, max_distance, distance_steps)
    hit_percentages = []
    for current_distance in distance_range:
        hits = 0
        current_source_pos = source_start_pos + np.array([current_distance, 0.0, 0.0])
        for i in range(particles_per_step):
            particle_dir = napravlenie()
            if peresechenie(current_source_pos, particle_dir, cube_center, cube_side):
                hits += 1
        hit_percent = 100 * hits / particles_per_step
        hit_percentages.append(hit_percent)
        print(f"Расстояние: {current_distance:.2f}, Попаданий: {hit_percent:.2f}%")

    plt.figure(figsize=(14, 10))
    plt.step(distance_range, hit_percentages, color="lime", linewidth=4)
    plt.xlabel('Расстояние от центра куба', fontsize=14)
    plt.ylabel('Процент попаданий, [%]', fontsize=14)
    plt.title('Зависимость попаданий частиц в куб от расстояния до источника', fontsize=16)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    print("Тесты для разных положений источника")
    tests()

    print("Источник в начале координат:")
    popadanie(np.array([0.5, 0.5, 0.0]))
