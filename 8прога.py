import numpy as np
import random
import math
import matplotlib.pyplot as plt

# Параметры детектора и источника
R = 5.0  # радиус цилиндра (см)
D = 10.0  # высота цилиндра (см)
X0, Y0, Z0 = 0, 0, -10  #(см)
z1 = -D / 2
z2 = D / 2
N_photons = 50000
N_channels = 1024
E_min = 0.05
E_max = 1.0
channel_width = (E_max - E_min) / N_channels
E_gamma = 0.662  # энергия Cs-137 (МэВ)
me_c2 = 0.511
Na_Z, I_Z = 11, 53

# Постоянные для расчета макроскопических сечений
NA = 6.022e23  # число Авогадро
density = 3.67  # плотность NaI (г/см^3)
M_NaI = 149.89  # г/моль
N_atoms = NA * density / M_NaI

# Массив для спектра
spectrum = [0] * N_channels

# ----------------- ФУНКЦИИ -----------------

def ray():
    cos_theta = random.uniform(-1, 1)
    sin_theta = math.sqrt(1 - cos_theta**2)
    phi = random.uniform(0, 2 * math.pi)
    return sin_theta * math.cos(phi), sin_theta * math.sin(phi), cos_theta

def sigmaPh(E, Z):
    return (5/4)*(6.651e-25 * 4 * math.sqrt(2) * (Z**5 / 137**4) * ((me_c2 / E) ** 3.5))

def sigmaK(E, Z):
    gamma = E / me_c2
    ln_term = math.log(2 * gamma + 1)
    bracket = (1 - (2 * (gamma + 1)) / gamma**2) * ln_term
    bracket += 0.5 + 4 / gamma - 1 / (2 * (2 * gamma + 1)**2)
    return 6.651e-25 * 3 * Z / (8 * gamma) * bracket

def sigma_macro(E):
    sigma_ph_Na = sigmaPh(E, Na_Z)
    sigma_ph_I = sigmaPh(E, I_Z)
    sigma_K_Na = sigmaK(E, Na_Z)
    sigma_K_I = sigmaK(E, I_Z)
    mu_ph = 5 / 4 * N_atoms * (sigma_ph_Na + sigma_ph_I)
    mu_K = N_atoms * ((Na_Z / 22) * sigma_K_Na + (I_Z / 127) * sigma_K_I)
    return mu_ph + mu_K, mu_ph, mu_K

def crossFlat(x0, y0, z0, u, z_face):
    ux, uy, uz = u
    if uz == 0:
        return None
    t = (z_face - z0) / uz
    if t <= 0:
        return None
    x = x0 + ux * t
    y = y0 + uy * t
    if x**2 + y**2 <= R**2:
        return t
    return None

def crossCil(x0, y0, z0, u):
    ux, uy, uz = u
    A = ux**2 + uy**2
    B = 2 * (x0 * ux + y0 * uy)
    C = x0**2 + y0**2 - R**2
    D = B**2 - 4 * A * C
    if D < 0 or A == 0:
        return []
    sqrtD = math.sqrt(D)
    t1 = (-B + sqrtD) / (2 * A)
    t2 = (-B - sqrtD) / (2 * A)
    valid = []
    for t in [t1, t2]:
        if t > 0:
            z = z0 + uz * t
            if z1 <= z <= z2:
                valid.append(t)
    return valid

def sample_interaction(E):
    mu_total, mu_ph, mu_K = sigma_macro(E)
    L = -math.log(random.random()) / mu_total
    if random.random() < mu_ph / mu_total:
        return "photo", L
    return "compton", L

def compton_scatter(E, u):
    while True:
        cos_theta = random.uniform(-1, 1)
        E_prime = E / (1 + (E / me_c2) * (1 - cos_theta))
        if E_prime > 0:
            break
    delta_E = E - E_prime

    phi = random.uniform(0, 2 * math.pi)
    sin_theta = math.sqrt(1 - cos_theta**2)
    ux, uy, uz = u
    if abs(ux) < 1e-6 and abs(uy) < 1e-6:
        a = (1, 0, 0)
    else:
        a = (-uy, ux, 0)
    a_len = math.sqrt(sum([i**2 for i in a]))
    a = tuple(i / a_len for i in a)
    b = (
        uy * a[2] - uz * a[1],
        uz * a[0] - ux * a[2],
        ux * a[1] - uy * a[0]
    )
    u_new = tuple(
        cos_theta * u[i] +
        sin_theta * math.cos(phi) * a[i] +
        sin_theta * math.sin(phi) * b[i]
        for i in range(3)
    )
    u_len = math.sqrt(sum(i**2 for i in u_new))
    u_new = tuple(i / u_len for i in u_new)
    return E_prime, delta_E, u_new

# ----------------- ОСНОВНОЙ ЦИКЛ -----------------

for _ in range(N_photons):
    u = ray()
    x0, y0, z0 = X0, Y0, Z0
    t_vals = []
    for z_face in [z1, z2]:
        t = crossFlat(x0, y0, z0, u, z_face)
        if t:
            t_vals.append(t)
    t_vals += crossCil(x0, y0, z0, u)
    t_vals = [t for t in t_vals if t > 0]
    if not t_vals:
        continue
    t_vals.sort()
    t_entry = t_vals[0]
    t_exit = t_vals[1] if len(t_vals) > 1 else t_vals[0]
    path_length = t_exit - t_entry
    x_entry = x0 + u[0] * t_entry
    y_entry = y0 +u[1] * t_entry
    z_entry = z0 + u[2] * t_entry
    E = E_gamma

    while True:
        interaction_type, L = sample_interaction(E)
        if L > path_length:
            break
        x_int = x_entry + u[0] * L
        y_int = y_entry + u[1] * L
        z_int = z_entry + u[2] * L
        if interaction_type == "photo":
            if E >= E_min:
                ch = int((E - E_min) / channel_width)
                if 0 <= ch < N_channels:
                    spectrum[ch] += 1
            break
        else:
            E, dE, u = compton_scatter(E, u)
            if dE > 0:
                E_bin = dE
                if E_bin >= E_min:
                    ch = int((E_bin - E_min) / channel_width)
                    if 0 <= ch < N_channels:
                        spectrum[ch] += 1
            x_entry, y_entry, z_entry = x_int, y_int, z_int
            t_vals = []
            for z_face in [z1, z2]:
                t = crossFlat(x_entry, y_entry, z_entry, u, z_face)
                if t:
                    t_vals.append(t)
            t_vals += crossCil(x_entry, y_entry, z_entry, u)
            t_vals = [t for t in t_vals if t > 0]
            if not t_vals:
                break
            t_vals.sort()
            path_length = t_vals[0]

energy_range = np.linspace(0.05, 1.0, 100)  # Энергия в МэВ

# Рассчитаем сечения для Na и I
sigma_ph_Na = [sigmaPh(E, Na_Z) for E in energy_range]  # Фотоэффект (Na)
sigma_ph_I = [sigmaPh(E, I_Z) for E in energy_range]    # Фотоэффект (I)
sigma_K_Na = [sigmaK(E, Na_Z) for E in energy_range]    # Комптон (Na)
sigma_K_I = [sigmaK(E, I_Z) for E in energy_range]      # Комптон (I)

# Настройка стиля графиков
plt.style.use('seaborn-v0_8')
plt.rcParams['figure.figsize'] = (12, 8)

# 1. График сечений фотоэффекта для Na и I
plt.figure()
plt.plot(energy_range * 1000, sigma_ph_Na, label='Na (Z=11)', linewidth=2)
plt.plot(energy_range * 1000, sigma_ph_I, label='I (Z=53)', linewidth=2)
plt.xlabel('Энергия (keV)', fontsize=12)
plt.ylabel('Сечение (см²/атом)', fontsize=12)
plt.title('Сечение фотоэффекта для Na и I', fontsize=14)
plt.yscale('log')
plt.grid(True, which="both", linestyle='--')
plt.legend(fontsize=12)
plt.tight_layout()

# 2. График сечений комптоновского рассеяния для Na и I
plt.figure()
plt.plot(energy_range * 1000, sigma_K_Na, label='Na (Z=11)', linewidth=2)
plt.plot(energy_range * 1000, sigma_K_I, label='I (Z=53)', linewidth=2)
plt.xlabel('Энергия (keV)', fontsize=12)
plt.ylabel('Сечение (см²/атом)', fontsize=12)
plt.title('Сечение комптоновского рассеяния для Na и I', fontsize=14)
plt.yscale('log')
plt.grid(True, which="both", linestyle='--')
plt.legend(fontsize=12)
plt.tight_layout()

# 3. График сечений для натрия (оба эффекта)
plt.figure()
plt.plot(energy_range * 1000, sigma_ph_Na, label='Фотоэффект (Na)', linewidth=2, linestyle='-')
plt.plot(energy_range * 1000, sigma_K_Na, label='Комптон (Na)', linewidth=2, linestyle='--')
plt.xlabel('Энергия (keV)', fontsize=12)
plt.ylabel('Сечение (см²/атом)', fontsize=12)
plt.title('Сравнение сечений для натрия (Z=11)', fontsize=14)
plt.yscale('log')
plt.grid(True, which="both", linestyle='--')
plt.legend(fontsize=12)
plt.tight_layout()

# 4. График сечений для йода (оба эффекта)
plt.figure()
plt.plot(energy_range * 1000, sigma_ph_I, label='Фотоэффект (I)', linewidth=2, linestyle='-')
plt.plot(energy_range * 1000, sigma_K_I, label='Комптон (I)', linewidth=2, linestyle='--')
plt.xlabel('Энергия (keV)', fontsize=12)
plt.ylabel('Сечение (см²/атом)', fontsize=12)
plt.title('Сравнение сечений для йода (Z=53)', fontsize=14)
plt.yscale('log')
plt.grid(True, which="both", linestyle='--')
plt.legend(fontsize=12)
plt.tight_layout()

# Показать все графики
plt.show()

# Оригинальный график спектра (оставляем без изменений)
energies = [E_min + (i + 0.5) * channel_width for i in range(N_channels)]
plt.plot([e * 1000 for e in energies], spectrum)
plt.xlabel('Энергия (keV)')
plt.ylabel('Число событий')
plt.title('Смоделированный спектр Cs-137 в детекторе NaI')
plt.grid(True)
plt.show()
