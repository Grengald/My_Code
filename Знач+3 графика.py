from scipy.io import readsav
import ctypes
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates

file_path = 'C:\\Users\\Grengald\\Desktop\\Котинг\\ИКИ РАН\\Sav bd\GRB160625B.sav'
data = readsav(file_path)

data_structure = {key: type(value) for key, value in data.items()}
data_structure

JD_Long_1 = data['jd_long_1']
CR_Long_1 = data['cr_long_1']
Short_Profile = data['short_profile']
CR_Long_1_Err = data['cr_long_1_err']
Short_Profile_Err = data['short_profile_err']
JD_Short = data['jd_short']

N_GRB = 0
GRB_Names = ['GRB160625B']
SJDs_Max = [5565.44649305]
T90s = [22.28]

GRB_Name = GRB_Names[N_GRB]
SJD_Max = SJDs_Max[N_GRB]
T90 = T90s[N_GRB]

DT = 120.0 / 86400.0
sjd_shift = 2452000.0
JD_Max = SJD_Max + sjd_shift


JD1 = sjd_shift + SJD_Max - 1.0 * DT
JD2 = sjd_shift + SJD_Max + 2.0 * DT

T_Start = datetime.strptime("2009-11-09 21:49:03", "%Y-%m-%d %H:%M:%S")
T_Stop = T_Start + timedelta(seconds=2 * 86400)


JD_BG_L_1 = JD_Max - (10.0 + T90 + T90 + 60.0) / 86400.0
JD_BG_L_2 = JD_Max - (10.0 + T90) / 86400.0
JD_BG_R_1 = JD_Max + (10.0 + 2.0 * T90) / 86400.0
JD_BG_R_2 = JD_Max + (10.0 + 2.0 * T90 + T90 + 60.0) / 86400.0

NN_BG_1_L = np.where((JD_Long_1 >= JD_BG_L_1) & (JD_Long_1 <= JD_BG_L_2))[0]
NN_BG_1_R = np.where((JD_Long_1 >= JD_BG_R_1) & (JD_Long_1 <= JD_BG_R_2))[0]
NN_BG_1 = np.concatenate((NN_BG_1_L, NN_BG_1_R))
Cut_JD_1 = (JD_Long_1 - JD_Long_1[0]) * 86400.0
Coeff_1 = np.polyfit(Cut_JD_1[NN_BG_1], CR_Long_1[NN_BG_1], 2)
BG_Long = np.polyval(Coeff_1, Cut_JD_1)

NN_BG_2_L = np.where((JD_Short >= JD_BG_L_1) & (JD_Short <= JD_BG_L_2))[0]
NN_BG_2_R = np.where((JD_Short >= JD_BG_R_1) & (JD_Short <= JD_BG_R_2))[0]
NN_BG_2 = np.concatenate((NN_BG_2_L, NN_BG_2_R))
Cut_JD_2 = (JD_Short - JD_Short[0]) * 86400.0
Coeff_2 = np.polyfit(Cut_JD_2[NN_BG_2], Short_Profile[NN_BG_2], 2)
BG_Short = np.polyval(Coeff_2, Cut_JD_2)

Signal_Short = Short_Profile - BG_Short
Signal_Long_1 = CR_Long_1 - BG_Long


HR = Signal_Short / Signal_Long_1
HR_Err = np.sqrt((Short_Profile_Err / Signal_Long_1)**2 + (Signal_Short * CR_Long_1_Err / CR_Long_1**2)**2)

CR_Long_1_Max = np.max(CR_Long_1)
N_CR_Long_1_Max = np.argmax(CR_Long_1)
print('Max CR long', CR_Long_1_Max, Short_Profile[N_CR_Long_1_Max])
print('Max HR long', HR[N_CR_Long_1_Max], HR_Err[N_CR_Long_1_Max])

Short_Profile_Max = np.max(Short_Profile)
N_Short_Profile_Max = np.argmax(Short_Profile)
print('Maximum signals(Short): ', Signal_Short[N_Short_Profile_Max], Signal_Long_1[N_Short_Profile_Max])
print('HR max bin(short):', HR[N_Short_Profile_Max], HR_Err[N_Short_Profile_Max])
########################################################################################################################################
Left_Short_Profile = Signal_Short[:N_Short_Profile_Max]
Right_Short_Profile = Signal_Short[(N_Short_Profile_Max+1):]
Left_Signal_Long_1 = Signal_Long_1[:N_Short_Profile_Max] 
Right_Signal_Long_1 = Signal_Long_1[(N_Short_Profile_Max+1):]

Left_Short_Profile_Err = Short_Profile_Err[:N_Short_Profile_Max] 
Right_Short_Profile_Err = Short_Profile_Err[(N_Short_Profile_Max+1):]
Left_CR_Long_1_Err = CR_Long_1_Err[:N_Short_Profile_Max] 
Right_CR_Long_1_Err = CR_Long_1_Err[(N_Short_Profile_Max+1):]

Left_CR_Long_1 = CR_Long_1[:N_Short_Profile_Max] 
Right_CR_Long_1 = CR_Long_1[(N_Short_Profile_Max+1):]

Left_JD_Long_1 = JD_Long_1[:N_Short_Profile_Max]
Right_JD_Long_1 = JD_Long_1[(N_Short_Profile_Max+1):]

Signal_Long_1_S = Signal_Long_1[N_Short_Profile_Max]
Signal_Short_S = Signal_Short[N_Short_Profile_Max]
Short_Profile_Err_S=Short_Profile_Err[N_Short_Profile_Max]
CR_Long_1_Err_S=CR_Long_1_Err[N_Short_Profile_Max]
CR_Long_1_S=CR_Long_1[N_Short_Profile_Max]

SHR = []
S_y = []
S_x = []
S_err_long = []
S_w = []
S_z = []
S_err_short = []

SHR_0 = HR[N_Short_Profile_Max] / HR_Err[N_Short_Profile_Max]
SHR.append(SHR_0)
S_y.append(Signal_Long_1[N_Short_Profile_Max])
S_x.append(JD_Long_1[N_Short_Profile_Max])
S_err_long.append(CR_Long_1_Err[N_Short_Profile_Max])
S_w.append(Signal_Short[N_Short_Profile_Max])
S_z.append(JD_Short[N_Short_Profile_Max])
S_err_short.append(Short_Profile_Err[N_Short_Profile_Max])

n = 1
m = 1
CR_Long_1_Err_S_1=CR_Long_1_Err[N_Short_Profile_Max] ** 2
Short_Profile_Err_S_1=Short_Profile_Err[N_Short_Profile_Max] ** 2
while True:
    if n < len(Left_Short_Profile) and (m >= len(Right_Short_Profile) or Left_Short_Profile[-n] > Right_Short_Profile[m-1]):
        Signal_Long_1_S += Left_Signal_Long_1[-n]
        Signal_Short_S += Left_Short_Profile[-n]
        CR_Long_1_S+= Left_CR_Long_1[-n]
        CR_Long_1_Err_S_1+=Left_CR_Long_1_Err[-n] ** 2
        Short_Profile_Err_S_1+=Left_Short_Profile_Err[-n] ** 2
        CR_Long_1_Err_S=np.sqrt(CR_Long_1_Err_S_1)
        Short_Profile_Err_S=np.sqrt(Short_Profile_Err_S_1)
        HR_Err_S = np.sqrt((Short_Profile_Err_S / Signal_Long_1_S)**2 + (Signal_Short_S * CR_Long_1_Err_S /  CR_Long_1_S**2)**2)
        JD_Long_1_S=Left_JD_Long_1[-n]
        
        S_y.append(Left_Signal_Long_1[-n])
        S_x.append(Left_JD_Long_1[-n])
        S_err_long.append(Left_CR_Long_1_Err[-n])
        S_w.append(Left_Short_Profile[-n])
        S_z.append(Left_JD_Long_1[-n])
        S_err_short.append(Left_Short_Profile_Err[-n])
        
        n += 1
    elif m < len(Right_Short_Profile):  # Добавляем проверку
        Signal_Long_1_S += Right_Signal_Long_1[m-1]
        Signal_Short_S += Right_Short_Profile[m-1]
        CR_Long_1_S+= Right_CR_Long_1[m-1]
        CR_Long_1_Err_S_1+=Right_CR_Long_1_Err[m-1] ** 2
        Short_Profile_Err_S_1+=Right_Short_Profile_Err[m-1] ** 2
        CR_Long_1_Err_S=np.sqrt(CR_Long_1_Err_S_1)
        Short_Profile_Err_S=np.sqrt(Short_Profile_Err_S_1)
        HR_Err_S = np.sqrt((Short_Profile_Err_S / Signal_Long_1_S)**2 + (Signal_Short_S * CR_Long_1_Err_S /  CR_Long_1_S**2)**2)
        JD_Long_1_S=Right_JD_Long_1[m-1]
        
        S_y.append(Right_Signal_Long_1[m-1])
        S_x.append(Right_JD_Long_1[m-1])
        S_err_long.append(Right_CR_Long_1_Err[m-1])
        S_w.append(Right_Short_Profile[m-1])
        S_z.append(Right_JD_Long_1[m-1])
        S_err_short.append(Right_Short_Profile_Err[m-1])
        
        m += 1
    else:
        break
    
    HR_S = Signal_Short_S / Signal_Long_1_S
    SHR_n = HR_S / HR_Err_S
    SHR.append(SHR_n)
    
    if Signal_Short_S < 3 * Short_Profile_Err_S:
        break

#print('Max SHR:', np.argmax(SHR))

##################################################################################################################################################
S_y_1=S_y[:np.argmax(SHR)+1]
S_x_1=S_x[:np.argmax(SHR)+1]
S_err_long_1=S_err_long[:np.argmax(SHR)+1]

S_w_1=S_w[:np.argmax(SHR)+1]
S_z_1=S_z[:np.argmax(SHR)+1]
S_err_short_1=S_err_short[:np.argmax(SHR)+1]
##################################################################################################################################################
y = np.array(S_y_1)
x = np.array(S_x_1)
err_long = np.array(S_err_long_1)

sort_order= np.argsort(x)
x_sorted = x[sort_order]
y_sorted = y[sort_order]
err_long_sorted = err_long[sort_order]

integral_1 = np.trapz(y_sorted, x_sorted)
dx = np.diff(x_sorted)
err_long_shortened = err_long_sorted[:-1]
error_integral_1 = np.sqrt(np.sum((err_long_shortened * dx)**2))
print("Integral Long Profile:", integral_1)
print('Error Integral Long Profile:',error_integral_1)
#print(x_sorted )

w = np.array(S_w_1)
#print(S_w_1)
z = np.array(S_z_1)
#print(S_z_1)
err_short = np.array(S_err_short_1)

w_sorted=w[sort_order]
z_sorted=z[sort_order]
err_short_sorted=err_short[sort_order]

integral_2 = np.trapz(w_sorted, z_sorted)
dz = np.diff(z_sorted)
err_short_shortened = err_short_sorted[:-1]
error_integral_2 = np.sqrt(np.sum((err_short_shortened * dz)**2))

print("Integral Short Profile:", integral_2)
print('Error Integral Short Profile:',error_integral_2)

HR_I = integral_2 / integral_1
HR_I_error = HR_I * np.sqrt((error_integral_1 / integral_1)**2 + (error_integral_2 / integral_2)**2)
print('HR:', HR_I)
print('Error HR:',HR_I_error)
print('Колличество бинов', np.argmax(SHR)+1)
print('Max Sign:', max(SHR))
##############################################################################################################################################################

import matplotlib.dates as mdates
from datetime import datetime, timedelta

# Функция для преобразования JD в datetime
def jd_to_datetime(jd):
    jd_epoch = 2400000.5  # JD эпоха для преобразования
    return datetime(1858, 11, 17) + timedelta(days=(jd - jd_epoch))

# Преобразование массива юлианских дат в datetime
JD_Long_1_datetime = [jd_to_datetime(jd) for jd in JD_Long_1]
JD_Short_datetime = [jd_to_datetime(jd) for jd in JD_Short]

# Построение трех графиков
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(14, 15), sharex=False)

# Синхронизация оси X для ax1 и ax2
ax2.sharex(ax1)

# Верхний график: длинный и короткий профиль в виде гистограммы с погрешностями
width = 0.0000115746  # ширина столбцов гистограммы

# Длинный профиль
ax1.bar(JD_Long_1_datetime, CR_Long_1, width=width, edgecolor='blue', facecolor='none', label='Long Profile')
ax1.errorbar(JD_Long_1_datetime, CR_Long_1, yerr=CR_Long_1_Err, fmt='.', color='blue')
ax1.plot(JD_Long_1_datetime, BG_Long, color='#969696', label='Background Long', linestyle='--')

# Короткий профиль
ax1.bar(JD_Short_datetime, Short_Profile*5, width=width, edgecolor='green', facecolor='none', label='Short Profile*5')
ax1.errorbar(JD_Short_datetime, Short_Profile*5, yerr=Short_Profile_Err*5, fmt='.', color='green')
ax1.plot(JD_Short_datetime, BG_Short, color='#FFD700', label='Background Short', linestyle='--')

# Бины с максимальной значимостью 
ax1.bar([jd_to_datetime(jd) for jd in S_z_1], np.array(S_w_1), width=width, edgecolor='red', facecolor='red', label='Max Significance')

# Вертикальная линия для максимальной счетной скорости
ax1.axvline(jd_to_datetime(JD_Short[N_Short_Profile_Max]), linestyle='--', color='red', label='Max Count Rate')
ax1.set_ylabel('Count Rate, CNT/S')
ax1.legend()
ax1.grid(True)

# Убрать метки по оси X на первом графике
ax1.set_xticklabels([])

# Средний график: отношение жесткости
NN_HR = np.where(Signal_Short >= Short_Profile_Err * 3)[0]

if len(NN_HR) == 0:
    raise ValueError

ax2.plot([JD_Long_1_datetime[i] for i in NN_HR], HR[NN_HR], '+', markersize=5, label='Harness Ratio')
ax2.errorbar([JD_Long_1_datetime[i] for i in NN_HR], HR[NN_HR], yerr=HR_Err[NN_HR], color='blue', fmt='+')
ax2.axvline(jd_to_datetime(JD_Long_1[N_CR_Long_1_Max]), linestyle='--', color='red')
ax2.axhline(0, linestyle='--', color='black')
ax2.set_ylabel('Harness Ratio')
ax2.set_ylim([-0.1, 0.2])
ax2.legend()
ax2.grid(True)

# Нижний график: значимость HR (SHR)
ax3.plot(range(1, len(SHR) + 1), SHR, marker='.')
ax3.set_xlabel('Количество бинов')
ax3.set_ylabel('Значимость HR (SHR)')
ax3.set_xlim(0, len(SHR) + 1)  # Пределы для количества бинов
ax3.set_ylim(0, 60)  # Пределы для значимости HR
ax3.grid(True)

# Форматирование оси X для отображения времени
formatter = mdates.DateFormatter('%H:%M:%S')  # формат времени ЧЧ:ММ:СС
ax1.xaxis.set_major_formatter(formatter)
ax2.xaxis.set_major_formatter(formatter)

# Добавляем отступы для оси X, чтобы избежать обрезки
ax1.set_xlim(min(JD_Long_1_datetime) - timedelta(seconds=0), max(JD_Long_1_datetime) + timedelta(seconds=0))
ax2.set_xlim(min(JD_Long_1_datetime) - timedelta(seconds=0), max(JD_Long_1_datetime) + timedelta(seconds=0))

# Принудительное обновление границ
ax1.autoscale_view()

# Общий заголовок для всей фигуры
plt.suptitle(GRB_Name)

# Показ графиков
plt.show()
