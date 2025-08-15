from scipy.io import readsav
import ctypes
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates


file_path = 'C:\\Users\\Grengald\\Desktop\\Котинг\\ИКИ РАН\\Sav bd\\GRB130427A.sav'
data = readsav(file_path)

long_time = data['jd_long_1']
long_time = np.array(long_time)
time = long_time - 2452000

N_GRB = 0
GRB_Names = ['GRB130427A']
SJDs_Max = [4409.83236026]
T90s = [6.46]

GRB_Name = GRB_Names[N_GRB]
SJD_Max = SJDs_Max[N_GRB]
T90 = T90s[N_GRB]

intensity_1 = data['cr_long_1']
intensity_2 = data['cr_long_1']
intensity_1 = list(intensity_1)
intensity_2 = list(intensity_2)

JD_Long_1 = data['jd_long_1']
CR_Long_1 = data['cr_long_1']
Short_Profile = data['short_profile']
CR_Long_1_Err = data['cr_long_1_err']
Short_Profile_Err = data['short_profile_err']
JD_Short = data['jd_short']

CR_long_1 = data['cr_long_1']
CR_long_1 = np.array(CR_long_1)

i=50
n=50
pre_line_1 = sum(intensity_1[:50])/50
while True:
    if i<len(intensity_1):
        if intensity_1[i] < 1.5*pre_line_1:
            pre_line_1=sum(intensity_1[:i])/i
            i+=1
        elif intensity_1[i] > 1.5*pre_line_1:
            del intensity_1[i]
            i=i
    else:
        break
print('Фон, начиная с лева:', pre_line_1)

intensity_2.reverse()
pre_line_2 = sum(intensity_2[:50])/50
while True:
    if n<len(intensity_2):
        if intensity_2[n] < 1.5*pre_line_2:
            pre_line_2=sum(intensity_2[:n])/n
            n+=1
        elif intensity_2[n] > 1.5*pre_line_2:
            del intensity_2[n]
            n=n
    else:
        break
print("Фон, начиная справа:", pre_line_2)

line_BG=(pre_line_1+pre_line_2)/2
print('Среднее значение фона:', line_BG)
clean_CR=CR_long_1-line_BG
total_clean_CR=sum(clean_CR)

sorted_clean_CR = np.sort(clean_CR)[::-1]
C=0
m=0
while True:
    if C<0.5*total_clean_CR:
        C=C+sorted_clean_CR[m]
        m=m+1
    else:
        print(m)
        print(C)
        break

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

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10), sharex=True)

width = 0.0000116

# Определяем эмиссионное время по самым ярким бинам
largest_indices = np.argpartition(CR_Long_1, -m)[-m:]
largest_indices = largest_indices[np.argsort(CR_Long_1[largest_indices])][::-1]
JD_Long_1_largest = JD_Long_1[largest_indices]
CR_Long_1_largest = CR_Long_1[largest_indices]

# Находим соответствующие точки в коротком профиле
emission_mask_long = (JD_Long_1 >= JD_Long_1_largest.min()) & (JD_Long_1 <= JD_Long_1_largest.max())
emission_mask_short = (JD_Short >= JD_Long_1_largest.min()) & (JD_Short <= JD_Long_1_largest.max())

# Суммируем сигналы на эмиссионном интервале
sum_signal_short = np.sum(Signal_Short[emission_mask_short])
sum_signal_long = np.sum(Signal_Long_1[emission_mask_long])

# Рассчитываем HR для всего эмиссионного интервала
HR_emission = sum_signal_short / sum_signal_long

# Ошибки для суммарных сигналов
sum_signal_short_err = np.sqrt(np.sum(Short_Profile_Err[emission_mask_short]**2))
sum_signal_long_err = np.sqrt(np.sum(CR_Long_1_Err[emission_mask_long]**2))

# Ошибка для HR
HR_emission_err = HR_emission * np.sqrt((sum_signal_short_err/sum_signal_short)**2 + 
                                      (sum_signal_long_err/sum_signal_long)**2)

print(f"\nСпектральная жёсткость для всего эмиссионного интервала:")
print(f"Длительность эмиссии: {m} секунд")
print(f"Сумма сигналов в коротком профиле: {sum_signal_short:.1f} ± {sum_signal_short_err:.1f}")
print(f"Сумма сигналов в длинном профиле: {sum_signal_long:.1f} ± {sum_signal_long_err:.1f}")
print(f"HR для эмиссионного интервала: {HR_emission:.3f} ± {HR_emission_err:.3f}")

# Верхний график
ax1.bar(JD_Long_1, CR_Long_1, width=width, edgecolor='blue', facecolor='none', label='Профиль с CsI')
ax1.errorbar(JD_Long_1, CR_Long_1, yerr=CR_Long_1_Err, fmt='.', color='blue')
ax1.bar(JD_Long_1_largest, CR_Long_1_largest, width=width, edgecolor='orange', facecolor='orange', label='Эмиссионное время')
ax1.plot(JD_Long_1, BG_Long, color='black', label='Уровень фона для CsI', linestyle='--')

ax1.bar(JD_Short, Short_Profile, width=width, edgecolor='green', facecolor='none', label='Профиль со стильбена')
ax1.errorbar(JD_Short, Short_Profile, yerr=Short_Profile_Err, fmt='.', color='green')
ax1.plot(JD_Short, BG_Short, color='purple', label='Уровень фона для стильбена', linestyle='--')

#ax1.axvline(JD_Long_1[N_Short_Profile_Max], linestyle='--', color='red', label='Max Count Rate')
ax1.set_ylabel('Темп счёта, отсчёты/с', fontsize=16)
ax1.legend()
ax1.grid(True)

# Нижний график
NN_HR = np.where(Signal_Short >= Short_Profile_Err * 3)[0]
if len(NN_HR) == 0:
    raise ValueError

ax2.errorbar(JD_Long_1[NN_HR], HR[NN_HR], yerr=HR_Err[NN_HR], fmt='.', color='blue', label=':Жёсткость')
ax2.axhline(HR_emission, linestyle='--', color='green', 
           label=f'HR эмиссионного интервала: {HR_emission:.2f}±{HR_emission_err:.2f}')
ax2.axvline(JD_Long_1[N_Short_Profile_Max], linestyle='--', color='red')
ax2.axhline(0, linestyle='--', color='black')
ax2.set_xlabel('UTC', fontsize=16)
ax2.set_ylabel('Harness Ratio')
ax2.set_ylim([-0.1, 0.25])
ax2.legend()
ax2.grid(True)

ax2.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))

xlim_min = JD_Long_1.min()
xlim_max = JD_Long_1.max()
ax1.set_xlim([xlim_min, xlim_max])

plt.suptitle(GRB_Name)
plt.show()
