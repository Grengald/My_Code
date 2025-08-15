from scipy.io import readsav
import ctypes
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates
import math

file_path = 'C:\\Users\\Grengald\\Desktop\\Котинг\\ИКИ РАН\\Sav bd\\GRB100324B.sav'
data = readsav(file_path)


data_structure = {key: type(value) for key, value in data.items()}
data_structure

JD_Long_1 = data['jd_long_1']
CR_Long_1 = data['cr_long_1']
Short_Profile = data['short_profile']
#CR_Long_1_Err = data['cr_long_1_err']
#Short_Profile_Err = data['short_profile_err']
JD_Short = data['jd_short']

N_GRB = 0
GRB_Names = ['GRB100324B']
SJDs_Max = [3279.6722242600]
T90s = [20]

GRB_Name = GRB_Names[N_GRB]
SJD_Max = SJDs_Max[N_GRB]
T90 = T90s[N_GRB]

DT = 120.0 / 86400.0
sjd_shift = 2452000.0
JD_Max = SJD_Max + sjd_shift
#print(SJD_Max)


JD1 = sjd_shift + SJD_Max - 1.0 * DT
JD2 = sjd_shift + SJD_Max + 2.0 * DT

T_Start = datetime.strptime("2009-11-09 21:49:03", "%Y-%m-%d %H:%M:%S")
T_Stop = T_Start + timedelta(seconds=2 * 86400)


JD_BG_L_1 = JD_Max - (10.0 + T90 + T90 + 60.0) / 86400.0
#print(JD_BG_L_1)
JD_BG_L_2 = JD_Max - (10.0 + T90) / 86400.0
#print(JD_BG_L_2)
JD_BG_R_1 = JD_Max + (10.0 + 2.0 * T90) / 86400.0
#print(JD_BG_R_1)
JD_BG_R_2 = JD_Max + (10.0 + 2.0 * T90 + T90 + 60.0) / 86400.0


#N_Short_Profile_Max = np.argmax(Short_Profile)
#N_Short_Profile = len(Short_Profile)
N_CR_Long_1_Max = np.argmax(CR_Long_1)
N_CR_Long_1 = len(CR_Long_1)
#Перебинировка(по максимуму короткого профиля)
##################################################################################################################
if CR_Long_1[N_CR_Long_1_Max-1]>CR_Long_1[N_CR_Long_1_Max+1]:
    if CR_Long_1[N_CR_Long_1_Max]%2==0:
        N_CR_Long_1_Rebin = N_CR_Long_1 // 2  
        n = 0
        a = 1
        m = 0
        Short_Profile_Rebin = [0] * N_CR_Long_1_Rebin
        CR_Long_1_Rebin = [0] * N_CR_Long_1_Rebin
        JD_Long_1_Rebin = [0] * N_CR_Long_1_Rebin
        JD_Short_Rebin = [0] * N_CR_Long_1_Rebin
        while a < N_CR_Long_1:
            Short_Profile_Rebin[m] = (Short_Profile[n] + Short_Profile[a]) / 2
            CR_Long_1_Rebin[m] = (CR_Long_1[n]+CR_Long_1[a])/2
            JD_Long_1_Rebin[m] = (JD_Long_1[n]+JD_Long_1[a])/2
            JD_Short_Rebin[m] = (JD_Short[n]+JD_Short[a])/2
            n += 2
            a = n + 1
            m += 1
        #print(Short_Profile_Rebin)

            
    elif CR_Long_1[N_CR_Long_1_Max]%2==1:
        N_CR_Long_1_Rebin = N_CR_Long_1 // 2 - 1  
        n = 1
        a = 2
        m = 0
        Short_Profile_Rebin = [0] * N_CR_Long_1_Rebin
        CR_Long_1_Rebin = [0] * N_CR_Long_1_Rebin
        JD_Long_1_Rebin = [0] * N_CR_Long_1_Rebin
        JD_Short_Rebin = [0] * N_CR_Long_1_Rebin
        while a < N_CR_Long_1:
            Short_Profile_Rebin[m] = (Short_Profile[n] + Short_Profile[a]) / 2  
            CR_Long_1_Rebin[m] = (CR_Long_1[n]+CR_Long_1[a])/2
            JD_Long_1_Rebin[m] = (JD_Long_1[n]+JD_Long_1[a])/2
            JD_Short_Rebin[m] = (JD_Short[n]+JD_Short[a])/2
            n += 2
            a = n + 1
            m += 1

elif CR_Long_1[N_CR_Long_1_Max-1]<CR_Long_1[N_CR_Long_1_Max+1]:
    if CR_Long_1[N_CR_Long_1_Max]%2==0:
        N_CR_Long_1_Rebin = N_CR_Long_1 // 2 - 1
        n = 1
        a = 2
        m = 0
        Short_Profile_Rebin = [0] * N_CR_Long_1_Rebin
        CR_Long_1_Rebin = [0] * N_CR_Long_1_Rebin
        JD_Long_1_Rebin = [0] * N_CR_Long_1_Rebin
        JD_Short_Rebin = [0] * N_CR_Long_1_Rebin
        while a < N_CR_Long_1:
            Short_Profile_Rebin[m] = (Short_Profile[n] + Short_Profile[a]) / 2
            CR_Long_1_Rebin[m] = (CR_Long_1[n]+CR_Long_1[a])/2
            JD_Long_1_Rebin[m] = (JD_Long_1[n]+JD_Long_1[a])/2
            JD_Short_Rebin[m] = (JD_Short[n]+JD_Short[a])/2
            n += 2
            a = n + 1
            m += 1
        #print(Short_Profile_Rebin)            
    elif CR_Long_1[N_CR_Long_1_Max]%2==1:
        N_CR_Long_1_Rebin = N_CR_Long_1 // 2 
        n = 0
        a = 1
        m = 0
        Short_Profile_Rebin = [0] * N_CR_Long_1_Rebin
        CR_Long_1_Rebin = [0] * N_CR_Long_1_Rebin
        JD_Long_1_Rebin = [0] * N_CR_Long_1_Rebin
        JD_Short_Rebin = [0] * N_CR_Long_1_Rebin
        while a < N_CR_Long_1:
            Short_Profile_Rebin[m] = (Short_Profile[n] + Short_Profile[a]) / 2
            CR_Long_1_Rebin[m] = (CR_Long_1[n]+CR_Long_1[a])/2
            JD_Long_1_Rebin[m] = (JD_Long_1[n]+JD_Long_1[a])/2
            JD_Short_Rebin[m] = (JD_Short[n]+JD_Short[a])/2
            n += 2
            a = n + 1
            m += 1
###################################################################################################################             
#Новые ошибки
CR_Long_1_Rebin_1 = np.array(CR_Long_1_Rebin)
CR_Long_1_Err_Rebin = np.sqrt(CR_Long_1_Rebin_1*8)/4
Short_Profile_Rebin_1 = np.array(Short_Profile_Rebin)
Short_Profile_Err_Rebin = np.sqrt(Short_Profile_Rebin_1/2)
###################################################################################################################

JD_Long_1_Rebin = np.array(JD_Long_1_Rebin)
JD_Short_Rebin = np.array(JD_Short_Rebin)


NN_BG_1_L = np.where((JD_Long_1_Rebin >= JD_BG_L_1) & (JD_Long_1_Rebin <= JD_BG_L_2))[0]
#print(NN_BG_1_L)
NN_BG_1_R = np.where((JD_Long_1_Rebin >= JD_BG_R_1) & (JD_Long_1_Rebin <= JD_BG_R_2))[0]
#print(NN_BG_1_R)
NN_BG_1 = np.concatenate((NN_BG_1_L, NN_BG_1_R))
#print(NN_BG_1)
Cut_JD_1 = (JD_Long_1_Rebin - JD_Long_1_Rebin[0]) * 86400.0
#print(len(Cut_JD_1[NN_BG_1]))
#print(Cut_JD_1[NN_BG_1])
Cut_JD_1 = np.array(Cut_JD_1)
CR_Long_1_Rebin = np.array(CR_Long_1_Rebin)


Coeff_1 = np.polyfit(Cut_JD_1[NN_BG_1], CR_Long_1_Rebin[NN_BG_1], 2)
BG_Long = np.polyval(Coeff_1, Cut_JD_1)

NN_BG_2_L = np.where((JD_Short_Rebin >= JD_BG_L_1) & (JD_Short_Rebin <= JD_BG_L_2))[0]
NN_BG_2_R = np.where((JD_Short_Rebin >= JD_BG_R_1) & (JD_Short_Rebin <= JD_BG_R_2))[0]
NN_BG_2 = np.concatenate((NN_BG_2_L, NN_BG_2_R))
Cut_JD_2 = (JD_Short_Rebin - JD_Short_Rebin[0]) * 86400.0
Cut_JD_2 = np.array(Cut_JD_2)
Short_Profile_Rebin = np.array(Short_Profile_Rebin)

Coeff_2 = np.polyfit(Cut_JD_2[NN_BG_2], Short_Profile_Rebin[NN_BG_2], 2)
BG_Short = np.polyval(Coeff_2, Cut_JD_2)



Signal_Short = Short_Profile_Rebin - BG_Short
Signal_Long_1 = CR_Long_1_Rebin - BG_Long

#print(Short_Profile_Err_Rebin)
#print(Short_Profile_Rebin)
#print(CR_Long_1_Err_Rebin)
#print(CR_Long_1_Rebin)


HR = Signal_Short / Signal_Long_1
HR_Err = np.sqrt((Short_Profile_Err_Rebin / Signal_Long_1)**2 + (Signal_Short * CR_Long_1_Err_Rebin / CR_Long_1_Rebin**2)**2)

CR_Long_1_Max = np.max(CR_Long_1_Rebin)
N_CR_Long_1_Max = np.argmax(CR_Long_1_Rebin)
print('Max CR long', CR_Long_1_Max, Short_Profile_Rebin[N_CR_Long_1_Max])
print('Max HR long', HR[N_CR_Long_1_Max], HR_Err[N_CR_Long_1_Max])

Short_Profile_Max = np.max(Short_Profile_Rebin)
#print(Short_Profile_Max)
N_Short_Profile_Max = np.argmax(Short_Profile_Rebin)
#print(N_Short_Profile_Max)
print('Maximum signals(Short): ', Signal_Short[N_Short_Profile_Max], Signal_Long_1[N_Short_Profile_Max])
print('HR max bin(short):', HR[N_Short_Profile_Max], HR_Err[N_Short_Profile_Max])

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10), sharex=True)

# Верхний график: длинный и короткий профиль в виде гистограммы с погрешностями
width = 0.0000232  # ширина столбцов гистограммы

# Длинный профиль
ax1.bar(JD_Long_1_Rebin, CR_Long_1_Rebin, width=width, color='blue', alpha=0.6, label='Long Profile')
ax1.errorbar(JD_Long_1_Rebin, CR_Long_1_Rebin, yerr=CR_Long_1_Err_Rebin, fmt='.', color='blue')
ax1.plot(JD_Long_1_Rebin, BG_Long, color='#969696', label='Background Long', linestyle='--')

# Короткий профиль
ax1.bar(JD_Short_Rebin, Short_Profile_Rebin*5, width=width, color='green', alpha=0.6, label='Short Profile x5')
ax1.errorbar(JD_Short_Rebin, Short_Profile_Rebin*5, yerr=Short_Profile_Err_Rebin*5, fmt='.', color='green')
ax1.plot(JD_Short_Rebin, BG_Short*5, color='#FFD700', label='Background Short x5', linestyle='--')

# Вертикальная линия для максимальной счетной скорости
ax1.axvline(JD_Long_1_Rebin[N_CR_Long_1_Max], linestyle='--', color='red', label='Max Count Rate')
ax1.set_ylabel('Count Rate, CNT/2S')
ax1.legend()
ax1.grid(True)

NN_HR = np.where(Signal_Short >= Short_Profile_Err_Rebin * 3)[0]

if len(NN_HR) == 0:
    raise ValueError

# Нижний график: отношение жесткости
ax2.plot(JD_Long_1_Rebin[NN_HR], HR[NN_HR], '+', markersize=5, label='Harness Ratio')
ax2.errorbar(JD_Long_1_Rebin[NN_HR], HR[NN_HR], yerr=HR_Err[NN_HR],  color='blue', fmt='+')
ax2.axvline(JD_Long_1_Rebin[N_CR_Long_1_Max], linestyle='--', color='red')
ax2.axhline(0, linestyle='--', color='black')
ax2.set_xlabel('UTC')
ax2.set_ylabel('Harness Ratio')
ax2.set_ylim([-0.1, 0.2]) 
ax2.legend()
ax2.grid(True)

# Настройка формата времени по оси X
ax2.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))


# Настройка пределов осей и заголовков
xlim_min = JD_Long_1.min()
xlim_max = JD_Long_1.max()
ax1.set_xlim([xlim_min, xlim_max])

# Общий заголовок для всей фигуры
plt.suptitle(GRB_Name)
plt.show()
