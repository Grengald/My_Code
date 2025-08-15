from scipy.io import readsav
import ctypes
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates

file_path = 'C:\\Users\\Grengald\\Desktop\\Котинг\\ИКИ РАН\\Sav bd\\GRB080319B.sav'
data = readsav(file_path)

JD_Long_1 = data['jd_long_1']
CR_Long_1 = data['cr_long_1']
Short_Profile = data['short_profile']
CR_Long_1_Err = data['cr_long_1_err']
Short_Profile_Err = data['short_profile_err']
JD_Short = data['jd_short']

N_GRB = 0
GRB_Names = ['GRB080319B']
SJDs_Max = [2544.7599999998]
T90s = [43.600]

GRB_Name = GRB_Names[N_GRB]
SJD_Max = SJDs_Max[N_GRB]
T90 = T90s[N_GRB]

DT = 120.0 / 86400.0
sjd_shift = 2452000.0
JD_Max = SJD_Max + sjd_shift

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

N_Short_Profile_Max = np.argmax(Signal_Short)

# Предварительно уберем ненужные переменные для ясности и сделаем проверку на выход за границы массива
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

Left_JD_Short =JD_Short[:N_Short_Profile_Max]
Right_JD_Short =JD_Short[(N_Short_Profile_Max+1):]

Signal_Long_1_S = Signal_Long_1[N_Short_Profile_Max]
print(Short_Profile_Err[N_Short_Profile_Max])
Signal_Short_S = Signal_Short[N_Short_Profile_Max]
Short_Profile_Err_S=Short_Profile_Err[N_Short_Profile_Max]
CR_Long_1_Err_S=CR_Long_1_Err[N_Short_Profile_Max]
CR_Long_1_S=CR_Long_1[N_Short_Profile_Max]

SHR = []
SHR_0 = HR[N_Short_Profile_Max] / HR_Err[N_Short_Profile_Max]
SHR.append(SHR_0)
n = 1
m = 1
CR_Long_1_Err_S_1=CR_Long_1_Err[N_Short_Profile_Max] ** 2
Short_Profile_Err_S_1=Short_Profile_Err[N_Short_Profile_Max] ** 2
while True:
    if n < len(Left_Short_Profile) and (m >= len(Right_Short_Profile) or Left_Short_Profile[-n] > Right_Short_Profile[m-1]):
        Signal_Long_1_S += Left_Signal_Long_1[-n]
        Signal_Short_S += Left_Short_Profile[-n]
        #print(Left_Short_Profile[-n])
        CR_Long_1_S+= Left_CR_Long_1[-n]
        CR_Long_1_Err_S_1+=Left_CR_Long_1_Err[-n] ** 2
        #print(Left_CR_Long_1_Err[-n])
        Short_Profile_Err_S_1+=Left_Short_Profile_Err[-n] ** 2
        print(Left_Short_Profile_Err[-n])
        CR_Long_1_Err_S=np.sqrt(CR_Long_1_Err_S_1)
        Short_Profile_Err_S=np.sqrt(Short_Profile_Err_S_1)
        HR_Err_S = np.sqrt((Short_Profile_Err_S / Signal_Long_1_S)**2 + (Signal_Short_S * CR_Long_1_Err_S /  CR_Long_1_S**2)**2)
        n += 1
    else:
        Signal_Long_1_S += Right_Signal_Long_1[m-1]
        Signal_Short_S += Right_Short_Profile[m-1]
        #print(Right_Short_Profile[m-1])
        CR_Long_1_S+= Right_CR_Long_1[m-1]
        CR_Long_1_Err_S_1+=Right_CR_Long_1_Err[m-1] ** 2
        #print(Right_CR_Long_1_Err[m-1])
        Short_Profile_Err_S_1+=Right_Short_Profile_Err[m-1] ** 2
        CR_Long_1_Err_S=np.sqrt(CR_Long_1_Err_S_1)
        Short_Profile_Err_S=np.sqrt(Short_Profile_Err_S_1)
        print(Right_Short_Profile_Err[m-1])
        HR_Err_S = np.sqrt((Short_Profile_Err_S / Signal_Long_1_S)**2 + (Signal_Short_S * CR_Long_1_Err_S /  CR_Long_1_S**2)**2)
        m += 1
    HR_S = Signal_Short_S / Signal_Long_1_S
    SHR_n = HR_S / HR_Err_S
    SHR.append(SHR_n)
    
    if Signal_Short_S< 3*Short_Profile_Err_S:
        break
    

# Построение графика с заданными пределами
plt.figure(figsize=(10, 6))
plt.plot(range(1, len(SHR) + 1), SHR, marker='.')
plt.xlabel('Количество бинов')
plt.ylabel('Значимость HR (SHR)')
plt.title(f'Зависимость значимости HR для {GRB_Name}')
plt.xlim(0, len(SHR)+1)  # Пределы для количества бинов
plt.ylim(0, 60)    # Пределы для значимости HR
plt.grid(True)    # Добавляем сетку для наглядности
plt.show()







