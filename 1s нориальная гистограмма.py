from scipy.io import readsav
import ctypes
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates

file_path = 'C:\\Users\\Grengald\\Desktop\\Котинг\\ИКИ РАН\\Sav bd\\GRB180720B.sav'
data = readsav(file_path)


data_structure = {key: type(value) for key, value in data.items()}
data_structure

JD_Long_1 = data['jd_long_1']
CR_Long_1 = data['cr_long_1']
Short_Profile = data['short_profile']
CR_Long_1_Err = data['cr_long_1_err']
Short_Profile_Err = data['short_profile_err']
JD_Short = data['jd_short']


#print("jd_long_1:")
#print(JD_Long_1)

#print("jd_short:")
#print(JD_Short)

#print("\ncr_long_1:")
print(len(CR_Long_1))

#print("\nshort_profile:")
#print(Short_Profile)

#print("\ncr_long_1_err:")
#print(CR_Long_1_Err)

#print("\nshort_profile_err:")
#print(Short_Profile_Err)


N_GRB = 0
GRB_Names = ['GRB180720B']
SJDs_Max = [6320.09707600]
T90s = [24.46]

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
#print(JD_BG_R_2)
#print("jd_long_1:")
#print(JD_Long_1)
#stop

NN_BG_1_L = np.where((JD_Long_1 >= JD_BG_L_1) & (JD_Long_1 <= JD_BG_L_2))[0]
#print(NN_BG_1_L)
NN_BG_1_R = np.where((JD_Long_1 >= JD_BG_R_1) & (JD_Long_1 <= JD_BG_R_2))[0]
#print(NN_BG_1_R)
NN_BG_1 = np.concatenate((NN_BG_1_L, NN_BG_1_R))
#print(NN_BG_1)
Cut_JD_1 = (JD_Long_1 - JD_Long_1[0]) * 86400.0
#print(len(Cut_JD_1[NN_BG_1]))
#print(Cut_JD_1[NN_BG_1])
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
#print(Short_Profile_Max)
N_Short_Profile_Max = np.argmax(Short_Profile)
#print(N_Short_Profile_Max)
print('Maximum signals(Short): ', Signal_Short[N_Short_Profile_Max], Signal_Long_1[N_Short_Profile_Max])
print('HR max bin(short):', HR[N_Short_Profile_Max], HR_Err[N_Short_Profile_Max])

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10), sharex=True)

# Define bin width
width = 0.0000116  # Adjust if needed

# Plot long profile as an outlined histogram
ax1.hist(JD_Long_1, bins=len(JD_Long_1), weights=CR_Long_1, histtype="step",
         edgecolor='blue', label='Профиль с CsI')

# Plot short profile as an outlined histogram, scaled by 5
ax1.hist(JD_Short, bins=len(JD_Short), weights=Short_Profile, histtype="step",
         edgecolor='green', label='Профиль со стильбена')

# Background curves
ax1.plot(JD_Long_1, BG_Long, color='black', label='Уровень фона для CsI', linestyle='--')
ax1.plot(JD_Short, BG_Short, color='purple', label='Уровень фона для стильбена', linestyle='--')

# Max count rate vertical line
#ax1.axvline(JD_Long_1[N_Short_Profile_Max], linestyle='--', color='red', label='Max Count Rate')

ax1.set_ylabel('Темп счёта, отсчёты/с', fontsize=16)
ax1.legend()
ax1.grid(True)
NN_HR = np.where(Signal_Short>=Short_Profile_Err*3)[0]

if len(NN_HR) == 0:
    raise ValueError
# Lower plot: Hardness Ratio
ax2.plot(JD_Long_1[NN_HR], HR[NN_HR], '+', markersize=5, label='Жёсткость')
ax2.errorbar(JD_Long_1[NN_HR], HR[NN_HR], yerr=HR_Err[NN_HR], color='blue', fmt='+')
#ax2.axvline(JD_Long_1[N_Short_Profile_Max], linestyle='--', color='red')
ax2.axhline(0, linestyle='--', color='black')
ax2.set_xlabel('UTC', fontsize=16)
ax2.set_ylabel('Жёсткость', fontsize=16)
ax2.set_ylim([-0.1, 0.25])
ax2.legend()
ax2.grid(True)

# Time formatting for x-axis
ax2.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))

# X-axis limits
ax1.set_xlim([JD_Long_1.min(), JD_Long_1.max()])

# Overall title
plt.suptitle(GRB_Name)
plt.show()
