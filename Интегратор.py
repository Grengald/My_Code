import numpy as np

y = np.array([1915.8326828512188, 1917.8189364750706, 1755.846400484134, 1855.805159143349, 889.8600877163623, 973.7913525171108, 745.7775160426701,
              367.8737451003882, 279.8873726362118, 847.9009703238331, 271.91453870853127])
x = np.array([2452555.7816361976, 2452555.781624623, 2452555.7816477716, 2452555.781613049, 2452555.7816593456, 2452555.781601475, 2452555.781589901,
              2452555.7816709196, 2452555.7816824936, 2452555.7816940676, 2452555.7817056417])
err_long = np.array([15.6205, 15.0996685, 16.0, 15.0996685, 15.362291, 14.6969385, 15.748015, 16.0, 15.748016, 15.874508, 15.491934])

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

w = np.array([157.23481104254626, 149.23655627282793, 77.23289446067635, 53.23813001109406, 51.230806460452904, 29.239532331006654, 9.240763232565719, 7.228547041875924,
              9.226116204945416, 3.223513949661381, 9.220740276023815])
z = np.array([2452555.7816361976, 2452555.781624623, 2452555.7816477716, 2452555.781613049, 2452555.7816593456, 2452555.781601475, 2452555.781589901,
              2452555.7816709196, 2452555.7816824936, 2452555.7816940676, 2452555.7817056417])
err_short = np.array([9.0, 8.774964, 6.4031243, 5.3851647, 5.2915025, 4.1231055, 2.6457512, 2.4494898, 2.6457512, 2.0, 2.6457512])

w_sorted=w[sort_order]
z_sorted=z[sort_order]
err_short_sorted=err_short[sort_order]

integral_2 = np.trapz(w_sorted, z_sorted)
dz = np.diff(z_sorted)
err_short_shortened = err_short_sorted[:-1]
error_integral_2 = np.sqrt(np.sum((err_short_shortened * dz)**2))
print("Integral Short Profile:", integral_2)
print('Error Integral Short Profile:',error_integral_2)

HR = integral_2 / integral_1
HR_error = HR * np.sqrt((error_integral_1 / integral_1)**2 + (error_integral_2 / integral_2)**2)
print('HR:', HR)
print('Error HR:',HR_error)



