import scipy.io
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates

file_dir = 'E:/Research/Work/202405_solar_storm/helicity/'
file_name = 'hf_13664.sav'
data = scipy.io.readsav(file_dir+file_name)

print(data.keys())

HOD= data['time']
dhm = data['dhm']
dhm_shear = data['dhm_shear']
dhm_emerge = data['dhm_emerge']
hm = data['hm']
hc= data['hc']
de = data['de']
e = data['e']
jz_total = data['jz_total']

start_time = datetime(2024, 5, 1, 0, 0, 0)
date_times = [start_time + timedelta(hours=float(hours)) for hours in HOD]

plt.figure(figsize=(6,10))

plt.subplot(6,1,1)
plt.plot(date_times, dhm_shear)
plt.xticks([])
plt.ylabel('dH/dt_braid')

plt.subplot(6,1,2)
plt.plot(date_times, dhm_emerge)
plt.xticks([])
plt.ylabel('dH/dt_emerge')

plt.subplot(6,1,3)
plt.plot(date_times, dhm)
plt.xticks([])
plt.ylabel('dH/dt [Wb m3/s]')

plt.subplot(6,1,4)
plt.plot(date_times, hm)
plt.xticks([])
plt.ylabel('H [Wb m3]')

plt.subplot(6,1,5)
plt.plot(date_times, de)
plt.xticks([])
plt.ylabel('dE/dt [erg/s]')

plt.subplot(6,1,6)
plt.plot(date_times, e)
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m/%d-%H'))
plt.gcf().autofmt_xdate()
plt.ylabel('E [erg]')

plt.show()

db