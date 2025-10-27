import numpy as np
import matplotlib.pyplot as plt
import sunpy.map

NOAA_num = 3664
SHARP_num = 11149
REC_date = '20240507'
REC_time = '000000'

file_dir = 'E:/Research/Data/SDO/HMI/SHARP/AR' + str(NOAA_num) + '/'
file_Bp = 'hmi.sharp_cea_720s.' + str(SHARP_num) + '.' + REC_date + '_' + REC_time + '_TAI.Bp.fits'
file_Bt = 'hmi.sharp_cea_720s.' + str(SHARP_num) + '.' + REC_date + '_' + REC_time + '_TAI.Bt.fits'
file_Br = 'hmi.sharp_cea_720s.' + str(SHARP_num) + '.' + REC_date + '_' + REC_time + '_TAI.Br.fits'

Br_map = sunpy.map.Map(file_dir + file_Br)

fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(projection=Br_map)
Br_map.plot(axes=ax, vmin=-1500, vmax=1500)
plt.xlabel('')
plt.ylabel('')

plt.show()