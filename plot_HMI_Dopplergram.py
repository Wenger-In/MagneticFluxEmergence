import numpy as np
import matplotlib.pyplot as plt
import sunpy.map

fdisk_dir = 'C:/Users/rzhuo/Desktop/'
fdisk_fn = 'hmi.V_720s.20240507_000000_TAI.3.Dopplergram.fits'
sharp_dir = 'E:/Research/Data/SDO/HMI/SHARP/Dopplergram/'
sharp_fn = 'hmi.sharp_cea_720s.11149.20240507_000000_TAI.Dopplergram.fits'

fdisk_map = sunpy.map.Map(fdisk_dir+fdisk_fn)
sharp_map = sunpy.map.Map(sharp_dir+sharp_fn)

fdisk_data = fdisk_map.data
sharp_data = sharp_map.data

plt.figure()
plt.subplot(121)
plt.pcolormesh(np.fliplr(np.flipud(fdisk_data)), cmap='bwr', vmin=1000, vmax=3000)
plt.colorbar()

plt.subplot(122)
plt.pcolormesh(sharp_data, cmap='bwr', vmin=1000, vmax=3000)
plt.colorbar()

plt.show()