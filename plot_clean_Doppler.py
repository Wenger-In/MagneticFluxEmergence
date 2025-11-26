import numpy as np
import matplotlib.pyplot as plt
import sunpy.map

Br_path = 'E:\Research\Data\SDO\HMI\SHARP\AR3664\hmi.sharp_cea_720s.11149.20240507_000000_TAI.Br.fits'
Dop_raw_path = 'C:/Users/rzhuo/Desktop/hmi.sharp_cea_720s.11149.20240507_000000_TAI.Dopplergram.fits'
Dop_res_path = 'C:/Users/rzhuo/Desktop/HMIDATA/data_analysisresidual_2024_000.npy'

Br_map = sunpy.map.Map(Br_path)
Br = Br_map.data
Dop_raw_map = sunpy.map.Map(Dop_raw_path)
Dop_raw = Dop_raw_map.data
Dop_res = np.load(Dop_res_path)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8)) 

im1 = ax1.pcolormesh(Dop_raw, cmap='PiYG', vmin=-3000, vmax=3000)
cbar1 = fig.colorbar(im1, ax=ax1, label='Doppler (m/s)')

levels = [-1000, -500, 500, 1000]
im2 = ax1.contour(Br, levels=levels, cmap='bwr')
cbar2 = fig.colorbar(im2, ax=ax1, label='Bz (G)', extend='both')
ax1.axis('equal')
ax1.set_title('Raw')

im3 = ax2.pcolormesh(Dop_res, cmap='PiYG', vmin=-1500, vmax=1500)
cbar3 = fig.colorbar(im3, ax=ax2, label='Doppler (m/s)')

levels = [-1000, -500, 500, 1000]
im4 = ax2.contour(Br, levels=levels, cmap='bwr')
cbar4 = fig.colorbar(im4, ax=ax2, label='Bz (G)', extend='both')
ax2.axis('equal')
ax2.set_title('Residual')

plt.suptitle('Dopplergram at 20240507_000000', fontsize=14)
plt.tight_layout()

plt.show()