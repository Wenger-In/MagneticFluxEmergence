import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
import sunpy.map
from sunpy.map.header_helper import make_heliographic_header
from astropy.io import fits

HARPNUM = '11426'
REC_date = '20240630'
REC_time = '000000'

## import data
Br_dir = 'E:/Research/Data/HMI/SHARP/AR3723/'
Br_file = 'hmi.sharp_cea_720s.' + HARPNUM + '.' + REC_date + '_' + REC_time + '_TAI.Br.fits'
Br_path = Br_dir + Br_file
Br_data, Br_header = fits.getdata(Br_path, header=True)

Br_map = sunpy.map.Map(Br_path)
fig = plt.figure(figsize=(18,5))
abs_max = np.nanmax(np.abs(Br_data))
ax = fig.add_subplot(projection=Br_map)
im = Br_map.plot(axes=ax, cmap='RdBu', vmin=-abs_max, vmax=abs_max)
cbar = fig.colorbar(im, ax=ax)
cbar.set_label('Br')

print('NOAA AR', Br_header['NOAA_AR'])
print('Carr Lon', Br_header['CRVAL1'])
print('Carr Lat', Br_header['CRVAL2'])

plt.show()

db