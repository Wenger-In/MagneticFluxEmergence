import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
import sunpy.map
from sunpy.map.header_helper import make_heliographic_header
from astropy.io import fits

HARPNUM = '11149' # NOAANUM AR3664
REC_date = '20240508'
REC_time = '120000'

## import data
Br_dir = 'E:/Research/Data/HMI/SHARP/'
Br_file = 'hmi.sharp_cea_720s.' + HARPNUM + '.' + REC_date + '_' + REC_time + '_TAI.Br.fits'
Br_path = Br_dir + Br_file
Br_data, Br_header = fits.getdata(Br_path, header=True)

Br_map = sunpy.map.Map(Br_path)
fig = plt.figure()
ax = fig.add_subplot(projection=Br_map)
Br_map.plot(axes=ax)

print(Br_header['NOAA_AR'])

plt.show()

db