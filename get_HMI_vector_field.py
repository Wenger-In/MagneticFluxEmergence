import os

os.environ['PYDEVD_WARN_SLOW_RESOLVE_TIMEOUT'] = '3.0'

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.io import fits
import sunpy.coordinates
from sunpy.map import Map
from datetime import datetime as dt_obj
import drms
from disambiguation import Basic, CoordinateTransform

recordset = 'hmi.B_720s[2024.05.07_00:00:00]'
method = 2  # Method used to disambiguate the data
save_dir = 'E:/Research/Data/HMI/VectorField/Bptr/'

basic_instance = Basic(recordset, method)
data = Basic.get_data(basic_instance)

# disambiguation
disambiguated_azimuth = Basic.perform_disambiguation(basic_instance, data[1], data[4])

# convert field to spherical coordinate components (B_phi, B_theta, radial B_r) on the CCD grid
coord_transform = CoordinateTransform(disambiguated_azimuth, data[2], data[3], data[0])
latlon, bptr = coord_transform.ccd(coord_transform)

# set filename
date_time_str = recordset.split('[')[1].split(']')[0]
date_time_format = '%Y.%m.%d_%H:%M:%S'
file_date_time = dt_obj.strptime(date_time_str, date_time_format)
filename_date_time = file_date_time.strftime('%Y%m%d%H%M%S')

# save file
filename = f"{filename_date_time}_bptr.dat"
bptr.tofile(os.path.join(save_dir, filename))

# plot Br
plt.figure()
abs_max = np.nanmax(np.abs(bptr[:,:,2]))
plt.pcolormesh(bptr[:,:,2], cmap='RdBu', vmin=-abs_max, vmax=abs_max)
plt.colorbar()
plt.show()