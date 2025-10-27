import numpy as np
import matplotlib.pyplot as plt
import sunpy.map

hmi_dir = 'E:/Research/Data/SDO/HMI/Magnetogram/720s/'
hmi_file = 'hmi.M_720s.20240507_000000_TAI.3.magnetogram.fits'
hmi_path = hmi_dir + hmi_file

hmi_map = sunpy.map.Map(hmi_path)

hmi_rotated = hmi_map.rotate(order=3)

fig = plt.figure()
ax = fig.add_subplot(projection=hmi_rotated)
hmi_rotated.plot(axes=ax, label=hmi_rotated.unit.to_string())
plt.colorbar()
plt.grid('off')

plt.show()