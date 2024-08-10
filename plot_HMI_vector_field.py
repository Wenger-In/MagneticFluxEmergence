import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
import sunpy.map
from sunpy.map.header_helper import make_heliographic_header
from astropy.io import fits

record_date = '20240508'
record_time = '120000'

## import data
field_dir = 'E:/Research/Data/HMI/Field/'
field_file = 'hmi.B_720s.' + record_date + '_' + record_time + '_TAI.field.fits'
field_path = field_dir + field_file
field_data, field_header = fits.getdata(field_path, header=True)
field_map = sunpy.map.Map(field_path)

Bptr_dir = 'E:/Research/Data/HMI/VectorField/Bptr/'
Bptr_file = record_date + record_time + '_bptr.dat'
Bptr = np.fromfile(Bptr_dir + Bptr_file).reshape(4096, 4096, 3)
Bp = np.squeeze(Bptr[:,:,0])
Bt = np.squeeze(Bptr[:,:,1])
Br = np.squeeze(Bptr[:,:,2])

## construct Bp Bt Br map
Bp_map = sunpy.map.Map(Bp, field_map.meta)
Bt_map = sunpy.map.Map(Bt, field_map.meta)
Br_map = sunpy.map.Map(Br, field_map.meta)

# plot_map = Br_map
# abs_max = np.nanmax(np.abs(plot_map.data))
# fig = plt.figure()
# ax = fig.add_subplot(projection=plot_map)
# im = plot_map.plot(axes=ax, vmin=-abs_max, vmax=abs_max)
# plt.colorbar(im, ax=ax, orientation='vertical', label='Br (G)')
# plt.show()

## convert to Carrington Coordinate System
shape = (1800,3600)
carr_header = make_heliographic_header(field_map.date, field_map.observer_coordinate, shape, frame='carrington')

field_map_rep = field_map.reproject_to(carr_header)
Bp_map_rep = Bp_map.reproject_to(carr_header)
Bt_map_rep = Bt_map.reproject_to(carr_header)
Br_map_rep = Br_map.reproject_to(carr_header)

field_rep = field_map_rep.data
Bp_rep = Bp_map_rep.data
Bt_rep = Bt_map_rep.data
Br_rep = Br_map_rep.data

## construct grid
x = np.linspace(0, 360, 3600)
y = np.linspace(-90, 90, 1800)
X, Y = np.meshgrid(x, y)
x_range = (150, 185)
y_range = (-30, -10)
x_indices = np.where((x >= x_range[0]) & (x <= x_range[1]))[0]
y_indices = np.where((y >= y_range[0]) & (y <= y_range[1]))[0]

X_sub = X[y_indices[:, None], x_indices]
Y_sub = Y[y_indices[:, None], x_indices]
field_sub = field_rep[y_indices[:, None], x_indices]
Bp_sub = Bp_rep[y_indices[:, None], x_indices]
Bt_sub = Bt_rep[y_indices[:, None], x_indices]
Br_sub = Br_rep[y_indices[:, None], x_indices]

# ## plot figure 1: total field
# fig, ax = plt.subplots(figsize=(12, 6))

# # total field
# abs_max = np.nanmax(np.abs(field_sub))
# c = ax.pcolormesh(X_sub, Y_sub, field_sub, shading='auto', cmap='Greys')
# fig.colorbar(c, ax=ax, label='Btot (G)')

# # axis settings
# ax.set_aspect('equal')
# ax.set_xlabel('Longitude (deg.)')
# ax.set_ylabel('Latitude (deg.)')
# ax.set_title('HMI total field at ' + record_date + record_time)
# plt.grid()

## plot figure 2: field component Bptr
fig, ax = plt.subplots(figsize=(12, 6))

# Br as background
abs_max = np.nanmax(np.abs(Br_sub))
c = ax.pcolormesh(X_sub, Y_sub, Br_sub, shading='auto', cmap='RdBu', vmin=-abs_max, vmax=abs_max)
fig.colorbar(c, ax=ax, label='Br (G)')

# Bp and Bt as quiver
arrow_interval = 5
X_arrow = X_sub[::arrow_interval, ::arrow_interval]
Y_arrow = Y_sub[::arrow_interval, ::arrow_interval]
Bp_arrow = Bp_sub[::arrow_interval, ::arrow_interval]*0.0005 # shrink Bt 2000 times
Bt_arrow = Bt_sub[::arrow_interval, ::arrow_interval]*0.0005 # shrink Bt 2000 times
q = ax.quiver(X_arrow, Y_arrow, Bp_arrow, Bt_arrow, angles='xy', scale_units='xy', scale=1, color='black')
qk = ax.quiverkey(q, 0.8, 0.08, 0.5, 'Bh 1000G', labelpos='E', coordinates='figure')

# axis settings
ax.set_aspect('equal')
ax.set_xlabel('Longitude (deg.)')
ax.set_ylabel('Latitude (deg.)')
ax.set_title('HMI vector field at ' + record_date + record_time)
plt.grid()
plt.show()

db