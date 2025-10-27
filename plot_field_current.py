import pyvista as pv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from mpl_toolkits.mplot3d import Axes3D
from scipy.ndimage import gaussian_filter, uniform_filter, median_filter

from field_analyse_utils import plot_arrows, plot_streamlines, plot_boundary_plane, \
    calc_field_curl, calc_curl_B_based_on_SFD, calc_field_divergence

# Reading data from saved files
Ve = 0.5
sparse_factor = 3
save_dir = 'E:/Research/Work/202405_solar_storm/reconstruction/sparse_factor=' + str(sparse_factor) + '/'
data_type = 1 # 0 for original field; 1 for resulted field
if data_type == 0:
    save_name = 'origin_field_data_' + str(Ve) + '.csv'
elif data_type == 1: 
    save_name = 'result_field_data_' + str(Ve) + '.csv'
data = pd.read_csv(save_dir+save_name)

# Extracting dataframes
xx = data['xx'].values
yy = data['yy'].values
zz = data['zz'].values

if data_type == 0:
    Bx = data['data_cube_Bx'].values
    By = data['data_cube_By'].values
    Bz = data['data_cube_Bz'].values
    if sparse_factor == 3:
        nx, ny, nz = 401, 171, 149
    elif sparse_factor == 5:
        nx, ny, nz = 241, 103, 90
    elif sparse_factor == 10: 
        nx, ny, nz = 121, 52, 46
    elif sparse_factor == 25:
        nx, ny, nz = 49, 21, 19
elif data_type == 1:
    Bx = data['Bx'].values
    By = data['By'].values
    Bz = data['Bz'].values  
    if sparse_factor == 3:
        nx, ny, nz = 257, 171, 149
    elif sparse_factor == 5:
        nx, ny, nz = 129, 103, 89
    elif sparse_factor == 10: 
        nx, ny, nz = 65, 52, 46
    elif sparse_factor == 25:
        nx, ny, nz = 33, 21, 19

# Reshaping saved data based on considering time range
xx = xx.reshape((nx, ny, nz))
yy = yy.reshape((nx, ny, nz))
zz = zz.reshape((nx, ny, nz)) # [km]
Bx = Bx.reshape((nx, ny, nz))
By = By.reshape((nx, ny, nz))
Bz = Bz.reshape((nx, ny, nz)) # [G]

Bx = np.nan_to_num(Bx, nan=0.0)
By = np.nan_to_num(By, nan=0.0)
Bz = np.nan_to_num(Bz, nan=0.0)

# curl_B_x_reg, curl_B_y_reg, curl_B_z_reg = calc_field_curl(Bx, By, Bz, xx, yy, zz)
curl_B_x_sfd, curl_B_y_sfd, curl_B_z_sfd = calc_curl_B_based_on_SFD(Bx, By, Bz, xx, yy, zz)

# curl_B_x, curl_B_y, curl_B_z = curl_B_x_reg, curl_B_y_reg, curl_B_z_reg
curl_B_x, curl_B_y, curl_B_z = curl_B_x_sfd, curl_B_y_sfd, curl_B_z_sfd

div_curl_B = calc_field_divergence(curl_B_x, curl_B_y, curl_B_z, xx, yy, zz)

mu0 = 4 * np.pi * 1e-7
Jx, Jy, Jz = curl_B_x / mu0, curl_B_y / mu0, curl_B_z / mu0 # [1e-7 A]
Jx, Jy, Jz = Jx / 10, Jy / 10, Jz / 10 # [uA]

# Smoothing data
sigma_gs = 1.0
Jx_gs = gaussian_filter(Jx, sigma=sigma_gs)
Jy_gs = gaussian_filter(Jy, sigma=sigma_gs)
Jz_gs = gaussian_filter(Jz, sigma=sigma_gs)

size_us = 3
Jx_us = uniform_filter(Jx, size=size_us)
Jy_us = uniform_filter(Jy, size=size_us)
Jz_us = uniform_filter(Jz, size=size_us)

sc_factor = 1e-6

plotter = plot_streamlines(xx, yy, zz, Jx_us, Jy_us, Jz_us, label='J (uA)', scale_factor=sc_factor, start_point_type=1, Bc=5e4)

plot_boundary_plane(xx, yy, zz, Bx, By, Bz, plotter)

plotter.add_axes()
plotter.show()

Fl_x = np.zeros_like(Bx)
Fl_y = np.zeros_like(By)
Fl_z = np.zeros_like(Bz)
for ix in range(nx):
    for iy in range(ny):
        for iz in range(nz):
            J_vec = [Jx[ix,iy,iz], Jy[ix,iy,iz], Jz[ix,iy,iz]]
            B_vec = [Bx[ix,iy,iz], By[ix,iy,iz], Bz[ix,iy,iz]]
            Fl_vec = np.cross(J_vec, B_vec)
            Fl_x[ix,iy,iz], Fl_y[ix,iy,iz], Fl_z[ix,iy,iz] = Fl_vec[0], Fl_vec[1], Fl_vec[2]

plotter = plot_arrows(xx, yy, zz, Fl_x, Fl_y, Fl_z, step=2, scale_factor=1e-5)

plot_boundary_plane(xx, yy, zz, Bx, By, Bz, plotter)

plotter.add_axes()
plotter.show()

Fl_z_zmean = np.mean(Fl_z, axis=2)

fig = plt.figure()
Fl_z_zmean_abs_max = np.max(np.abs(Fl_z_zmean))
plt.pcolormesh(xx[:,:,0]/1e3, yy[:,:,0]/1e3, Fl_z_zmean, cmap='bwr', \
               vmin=-Fl_z_zmean_abs_max, vmax=Fl_z_zmean_abs_max)
plt.colorbar()
plt.axis('equal')
plt.xlabel('X [Mm]')
plt.ylabel('Y [Mm]')
plt.title('Fz averaged along z-axis')
fig.tight_layout()
plt.show()

db