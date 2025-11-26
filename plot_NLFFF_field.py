import pyvista as pv
import numpy as np
import pandas as pd
import pyfits as py
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from mpl_toolkits.mplot3d import Axes3D

from utils_field_visualize import plot_arrows, plot_streamlines, plot_boundary_plane

########## Importing magnetic fields under the photosphere

# Reading data from saved files
Ve = 0.5
sparse_factor_under = 3
save_under_dir = 'E:/Research/Work/202405_solar_storm/reconstruction/sparse_factor=' + str(sparse_factor_under) + '/'
data_type = 1 # 0 for original field; 1 for resulted field
if data_type == 0:
    save_under_name = 'origin_field_data_' + str(Ve) + '.csv'
elif data_type == 1: 
    save_under_name = 'result_field_data_' + str(Ve) + '.csv'
data = pd.read_csv(save_under_dir+save_under_name)

# Extracting dataframes
xx_under = data['xx'].values
yy_under = data['yy'].values
zz_under = data['zz'].values

if data_type == 0:
    Bx_under = data['data_cube_Bx'].values
    By_under = data['data_cube_By'].values
    Bz_under = data['data_cube_Bz'].values
    if sparse_factor_under == 3:
        nx, ny, nz = 401, 171, 149
    elif sparse_factor_under == 5:
        nx, ny, nz = 241, 103, 90
    elif sparse_factor_under == 10: 
        nx, ny, nz = 121, 52, 46
    elif sparse_factor_under == 25:
        nx, ny, nz = 49, 21, 19
elif data_type == 1:
    Bx_under = data['Bx'].values
    By_under = data['By'].values
    Bz_under = data['Bz'].values  
    if sparse_factor_under == 3:
        nx, ny, nz = 257, 171, 149
    elif sparse_factor_under == 5:
        nx, ny, nz = 129, 103, 89
    elif sparse_factor_under == 10: 
        nx, ny, nz = 65, 52, 46
    elif sparse_factor_under == 25:
        nx, ny, nz = 33, 21, 19

# Reshaping saved data based on considering time range
xx_under = xx_under.reshape((nx, ny, nz))
yy_under = yy_under.reshape((nx, ny, nz))
zz_under = zz_under.reshape((nx, ny, nz))
Bx_under = Bx_under.reshape((nx, ny, nz))
By_under = By_under.reshape((nx, ny, nz))
Bz_under = Bz_under.reshape((nx, ny, nz))

Bx_under = np.nan_to_num(Bx_under, nan=0.0)
By_under = np.nan_to_num(By_under, nan=0.0)
Bz_under = np.nan_to_num(Bz_under, nan=0.0)

########## Importing magnetic fields upon the photosphere
save_upper_dir = 'E:/Research/Program/else/NLFFF_opti/sparse_factor=10/'
save_upper_name = 'nonlinearfff.fits'

hdu = py.open(save_upper_dir + save_upper_name)
data = hdu[0].data
bx_upper = data[0,:,:,:]
by_upper = data[1,:,:,:]
bz_upper = data[2,:,:,:]

bx_upper = np.nan_to_num(bx_upper, nan=0.0)
by_upper = np.nan_to_num(by_upper, nan=0.0)
bz_upper = np.nan_to_num(bz_upper, nan=0.0)

m, n, k = bx_upper.shape

sparse_factor_upper = 10
dl_pix = 364.6 # [km]
dx, dy, dz = dl_pix * sparse_factor_upper, dl_pix * sparse_factor_upper, dl_pix * sparse_factor_upper

x_arr_upper = np.linspace(0, dx*(m-1), m)
y_arr_upper = np.linspace(0, dy*(n-1), n)
z_arr_upper = np.linspace(0, dz*(k-1), k)

yy_upper, xx_upper, zz_upper = np.meshgrid(y_arr_upper, x_arr_upper, z_arr_upper)

########## Plotting the fields upon the photosphere 
sc_factor = 1e1
plotter = plot_streamlines(xx_upper, yy_upper, zz_upper, bx_upper, by_upper, bz_upper, \
    label='|B|', scale_factor=sc_factor, start_point_type=2, Bc=0.3, thres_factor=10000, clim_factor=3, plot_box=False)

########## Plotting the boundaey fields 
sc_factor = 1e-3
plot_boundary_plane(xx_under, yy_under, zz_under, Bx_under, By_under, Bz_under, plotter, scale_factor=sc_factor, which_boundary='top')

plotter.add_axes()
plotter.show()

db