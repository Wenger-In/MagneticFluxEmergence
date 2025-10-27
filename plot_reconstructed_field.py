import pyvista as pv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from mpl_toolkits.mplot3d import Axes3D

from field_analyse_utils import plot_arrows, plot_streamlines, plot_boundary_plane

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
zz = zz.reshape((nx, ny, nz))
Bx = Bx.reshape((nx, ny, nz))
By = By.reshape((nx, ny, nz))
Bz = Bz.reshape((nx, ny, nz))

Bx = np.nan_to_num(Bx, nan=0.0)
By = np.nan_to_num(By, nan=0.0)
Bz = np.nan_to_num(Bz, nan=0.0)

# plotter = plot_arrows(xx, yy, zz, Bx, By, Bz, step=3, scale_factor=1)

sc_factor = 1e-3

plotter = plot_streamlines(xx, yy, zz, Bx, By, Bz, label='Bz', \
    scale_factor=sc_factor, start_point_type=1, Bc=2000, thres_factor=50, clim_factor=1, plot_box=False)

plot_boundary_plane(xx, yy, zz, Bx, By, Bz, plotter, scale_factor=sc_factor)

plotter.camera_position = [(np.max(xx)/4*3, -(np.max(xx)*5), 1), (np.max(xx)/4*3, np.max(yy)/2, np.min(zz)), (0, 0, 1)]

plotter.add_axes()
plotter.show()

db