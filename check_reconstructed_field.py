import pyvista as pv
import numpy as np
import pandas as pd
import pyfits as py
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from mpl_toolkits.mplot3d import Axes3D

# Reading data from saved files
Ve = 0.5
sparse_factor_under = 3
save_dir = 'E:/Research/Work/202405_solar_storm/reconstruction/sparse_factor=' + str(sparse_factor_under) + '/'

origin_field_name = 'origin_field_data_' + str(Ve) + '.csv'
cut_field_name = 'cut_field_data_' + str(Ve) + '.csv'
result_field_name = 'result_field_data_' + str(Ve) + '.csv'

origin_field = pd.read_csv(save_dir + origin_field_name)
cut_field = pd.read_csv(save_dir + cut_field_name)
result_field = pd.read_csv(save_dir + result_field_name)

# Extracting dataframes
xx_origin = origin_field['xx'].values
yy_origin = origin_field['yy'].values
zz_origin = origin_field['zz'].values
Bx_origin = origin_field['data_cube_Bx'].values
By_origin = origin_field['data_cube_By'].values
Bz_origin = origin_field['data_cube_Bz'].values

xx_cut = cut_field['xx'].values
yy_cut = cut_field['yy'].values
zz_cut = cut_field['zz'].values
Bx_cut = cut_field['Bx_cut'].values
By_cut = cut_field['By_cut'].values
Bz_cut = cut_field['Bz_cut'].values

xx_result = result_field['xx'].values
yy_result = result_field['yy'].values
zz_result = result_field['zz'].values
Bx_result = result_field['Bx'].values
By_result = result_field['By'].values
Bz_result = result_field['Bz'].values

if sparse_factor_under == 3:
    nx_origin, ny_origin, nz_origin = 401, 171, 149
    nx_cut, ny_cut, nz_cut = 257, 171, 149
    nx_result, ny_result, nz_result = 257, 171, 149
elif sparse_factor_under == 5:
    nx_origin, ny_origin, nz_origin = 241, 103, 90
    nx_cut, ny_cut, nz_cut = 129, 103, 89
    nx_result, ny_result, nz_result = 129, 103, 89
elif sparse_factor_under == 10: 
    nx_origin, ny_origin, nz_origin = 121, 52, 46
    nx_cut, ny_cut, nz_cut = 65, 52, 46
    nx_result, ny_result, nz_result = 65, 52, 46
elif sparse_factor_under == 25:
    nx_origin, ny_origin, nz_origin = 49, 21, 19
    nx_cut, ny_cut, nz_cut = 33, 21, 19
    nx_result, ny_result, nz_result = 33, 21, 19

# Reshaping saved data based on considering time range
xx_origin = xx_origin.reshape((nx_origin, ny_origin, nz_origin))
yy_origin = yy_origin.reshape((nx_origin, ny_origin, nz_origin))
zz_origin = zz_origin.reshape((nx_origin, ny_origin, nz_origin))
Bx_origin = Bx_origin.reshape((nx_origin, ny_origin, nz_origin))
By_origin = By_origin.reshape((nx_origin, ny_origin, nz_origin))
Bz_origin = Bz_origin.reshape((nx_origin, ny_origin, nz_origin))

xx_cut = xx_cut.reshape((nx_cut, ny_cut, nz_cut))
yy_cut = yy_cut.reshape((nx_cut, ny_cut, nz_cut))
zz_cut = zz_cut.reshape((nx_cut, ny_cut, nz_cut))
Bx_cut = Bx_cut.reshape((nx_cut, ny_cut, nz_cut))
By_cut = By_cut.reshape((nx_cut, ny_cut, nz_cut))
Bz_cut = Bz_cut.reshape((nx_cut, ny_cut, nz_cut))

xx_result = xx_result.reshape((nx_result, ny_result, nz_result))
yy_result = yy_result.reshape((nx_result, ny_result, nz_result))
zz_result = zz_result.reshape((nx_result, ny_result, nz_result))
Bx_result = Bx_result.reshape((nx_result, ny_result, nz_result))
By_result = By_result.reshape((nx_result, ny_result, nz_result))
Bz_result = Bz_result.reshape((nx_result, ny_result, nz_result))

def plot_slice_field(Bn_origin, Bn_cut, Bn_result, normal, index):
    if normal == 'x':
        Bn_origin_slice = Bn_origin[index,:,:]
        Bn_cut_slice = Bn_cut[index,:,:]
        Bn_result_slice = Bn_result[index,:,:]
        hh_origin, vv_origin = yy_origin[index,:,:], zz_origin[index,:,:]
        hh_cut, vv_cut = yy_cut[index,:,:], zz_cut[index,:,:]
        hh_result, vv_result = yy_result[index,:,:], zz_result[index,:,:]
    elif normal == 'y':
        Bn_origin_slice = Bn_origin[:,index,:]
        Bn_cut_slice = Bn_cut[:,index,:]
        Bn_result_slice = Bn_result[:,index,:]
        hh_origin, vv_origin = xx_origin[:,index,:], zz_origin[:,index,:]
        hh_cut, vv_cut = xx_cut[:,index,:], zz_cut[:,index,:]
        hh_result, vv_result = xx_result[:,index,:], zz_result[:,index,:]
    elif normal == 'z':
        Bn_origin_slice = Bn_origin[:,:,index]
        Bn_cut_slice = Bn_cut[:,:,index]
        Bn_result_slice = Bn_result[:,:,index]
        hh_origin, vv_origin = xx_origin[:,:,index], yy_origin[:,:,index]
        hh_cut, vv_cut = xx_cut[:,:,index], yy_cut[:,:,index]
        hh_result, vv_result = xx_result[:,:,index], yy_result[:,:,index]
    
    max_abs_origin = np.nanmax(np.abs(Bn_origin_slice))
    max_abs_cut = np.nanmax(np.abs(Bn_cut_slice))
    max_abs_result = np.nanmax(np.abs(Bn_result_slice))
    
    if normal == 'x':
        Bn_origin_slice = Bn_origin[index,:,:]
        Bn_cut_slice = Bn_cut[index,:,:]
        Bn_result_slice = Bn_result[index,:,:]
        suptitle = f'Bx(G) at x = {xx_origin[index, 0, 0]/1e3:.2f} Mm (ix={index})'
        xlabel, ylabel = 'y (Mm)', 'z (Mm)'
    elif normal == 'y':
        Bn_origin_slice = Bn_origin[:,index,:]
        Bn_cut_slice = Bn_cut[:,index,:]
        Bn_result_slice = Bn_result[:,index,:]
        suptitle = f'By(G) at y = {yy_origin[0, index, 0]/1e3:.2f} Mm (iy={index})'
        xlabel, ylabel = 'x (Mm)', 'z (Mm)'
    elif normal == 'z':
        Bn_origin_slice = Bn_origin[:,:,index]
        Bn_cut_slice = Bn_cut[:,:,index]
        Bn_result_slice = Bn_result[:,:,index]
        suptitle = f'Bz(G) at z = {zz_origin[0, 0, index]/1e3:.2f} Mm (iz={index})'
        xlabel, ylabel = 'x (Mm)', 'y (Mm)'
    
    def plot_subfigure(ax, hh_grid, vv_grid, Bn_slice, max_abs_Bn, title):
        c = ax.pcolormesh(hh_grid / 1e3, vv_grid / 1e3, Bn_slice, cmap='bwr', vmin=-max_abs_Bn, vmax=max_abs_Bn)
        plt.colorbar(c, ax=ax, orientation='horizontal')
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_aspect('equal', adjustable='box')
        ax.grid()
        return
    
    fig, axs = plt.subplots(1, 3, figsize=(15, 4))
    plot_subfigure(axs[0], hh_origin, vv_origin, Bn_origin_slice, max_abs_origin, 'origin')
    plot_subfigure(axs[1], hh_cut, vv_cut, Bn_cut_slice, max_abs_cut, 'before cleaning')
    plot_subfigure(axs[2], hh_result, vv_result, Bn_result_slice, max_abs_result, 'after cleaning')
    plt.suptitle(suptitle)
    
    return

# Plotting slices
z_indices = [-90,-85,-80,-75,-70]
for iz in z_indices:
    plot_slice_field(Bz_origin, Bz_cut, Bz_result, 'z', iz)
    plt.tight_layout()

plt.show()