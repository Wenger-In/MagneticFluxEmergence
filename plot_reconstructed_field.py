import pyvista as pv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from mpl_toolkits.mplot3d import Axes3D

# Reading data from saved files
Ve = 0.3
sparse_factor = 10
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
        nx, ny, nz = 171, 171, 149
    elif sparse_factor == 5:
        nx, ny, nz = 129, 103, 89
    elif sparse_factor == 10: 
        nx, ny, nz = 52, 52, 46
    elif sparse_factor == 25:
        nx, ny, nz = 21, 21, 19

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

def plot_streamlines(xx, yy, zz, Bx, By, Bz, scale_factor=1e-3):
    '''
    This function is used to plot the streamlines of a 3D vector field by pyvista.
    
    INPUT:
        - xx, yy, zz: (m,n,k) array
            Three-dimensional grid
        - Bx, By, Bz: (m,n,k) array
            Three components of the vector magnetic field defined on a 3D grid
        - scale_factor: float
            Scaling factor of the vectors, default to 1e-3
    
    OUTPUT:
        - plotter: object
            Pyvista object used for plotting
    '''
    # Scaling the vector fields
    Bx = Bx * scale_factor
    By = By * scale_factor
    Bz = Bz * scale_factor
    
    # Creating pyvista grid
    grid = pv.StructuredGrid(xx, yy, zz)
    grid['Bx'] = Bx.ravel(order='F')
    grid['By'] = By.ravel(order='F')
    grid['Bz'] = Bz.ravel(order='F')
    grid['vectors'] = np.stack((Bx, By, Bz), axis=-1).reshape(-1, 3, order='F')

    # Selecting start points
    xx_range = np.max(xx) - np.min(xx)
    yy_range = np.max(yy) - np.min(yy)
    zz_range = np.max(zz) - np.min(zz)
    start_x = np.linspace(np.min(xx)+xx_range/10, np.max(xx)-xx_range/10, 10)
    start_y = np.linspace(np.min(yy)+yy_range/10, np.max(yy)-yy_range/10, 5)
    start_z = [np.min(zz)+zz_range/10, np.max(zz)-zz_range/10]
    start_points = np.array([[x, y, z] for x in start_x for y in start_y for z in start_z])

    # Plotting the streamlines
    plotter = pv.Plotter()
    for point in start_points:
        streamline = grid.streamlines(start_position=point, integration_direction='both', max_steps=10000000)
        plotter.add_mesh(streamline, line_width=4)

    # Plotting outline box
    outline = grid.outline()
    plotter.add_mesh(outline, color='black', line_width=2)
    
    return plotter

def plot_arrows(xx, yy, zz, Bx, By, Bz, step=5, scale_factor=1):
    '''
    This function is used to plot the arrows of a 3D vector field by pyvista.
    
    INPUT:
        - xx, yy, zz: (m,n,k) array
            Three-dimensional grid
        - Bx, By, Bz: (m,n,k) array
            Three components of the vector magnetic field defined on a 3D grid
        - step: integer
            Sample step for plotting, default to 5
        - scale_factor: float
            Scaling factor of the arrows, default to 1
    
    OUTPUT:
        - plotter: object
            Pyvista object used for plotting
    '''
    # Resampling grid and data
    xx_step = xx[::step, ::step, ::step]
    yy_step = yy[::step, ::step, ::step]
    zz_step = zz[::step, ::step, ::step]
    Bx_step = Bx[::step, ::step, ::step]
    By_step = By[::step, ::step, ::step]
    Bz_step = Bz[::step, ::step, ::step]

    # Creating pyVista grids
    vectors = np.stack((Bx_step.ravel(), By_step.ravel(), Bz_step.ravel()), axis=1)
    points = np.stack((xx_step.ravel(), yy_step.ravel(), zz_step.ravel()), axis=1)

    # Creating arrows
    grid = pv.PolyData(points)
    grid['vectors'] = vectors
    arrows = grid.glyph(orient='vectors', scale=scale_factor, factor=8)

    # Plotting the arrows
    plotter = pv.Plotter()
    plotter.add_mesh(arrows)

    # Plotting outline box
    outline = grid.outline()
    plotter.add_mesh(outline, color='black', line_width=2)
    
    return plotter

def plot_boundary_plane(xx, yy, zz, Bx, By, Bz, plotter):
    '''
    This function is used to plot the arrows of a 3D vector field by pyvista.
    
    INPUT:
        - xx, yy, zz: (m,n,k) array
            Three-dimensional grid
        - Bx, By, Bz: (m,n,k) array
            Three components of the vector magnetic field defined on a 3D grid
        - plotter: object
            Pyvista object used for plotting
    '''
    # Defining the bottom/top plane
    Bz_bot = Bz[:, :, -2].T
    Bz_top = Bz[:, :, 1].T
    Bz_top_abs_max = np.max(np.abs(Bz_top))
    Bz_bot_abs_max = np.max(np.abs(Bz_bot))
    
    # Creating pyvista grid for the bottom plane
    z_min = np.min(zz)
    plane_bot = pv.StructuredGrid(xx[:, :, -2], yy[:, :, -2], np.full_like(xx[:, :, -2], z_min))
    plane_bot['Bz'] = Bz_bot.flatten()
    
    # Creating pyvista grid for the top plane
    z_max = np.max(zz)
    plane_top = pv.StructuredGrid(xx[:, :, 1], yy[:, :, 1], np.full_like(xx[:, :, 1], z_max))
    plane_top['Bz'] = Bz_top.flatten()
    
    # Plotting the bottom/top plane
    plotter.add_mesh(plane_bot, scalars='Bz', cmap='RdBu', clim=[-Bz_bot_abs_max, Bz_bot_abs_max], opacity=0.9)
    plotter.add_mesh(plane_top, scalars='Bz', cmap='RdBu', clim=[-Bz_top_abs_max, Bz_top_abs_max], opacity=0.9)

# plotter = plot_arrows(xx, yy, zz, Bx, By, Bz, step=3, scale_factor=1)

plotter = plot_streamlines(xx, yy, zz, Bx, By, Bz, scale_factor=1e-3)

plot_boundary_plane(xx, yy, zz, Bx, By, Bz, plotter)

plotter.add_axes()
plotter.show()

db