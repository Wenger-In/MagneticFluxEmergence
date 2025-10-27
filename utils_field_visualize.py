import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyvista as pv

def select_start_indices(F, xx, yy, zz, z_ind, thres, save_dir):
    """
    获取三维数组指定z剖面内数值超过阈值的点,并返回其坐标
    
    参数:
        F: 三维数组
        xx, yy, zz: 对应的网格坐标
        z_ind: z方向的剖面序号
        thres: 阈值
    """
    slice_data = F[:, :, z_ind]
    
    x_indices, y_indices = np.where(np.abs(slice_data) > thres)

    indices_data = []
    
    for x_ind, y_ind in zip(x_indices, y_indices):
        indices_data.append([x_ind, y_ind, z_ind, F[x_ind, y_ind, z_ind]])
    
    df = pd.DataFrame(indices_data, columns=['x_ind', 'y_ind', 'z_ind', 'field'])
    df.to_csv(save_dir+'start_indices.csv', index=False)

def load_start_points(xx, yy, zz, save_dir, save_name='start_indices.csv'):
    """
    从保存的索引CSV文件中加载点，并生成对应的坐标点列表
    
    参数:
        xx, yy, zz: 原始网格坐标（三维数组）
        csv_path: 索引CSV文件路径，默认为'start_indices.csv'
        
    返回:
        start_points: 包含坐标的numpy数组列表，与select_start_points返回格式一致
    """
    df = pd.read_csv(save_dir + save_name)
    start_points = []
    
    for _, row in df.iterrows():
        x_ind = int(row['x_ind'])
        y_ind = int(row['y_ind'])
        z_ind = int(row['z_ind'])
        
        x_coord = xx[x_ind, y_ind, z_ind]
        y_coord = yy[x_ind, y_ind, z_ind]
        z_coord = zz[x_ind, y_ind, z_ind]
        
        start_points.append(np.array([x_coord, y_coord, z_coord]))
    
    return start_points

def plot_streamlines(xx, yy, zz, Bx, By, Bz, label='|B|', scale_factor=1e-3, start_point_type=1, start_points=[], \
    Bc=1e3, thres_factor=20, clim_factor=1, plot_box=True, plotter_existed=False):
    '''
    Plot the streamlines of a 3D vector field by pyvista.
    
    INPUT:
        - xx, yy, zz: (m,n,k) array
            Three-dimensional grid
        - Bx, By, Bz: (m,n,k) array
            Three components of the vector magnetic field defined on a 3D grid
        - scale_factor: float
            Scaling factor of the vectors, default to 1e-3
        - start_point_type: integer
            0 for even distributed; 1 for magnitude in planes; 2 for magnitude in volume. Default to 1
        - Bc: float
            Field magnitude criterion for selecting start points      
    
    OUTPUT:
        - plotter: object
            Pyvista object used for plotting
    '''
    # Scaling the vector fields
    Bmod = np.sqrt(Bx**2 + By**2 + Bz**2)
    Bx = Bx * scale_factor
    By = By * scale_factor
    Bz = Bz * scale_factor
    m, n, k = Bmod.shape
    
    B_mod = np.sqrt(Bx**2 + By**2 + Bz**2)
    thres = np.max(B_mod) / thres_factor
    mask = B_mod < thres
    Bx[mask] = 0
    By[mask] = 0
    Bz[mask] = 0
    
    # Creating pyvista grid
    grid = pv.StructuredGrid(xx, yy, zz)
    grid['Bx'] = Bx.ravel(order='F')
    grid['By'] = By.ravel(order='F')
    grid['Bz'] = Bz.ravel(order='F')
    grid['|B|'] = np.stack((Bx, By, Bz), axis=-1).reshape(-1, 3, order='F')

    # Selecting start points
    xx_range = np.max(xx) - np.min(xx)
    yy_range = np.max(yy) - np.min(yy)
    zz_range = np.max(zz) - np.min(zz)
    if start_point_type == 0: # even distributed
        start_x = np.linspace(np.min(xx)+xx_range/10, np.max(xx)-xx_range/10, 12)
        start_y = np.linspace(np.min(yy)+yy_range/10, np.max(yy)-yy_range/10, 12)
        # start_z = [np.min(zz)+zz_range/10, np.max(zz)-zz_range/10]
        start_z = [np.min(zz)+zz_range/10]
        start_points = np.array([[x, y, z] for x in start_x for y in start_y for z in start_z])
    elif start_point_type == 1: # by field magnitude in planes
        start_points = []
        iz_top, iz_bot = 5, -5
        for ix in range(0, m, 2):
            for iy in range(0, n, 2):
                if Bmod[ix, iy, iz_top] > Bc:
                    start_points.append(np.array([xx[ix, iy, iz_top], yy[ix, iy, iz_top], zz[ix, iy, iz_top]]))
                if Bmod[ix, iy, iz_bot] > Bc:
                    start_points.append(np.array([xx[ix, iy, iz_bot], yy[ix, iy, iz_bot], zz[ix, iy, iz_bot]]))
        start_points = np.array(start_points)
    elif start_point_type == 2: # by field magnitude in volume
        start_points = []
        for ix in range(0, m, 2):
            for iy in range(0, n, 2):
                for iz in range(0, k, 2):
                    if Bmod[ix, iy, iz] > Bc:
                        start_points.append(np.array([xx[ix, iy, iz], yy[ix, iy, iz], zz[ix, iy, iz]]))
        start_points = np.array(start_points)
    elif start_point_type == 3:
        start_points = np.array(start_points)

    # Plotting the streamlines
    if plotter_existed:
        plotter = plotter_existed
    else:
        plotter = pv.Plotter()
    edge_thres = 1e3
    for point in start_points:
        try:
            streamline = grid.streamlines(start_position=point, integration_direction='both', max_steps=20000000)
            if streamline.points.size > 0:
                end_point = streamline.points[-1]
                if (np.abs(end_point[0] - np.min(xx)) < edge_thres or np.abs(end_point[0] - np.max(xx)) < edge_thres or \
                    np.abs(end_point[1] - np.min(yy)) < edge_thres or np.abs(end_point[1] - np.max(yy)) < edge_thres or \
                    np.abs(end_point[2] - np.min(zz)) < edge_thres or np.abs(end_point[2] - np.max(zz)) < edge_thres):
                    if label == '|B|':
                        plotter.add_mesh(streamline, line_width=2, scalars='|B|', cmap='jet', 
                                        clim=[-np.max(np.abs(Bz))/clim_factor, np.max(np.abs(Bz))/clim_factor], \
                                        scalar_bar_args={'position_x': 0, 'position_y': 0.9, 'width': 0.5, 'height': 0.08})
                    elif label == 'Bz': 
                        plotter.add_mesh(streamline, line_width=2, scalars='Bz', cmap='PiYG', 
                                        clim=[-np.max(np.abs(Bz))/clim_factor, np.max(np.abs(Bz))/clim_factor], \
                                        scalar_bar_args={'position_x': 0, 'position_y': 0.9, 'width': 0.5, 'height': 0.08})
                    # plotter.add_mesh(streamline.tube(radius=0.2), lighting=False)
        except ValueError as e:
            print(f"Encountered an error for start point {point}: {e}")
            continue 

    # Plotting outline box
    if plot_box:
        outline = grid.outline()
        plotter.add_mesh(outline, color='black', line_width=2)
    
    return plotter

def plot_arrows(xx, yy, zz, Bx, By, Bz, step=5, scale_factor=1, offset = 5000):
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
    Bx_step = Bx[::step, ::step, ::step] * scale_factor
    By_step = By[::step, ::step, ::step] * scale_factor
    Bz_step = Bz[::step, ::step, ::step] * scale_factor

    # Creating pyVista grids
    vectors = np.stack((Bx_step.ravel(), By_step.ravel(), Bz_step.ravel()), axis=1)
    points = np.stack((xx_step.ravel(), yy_step.ravel(), zz_step.ravel()), axis=1)

    # Creating arrows
    grid = pv.PolyData(points)
    grid['vectors'] = vectors
    Bz_step_offset = (Bz_step + offset) / np.sqrt(3)
    grid['Bz_value'] = np.stack((Bz_step_offset.ravel(), Bz_step_offset.ravel(), Bz_step_offset.ravel()), axis=1)
    arrows = grid.glyph(orient='vectors', scale=1, factor=8)
        
    # Plotting the arrows
    plotter = pv.Plotter()
    Bz_step_abs_max = np.nanmax(np.abs(Bz_step))
    plotter.add_mesh(arrows, scalars='Bz_value', cmap='PiYG', clim=[offset-Bz_step_abs_max, offset+Bz_step_abs_max], \
        scalar_bar_args={'title': 'Fz', 'position_x': 0, 'position_y': 0.9, 'width': 0.5, 'height': 0.08})
    
    # Plotting outline box
    outline = grid.outline()
    plotter.add_mesh(outline, color='black', line_width=2)
    
    return plotter

def plot_boundary_plane(xx, yy, zz, Bx, By, Bz, plotter, scale_factor=1e-3, which_boundary='both'):
    '''
    This function is used to plot the arrows of a 3D vector field by pyvista.
    
    INPUT:
        - xx, yy, zz: (m,n,k) array
            Three-dimensional grid
        - Bx, By, Bz: (m,n,k) array
            Three components of the vector magnetic field defined on a 3D grid
        - plotter: object
            Pyvista object used for plotting
        - scale_factor: float
            Scaling factor of the vectors, default to 1e-3
    '''
    # Scaling the vector fields
    Bx = Bx * scale_factor
    By = By * scale_factor
    Bz = Bz * scale_factor
    
    # Defining the bottom/top plane
    Bz_bot = Bz[:, :, -2].T
    Bz_top = Bz[:, :, 1].T
    Bz_top_abs_max = np.max(np.abs(Bz_top))
    Bz_bot_abs_max = np.max(np.abs(Bz_bot))
    
    # Creating pyvista grid for the bottom plane
    z_min = np.min(zz)
    plane_bot = pv.StructuredGrid(xx[:, :, -2], yy[:, :, -2], np.full_like(xx[:, :, -2], z_min))
    plane_bot['Bz (kG)'] = Bz_bot.flatten()
    
    # Creating pyvista grid for the top plane
    z_max = np.max(zz)
    plane_top = pv.StructuredGrid(xx[:, :, 1], yy[:, :, 1], np.full_like(xx[:, :, 1], z_max))
    plane_top['Bz (kG)'] = Bz_top.flatten()
    
    Bz_abs_clim = np.min([Bz_bot_abs_max, Bz_top_abs_max])
    
    plot_top = True
    plot_bottom = True
    if which_boundary == 'top':
        plot_bottom = False
    elif which_boundary == 'bottom':
        plot_top = False
    
    # Plotting the bottom/top plane
    if plot_top:
        plotter.add_mesh(plane_top, scalars='Bz (kG)', cmap='bwr', clim=[-Bz_abs_clim, Bz_abs_clim], opacity=0.9, \
            scalar_bar_args={'title': 'Bz (kG)', 'position_x': 0, 'position_y': 0.8, 'width': 0.5, 'height': 0.08})
    if plot_bottom:
        plotter.add_mesh(plane_bot, scalars='Bz (kG)', cmap='bwr', clim=[-Bz_abs_clim, Bz_abs_clim], opacity=0.9, \
            scalar_bar_args={'title': 'Bz (kG)', 'position_x': 0, 'position_y': 0.8, 'width': 0.5, 'height': 0.08})
    
    return
