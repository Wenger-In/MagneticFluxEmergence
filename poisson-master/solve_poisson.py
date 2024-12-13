from grids import Domain, Grid
from poisson import SimpleSolver, MultiGridSolver
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os

utils_path = os.path.dirname('E:/Research/Program/MagneticFluxEmergence/')
sys.path.append(utils_path)
from field_analyse_utils import calc_field_gradient, calc_field_divergence

# Reading data from saved files
Ve = 0.3
sparse_factor = 10
cal_done = 0
data_dir = 'E:/Research/Work/202405_solar_storm/reconstruction/sparse_factor=' + str(sparse_factor) + '/'
data_name = 'origin_field_data_' + str(Ve) + '.csv'
save_dir = 'E:/Research/Work/202405_solar_storm/reconstruction/sparse_factor=' + str(sparse_factor) + '/'
save_name = 'Poisson_solution_' + str(Ve) + '.csv'
data = pd.read_csv(data_dir+data_name)

# Extracting dataframes
xx = data['xx'].values
yy = data['yy'].values
zz = data['zz'].values
Bx = data['data_cube_Bx'].values
By = data['data_cube_By'].values
Bz = data['data_cube_Bz'].values
div_B = data['div_B'].values

# Reshaping saved data based on considering time 
if sparse_factor == 3:
    nx, ny, nz = 401, 171, 149
elif sparse_factor == 5: 
    nx, ny, nz = 241, 103, 90
elif sparse_factor == 10: 
    nx, ny, nz = 121, 52, 46
elif sparse_factor == 25:
    nx, ny, nz = 49, 21, 19

xx = xx.reshape((nx, ny, nz))
yy = yy.reshape((nx, ny, nz))
zz = zz.reshape((nx, ny, nz))
Bx = Bx.reshape((nx, ny, nz))
By = By.reshape((nx, ny, nz))
Bz = Bz.reshape((nx, ny, nz))
div_B = div_B.reshape((nx, ny, nz))

if sparse_factor == 3: 
    nx_cut, nz_cut = 171, 149
elif sparse_factor == 5: 
    nx_cut, nz_cut = 129, 89
elif sparse_factor == 10: 
    nx_cut, nz_cut = 52, 46
elif sparse_factor == 25:
    nx_cut, nz_cut = 21, 19

nx_beg = (nx-nx_cut) // 2
nx_end = nx_beg + nx_cut
nz_beg = (nz-nz_cut) // 2
nz_end = nz_beg + nz_cut

xx_cut = xx[nx_beg:nx_end, :, nz_beg:nz_end]
yy_cut = yy[nx_beg:nx_end, :, nz_beg:nz_end]
zz_cut = zz[nx_beg:nx_end, :, nz_beg:nz_end]
Bx_cut = Bx[nx_beg:nx_end, :, nz_beg:nz_end]
By_cut = By[nx_beg:nx_end, :, nz_beg:nz_end]
Bz_cut = Bz[nx_beg:nx_end, :, nz_beg:nz_end]
div_B_cut = div_B[nx_beg:nx_end, :, nz_beg:nz_end]

x_center = (np.max(xx_cut) - np.min(xx_cut))//2
y_center = (np.max(yy_cut) - np.min(yy_cut))//2
z_center = (np.min(zz_cut) - np.max(zz_cut))//2
x_edge = np.max(xx_cut) - np.min(xx_cut)
y_edge = np.max(yy_cut) - np.min(yy_cut)
z_edge = np.min(zz_cut) - np.max(zz_cut)

# Making the grid on which to solve the Poisson equation
domain = Domain(center=(x_center, y_center, z_center), edges=(x_edge, y_edge, z_edge))
grid = Grid(domain, shape=(nx_cut, ny, nz_cut))

# Defining the boundary condiction
def g(x, y, z):
    """ This function is used to produce the boundary conditions """
    return 0

# Preparing the boundary conditions
bc = {}
for index in grid.boundary:
    x, y, z = grid.loc(index)
    bc[index] = g(x, y, z)

rhs = grid.field_from_array(div_B_cut)

if cal_done == 0:
    # Solving Poisson Equation
    solver = MultiGridSolver(rhs, bc, atol=1.0E-11)
    solver.solve()
    u = solver.solution()
    psi = u.values
elif cal_done == 1:
    # Importing the calculated results
    result = pd.read_csv(save_dir+save_name)
    psi = result['psi'].values
    psi = psi.reshape((nx_cut, ny, nz_cut))

# Checking convergence
Bx_res, By_res, Bz_res = calc_field_gradient(psi, xx_cut, yy_cut, zz_cut)

Bx = Bx_cut - Bx_res
By = By_cut - By_res
Bz = Bz_cut - Bz_res
div_B_res = calc_field_divergence(Bx, By, Bz, xx_cut, yy_cut, zz_cut)



# Plotting examination
plt.figure()
plot_xgrid = xx_cut[:,:,10]/1e3
plot_ygrid = yy_cut[:,:,10]/1e3
xmin, xmax = np.min(plot_xgrid), np.max(plot_xgrid)
ymin, ymax = np.min(plot_ygrid), np.max(plot_ygrid)

plt.subplot(2,2,1)
plot_data1 = div_B_cut[:,:,10]
div_B_max = np.nanmax(np.abs(plot_data1))
plt.pcolormesh(plot_xgrid, plot_ygrid, plot_data1, cmap='RdBu', vmin=-div_B_max, vmax=div_B_max)
plt.title('divergence before')
plt.colorbar()
plt.xlabel('X [Mm]')
plt.ylabel('Y [Mm]')
plt.axis('equal')

plt.subplot(2,2,2)
plot_data2 = div_B_res[:,:,10]
div_B_max = np.nanmax(np.abs(plot_data2))
plt.pcolormesh(plot_xgrid, plot_ygrid, plot_data2, cmap='RdBu', vmin=-div_B_max, vmax=div_B_max)
plt.title('divergence now')
plt.colorbar()
plt.xlabel('X [Mm]')
plt.ylabel('Y [Mm]')
plt.axis('equal')

plt.subplot(2,2,3)
plot_data3 = plot_data1 - plot_data2
div_B_max = np.nanmax(np.abs(plot_data3))
plt.pcolormesh(plot_xgrid, plot_ygrid, plot_data3, cmap='RdBu', vmin=-div_B_max, vmax=div_B_max)
plt.title('before - now')
plt.colorbar()
plt.xlabel('X [Mm]')
plt.ylabel('Y [Mm]')
plt.axis('equal')

plt.show()

if cal_done == 0:
    # Saving interpolated field data
    xx_cut_flat, yy_cut_flat, zz_cut_flat = xx_cut.flatten(), yy_cut.flatten(), zz_cut.flatten()
    psi_flat = psi.flatten()

    # Organizing dataframes 
    data_save = pd.DataFrame({
        'xx_cut': xx_cut_flat,
        'yy_cut': yy_cut_flat,
        'zz_cut': zz_cut_flat,
        'psi': psi_flat
    })

    # Saving data to csv files
    data_save.to_csv(save_dir + save_name, index=False)

db