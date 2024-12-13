import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
from matplotlib.colors import SymLogNorm
import sunpy.map
from sunpy.map.header_helper import make_heliographic_header
from scipy.sparse import diags, csr_matrix
from scipy.sparse.linalg import LinearOperator, cg

def plot_vector_field(REC_time, X, Y, Fx, Fy, Fz, label_xy, label_z, rule=1, scale_xy=1, arrow_interval=5):
    '''
    This function is used to plot the magnetic field or the velocity field.
    '''
    fig, ax = plt.subplots(figsize=(12, 6))

    # Fz as background
    abs_max = np.nanmax(np.abs(Fz))
    c = ax.pcolormesh(X, Y, Fz, shading='auto', cmap='RdBu', vmin=-abs_max, vmax=abs_max)
    fig.colorbar(c, ax=ax, label=label_z)

    # Fx & Fy as quivers
    X_arrow = X[::arrow_interval, ::arrow_interval]
    Y_arrow = Y[::arrow_interval, ::arrow_interval]
    Fx_arrow = Fx[::arrow_interval, ::arrow_interval]
    Fy_arrow = Fy[::arrow_interval, ::arrow_interval] 
    q = ax.quiver(X_arrow, Y_arrow, Fx_arrow, Fy_arrow, angles='xy', scale_units='xy', scale=scale_xy, color='black')
    qk = ax.quiverkey(q, 0.8, 1.1, rule, label_xy + '=' + str(rule), labelpos='E')

    # axis settings
    ax.set_aspect('equal')
    ax.set_xlabel('Pixel')
    ax.set_ylabel('Pixel')
    ax.set_title('Photosperic field at ' + REC_time)
    plt.grid()

def plot_scalar_field(REC_time, X, Y, F, label, use_symlognorm=True, linthresh=10, linscale=1):
    '''
    This function is used to plot the scalar field.
    '''
    fig, ax = plt.subplots(figsize=(12, 6))

    # check scaling method and plot field
    abs_max = np.nanmax(np.abs(F))
    if use_symlognorm:
        norm = SymLogNorm(linthresh, linscale, vmin=-abs_max, vmax=abs_max)
        c = ax.pcolormesh(X, Y, F, shading='auto', cmap='RdBu', norm=norm)
    else:
        c = ax.pcolormesh(X, Y, F, shading='auto', cmap='RdBu', vmin=-abs_max, vmax=abs_max)
    
    # plot scalar field
    fig.colorbar(c, ax=ax, label=label)

    # axis settings
    ax.set_aspect('equal')
    ax.set_xlabel('Pixel')
    ax.set_ylabel('Pixel')
    ax.set_title(label + ' at ' + REC_time)
    plt.grid()

def calc_field_gradient(F, xx, yy, zz):
    '''
    This function is used to calculate the gradient of a 3D scalar field.
    
    INPUT:
    - F: (m,n,k) array
        Scalar field defined on a 3D grid
    - xx, yy, zz: (m,n,k) array
        Three-dimensional grid coordinates
    
    OUTPUT:
    - dF_dx, dF_dy, dF_dz : (m,n,k) or (m-2,n-2,k-2) array
        Gradient of the field, corresponding to three directions
    '''
    
    m, n, k = F.shape
    F = np.array(F)
    
    Delta_x = xx[1, 0, 0] - xx[0, 0, 0]
    Delta_y = yy[0, 1, 0] - yy[0, 0, 0]
    Delta_z = zz[0, 0, 1] - zz[0, 0, 0]
    
    dF_dx = np.diff(F, axis=0) / Delta_x
    nan_x = np.full((1, n, k), np.nan)
    dF_dx = np.concatenate((dF_dx, nan_x), axis=0)
    
    dF_dy = np.diff(F, axis=1) / Delta_y
    nan_y = np.full((m, 1, k), np.nan)
    dF_dy = np.concatenate((dF_dy, nan_y), axis=1)
    
    dF_dz = np.diff(F, axis=2) / np.abs(Delta_z)
    nan_z = np.full((m, n, 1), np.nan)
    dF_dz = np.concatenate((dF_dz, nan_z), axis=2)

    return dF_dx, dF_dy, dF_dz

def calc_field_divergence(Fx, Fy, Fz, xx, yy, zz):
    '''
    This function is used to calculate the divergence of a 3D vector field.
    
    INPUT:
    - Fx, Fy, Fz: (m,n,k) array
        Three components of the vector field defined on a 3D grid
    - xx, yy, zz: (m,n,k) array
        Three-dimensional grid
    
    OUTPUT:
    - div_F: (m,n,k) or (m-2,n-2,k-2) array
        Divergence of the field
    '''

    m, n, k = Fx.shape
    
    Delta_x = xx[1, 0, 0] - xx[0, 0, 0]
    Delta_y = yy[0, 1, 0] - yy[0, 0, 0]
    Delta_z = zz[0, 0, 1] - zz[0, 0, 0]
    
    dFx_dx = np.diff(Fx, axis=0) / Delta_x
    nan_x = np.full((1, n, k), np.nan)
    dFx_dx = np.concatenate((nan_x, dFx_dx), axis=0)
    
    dFy_dy = np.diff(Fy, axis=1) / Delta_y
    nan_y = np.full((m, 1, k), np.nan)
    dFy_dy = np.concatenate((nan_y, dFy_dy), axis=1)
    
    dFz_dz = np.diff(Fz, axis=2) / np.abs(Delta_z)
    nan_z = np.full((m, n, 1), np.nan)
    dFz_dz = np.concatenate((nan_z, dFz_dz), axis=2)
    
    # Calculating divergence
    div_F = dFx_dx + dFy_dy + dFz_dz

    return div_F

def Poisson_solver(source, dx, dy, dz):
    '''
    This function is used to solve 3D Poisson Equation based on Finite Difference method, 
        it is designed for the gird spacing uniform in X & Y axes, and nonuniform in Z axis. 
    The corresponding Poisson Equation is: {nabla}^2 {psi} = source

    INPUT:
    - source: (m,n,k) array
        Source term in Poisson Equ.
    - dx: float
        Uniform step spacing in X axis
    - dy: float
        Uniform step spacing in Y axis
    - dz: (k-1) array
        Non-uniform step spacing in Z axis

    OUTPUT:
    - psi: (m, n, k) array
        Solution of Poisson Equ.
    '''
    m, n, k = source.shape
    
    # # Constructing the finite difference matrix for Laplacian operator
    # # for X axis (uniform spacing)
    # D2x = diags([1, -2, 1], [-1, 0, 1], shape=(m, m)) / (dx ** 2)
    # # for Y axis (uniform spacing)
    # D2y = diags([1, -2, 1], [-1, 0, 1], shape=(n, n)) / (dy ** 2)
    # # for Z axis (non-uniform spacing)
    # D2z = np.zeros((k, k))
    # for i in range(1, k-1):
    #     D2z[i, i-1] = 1 / dz[i-1]**2  # forward difference
    #     D2z[i, i] = -2 / dz[i-1]**2 - 2 / dz[i]**2  # central difference
    #     D2z[i, i+1] = 1 / dz[i]**2  # backward difference
    
    # Constructing the finite difference matrix for Laplacian operator as sparse matrices
    # for X axis (uniform spacing)
    D2x = diags([1, -2, 1], [-1, 0, 1], shape=(m, m), format='csr') / (dx ** 2)
    
    # for Y axis (uniform spacing)
    D2y = diags([1, -2, 1], [-1, 0, 1], shape=(n, n), format='csr') / (dy ** 2)
    
    # for Z axis (non-uniform spacing)
    D2z = csr_matrix((k, k))
    diags_values = np.zeros((3, k))
    diags_values[0, 1:-1] = 1 / dz[:-1]**2  # forward difference
    diags_values[1, 1:-1] = -2 / dz[:-1]**2 - 2 / dz[1:]**2  # central difference
    diags_values[2, 1:-1] = 1 / dz[1:]**2  # backward difference
    D2z = diags(diags_values, offsets=[-1, 0, 1], shape=(k, k), format='csr')

    # Define the linear operator for the Laplacian
    def apply_laplacian(phi_flat):
        phi = phi_flat.reshape(m, n, k)
        phi_x = np.zeros_like(phi)
        phi_y = np.zeros_like(phi)
        phi_z = np.zeros_like(phi)

        # Applying D2x
        for j in range(n):
            phi_x[:, j, :] = D2x @ phi[:, j, :]
        # Applying D2y
        for i in range(m):
            phi_y[i, :, :] = D2y @ phi[i, :, :]
        # Applying D2z
        for j in range(n):
            for i in range(m):
                phi_z[i, j, :] = D2z @ phi[i, j, :]
        
        # # Applying Neumann boundary conditions
        # # for X direction
        # phi_x[0, :, :] = phi_x[1, :, :]
        # phi_x[-1, :, :] = phi_x[-2, :, :]
        # # for Y direction
        # phi_y[:, 0, :] = phi_y[:, 1, :]
        # phi_y[:, -1, :] = phi_y[:, -2, :]
        # # for Z direction
        # phi_z[:, :, 0] = phi_z[:, :, 1]
        # phi_z[:, :, -1] = phi_z[:, :, -2]
        
        return (phi_x + phi_y + phi_z).ravel()

    # Solving the linear system using Conjugate Gradient method
    source_flat = source.ravel()
    laplacian_operator = LinearOperator((m*n*k, m*n*k), apply_laplacian)
    psi_flat, info = cg(laplacian_operator, source_flat, tol=1e-10, maxiter=100*m*n*k)
    psi = psi_flat.reshape(m, n, k)
    
    return psi
