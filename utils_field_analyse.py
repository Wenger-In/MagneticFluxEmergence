import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LogNorm, SymLogNorm

def plot_2d_vector_field_quiver(fig, ax, REC_time, label, unit, 
                                X, Y, Fx, Fy, 
                                arrow_interval=5, rule=100, scale=1, 
                                color='k', width=0.01, headwidth=5, headlength=7, 
                                key_x=0.8, key_y=1.1, 
                                fontsize=10,
                                show_key=True, show_grid=False):
    '''
    This function is used to plot the 2D vector field as quivers.
    '''
    # Plot 2D vector field
    X_arrow = X[::arrow_interval, ::arrow_interval]
    Y_arrow = Y[::arrow_interval, ::arrow_interval]
    Fx_arrow = Fx[::arrow_interval, ::arrow_interval]
    Fy_arrow = Fy[::arrow_interval, ::arrow_interval] 
    q = ax.quiver(X_arrow, Y_arrow, Fx_arrow, Fy_arrow, 
                  angles='xy', scale_units='xy', scale=scale,
                  color=color, width=width, 
                  headwidth=headwidth, headlength=headlength)
    
    # Show quiver key
    if show_key:
        qk = ax.quiverkey(q, key_x, key_y, rule, label=f'{label}={rule}{unit}', labelpos='E')

    # Axis settings
    ax.set_aspect('equal')
    ax.set_xlabel('pixel', fontsize=fontsize)
    ax.set_ylabel('pixel', fontsize=fontsize)
    ax.set_title(f'{REC_time}', fontsize=fontsize)
    
    # Show grids
    if show_grid:
        ax.grid(alpha=0.3)

def plot_2d_vector_field_streamline(fig, ax, REC_time, label, 
                                    X, Y, Fx, Fy, 
                                    color, cmap,
                                    density=1, linewidth=1, fontsize=10,
                                    arrowsize=1.5, arrowstyle='->', 
                                    show_colorbar=True, show_grid=False):
    s = ax.streamplot(X, Y, Fx, Fy, 
                      density=density, linewidth=linewidth, 
                      color=color, cmap=cmap, 
                      arrowsize=arrowsize, arrowstyle=arrowstyle)

    # Show colorbar
    if show_colorbar:
        fig.colorbar(s.lines, ax=ax, label=label)
    
    # Axis settings
    ax.set_aspect('equal')
    ax.set_xlabel('pixel', fontsize=fontsize)
    ax.set_ylabel('pixel', fontsize=fontsize)
    ax.set_title(f'{REC_time}', fontsize=fontsize)
    
    # Show grids
    if show_grid:
        ax.grid(alpha=0.3)

def plot_scalar_field_pcolor(fig, ax, REC_time, label, 
                             X, Y, F, 
                             cmap, 
                             use_symlognorm=False, linthresh=10, linscale=1,
                             fontsize=10,
                             show_colorbar=True, show_grid=False):
    '''
    This function is used to plot the scalar field as 'pcolor'.
    '''
    # Plot scalar field
    abs_max = np.nanmax(np.abs(F))
    if use_symlognorm:
        norm = SymLogNorm(linthresh, linscale, vmin=-abs_max, vmax=abs_max)
    else:
        norm = Normalize(vmin=-abs_max, vmax=abs_max)
    c = ax.pcolormesh(X, Y, F, shading='auto', cmap=cmap, norm=norm)
    
    # Show colorbar
    if show_colorbar:
        fig.colorbar(c, ax=ax, label=label)
    
    # Axis settings
    ax.set_aspect('equal')
    ax.set_xlabel('pixel', fontsize=fontsize)
    ax.set_ylabel('pixel', fontsize=fontsize)
    ax.set_title(f'{label} at {REC_time}', fontsize=fontsize)
    
    # Show grids
    if show_grid:
        ax.grid(alpha=0.3)
    
def plot_scalar_field_contour(fig, ax, REC_time, label, 
                              X, Y, F, 
                              levels, cmap, 
                              linestyle='-', linewidth=2, fontsize=10,
                              show_colorbar=False, show_grid=False):
    '''
    This function is used to plot the scalar field as 'contour'.
    '''

    # Plot scalar field
    c = ax.contour(X, Y, F, 
                   levels=levels, 
                   cmap=cmap,
                   linestyles=linestyle,
                   linewidths=linewidth)
    
    # Set contour value
    ax.clabel(c, inline=True, fontsize=fontsize, fmt='%.1f')
    
    # Show colorbar
    if show_colorbar:
        fig.colorbar(c, ax=ax, label=label)
    
    # Axis settings
    ax.set_aspect('equal')
    ax.set_xlabel('pixel', fontsize=fontsize)
    ax.set_ylabel('pixel', fontsize=fontsize)
    ax.set_title(f'{label} at {REC_time}', fontsize=fontsize)
    
    # Show grids
    if show_grid:
        ax.grid(alpha=0.3)

def plot_distribution_hist2d(fig, ax, REC_time, xlabel, ylabel,
                             x, y, 
                             bins=50, cmap='jet', 
                             fontsize=10, 
                             uselog=True, mask_counts=None, 
                             show_colorbar=True, show_grid=True):
    if uselog:
        norm = LogNorm()
    else:
        norm = Normalize()
    
    # Plot distribution
    x = x.flatten()
    y = y.flatten()
    range = [[np.nanmin(x), np.nanmax(x)], [np.nanmin(y), np.nanmax(y)]]
    counts, xedges, yedges, ht = ax.hist2d(x, y, bins=bins, range=range, cmap=cmap, norm=norm)
    
    if mask_counts is not None:
        mask = counts < mask_counts
        cmap_obj = plt.get_cmap(cmap)
        cmap_obj.set_bad(alpha=0)
        counts_masked = np.ma.masked_where(mask, counts)
        
        ax.clear()
        ht = ax.pcolormesh(xedges, yedges, counts_masked.T,
                           cmap=cmap_obj, norm=norm, shading='flat')
    
    # Axis settings
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.set_title(f'{xlabel}-{ylabel} at {REC_time}', fontsize=fontsize)
    
    # Show colorbar and grids
    if show_colorbar:
        fig.colorbar(ht, ax=ax, label='counts')
    if show_grid:
        ax.grid(alpha=0.3)
    
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

def calc_field_divergence(Fx, Fy, Fz, 
                          xx, yy, zz):
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

def calc_field_curl(Fx, Fy, Fz, 
                    xx, yy, zz):
    '''
    This function is used to calculate the curl of a 3D vector field.

    INPUT:
    - Fx, Fy, Fz: (m,n,k) array
        Three components of the vector field defined on a 3D grid
    - xx, yy, zz: (m,n,k) array
        Three-dimensional grid

    OUTPUT:
    - curl_F_x, curl_F_y, curl_F_z: (m,n,k) array
        Curl of the field, with the last dimension representing the x, y, z components of the curl vector respectively
    '''
    m, n, k = Fx.shape

    Delta_x = xx[1, 0, 0] - xx[0, 0, 0]
    Delta_y = yy[0, 1, 0] - yy[0, 0, 0]
    Delta_z = zz[0, 0, 1] - zz[0, 0, 0]

    # Calculate the partial derivatives for the curl components
    dFz_dy = np.diff(Fz, axis=1) / Delta_y
    nan_y_1 = np.full((m, 1, k), np.nan)
    dFz_dy = np.concatenate((nan_y_1, dFz_dy), axis=1)
    dFy_dz = np.diff(Fy, axis=2) / Delta_z
    nan_z_1 = np.full((m, n, 1), np.nan)
    dFy_dz = np.concatenate((nan_z_1, dFy_dz), axis=2)

    dFx_dz = np.diff(Fx, axis=2) / Delta_z
    nan_z_2 = np.full((m, n, 1), np.nan)
    dFx_dz = np.concatenate((nan_z_2, dFx_dz), axis=2)
    dFz_dx = np.diff(Fz, axis=0) / Delta_x
    nan_x_1 = np.full((1, n, k), np.nan)
    dFz_dx = np.concatenate((nan_x_1, dFz_dx), axis=0)

    dFy_dx = np.diff(Fy, axis=0) / Delta_x
    nan_x_2 = np.full((1, n, k), np.nan)
    dFy_dx = np.concatenate((nan_x_2, dFy_dx), axis=0)
    dFx_dy = np.diff(Fx, axis=1) / Delta_y
    nan_y_2 = np.full((m, 1, k), np.nan)
    dFx_dy = np.concatenate((nan_y_2, dFx_dy), axis=1)

    # Calculate the curl components
    curl_F_x = dFz_dy - dFy_dz
    curl_F_y = dFx_dz - dFz_dx
    curl_F_z = dFy_dx - dFx_dy
    
    return curl_F_x, curl_F_y, curl_F_z

def calc_curl_B_based_on_SFD(Bx, By, Bz, 
                             xx, yy, zz):
    """
    This function is used to calculate the curl of magnetic fields using Staggered-grid Finite-Difference approach.

    INPUT:
    - Bx, By, Bz: (m,n,k) array
        Magnetic field components defined on regular grid points.
    - xx, yy, zz: (m,n,k) array
        Three-dimensional regular grid coordinates.

    OUTPUT:
    - curl_B_x, curl_B_y, curl_B_z: (m,n,k) array
        Curl of the magnetic field defined on regular grid points.
    """
    m, n, k = Bx.shape
    
    Delta_x = xx[1, 0, 0] - xx[0, 0, 0]
    Delta_y = yy[0, 1, 0] - yy[0, 0, 0]
    Delta_z = zz[0, 0, 1] - zz[0, 0, 0]

    # Calculate the components of magnetic field on staggered grid (ON AREA CENTER)
    Bx_sfd = (Bx[:, 1:, 1:] + Bx[:, 1:, :-1] + Bx[:, :-1, 1:] + Bx[:, :-1, :-1]) / 4 # (m,n-1,k-1)
    By_sfd = (By[1:, :, 1:] + By[1:, :, :-1] + By[:-1, :, 1:] + By[:-1, :, :-1]) / 4 # (m-1,n,k-1)
    Bz_sfd = (Bz[1:, 1:, :] + Bz[1:, :-1, :] + Bz[:-1, 1:, :] + Bz[:-1, :-1, :]) / 4 # (m-1,n-1,k)
    # Bx_sfd = Bx[:, :-1, :-1] # (m,n-1,k-1)
    # By_sfd = By[:-1, :, :-1] # (m-1,n,k-1)
    # Bz_sfd = Bz[:-1, :-1, :] # (m-1,n-1,k)
    
    # Calculate the components of curl B on staggered grid (ON EDGE CENTER)
    curl_B_x_sfd = (Bz_sfd[:, 1:, 1:-1] - Bz_sfd[:, :-1, 1:-1]) / Delta_y \
        - (By_sfd[:, 1:-1, 1:] - By_sfd[:, 1:-1, :-1]) / Delta_z  # (m-1,n-2,k-2)
    curl_B_y_sfd = (Bx_sfd[1:-1, :, 1:] - Bx_sfd[1:-1, :, :-1]) / Delta_z \
        - (Bz_sfd[1:, :, 1:-1] - Bz_sfd[:-1, :, 1:-1]) / Delta_x  # (m-2,n-1,k-2)
    curl_B_z_sfd = (By_sfd[1:, 1:-1, :] - By_sfd[:-1, 1:-1, :]) / Delta_x \
        - (Bx_sfd[1:-1, 1:, :] - Bx_sfd[1:-1, :-1, :]) / Delta_y  # (m-2,n-2,k-1)
    
    # Calculate the components of curl B on regular grid points
    curl_B_x, curl_B_y, curl_B_z = np.zeros_like(Bx), np.zeros_like(By), np.zeros_like(Bz)
    curl_B_x[1:-1, 1:-1, 1:-1] = (curl_B_x_sfd[1:, :, :] + curl_B_x_sfd[:-1, :, :]) / 2
    curl_B_y[1:-1, 1:-1, 1:-1] = (curl_B_y_sfd[:, 1:, :] + curl_B_y_sfd[:, :-1, :]) / 2
    curl_B_z[1:-1, 1:-1, 1:-1] = (curl_B_z_sfd[:, :, 1:] + curl_B_z_sfd[:, :, :-1]) / 2
    # curl_B_x[1:-1, 1:-1, 1:-1] = curl_B_x_sfd[1:, :, :]
    # curl_B_y[1:-1, 1:-1, 1:-1] = curl_B_y_sfd[:, 1:, :]
    # curl_B_z[1:-1, 1:-1, 1:-1] = curl_B_z_sfd[:, :, 1:]

    return curl_B_x, curl_B_y, curl_B_z
