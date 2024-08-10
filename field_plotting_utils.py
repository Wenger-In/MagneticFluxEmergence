import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
import sunpy.map
from sunpy.map.header_helper import make_heliographic_header

def plot_field(REC_time, X, Y, Fx, Fy, Fz, label_xy, label_z, scale_xy=1, arrow_interval=5):
    '''
    This function is used to plot the magnetic field or the velocity field
    '''
    fig, ax = plt.subplots(figsize=(12, 6))

    # Fz as background
    abs_max = np.nanmax(np.abs(Fz))
    c = ax.pcolormesh(X, Y, Fz, shading='auto', cmap='RdBu', vmin=-abs_max, vmax=abs_max)
    fig.colorbar(c, ax=ax, label=label_z)

    # Fx & Fy as quivers
    X_arrow = X[::arrow_interval, ::arrow_interval]
    Y_arrow = Y[::arrow_interval, ::arrow_interval]
    Fx_arrow = Fx[::arrow_interval, ::arrow_interval] * scale_xy # shrink Fx & Fy by 1/scale_xy
    Fy_arrow = Fy[::arrow_interval, ::arrow_interval] * scale_xy
    q = ax.quiver(X_arrow, Y_arrow, Fx_arrow, Fy_arrow, angles='xy', scale_units='xy', scale=1, color='black')
    qk = ax.quiverkey(q, 0.8, 0.08, 1, label_xy, labelpos='E', coordinates='figure')

    # axis settings
    ax.set_aspect('equal')
    ax.set_xlabel('Pixel')
    ax.set_ylabel('Pixel')
    ax.set_title('Photosperic field at ' + REC_time)
    plt.grid()
    plt.show()
