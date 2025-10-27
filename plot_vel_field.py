import numpy as np
import matplotlib.pyplot as plt
import os
from utils_field_analyse import plot_2d_vector_field_quiver, plot_scalar_field_pcolor, plot_scalar_field_contour

def loadField(target, fn):
    file_path = os.path.join(target, fn)
    fields = np.load(file_path)
    
    Bx = fields['Bx']
    By = fields['By']
    Bz = fields['Bz']
    Vx = fields['Vx']
    Vy = fields['Vy']
    Vz = fields['Vz']
    
    fields.close()
    
    return Bx, By, Bz, Vx, Vy, Vz

def plot_vel_histogram(Vx, Vy, Vz, OBS_time_cleaned):
    plt.figure()
    plt.hist(Vx.flatten(), bins=100, range=(-2, 2), fill=False, edgecolor='r', label='Vx')
    plt.hist(Vy.flatten(), bins=100, range=(-2, 2), fill=False, edgecolor='g', label='Vy')
    plt.hist(Vz.flatten(), bins=100, range=(-2, 2), fill=False, edgecolor='b', label='Vz')
    plt.xlabel('plasma velocity esimated by DAVE4VM (km/s)')
    plt.ylabel('counts')
    plt.yscale('log')
    plt.legend()
    plt.title(f'plasma velocity at {OBS_time_cleaned}')
    plt.show()

data_dir = 'E:/Research/Work/202405_solar_storm/HMI_field/FieldData/'
save_dir = 'E:/Research/Work/202405_solar_storm/HMI_field/'
save_or_not = 0

for fn in os.listdir(data_dir):
    if fn.endswith('.npz'):
        Bx, By, Bz, Vx, Vy, Vz = loadField(data_dir, fn)
        OBS_time_cleaned = fn[7:-5]
        
        # Getting the Dimension and Reference Pixel of the Map
        AR_shape = Bx.shape
        AR_center = [(AR_shape[1]-1)/2, (AR_shape[0]-1)/2]

        # Constructing the grid
        x_indices, y_indices = np.arange(AR_shape[1]), np.arange(AR_shape[0])
        x_coords = (x_indices - AR_center[0]) * 1
        y_coords = (y_indices - AR_center[1]) * 1
        X, Y = np.meshgrid(x_coords, y_coords)
        
        # Plotting velocity histogram
        plot_vel_histogram(Vx, Vy, Vz, OBS_time_cleaned)
        if save_or_not == 1:
            plt.savefig(save_dir+'VzVh/'+'Vx-Vy-Vz-'+OBS_time_cleaned+'.png', format='png')
        
        # Plotting magnetic vector field
        fig, ax = plt.subplots(figsize=(10, 4))
        plot_scalar_field_pcolor(fig, ax, OBS_time_cleaned, 'Bz(G)',
                                X, Y, Bz, 
                                cmap='bwr', 
                                show_grid=True)
        plot_2d_vector_field_quiver(fig, ax, OBS_time_cleaned, 'Bh', 'G',
                                    X, Y, Bx, By, 
                                    arrow_interval=10, rule=1000, scale=50, 
                                    color='k', width=0.002, 
                                    key_x=0.8, key_y=1.1)
        if save_or_not == 1:
            plt.savefig(save_dir+'BzBh/'+'Bz-Bh-'+OBS_time_cleaned+'.png', format='png')

        # Plotting velocity vector field
        fig, ax = plt.subplots(figsize=(10, 4))
        plot_scalar_field_pcolor(fig, ax, OBS_time_cleaned, 'Vz(km/s)',
                                X, Y, Vz, 
                                cmap='PiYG', 
                                show_grid=True)
        plot_2d_vector_field_quiver(fig, ax, OBS_time_cleaned, 'Vh', 'km/s',
                                    X, Y, Vx, Vy,  
                                    arrow_interval=10, rule=1, scale=0.03, 
                                    color='k', width=0.002, 
                                    key_x=0.8, key_y=1.05)
        if save_or_not == 1:
            plt.savefig(save_dir+'VzVh/'+'Vz-Vh-'+OBS_time_cleaned+'.png', format='png')

        # Plotting horizontal velocity field over vertical magnetic field (Vx, Vy, Bz)
        fig, ax = plt.subplots(figsize=(10, 4))
        plot_scalar_field_pcolor(fig, ax, OBS_time_cleaned, 'Bz(G)',
                                X, Y, Bz, 
                                cmap='bwr', 
                                show_grid=True)
        plot_2d_vector_field_quiver(fig, ax, OBS_time_cleaned, 'Vh', 'km/s',
                                    X, Y, Vx, Vy,  
                                    arrow_interval=10, rule=1, scale=0.03, 
                                    color='k', width=0.002, 
                                    key_x=0.8, key_y=1.05)
        if save_or_not == 1:
            plt.savefig(save_dir+'BzVh/'+'Bz-Vh-'+OBS_time_cleaned+'.png', format='png')

        plt.show()