import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import convolve2d
import os
from utils_field_analyse import *

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

def smooth_with_convolution(A, kernel_size=3, mode='same'):
    kernel = np.ones((kernel_size, kernel_size)) / (kernel_size ** 2)
    A_smoothed = convolve2d(A, kernel, mode=mode, boundary='symm')
    
    return A_smoothed

def calc_curl_2d(Ax, Ay):
    '''
    Calculate the curl of 2D field
    
    Input:
    - Ax: x-component of the field 
    - Ay: y-component of the field
    
    Output:
    - curl_A: curl of the field
    '''
    dAydx = np.gradient(Ay, axis=1)
    dAxdy = np.gradient(Ax, axis=0)
    curl_A = dAydx - dAxdy
    
    return curl_A

data_dir = 'E:/Research/Work/202405_solar_storm/HMI_field/FieldData/'
save_dir = 'E:/Research/Work/202405_solar_storm/HMI_field/'
save_or_not = 0

for fn in os.listdir(data_dir):
    if fn.endswith('.npz'):
        Bx, By, Bz, Vx, Vy, Vz = loadField(data_dir, fn)
        OBS_time_cleaned = fn[9:-4]
        
        # Getting the Dimension and Reference Pixel of the Map
        AR_shape = Bx.shape
        AR_center = [(AR_shape[1]-1)/2, (AR_shape[0]-1)/2]

        # Constructing the grid
        x_indices, y_indices = np.arange(AR_shape[1]), np.arange(AR_shape[0])
        x_coords = (x_indices - AR_center[0]) * 1
        y_coords = (y_indices - AR_center[1]) * 1
        X, Y = np.meshgrid(x_coords, y_coords)
        
        # Calculating curl of fields
        kernel_size = 31
        U_sm = smooth_with_convolution(Vx, kernel_size)
        V_sm = smooth_with_convolution(Vy, kernel_size)
        bx_sm = smooth_with_convolution(Bx, kernel_size)
        by_sm = smooth_with_convolution(By, kernel_size)
        bz_sm = smooth_with_convolution(Bz, kernel_size)

        curl_Vh = calc_curl_2d(U_sm, V_sm)
        curl_Bh = calc_curl_2d(bx_sm, by_sm)
        alpha_t = curl_Vh * curl_Bh * bz_sm
        
        # Plotting the distribution
        fig, ax = plt.subplots(figsize=(10, 4))
        plot_2d_vector_field_quiver(fig, ax, OBS_time_cleaned, 'Bh', 'G',
                                    X, Y, Bx, By, 
                                    arrow_interval=10, rule=1000, scale=50, 
                                    color='r', width=0.002, 
                                    key_x=0.8, key_y=1.1)
        plot_2d_vector_field_quiver(fig, ax, OBS_time_cleaned, 'Vh', 'km/s', 
                                    X, Y, Vx, Vy,  
                                    arrow_interval=10, rule=1, scale=0.03, 
                                    color='b', width=0.002, 
                                    key_x=0.8, key_y=1.05)
        if save_or_not == 1:
            plt.savefig(save_dir+'winding/'+'curlV-curlB-'+OBS_time_cleaned+'.png', format='png')

        from utils_field_analyse import plot_2d_vector_field_streamline

        fig, ax = plt.subplots(figsize=(10, 4))
        plot_2d_vector_field_streamline(fig, ax, OBS_time_cleaned, 'Bh(G)', 
                                        X, Y, bx_sm, by_sm, 
                                        color=np.sqrt(bx_sm**2+by_sm**2), cmap='Reds')
        plot_2d_vector_field_streamline(fig, ax, OBS_time_cleaned, 'Vh(km/s)', 
                                        X, Y, U_sm, V_sm, 
                                        color=np.sqrt(Vx**2+Vy**2), cmap='Blues')
        if save_or_not == 1:
            plt.savefig(save_dir+'winding/'+'streamV-streamB-'+OBS_time_cleaned+'.png', format='png')

        fig, ax = plt.subplots(figsize=(10, 4))
        plot_scalar_field_pcolor(fig, ax, OBS_time_cleaned, '[curl(Vh) curl(Bh) Bz]', 
                                X, Y, alpha_t, 
                                cmap='PiYG', 
                                use_symlognorm=True)
        plot_scalar_field_contour(fig, ax, OBS_time_cleaned, 'Bz', 
                                X, Y, bz_sm, 
                                levels=[-1000, -500, 500, 1000], cmap='bwr', 
                                linestyle='-', linewidth=1.5)
        if save_or_not == 1:
            plt.savefig(save_dir+'winding/'+'curlV-curlB-Bz-'+OBS_time_cleaned+'.png', format='png')

        fig, ax = plt.subplots(figsize=(10, 4))
        plot_distribution_hist2d(fig, ax, OBS_time_cleaned, 'Bz(G)', '[curl(Vh) curl(Bh) Bz]',
                                bz_sm, alpha_t, mask_counts=1e1)
        if save_or_not == 1:
            plt.savefig(save_dir+'winding/'+'Bz-curlV-curlB-Bz-'+OBS_time_cleaned+'.png', format='png')
        
        plt.show()
        db