import numpy as np
import matplotlib.pyplot as plt
import sunpy.map
import os
import re

NOAA_num = 3664
SHARP_num = 11149

Br_dir = f'E:/Research/Data/SDO/HMI/SHARP/AR{NOAA_num}/'
Dop_dir = 'E:/Research/Data/SDO/HMI/SHARP/Dopplergram/'
save_dir = 'E:/Research/Work/Flux_emergence_of_bipolars_by_SDO/HMI_field/Doppler/distrib/'

if not os.path.exists(save_dir):
    os.makedirs(save_dir)

# (hmi.sharp_cea_720s.SHARP_num.yyyymmdd_hhMMSS_TAI.xxx.fits)
pattern = re.compile(
    rf'hmi\.sharp_cea_720s\.{SHARP_num}\.(\d{{8}})_(\d{{6}})_TAI\.[A-Za-z]+\.fits'
)

dop_files = [f for f in os.listdir(Dop_dir) if f.endswith('.fits')]

for dop_fn in dop_files:
    match = pattern.match(dop_fn)
    if not match:
        continue
    
    REC_date = match.group(1)
    REC_time = match.group(2)
    
    br_fn = f'hmi.sharp_cea_720s.{SHARP_num}.{REC_date}_{REC_time}_TAI.Br.fits'
    
    if not os.path.exists(os.path.join(Br_dir, br_fn)):
        print(f"{br_fn} is not existing")
        continue
    
    Br_map = sunpy.map.Map(os.path.join(Br_dir, br_fn))
    
    Dop_map = sunpy.map.Map(os.path.join(Dop_dir, dop_fn))
    
    Br = Br_map.data
    Bz = Br
    Dop = Dop_map.data
    Dop = Dop - np.nanmean(Dop[200:250,:50])
    Vlos = Dop
    
    fig, ax = plt.subplots(figsize=(8, 3))
    
    im1 = ax.pcolormesh(Vlos, cmap='PiYG', vmin=-1200, vmax=1200)
    cbar1 = fig.colorbar(im1, ax=ax, label=r'$V_{LOS}$ (m/s)')
    
    rect_x = [300, 350, 350, 300, 300]
    rect_y = [0, 0, 50, 50, 0]
    ax.plot(rect_x, rect_y, 'k-', linewidth=2) 
    
    levels = [-1000, -500, 500, 1000]
    im2 = ax.contour(Br, levels=levels, cmap='bwr')
    cbar2 = fig.colorbar(im2, ax=ax, label='Bz (G)', extend='both')
    
    ax.axis('equal')
    plt.title(f'{REC_date}_{REC_time}', fontsize=14)
    plt.tight_layout()
    
    save_filename = f'{REC_date}_{REC_time}_Doppler.png'
    save_path = os.path.join(save_dir, save_filename)
    
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    plt.close()
    # plt.show()
    # db