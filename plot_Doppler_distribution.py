import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import sunpy.map
import os
import re

NOAA_num = 3664
SHARP_num = 11149

# 定义文件夹路径
Br_dir = f'E:/Research/Data/SDO/HMI/SHARP/AR{NOAA_num}/'
Dop_dir = 'E:/Research/Data/SDO/HMI/SHARP/Dopplergram/'
# 图片保存路径
save_dir = 'E:/Research/Work/Flux_emergence_of_bipolars_by_SDO/HMI_field/Doppler/hist/'

# 创建保存文件夹（若不存在）
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

# 匹配文件名格式的正则表达式
pattern = re.compile(
    rf'hmi\.sharp_cea_720s\.{SHARP_num}\.(\d{{8}})_(\d{{6}})_TAI\.Dopplergram\.fits'
)

# 获取Doppler文件夹下所有fits文件
dop_files = [f for f in os.listdir(Dop_dir) if f.endswith('.fits')]

# 遍历所有Doppler文件
for dop_fn in dop_files:
    # 匹配文件名提取日期和时间
    match = pattern.match(dop_fn)
    if not match:
        continue
    
    REC_date = match.group(1)
    REC_time = match.group(2)
    
    br_fn = f'hmi.sharp_cea_720s.{SHARP_num}.{REC_date}_{REC_time}_TAI.Br.fits'
    
    if not os.path.exists(os.path.join(Br_dir, br_fn)):
        print(f"{br_fn} is not existing")
        continue
    
    # 读取Doppler数据
    try:
        Dop_map = sunpy.map.Map(os.path.join(Dop_dir, dop_fn))
    except Exception as e:
        print(f"读取文件 {dop_fn} 出错: {e}，已跳过")
        continue
    
    # 提取Doppler速度数据
    Dop = Dop_map.data
    Dop = Dop - np.nanmean(Dop[200:250,:50])
    Vlos = Dop
    Vlos_flat = Vlos.flatten()
    
    # 绘制直方图
    plt.figure(figsize=(10, 6), facecolor='white')
    
    n, bins, patches = plt.hist(
        Vlos_flat,
        bins=48,
        range=[-1200, 1200],
        density=True,
        edgecolor='white',
        linewidth=0.8,
        alpha=0.9
    )
    
    # 计算并绘制平均值线
    mean_vlos = np.mean(Vlos_flat)
    plt.axvline(mean_vlos, color='k', linestyle='--', linewidth=2, 
                label=f'Avg. {mean_vlos:.2f} m/s')
    
    # 设置图表属性
    plt.title(f'{REC_date}_{REC_time}', fontsize=14)
    plt.xlabel(r'$V_{los}$ (m/s)', fontsize=12)
    plt.ylabel('Probability Density', fontsize=12)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.grid()
    plt.legend(fontsize=10, loc='upper right')
    plt.xlim([-5000, 5000])
    plt.tight_layout()
    
    # 构建保存路径和文件名
    save_filename = f'{REC_date}_{REC_time}_Doppler_hist.png'
    save_path = os.path.join(save_dir, save_filename)
    
    # 保存图片
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"直方图已保存至: {save_path}")
    
    plt.close()  # 关闭图像释放内存