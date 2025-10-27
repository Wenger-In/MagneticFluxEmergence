import numpy as np
import matplotlib.pyplot as plt

from utils_field_analyse import plot_distribution_hist2d
        
# 生成测试数据
def generate_test_data(n_samples=100000):
    """
    生成二维测试数据，包含密集区域和稀疏区域
    
    返回:
    x, y: 两个一维数组，模拟待分析的数据
    """
    # 设置随机种子，保证结果可复现
    np.random.seed(42)
    
    # 生成三个主要密集区域（高计数）
    # 区域1: 正态分布
    x1 = np.random.normal(loc=0, scale=1.0, size=int(n_samples*0.4))
    y1 = np.random.normal(loc=0, scale=1.0, size=int(n_samples*0.4))
    
    # 区域2: 另一个正态分布簇
    x2 = np.random.normal(loc=5, scale=0.8, size=int(n_samples*0.3))
    y2 = np.random.normal(loc=3, scale=0.8, size=int(n_samples*0.3))
    
    # 区域3: 均匀分布的密集区域
    x3 = np.random.uniform(low=-4, high=-2, size=int(n_samples*0.2))
    y3 = np.random.uniform(low=2, high=5, size=int(n_samples*0.2))
    
    # 生成稀疏分布的背景噪声（低计数）
    x_noise = np.random.uniform(low=-6, high=8, size=int(n_samples*0.1))
    y_noise = np.random.uniform(low=-3, high=7, size=int(n_samples*0.1))
    
    # 合并所有数据
    x = np.concatenate([x1, x2, x3, x_noise])
    y = np.concatenate([y1, y2, y3, y_noise])
    
    return x, y

# 测试绘图函数
if __name__ == "__main__":
    # 生成测试数据
    x_data, y_data = generate_test_data()
    
    # 创建图形和轴
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # 导入之前的绘图函数（确保在同一文件或已正确导入）
    # 这里假设plot_distribution_hist2d函数已定义
    
    # 绘制原始直方图（无掩蔽）
    plot_distribution_hist2d(
        fig, ax1, "Test Data", "X Value", "Y Value",
        x_data, y_data,
        bins=50, cmap='jet',
        uselog=True, mask_counts=None,  # 不使用掩蔽
        show_colorbar=True, show_grid=True
    )
    ax1.set_title("原始直方图（含低计数区域）")
    
    # 绘制带掩蔽的直方图（掩蔽计数<100的区域）
    plot_distribution_hist2d(
        fig, ax2, "Test Data", "X Value", "Y Value",
        x_data, y_data,
        bins=50, cmap='jet',
        uselog=True, mask_counts=100,  # 掩蔽低计数区域
        show_colorbar=True, show_grid=True
    )
    ax2.set_title("掩蔽低计数区域（counts < 100）")
    
    plt.tight_layout()
    plt.show()
