import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
from scipy import stats

# 全局启用LaTeX解析，确保\phi等LaTeX语法生效
plt.rcParams['text.usetex'] = True  # 核心配置1：允许解析LaTeX
plt.rcParams['font.family'] = 'serif'

df = pd.read_csv('phi_123.csv')

Phi_sub1 = df['Phi_sub1'].to_numpy()
Phi_sub2 = df['Phi_sub2'].to_numpy()
Phi_sub3 = df['Phi_sub3'].to_numpy()

start_index = 600
end_index = -1
Phi_sub1_slice = Phi_sub1[start_index:end_index]
Phi_sub2_slice = Phi_sub2[start_index:end_index]
Phi_sub3_slice = Phi_sub3[start_index:end_index]

# 提取两组投影点（仅索引600之后的部分）
# 1. X-Y投影点 (x=Phi_sub1, y=Phi_sub2)
xy_x = Phi_sub1_slice
xy_y = Phi_sub2_slice

# 2. Y-Z投影点 (x=Phi_sub2, y=Phi_sub3)
yz_x = Phi_sub2_slice
yz_y = Phi_sub3_slice

# 对X-Y投影点进行线性拟合
xy_slope, xy_intercept, xy_r_value, xy_p_value, xy_std_err = stats.linregress(xy_x, xy_y)

# 对Y-Z投影点进行线性拟合
yz_slope, yz_intercept, yz_r_value, yz_p_value, yz_std_err = stats.linregress(yz_x, yz_y)

# 输出拟合结果
print("X-Y投影点（索引600之后）线性拟合结果：")
print(f"斜率: {xy_slope:.6f}")
print(f"截距: {xy_intercept:.6f}")
print(f"相关系数: {xy_r_value:.6f}")
print(f"相关系数平方（决定系数）: {xy_r_value**2:.6f}\n")

print("Y-Z投影点（索引600之后）线性拟合结果：")
print(f"斜率: {yz_slope:.6f}")
print(f"截距: {yz_intercept:.6f}")
print(f"相关系数: {yz_r_value:.6f}")
print(f"相关系数平方（决定系数）: {yz_r_value**2:.6f}")

points = np.column_stack((Phi_sub1, Phi_sub2, Phi_sub3))
cloud = pv.PolyData(points)
cloud['index'] = np.arange(len(points)) / 5  # 添加索引用于着色

# 创建投影点
xy_proj_points = np.column_stack((Phi_sub1, Phi_sub2, np.zeros_like(Phi_sub3)))
yz_proj_points = np.column_stack((np.zeros_like(Phi_sub1), Phi_sub2, Phi_sub3))

xy_proj = pv.PolyData(xy_proj_points)
yz_proj = pv.PolyData(yz_proj_points)
xy_proj['index'] = np.arange(len(points))
yz_proj['index'] = np.arange(len(points))

# 创建绘图器，不使用bgcolor参数
plotter = pv.Plotter(window_size=[1200, 1000], border=True)
# 使用set_background方法设置白色背景（兼容旧版本PyVista）
plotter.set_background('white')

# 绘制三维点，设置colorbar为竖直并放在左侧
plotter.add_mesh(cloud, point_size=5, cmap='jet', 
                 show_scalar_bar=True, 
                 scalar_bar_args={
                    #  'title': 'Hours from 2024-05-02 12:00:00',
                     'title': '',
                     'color': 'k', 
                     'vertical': True,  # 竖直colorbar
                     'position_x': 0.02,  # 左侧位置
                     'position_y': 0.1,   # 垂直方向位置
                     'width': 0.06,       # colorbar宽度
                     'height': 0.8        # colorbar高度（占窗口高度80%）
                 })

# -------------------------- 投影点按索引分色 --------------------------
index_mask = np.arange(len(points)) <= start_index
# X-Y投影：600及之前→灰色，之后→黑色
xy_proj_gray = pv.PolyData(xy_proj_points[index_mask])
xy_proj_black = pv.PolyData(xy_proj_points[~index_mask])
plotter.add_mesh(xy_proj_gray, point_size=5, color='gray', opacity=0.7)
plotter.add_mesh(xy_proj_black, point_size=5, color='k', opacity=0.7)

# Y-Z投影：600及之前→灰色，之后→黑色
yz_proj_gray = pv.PolyData(yz_proj_points[index_mask])
yz_proj_black = pv.PolyData(yz_proj_points[~index_mask])
plotter.add_mesh(yz_proj_gray, point_size=5, color='gray', opacity=0.7)
plotter.add_mesh(yz_proj_black, point_size=5, color='k', opacity=0.7)
# ----------------------------------------------------------------------

# 确定坐标轴范围（取最大值的1.2倍作为范围）
max_x = np.max(Phi_sub1) * 1.2
max_y = np.max(Phi_sub2) * 1.2
max_z = np.max(Phi_sub3) * 1.2

# -------------------------- 坐标轴刻度 + label --------------------------
# 生成刻度位置（0到最大值，分5段）
x_ticks = np.array([0, 1e9])
y_ticks = np.array([2e9, 4e9])
z_ticks = np.array([2e9, 4e9])

# 绘制X轴刻度（红色，与坐标轴同色）
for tick in x_ticks:
    tick_line = pv.Line(pointa=(tick, 0, 0), pointb=(tick, max_y*0.05, 0))
    plotter.add_mesh(tick_line, color='red', line_width=3)
    # plotter.add_point_labels([(tick, -max_y*0.05, 0)], [f'{tick:.1e}'],
    #                          text_color='red', font_size=15)
# X轴label
# plotter.add_point_labels([(max_x*0.9, -max_y*0.12, 0)], [r'$\Phi_A$'], 
#                          text_color='red', font_size=20)

# 绘制Y轴刻度（绿色，与坐标轴同色）
for tick in y_ticks:
    tick_line = pv.Line(pointa=(0, tick, 0), pointb=(0, tick, max_z*0.05))
    plotter.add_mesh(tick_line, color='green', line_width=3)
    # plotter.add_point_labels([(-max_x*0.05, tick, 0)], [f'{tick:.1e}'],
    #                          text_color='green', font_size=15)
# Y轴label
# plotter.add_point_labels([(-max_x*0.12, max_y*0.9, 0)], [r'$\Phi_B$'], 
#                          text_color='green', font_size=20)

# 绘制Z轴刻度（蓝色，与坐标轴同色）
for tick in z_ticks:
    tick_line = pv.Line(pointa=(0, 0, tick), pointb=(0, max_y*0.05, tick))
    plotter.add_mesh(tick_line, color='blue', line_width=3)
#     plotter.add_point_labels([(-max_x*0.05, 0, tick)], [f'{tick:.1e}'],
#                              text_color='blue', font_size=15)
# Z轴label
# plotter.add_point_labels([(-max_x*0.12, 0, max_z*0.9)], [r'$\Phi_C$'], 
#                          text_color='blue', font_size=20)
# ----------------------------------------------------------------------

# 绘制仅含正方向的坐标轴
x_axis = pv.Line(pointa=(0, 0, 0), pointb=(max_x, 0, 0))
y_axis = pv.Line(pointa=(0, 0, 0), pointb=(0, max_y, 0))
z_axis = pv.Line(pointa=(0, 0, 0), pointb=(0, 0, max_z))

plotter.add_mesh(x_axis, color='red', line_width=5)
plotter.add_mesh(y_axis, color='green', line_width=5)
plotter.add_mesh(z_axis, color='blue', line_width=5)

# -------------------------- 拟合直线生成+末端标注 --------------------------
def create_reg_line_with_label(slope, intercept, x_min, x_max, fixed_axis, fixed_val, line_color, label_color, r_value):
    """生成回归直线，并在末端分两行标注斜率和相关系数"""
    x1, x2 = x_min, x_max
    y1, y2 = slope * x1 + intercept, slope * x2 + intercept
    
    if fixed_axis == 'z':
        line = pv.Line(pointa=(x1, y1, fixed_val), pointb=(x2, y2, fixed_val))
        label_pos = (x2 + (x_max-x_min)*0.05, y2 + (max_y)*0.03, fixed_val)
    else:
        line = pv.Line(pointa=(fixed_val, x1, y1), pointb=(fixed_val, x2, y2))
        label_pos = (fixed_val + (max_x)*0.03, x2, y2 + (max_z)*0.05)
    
    label_text = [f'slope: {slope:.4f}\nr: {np.abs(r_value):.4f}']
    # plotter.add_point_labels([label_pos], label_text, 
    #                          text_color=label_color, font_size=16, bold=True)
    
    return line

# X-Y平面回归直线
xy_reg_line = create_reg_line_with_label(xy_slope, xy_intercept, np.min(xy_x), np.max(xy_x), 
                                         'z', 0, line_color='red', label_color='red', r_value=xy_r_value)
plotter.add_mesh(xy_reg_line, color='red', line_width=4)

# Y-Z平面回归直线
yz_reg_line = create_reg_line_with_label(yz_slope, yz_intercept, np.min(yz_x), np.max(yz_x), 
                                         'x', 0, line_color='blue', label_color='blue', r_value=yz_r_value)
plotter.add_mesh(yz_reg_line, color='blue', line_width=4)
# ----------------------------------------------------------------------

# 创建第一象限的X-Y平面（z=0）
xy_plane = pv.Plane(center=(max_x/2, max_y/2, 0), i_size=max_x, j_size=max_y)
plotter.add_mesh(xy_plane, color='green', opacity=0.1, show_edges=True, edge_color='darkgray', line_width=2)

# 创建第一象限的Y-Z平面（x=0）
yz_plane = pv.Plane(center=(max_z/2, max_y/2, 0), i_size=max_z, j_size=max_y)
yz_plane.rotate_y(-90, inplace=True)  # 旋转到Y-Z平面
plotter.add_mesh(yz_plane, color='green', opacity=0.1, show_edges=True, edge_color='darkgray', line_width=2)

# 设置相机位置和聚焦点
plotter.set_position((max_x*4.5, -max_y*0.5, max_z*1.5))
plotter.set_focus((max_x/2, max_y/2, max_z/2))
plotter.reset_camera_clipping_range()

# 显示交互窗口
plotter.show()
    