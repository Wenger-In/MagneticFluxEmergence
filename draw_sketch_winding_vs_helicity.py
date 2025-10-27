from utils_sketch_draw import *

save_or_not = 1
wind_flag = 1 # [0] weak winding;  [1] strong winding
heli_flag = 0 # [0] weak helicity; [1] strong helicity

####################################################################################################
# Setting axis parameters
axis_center = (0, 0, 0)
axis_radius = 14

axis_angle_range = (-np.pi/5, np.pi)
axis_num_points = 500
axis_color = 'y'
axis_opacity = 0.5
axis_point_size = 10

# Setting horizontal plane parameters
plane_center = (0, 0, 0)
plane_length = 60
plane_width = 50
plane_color = 'gray'
plane_opacity = 1

# Setting helix parameters
helix_radius = 3

if wind_flag == 0:
    helix_turns = 3
elif wind_flag == 1:
    helix_turns = 5

if heli_flag == 0:
    helix_line_width = 6
elif heli_flag == 1:
    helix_line_width = 15

helix1_color = 'r'
helix2_color = 'g'
helix3_color = 'b'
helix1_final_angle = np.pi / 3
helix2_final_angle = np.pi
helix3_final_angle = -np.pi / 3

# Setting torus parameters
torus_radius = 2.95
torus_color = 'c'
torus_opacity = 0.5

####################################################################################################
# Creating the axis of flux tube
axis, axis_points = create_axis(
    center=axis_center,
    radius=axis_radius,
    angle_range=axis_angle_range,
    num_points=axis_num_points,
    point_size=axis_point_size
)

# Creating the horizontal plane
plane = create_horizontal_plane(
    center=plane_center,
    length=plane_length,
    width=plane_width
)

# Creating the tube
helix1 = create_helix(
    axis_points=axis_points,
    helix_radius=helix_radius,
    turns=helix_turns,
    final_angle=helix1_final_angle
)
helix2 = create_helix(
    axis_points=axis_points,
    helix_radius=helix_radius,
    turns=helix_turns,
    final_angle=helix2_final_angle
)
helix3 = create_helix(
    axis_points=axis_points,
    helix_radius=helix_radius,
    turns=helix_turns,
    final_angle=helix3_final_angle
)
torus = create_torus(
    axis_points=axis_points,
    torus_radius=torus_radius
)

# Plotting
pl = pv.Plotter(off_screen=True)
pl.set_background('white')

pl.add_mesh(plane, color=plane_color, opacity=plane_opacity, label='solar surface')
# pl.add_mesh(axis, color=axis_color, opacity=axis_opacity, point_size=axis_point_size, label='axis')

pl.add_mesh(torus, color=torus_color, opacity=torus_opacity, label='flux tube')
pl.add_mesh(helix1, color=helix1_color, line_width=helix_line_width, label='field line 1')
pl.add_mesh(helix2, color=helix2_color, line_width=helix_line_width, label='field line 2')
pl.add_mesh(helix3, color=helix3_color, line_width=helix_line_width, label='field line 3')

pl.camera.position = (0, 60, 24)
pl.camera.focal_point = (0, 0, 10)

pl.add_legend()
pl.legend.SetBackgroundColor(0.75, 0.75, 0.75)

# pl.add_axes()

####################################################################################################
# Saving figures
save_dir = 'E:/Research/Work/202405_solar_storm/winding/sketch/winding_vs_helicity/'

if save_or_not == 1:
    if wind_flag == 0:
        save_name1 = 'windweak'
    elif wind_flag == 1:
        save_name1 = 'windstro'
    if heli_flag == 0:
        save_name2 = 'heliweak'
    elif heli_flag == 1:
        save_name2 = 'helistro'
    pl.screenshot(save_dir + save_name1 + '_' +  save_name2 + '.png')
    pl.export_html(save_dir + save_name1 + '_' + save_name2 + '.html')