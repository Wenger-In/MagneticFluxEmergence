from utils_sketch_draw import *

save_or_not = 1
plot_scenario = 'em' # [ex] pre-existing; [em] emerging; [bd] braiding

####################################################################################################
# Setting axis parameters
if plot_scenario == 'em':
    axis_center = (2, 0, 0)
    axis_radius = 16
else:
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
helix_turns = 3
helix_line_width = 6

helix1_final_angle = np.pi / 3
helix2_final_angle = np.pi
helix3_final_angle = -np.pi / 3

# Setting torus parameters
torus_radius = 2.95

# Setting colors
helix1_color_ex = (0.5, 0, 0)
helix1_color_em = 'r'
helix2_color_ex = (0, 0.3, 0)
helix2_color_em = 'g'
helix3_color_ex = (0, 0, 0.5)
helix3_color_em = 'b'
torus_color_ex = 'c'
torus_color_em = 'y'
torus_opacity_ex = 0.5
torus_opacity_em = 0.2

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

# Defining axis segment for plot
if plot_scenario == 'ex':
    axis_ex = axis_points
    helix_turns_ex = 3
elif plot_scenario == 'em':
    axis_ex = axis_points[int(axis_num_points/4):]
    axis_em = axis_points[:int(axis_num_points/4)]
    helix_turns_ex = 3
    helix_turns_em = 1
    axis_center = (-2, 0, 0)
    axis_radius = 16
elif plot_scenario == 'bd':
    axis_ex = axis_points
    helix_turns_ex = 2.5

####################################################################################################
# Creating existing part of the tube
helix1_ex = create_helix(
    axis_points=axis_ex,
    helix_radius=helix_radius,
    turns=helix_turns_ex,
    final_angle=helix1_final_angle
)
helix2_ex = create_helix(
    axis_points=axis_ex,
    helix_radius=helix_radius,
    turns=helix_turns_ex,
    final_angle=helix2_final_angle
)
helix3_ex = create_helix(
    axis_points=axis_ex,
    helix_radius=helix_radius,
    turns=helix_turns_ex,
    final_angle=helix3_final_angle
)
torus_ex = create_torus(
    axis_points=axis_points,
    torus_radius=torus_radius
)

# Plotting
pl = pv.Plotter(off_screen=True)
pl.set_background('white')

pl.add_mesh(plane, color=plane_color, opacity=plane_opacity, label='solar surface')
# pl.add_mesh(axis, color=axis_color, opacity=axis_opacity, point_size=axis_point_size, label='axis')

pl.add_mesh(torus_ex, color=torus_color_ex, opacity=torus_opacity_ex, label='pre-existing flux tube')
pl.add_mesh(helix1_ex, color=helix1_color_ex, line_width=helix_line_width, label='field line 1')
pl.add_mesh(helix2_ex, color=helix2_color_ex, line_width=helix_line_width, label='field line 2')
pl.add_mesh(helix3_ex, color=helix3_color_ex, line_width=helix_line_width, label='field line 3')

####################################################################################################
# Creating emerging part of the tube
if plot_scenario == 'em':
    helix1_em = create_helix(
        axis_points=axis_em,
        helix_radius=helix_radius,
        turns=helix_turns_em,
        final_angle=helix1_final_angle
    )
    helix2_em = create_helix(
        axis_points=axis_em,
        helix_radius=helix_radius,
        turns=helix_turns_em,
        final_angle=helix2_final_angle
    )
    helix3_em = create_helix(
        axis_points=axis_em,
        helix_radius=helix_radius,
        turns=helix_turns_em,
        final_angle=helix3_final_angle
    )
    torus_em = create_torus(
        axis_points=axis_em,
        torus_radius=torus_radius
    )
    
    pl.add_mesh(torus_em, color=torus_color_em, opacity=torus_opacity_em, label='emerging flux tube')
    pl.add_mesh(helix1_em, color=helix1_color_em, line_width=helix_line_width, label='field line 1\'')
    pl.add_mesh(helix2_em, color=helix2_color_em, line_width=helix_line_width, label='field line 2\'')
    pl.add_mesh(helix3_em, color=helix3_color_em, line_width=helix_line_width, label='field line 3\'')

pl.camera.position = (0, 60, 24)
pl.camera.focal_point = (0, 0, 10)

pl.add_legend()
pl.legend.SetBackgroundColor(0.75, 0.75, 0.75)

# pl.add_axes()

####################################################################################################
# Saving figures
save_dir = 'E:/Research/Work/202405_solar_storm/winding/sketch/braid_vs_emerge/'

if save_or_not == 1:
    if plot_scenario == 'ex':
        save_name = 'pre-existing'
    elif plot_scenario == 'em':
        save_name = 'emerging'
    elif plot_scenario == 'bd':
        save_name = 'braiding'
    pl.screenshot(save_dir + save_name + '.png')
    pl.export_html(save_dir + save_name + '.html')

# pl.show(full_screen=True)