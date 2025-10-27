import numpy as np
import pyvista as pv
from numpy.linalg import norm

def create_axis(center, radius, angle_range, num_points, point_size=10):
    '''
    Create an arc in X-Z plane as the axis of the helix, represented as points
    
    Input:
    - center:       arc center
    - radius:       arc radius
    - angle_range:  the angle range of the arc, in [rad.]
    - num_points:   number of points to show the arc
    - point_size:   size of the points
    
    Output: 
    - axis:         Pyvista point cloud object for the axis
    - axis_points:  coordinates of the axis points
    '''
    # Creating angles of arc points
    angle_beg, angle_end = angle_range
    theta = np.linspace(angle_beg, angle_end, num_points)
    
    # Calculating point coordinate along the axis
    coord_x = center[0] + radius * np.cos(theta)
    coord_y = np.full_like(theta, center[1])
    coord_z = center[2] + radius * np.sin(theta)
    axis_points = np.column_stack((coord_x, coord_y, coord_z))
    
    # Creating Pyvista point cloud object
    axis = pv.PolyData(axis_points)
    axis['point_size'] = np.full(num_points, point_size)  # Set point size attribute
    
    return axis, axis_points

def create_helix(axis_points, helix_radius, turns, final_angle):
    '''
    Create the helix around the axis with specified final angle
    
    Input:
    - axis_points:  coordinates of the axis points
    - helix_radius: helix radius
    - turns:        helix turns throughout the axis
    - final_angle:  angle of the end point at the end, in [rad.]
    
    Output:
    - helix:        Pyvista object for the helix
    '''
    num_points = len(axis_points)
    helix_points = []
    
    # Calculate total angle change (based on turns) plus the final angle
    total_angle_change = 2 * np.pi * turns
    start_angle = final_angle - total_angle_change
    
    for ip in range(num_points):
        # Extracting point on the axis
        p_axis = axis_points[ip]
        
        # Calculating the tangential direction
        if ip == 0: 
            tangent = axis_points[ip+1] - p_axis
        elif ip ==  num_points - 1:
            tangent = p_axis - axis_points[ip-1]
        else:
            tangent = axis_points[ip+1] - axis_points[ip-1]
        # Normalizing the tangential direction
        tangent = tangent / norm(tangent)
        
        # Calculating two orthogonal directions in the normal plane
        ref = np.array([0, 1, 0]) if abs(tangent[1]) < 0.9 else np.array([1, 0, 0])
        # The first normal direction
        normal1 = np.cross(tangent, ref)
        normal1 = normal1 / norm(normal1)
        # The second normal direction
        normal2 = np.cross(tangent, normal1)
        
        # Calculating angles of helix based on final angle
        # 从起点角度线性过渡到最终角度
        angle = start_angle + total_angle_change * (ip / num_points)
        
        # Calculating point coordinate along the helix
        helix_point = p_axis + helix_radius * (np.cos(angle) * normal1 + np.sin(angle) * normal2)
        helix_points.append(helix_point)
    
    # Creating Pyvista line object
    helix = pv.lines_from_points(np.array(helix_points))
    
    return helix

def create_torus(axis_points, torus_radius, num_segments=50):
    '''
    Create a torus with the axis as its center line
    
    Input:
    - axis_points:  coordinates of the axis points
    - torus_radius: radius of the torus cross-section
    - num_segments: number of segments for the torus cross-section
    
    Output:
    - torus:        Pyvista mesh object for the torus
    '''
    # Creating angles in the torus cross-section
    phi = np.linspace(0, 2 * np.pi, num_segments)
    torus_points = []
    
    # Creating points for each cross-section along the axis
    for ip in range(len(axis_points)):
        p_axis = axis_points[ip]
        
        # Calculating tangential direction of the axis
        if ip == 0:
            tangent = axis_points[ip+1] - p_axis
        elif ip == len(axis_points) - 1:
            tangent = p_axis - axis_points[ip-1]
        else:
            tangent = axis_points[ip+1] - axis_points[ip-1]
        tangent = tangent / norm(tangent)
        
        # Calculating two orthogonal directions in the normal plane
        ref = np.array([0, 1, 0]) if abs(tangent[1]) < 0.9 else np.array([1, 0, 0])
        normal1 = np.cross(tangent, ref)
        normal1 = normal1 / norm(normal1)
        normal2 = np.cross(tangent, normal1)
        
        # Creating points for current cross-section
        for angle in phi:
            torus_point = p_axis + torus_radius * (np.cos(angle) * normal1 + np.sin(angle) * normal2)
            torus_points.append(torus_point)
    
    # Creating Pyvista mesh object
    torus_poly = pv.PolyData(np.array(torus_points))
    torus = torus_poly.delaunay_3d(alpha=torus_radius * 1.2)
    
    return torus

def create_horizontal_plane(center, length, width):
    '''
    Create a horizontal rectangular plane (Z = constant)
    
    Input:
    - center:       plane center coordinates
    - length:       plane length in X direction
    - width:        plane width in Y direction
    
    Output:
    - plane:        Pyvista object for the plane
    '''
    # Calculating half dimensions
    half_len = length / 2
    half_wid = width / 2
    
    # Defining four vertex of the rectangle
    points = [
        [center[0] - half_len, center[1] - half_wid, center[2]],
        [center[0] + half_len, center[1] - half_wid, center[2]],
        [center[0] + half_len, center[1] + half_wid, center[2]],
        [center[0] - half_len, center[1] + half_wid, center[2]]
    ]
    # The face order
    faces = [4, 0, 1, 2, 3]
    
    # Creating PyVista plane object
    plane = pv.PolyData(points, faces)
    
    return plane