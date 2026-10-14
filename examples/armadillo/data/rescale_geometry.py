import trimesh
import numpy as np

def rescale_mesh_to_bounding_box(mesh, xmin, xmax):
    '''
    Rescales a mesh to fit within a bounding box defined by xmin and xmax for each axis.
    
    INPUT
        mesh: The input trimesh.Trimesh object
        xmin: Minimum bound for each axis as a tuple (xmin_x, xmin_y, xmin_z)
        xmax: Maximum bound for each axis as a tuple (xmax_x, xmax_y, xmax_z)
    
    OUTPUT
        Trimesh object rescaled to the specified bounding box
    '''
    # Calculate the current bounding box of the mesh
    current_min, current_max = mesh.bounds
    current_size = current_max - current_min

    # Target size based on specified xmin and xmax
    target_size = np.array(xmax) - np.array(xmin)

    # Scaling factor for each axis to fit the mesh within the target bounding box
    scale_factors = target_size / current_size
    scale = np.min(scale_factors)

    # Scale the mesh vertices
    mesh.vertices = np.array(mesh.vertices) * scale

    # Calculate the new bounding box after scaling
    new_min, _ = mesh.bounds

    # Translate the mesh to align it with the specified bounding box
    translation = np.array(xmin) - new_min
    mesh.vertices = np.array(mesh.vertices) + translation

    return mesh

# Load the OBJ file
mesh = trimesh.load("armadillo.obj")
xmin = (12, 0.5, 0.5)     # Minimum coordinates for the bounding box
xmax = (16, 7.5, 7.5)     # Maximum coordinates for the bounding box

# Rescale the mesh & output
rescaled_mesh = rescale_mesh_to_bounding_box(mesh, xmin, xmax)
rescaled_mesh.fix_normals()
rescaled_mesh.export("armadillo_withnormals.obj",include_normals=True)
