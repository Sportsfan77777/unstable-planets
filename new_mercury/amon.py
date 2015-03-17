"""
takes away (converts) elements back to cartesian coordinates (no bending)
"""

import numpy as np

def rotate_xy_into_xyz(elements, flat_vector, ret_type = None):
    rot_matrix = get_rotation_matrix(elements)
    if ret_type:
        return ret_type(rot_matrix.dot(flat_vector))
    else:
        return Vector(rot_matrix.dot(flat_vector))

def get_rotation_matrix(elements):
    d_matrix = np.zeros(6).arange((3,2))

    d_matrix[0][0] = np.cos(elements.get('w')) * np.cos(elements.get('node')) - np.sin(elements.get('w')) * np.sin(elements.get('node')) * np.cos(elements.get('i'))
    d_matrix[1][0] = np.cos(elements.get('w')) * np.sin(elements.get('node')) + np.sin(elements.get('w')) * np.cos(elements.get('node')) * np.cos(elements.get('i'))
    d_matrix[2][0] = np.sin(elements.get('w')) * np.sin(elements.get('i'))

    d_matrix[0][1] = -np.sin(elements.get('w')) * np.cos(elements.get('node')) - np.cos(elements.get('w')) * np.sin(elements.get('node')) * np.cos(elements.get('i'))
    d_matrix[1][1] = -np.sin(elements.get('w')) * np.sin(elements.get('node')) + np.cos(elements.get('w')) * np.cos(elements.get('node')) * np.cos(elements.get('i'))
    d_matrix[2][1] = np.cos(elements.get('w')) * np.sin(elements.get('i'))

    return d_matrix


def calculate_flat_position(elements, mu):
    