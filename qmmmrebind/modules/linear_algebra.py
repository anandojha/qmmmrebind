"""
Commonly used linear algebra capabilities.
"""

import numpy as np

# TODO: rename?
def unit_vector_N(u_BC, u_AB):

    """
    Calculates unit normal vector perpendicular to plane ABC.

    Parameters
    ----------
    u_BC : (.. , 1, 3) array
        Unit vector from atom B to atom C.

    u_AB : (..., 1, 3) array
        Unit vector from atom A to atom B.

    Returns
    -------
    u_N : (..., 1, 3) array
        Unit normal vector perpendicular to plane ABC.

    Examples
    --------
    >>> u_BC = [0.34040355, 0.62192853, 0.27011169]
    >>> u_AB = [0.28276792, 0.34232697, 0.02370306]
    >>> unit_vector_N(u_BC, u_AB)
    array([-0.65161629,  0.5726879 , -0.49741811])
    """
    cross_product = np.cross(u_BC, u_AB)
    norm_u_N = np.linalg.norm(cross_product)
    u_N = cross_product / norm_u_N
    return u_N

def u_PA_from_angles(atom_A, atom_B, atom_C, coords):

    """
    Returns the vector in the plane A,B,C and perpendicular to AB.

    Parameters
    ----------
    atom_A : int
        Index of atom A (left, starting from 0).

    atom_B : int
        Index of atom B (center, starting from 0).

    atom_C : int
        Index of atom C (right, starting from 0).

    coords : (..., N, 3) array
        An array which contains the coordinates of all
        the N atoms.

    """
    diff_AB = coords[atom_B, :] - coords[atom_A, :]
    norm_diff_AB = np.linalg.norm(diff_AB)
    u_AB = diff_AB / norm_diff_AB
    diff_CB = coords[atom_B, :] - coords[atom_C, :]
    norm_diff_CB = np.linalg.norm(diff_CB)
    u_CB = diff_CB / norm_diff_CB
    u_N = unit_vector_N(u_CB, u_AB)
    u_PA = np.cross(u_N, u_AB)
    norm_PA = np.linalg.norm(u_PA)
    u_PA = u_PA / norm_PA
    return u_PA

def dot_product(u_PA, eig_AB):

    """
    Returns the dot product of two vectors.

    Parameters
    ----------
    u_PA : (..., 1, 3) array
        Unit vector perpendicular to AB and in the
        plane of A, B, C.

    eig_AB : (..., 3, 3) array
        Eigenvectors of the hessian matrix for
        the bond AB.

    """
    x = 0
    for i in range(0, 3):
        x = x + u_PA[i] * eig_AB[i].conjugate()
    return x

def u_PA_from_angles(atom_A, atom_B, atom_C, coords):

    """
    Returns the vector in the plane A,B,C and perpendicular to AB.

    Parameters
    ----------
    atom_A : int
        Index of atom A (left, starting from 0).

    atom_B : int
        Index of atom B (center, starting from 0).

    atom_C : int
        Index of atom C (right, starting from 0).

    coords : (..., N, 3) array
        An array containing the coordinates of all the N atoms.

    Returns
    -------
    u_PA : (..., 1, 3) array
        Unit vector perpendicular to AB and in the plane of A, B, C.

    """
    diff_AB = coords[atom_B, :] - coords[atom_A, :]
    norm_diff_AB = np.linalg.norm(diff_AB)
    u_AB = diff_AB / norm_diff_AB
    diff_CB = coords[atom_B, :] - coords[atom_C, :]
    norm_diff_CB = np.linalg.norm(diff_CB)
    u_CB = diff_CB / norm_diff_CB
    u_N = unit_vector_N(u_CB, u_AB)
    u_PA = np.cross(u_N, u_AB)
    norm_PA = np.linalg.norm(u_PA)
    u_PA = u_PA / norm_PA
    return u_PA

