"""

"""

import math

import numpy as np

import modules.linear_algebra as linear_algebra

def force_angle_constant(
    atom_A,
    atom_B,
    atom_C,
    bond_lengths,
    eigenvalues,
    eigenvectors,
    coords,
    scaling_1,
    scaling_2,
):

    """
    Calculates force constant according to Equation 14 of
    Seminario calculation paper; returns angle (in kcal/mol/rad^2)
    and equilibrium angle (in degrees).

    Parameters
    ----------
    atom_A : int
        Index of atom A (left, starting from 0).

    atom_B : int
        Index of atom B (center, starting from 0).

    atom_C : int
        Index of atom C (right, starting from 0).

    bond_lengths : (N, N) array
        An N * N array containing the bond lengths for
        all the possible pairs of atoms.

    eigenvalues : (N, N, 3) array
        A numpy array of shape (N, N, 3) containing
        eigenvalues of the hessian matrix, where N
        is the total number of atoms.

    eigenvectors : (3, 3, N, N) array
        A numpy array of shape (3, 3, N, N) containing
        eigenvectors of the hessian matrix.

    coords : (N, 3) array
        A numpy array of shape (N, 3) having the X, Y and Z
        coordinates of all N atoms.

    scaling_1 : float
        Factor to scale the projections of eigenvalues for AB.

    scaling_2 : float
        Factor to scale the projections of eigenvalues for BC.

    Returns
    -------
    k_theta : float
        Force angle constant calculated using modified
        seminario method.

    k_0 : float
        Equilibrium angle between AB and BC.

    """
    # Vectors along bonds calculated
    diff_AB = coords[atom_B, :] - coords[atom_A, :]
    norm_diff_AB = np.linalg.norm(diff_AB)
    u_AB = diff_AB / norm_diff_AB
    diff_CB = coords[atom_B, :] - coords[atom_C, :]
    norm_diff_CB = np.linalg.norm(diff_CB)
    u_CB = diff_CB / norm_diff_CB
    # Bond lengths and eigenvalues found
    bond_length_AB = bond_lengths[atom_A, atom_B]
    eigenvalues_AB = eigenvalues[atom_A, atom_B, :]
    eigenvectors_AB = eigenvectors[0:3, 0:3, atom_A, atom_B]
    bond_length_BC = bond_lengths[atom_B, atom_C]
    eigenvalues_CB = eigenvalues[atom_C, atom_B, :]
    eigenvectors_CB = eigenvectors[0:3, 0:3, atom_C, atom_B]
    # Normal vector to angle plane found
    u_N = linear_algebra.unit_vector_N(u_CB, u_AB)
    u_PA = np.cross(u_N, u_AB)
    norm_u_PA = np.linalg.norm(u_PA)
    u_PA = u_PA / norm_u_PA
    u_PC = np.cross(u_CB, u_N)
    norm_u_PC = np.linalg.norm(u_PC)
    u_PC = u_PC / norm_u_PC
    sum_first = 0
    sum_second = 0
    # Projections of eigenvalues
    for i in range(0, 3):
        eig_AB_i = eigenvectors_AB[:, i]
        eig_BC_i = eigenvectors_CB[:, i]
        sum_first = sum_first + (
            eigenvalues_AB[i] * abs(linear_algebra.dot_product(u_PA, eig_AB_i))
        )
        sum_second = sum_second + (
            eigenvalues_CB[i] * abs(linear_algebra.dot_product(u_PC, eig_BC_i))
        )
    # Scaling due to additional angles - Modified Seminario Part
    sum_first = sum_first / scaling_1
    sum_second = sum_second / scaling_2
    # Added as two springs in series
    k_theta = (1 / ((bond_length_AB ** 2) * sum_first)) + (
        1 / ((bond_length_BC ** 2) * sum_second)
    )
    k_theta = 1 / k_theta
    k_theta = -k_theta  # Change to OPLS form
    k_theta = abs(k_theta * 0.5)  # Change to OPLS form
    # Equilibrium Angle
    theta_0 = math.degrees(math.acos(np.dot(u_AB, u_CB)))
    # If the vectors u_CB and u_AB are linearly dependent u_N cannot be defined.
    # This case is dealt with here :
    if abs(sum((u_CB) - (u_AB))) < 0.01 or (
        abs(sum((u_CB) - (u_AB))) > 1.99 and abs(sum((u_CB) - (u_AB))) < 2.01
    ):
        scaling_1 = 1
        scaling_2 = 1
        [k_theta, theta_0] = force_angle_constant_special_case(
            atom_A,
            atom_B,
            atom_C,
            bond_lengths,
            eigenvalues,
            eigenvectors,
            coords,
            scaling_1,
            scaling_2,
        )
    return k_theta, theta_0

def force_angle_constant_special_case(
    atom_A,
    atom_B,
    atom_C,
    bond_lengths,
    eigenvalues,
    eigenvectors,
    coords,
    scaling_1,
    scaling_2,
):

    """
    Calculates force constant according to Equation 14
    of Seminario calculation paper when the vectors
    u_CB and u_AB are linearly dependent and u_N cannot
    be defined. It instead takes samples of u_N across a
    unit sphere for the calculation; returns angle
    (in kcal/mol/rad^2) and equilibrium angle in degrees.

    Parameters
    ----------
    atom_A : int
        Index of atom A (left, starting from 0).

    atom_B : int
        Index of atom B (center, starting from 0).

    atom_C : int
        Index of atom C (right, starting from 0).

    bond_lengths : (N, N) array
        An N * N array containing the bond lengths for
        all the possible pairs of atoms.

    eigenvalues : (N, N, 3) array
        A numpy array of shape (N, N, 3) containing
        eigenvalues of the  hessian matrix, where N
        is the total number of atoms.

    eigenvectors : (3, 3, N, N) array
        A numpy array of shape (3, 3, N, N) containing
        eigenvectors of the hessian matrix.

    coords : (N, 3) array
        A numpy array of shape (N, 3) having the X, Y,
        and Z coordinates of all N atoms.

    scaling_1 : float
        Factor to scale the projections of eigenvalues for AB.

    scaling_2 : float
        Factor to scale the projections of eigenvalues for BC.

    Returns
    -------
    k_theta : float
        Force angle constant calculated using modified
        seminario method.
    k_0 : float
        Equilibrium angle between AB and BC.

    """
    # Vectors along bonds calculated
    diff_AB = coords[atom_B, :] - coords[atom_A, :]
    norm_diff_AB = np.linalg.norm(diff_AB)
    u_AB = diff_AB / norm_diff_AB
    diff_CB = coords[atom_B, :] - coords[atom_C, :]
    norm_diff_CB = np.linalg.norm(diff_CB)
    u_CB = diff_CB / norm_diff_CB
    # Bond lengths and eigenvalues found
    bond_length_AB = bond_lengths[atom_A, atom_B]
    eigenvalues_AB = eigenvalues[atom_A, atom_B, :]
    eigenvectors_AB = eigenvectors[0:3, 0:3, atom_A, atom_B]
    bond_length_BC = bond_lengths[atom_B, atom_C]
    eigenvalues_CB = eigenvalues[atom_C, atom_B, :]
    eigenvectors_CB = eigenvectors[0:3, 0:3, atom_C, atom_B]
    k_theta_array = np.zeros((180, 360))
    # Find force constant with varying u_N (with vector uniformly
    # sampled across a sphere)
    for theta in range(0, 180):
        for phi in range(0, 360):
            r = 1
            u_N = [
                r
                * math.sin(math.radians(theta))
                * math.cos(math.radians(theta)),
                r
                * math.sin(math.radians(theta))
                * math.sin(math.radians(theta)),
                r * math.cos(math.radians(theta)),
            ]
            u_PA = np.cross(u_N, u_AB)
            u_PA = u_PA / np.linalg.norm(u_PA)
            u_PC = np.cross(u_CB, u_N)
            u_PC = u_PC / np.linalg.norm(u_PC)
            sum_first = 0
            sum_second = 0
            # Projections of eigenvalues
            for i in range(0, 3):
                eig_AB_i = eigenvectors_AB[:, i]
                eig_BC_i = eigenvectors_CB[:, i]
                sum_first = sum_first + (
                    eigenvalues_AB[i] * abs(linear_algebra.dot_product(u_PA, eig_AB_i))
                )
                sum_second = sum_second + (
                    eigenvalues_CB[i] * abs(linear_algebra.dot_product(u_PC, eig_BC_i))
                )
            # Added as two springs in series
            k_theta_ij = (1 / ((bond_length_AB ** 2) * sum_first)) + (
                1 / ((bond_length_BC ** 2) * sum_second)
            )
            k_theta_ij = 1 / k_theta_ij
            k_theta_ij = -k_theta_ij  # Change to OPLS form
            k_theta_ij = abs(k_theta_ij * 0.5)  # Change to OPLS form
            k_theta_array[theta, phi] = k_theta_ij
    # Removes cases where u_N was linearly dependent of u_CB or u_AB.
    # Force constant used is taken as the mean.
    k_theta = np.mean(np.mean(k_theta_array))
    # Equilibrium Angle independent of u_N
    theta_0 = math.degrees(math.cos(np.dot(u_AB, u_CB)))
    return k_theta, theta_0


def force_constant_bond(atom_A, atom_B, eigenvalues, eigenvectors, coords):

    """
    Calculates the bond force constant for the bonds in the
    molecule according to equation 10 of seminario paper,
    given the bond atoms' indices and the corresponding
    eigenvalues, eigenvectors and coordinates matrices.

    Parameters
    ----------
    atom_A : int
        Index of Atom A.

    atom_B : int
        Index of Atom B.

    eigenvalues : (N, N, 3) array
        A numpy array of shape (N, N, 3) containing eigenvalues
        of the hessian matrix, where N is the total number
        of atoms.

    eigenvectors : (3, 3, N, N) array
        A numpy array of shape (3, 3, N, N) containing the
        eigenvectors of the hessian matrix.

    coords : (N, 3) array
        A numpy array of shape (N, 3) having the X, Y, and
        Z  coordinates of all N atoms.

    Returns
    --------
    k_AB : float
        Bond Force Constant value for the bond with atoms A and B.

    """
    # Eigenvalues and eigenvectors calculated
    eigenvalues_AB = eigenvalues[atom_A, atom_B, :]
    eigenvectors_AB = eigenvectors[:, :, atom_A, atom_B]
    # Vector along bond
    diff_AB = np.array(coords[atom_B, :]) - np.array(coords[atom_A, :])
    norm_diff_AB = np.linalg.norm(diff_AB)
    unit_vectors_AB = diff_AB / norm_diff_AB
    k_AB = 0
    # Projections of eigenvalues
    for i in range(0, 3):
        my_dot_product = abs(np.dot(unit_vectors_AB, eigenvectors_AB[:, i]))
        k_AB = k_AB + (eigenvalues_AB[i] * my_dot_product)
    k_AB = -k_AB * 0.5  # Convert to OPLS form
    return k_AB

