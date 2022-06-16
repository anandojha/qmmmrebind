"""

"""

import math

import numpy as np

def torsiondrive_input_to_xyz(psi_input_file, xyz_file):

    """
    Returns an xyz file from a torsiondrive formatted
    input file.

    Parameters
    ----------
    psi_input_file : str
        Input file for the psi4 QM engine.

    xyz_file : str
        XYZ format file to write the coords of the system.

    """
    with open(psi_input_file, "r") as f:
        lines = f.readlines()
    for i in range(len(lines)):
        if "molecule {" in lines[i]:
            to_begin = int(i)
        if "set {" in lines[i]:
            to_end = int(i)
    xyz_lines = lines[to_begin + 2 : to_end - 1]
    with open(xyz_file, "w") as f:
        f.write(str(len(xyz_lines)) + "\n")
        f.write(xyz_file + "\n")
        for i in xyz_lines:
            f.write(i)

def dihedral_energy(x, k1, k2, k3, k4=0):
    """
    Expression for the dihedral energy.
    """
    energy_1 = k1 * (1 + np.cos(1 * x * 0.01745))
    energy_2 = k2 * (1 - np.cos(2 * x * 0.01745))
    energy_3 = k3 * (1 + np.cos(3 * x * 0.01745))
    energy_4 = k4 * (1 - np.cos(4 * x * 0.01745))
    my_dihedral_energy = energy_1 + energy_2 + energy_3 + energy_4
    return my_dihedral_energy

def error_function(delta_qm, delta_mm):
    """
    Root Mean Squared Error.
    """
    squared_error = np.square(np.subtract(delta_qm, delta_mm))
    mean_squared_error = squared_error.mean()
    root_mean_squared_error = math.sqrt(mean_squared_error)
    return root_mean_squared_error

# TODO: remove? not used anywhere?
def error_function_boltzmann(delta_qm, delta_mm, T):
    """
    Boltzmann Root Mean Squared Error.
    """
    kb = 3.297623483 * 10 ** (-24)  # in cal/K
    delta_qm_boltzmann_weighted = [np.exp(-i / (kb * T)) for i in delta_qm]
    squared_error = (
        np.square(np.subtract(delta_qm, delta_mm))
        * delta_qm_boltzmann_weighted
    )
    mean_squared_error = squared_error.mean()
    root_mean_squared_error = math.sqrt(mean_squared_error)
    return root_mean_squared_error
