"""

"""

from operator import itemgetter
import itertools
import math
import re

import numpy as np
import pandas as pd
from biopandas.pdb import PandasPdb

import modules.constants as const
import modules.base as base
import modules.modified_seminario as modified_seminario
import modules.linear_algebra as linear_algebra

class ParameterizeGuest:

    """
    A class used to obtain force field parameters for the ligand (bond,
    angle and charge parameters) from QM calculations.

    This class contain methods to process the output files of the
    Gaussian QM output files (.chk, .fchk and .log files). Methods
    in the class extract the unprocessed hessian matrix from the
    Gaussian QM calculations, processes it and uses the Modified
    Seminario Method to ontain the bond and angle parameters. The
    class also extracts the QM charges from the log file.

    ...

    Attributes
    ----------
    xyz_file: str, optional
        XYZ file for ligand coordinates obtained from its corresponding
        formatted checkpoint file.

    coordinate_file: str, optional
        Text file containing the ligand coordinates (extracted
        from the formatted checkpoint file).

    unprocessed_hessian_file: str, optional
        Unprocessed hessian matrix of the ligand obtained from the
        formatted checkpoint file.

    bond_list_file: str, optional
        Text file containing the bond information of the ligand extracted
        from the log file.

    angle_list_file: str, optional
        Text file containing the angle information of the ligand extracted
        from the log file.

    hessian_file: str, optional
        Processed hessian matrix of the ligand.

    atom_names_file: str, optional
        Text file containing the list of atom names from the fchk file.

    bond_parameter_file: str, optional
        Text file containing the bond parameters for the ligand obtained
        using the Modified Seminario method.

    angle_parameter_file: str, optional
        Text file containing the angle parameters of the ligand obtained
        using the Modified Seminario method..

    charge_parameter_file: str, optional
        Text file containing the QM charges of the ligand.

    guest_pdb: str, optional
        Ligand PDB file with atom numbers beginning from 1.

    proper_dihedral_file: str, optional
        A text file containing proper dihedral angles of the ligand.

    functional: str, optional
        Exchange/Correlation or hybrid functional to use in the Gaussian
        QM calculation.

    basis_set: str, optional
        Basis set to use for the Gaussian QM calculation.

    """

    def __init__(
        self,
        xyz_file="guest_coords.xyz",
        coordinate_file="guest_coordinates.txt",
        unprocessed_hessian_file="guest_unprocessed_hessian.txt",
        bond_list_file="guest_bond_list.txt",
        angle_list_file="guest_angle_list.txt",
        hessian_file="guest_hessian.txt",
        atom_names_file="guest_atom_names.txt",
        bond_parameter_file="guest_bonds.txt",
        angle_parameter_file="guest_angles.txt",
        charge_parameter_file="guest_qm_surround_charges.txt",
        guest_pdb="guest_init_II.pdb",
        proper_dihedral_file="proper_dihedrals.txt",
        functional="B3LYP",
        basis_set="6-31G",
    ):

        self.xyz_file = xyz_file
        self.coordinate_file = coordinate_file
        self.unprocessed_hessian_file = unprocessed_hessian_file
        self.bond_list_file = bond_list_file
        self.angle_list_file = angle_list_file
        self.hessian_file = hessian_file
        self.atom_names_file = atom_names_file
        self.bond_parameter_file = bond_parameter_file
        self.angle_parameter_file = angle_parameter_file
        self.charge_parameter_file = charge_parameter_file
        self.guest_pdb = guest_pdb
        self.proper_dihedral_file = proper_dihedral_file
        self.functional = functional
        self.basis_set = basis_set

    def get_xyz(self):
        """
        Saves XYZ file from the formatted checkpoint file.
        """
        fchk_file = self.guest_pdb[:-4] + ".fchk"
        with open(fchk_file, "r") as f:
            lines = f.readlines()
        for i in range(len(lines)):
            if "Current cartesian coordinates" in lines[i]:
                no_coordinates = re.findall(r"\d+|\d+.\d+", lines[i])
                no_coordinates = int(no_coordinates[0])
                to_begin = int(i)
                
        cartesian_coords = lines[
            to_begin + 1 : to_begin + 1 + int(math.ceil(no_coordinates / 5))
        ]
        cartesian_list = []
        for i in range(len(cartesian_coords)):
            cartesian_list.append(cartesian_coords[i].strip().split())
        
        coordinates_list = [
            item for sublist in cartesian_list for item in sublist
        ]
        # Converted from Atomic units (Bohrs) to Angstroms
        list_coords = [float(x) * const.BOHRS_PER_ANGSTROM for x in coordinates_list]
        for i in range(len(lines)):
            if "Atomic numbers" in lines[i]:
                to_begin = int(i)
            if "Nuclear charges" in lines[i]:
                to_end = int(i)
        atomic_number_strings = lines[to_begin + 1 : to_end]
        atom_numbers_nested = []
        for i in range(len(atomic_number_strings)):
            atom_numbers_nested.append(atomic_number_strings[i].strip().split())
        numbers = [item for sublist in atom_numbers_nested for item in sublist]
        N = int(no_coordinates / 3)
        # Opens the new xyz file
        with open(self.xyz_file, "w") as file:
            file.write(str(N) + "\n \n")
            coords = np.zeros((N, 3))
            n = 0
            names = []
            # Gives name for atomic number
            for x in range(0, len(numbers)):
                names.append(const.element_list[int(numbers[x]) - 1][1])
            # Print coordinates to new input_coords.xyz file
            for i in range(0, N):
                for j in range(0, 3):
                    coords[i][j] = list_coords[n]
                    n = n + 1
                file.write(
                    names[i]
                    + str(round(coords[i][0], 3))
                    + " "
                    + str(round(coords[i][1], 3))
                    + " "
                    + str(round(coords[i][2], 3))
                    + "\n"
                )
            
        np.savetxt(self.coordinate_file, coords, fmt="%s")

    def get_unprocessed_hessian(self):
        """
        Saves a text file of the unprocessed hessian matrix from the
        formatted checkpoint file.
        """
        fchk_file = self.guest_pdb[:-4] + ".fchk"
        with open(fchk_file, "r") as f:
            lines = f.readlines()
        for i in range(len(lines)):
            if "Cartesian Force Constants" in lines[i]:
                no_hessian = re.findall(r"\d+|\d+.\d+", lines[i])
                no_hessian = int(no_hessian[0])
                to_begin = int(i)
        hessian = lines[
            to_begin + 1 : to_begin + 1 + int(math.ceil(no_hessian / 5))
        ]
        hessian_list = []
        for i in range(len(hessian)):
            hessian_list.append(hessian[i].strip().split())
        unprocessed_Hessian = [
            item for sublist in hessian_list for item in sublist
        ]
        np.savetxt(
            self.unprocessed_hessian_file, unprocessed_Hessian, fmt="%s",
        )

    def get_bond_angles(self):
        """
        Saves a text file containing bonds and angles from the gaussian
        log file.
        """
        log_file = self.guest_pdb[:-4] + ".log"
        with open(log_file, "r") as fid:
            tline = fid.readline()
            bond_list = []
            angle_list = []
            tmp = "R"  # States if bond or angle
            # Finds the bond and angles from the .log file
            while tline:
                tline = fid.readline()
                # Line starts at point when bond and angle list occurs
                if (
                    len(tline) > 80
                    and tline[0:81].strip()
                    == "! Name  Definition              Value          Derivative Info.                !"
                ):
                    tline = fid.readline()
                    tline = fid.readline()
                    # Stops when all bond and angles recorded
                    while (tmp[0] == "R") or (tmp[0] == "A"):
                        line = tline.split()
                        tmp = line[1]
                        # Bond or angles listed as string
                        list_terms = line[2][2:-1]
                        # Bond List
                        if tmp[0] == "R":
                            x = list_terms.split(",")
                            # Subtraction due to python array indexing at 0
                            x = [(int(i) - 1) for i in x]
                            bond_list.append(x)
                            # Angle List
                        if tmp[0] == "A":
                            x = list_terms.split(",")
                            # Subtraction due to python array indexing at 0
                            x = [(int(i) - 1) for i in x]
                            angle_list.append(x)
                        tline = fid.readline()
                    # Leave loop
                    tline = -1
            np.savetxt(self.bond_list_file, bond_list, fmt="%s")
            np.savetxt(self.angle_list_file, angle_list, fmt="%s")
        
    def get_hessian(self):
        """
        Extracts hessian matrix from the unprocessed hessian matrix
        and saves into a new file.
        """
        unprocessed_Hessian = np.loadtxt(self.unprocessed_hessian_file)
        fchk_file = self.guest_pdb[:-4] + ".fchk"
        with open(fchk_file, "r") as f:
            lines = f.readlines()
        for i in range(len(lines)):
            if "Current cartesian coordinates" in lines[i]:
                no_coordinates = re.findall(r"\d+|\d+.\d+", lines[i])
                no_coordinates = int(no_coordinates[0])
        
        N = int(no_coordinates / 3)
        length_hessian = 3 * N
        hessian = np.zeros((length_hessian, length_hessian))
        m = 0
        # Write the hessian in a 2D array format
        for i in range(0, length_hessian):
            for j in range(0, (i + 1)):
                hessian[i][j] = unprocessed_Hessian[m]
                hessian[j][i] = unprocessed_Hessian[m]
                m = m + 1
        hessian = (hessian * const.HARTREE_PER_KCAL_MOL) / (
            const.BOHRS_PER_ANGSTROM ** 2
        )  # Change from Hartree/bohr to kcal/mol/ang
        np.savetxt(self.hessian_file, hessian, fmt="%s")

    def get_atom_names(self):
        """
        Saves a list of atom names from the formatted checkpoint file.
        """
        fchk_file = self.guest_pdb[:-4] + ".fchk"
        with open(fchk_file, "r") as f:
            lines = f.readlines()
        for i in range(len(lines)):
            if "Atomic numbers" in lines[i]:
                to_begin = int(i)
            if "Nuclear charges" in lines[i]:
                to_end = int(i)
        atomic_numbers = lines[to_begin + 1 : to_end]
        atom_numbers = []
        for i in range(len(atomic_numbers)):
            atom_numbers.append(atomic_numbers[i].strip().split())
        numbers = [item for sublist in atom_numbers for item in sublist]
        names = []
        # Gives name for atomic number
        for x in range(0, len(numbers)):
            names.append(const.element_list[int(numbers[x]) - 1][1])
        atom_names = []
        for i in range(0, len(names)):
            atom_names.append(names[i].strip() + str(i + 1))
        np.savetxt(self.atom_names_file, atom_names, fmt="%s")

    def get_bond_angle_params(self):
        """
        Saves the bond and angle parameter files obtained from
        the formatted checkpoint file.
        """
        fchk_file = self.guest_pdb[:-4] + ".fchk"
        with open(fchk_file, "r") as f:
            lines = f.readlines()
        for i in range(len(lines)):
            if "Current cartesian coordinates" in lines[i]:
                no_coordinates = re.findall(r"\d+|\d+.\d+", lines[i])
                no_coordinates = int(no_coordinates[0])
        N = int(no_coordinates / 3)
        coords = np.loadtxt(self.coordinate_file)
        hessian = np.loadtxt(self.hessian_file)
        bond_list = np.loadtxt(self.bond_list_file, dtype=int)
        atom_names = np.loadtxt(self.atom_names_file, dtype=str)
        # Find bond lengths
        bond_lengths = np.zeros((N, N))
        for i in range(0, N):
            for j in range(0, N):
                diff_i_j = np.array(coords[i, :]) - np.array(coords[j, :])
                bond_lengths[i][j] = np.linalg.norm(diff_i_j)
        eigenvectors = np.empty((3, 3, N, N), dtype=complex)
        eigenvalues = np.empty((N, N, 3), dtype=complex)
        partial_hessian = np.zeros((3, 3))
        for i in range(0, N):
            for j in range(0, N):
                partial_hessian = hessian[
                    (i * 3) : ((i + 1) * 3), (j * 3) : ((j + 1) * 3)
                ]
                [a, b] = np.linalg.eig(partial_hessian)
                eigenvalues[i, j, :] = a
                eigenvectors[:, :, i, j] = b
        # Modified Seminario method to find the bond parameters and
        # print them to file
        file_bond = open(self.bond_parameter_file, "w")
        k_b = np.zeros(len(bond_list))
        bond_length_list = np.zeros(len(bond_list))
        unique_values_bonds = []  # Used to find average values
        for i in range(0, len(bond_list)):
            AB = modified_seminario.force_constant_bond(
                bond_list[i][0],
                bond_list[i][1],
                eigenvalues,
                eigenvectors,
                coords,
            )
            BA = modified_seminario.force_constant_bond(
                bond_list[i][1],
                bond_list[i][0],
                eigenvalues,
                eigenvectors,
                coords,
            )
            # Order of bonds sometimes causes slight differences,
            # find the mean
            k_b[i] = np.real((AB + BA) / 2)
            # Vibrational_scaling takes into account DFT deficities /
            # anharmocity
            vibrational_scaling = const.get_vibrational_scaling(
                functional=self.functional, basis_set=self.basis_set
            )
            vibrational_scaling_squared = vibrational_scaling ** 2
            k_b[i] = k_b[i] * vibrational_scaling_squared
            bond_length_list[i] = bond_lengths[bond_list[i][0]][
                bond_list[i][1]
            ]
            file_bond.write(
                atom_names[bond_list[i][0]]
                + "-"
                + atom_names[bond_list[i][1]]
                + "  "
            )
            file_bond.write(
                str("%#.5g" % k_b[i])
                + "   "
                + str("%#.4g" % bond_length_list[i])
                + "   "
                + str(bond_list[i][0] + 1)
                + "   "
                + str(bond_list[i][1] + 1)
            )
            file_bond.write("\n")
            unique_values_bonds.append(
                [
                    atom_names[bond_list[i][0]],
                    atom_names[bond_list[i][1]],
                    k_b[i],
                    bond_length_list[i],
                    1,
                ]
            )
        file_bond.close()
        angle_list = np.loadtxt(self.angle_list_file, dtype=int)
        # Modified Seminario method to find the angle parameters
        # and print them to file
        file_angle = open(self.angle_parameter_file, "w")
        k_theta = np.zeros(len(angle_list))
        theta_0 = np.zeros(len(angle_list))
        unique_values_angles = []  # Used to find average values
        # Modified Seminario part goes here ...
        # Connectivity information for Modified Seminario Method
        central_atoms_angles = []
        # A structure is created with the index giving the central
        # atom of the angle,
        # an array then lists the angles with that central atom.
        # i.e. central_atoms_angles{3} contains an array of angles
        # with central atom 3
        for i in range(0, len(coords)):
            central_atoms_angles.append([])
            for j in range(0, len(angle_list)):
                if i == angle_list[j][1]:
                    # For angle ABC, atoms A C are written to array
                    AC_array = [angle_list[j][0], angle_list[j][2], j]
                    central_atoms_angles[i].append(AC_array)
                    # For angle ABC, atoms C A are written to array
                    CA_array = [angle_list[j][2], angle_list[j][0], j]
                    central_atoms_angles[i].append(CA_array)
        # Sort rows by atom number
        for i in range(0, len(coords)):
            central_atoms_angles[i] = sorted(
                central_atoms_angles[i], key=itemgetter(0)
            )
        # Find normals u_PA for each angle
        unit_PA_all_angles = []
        for i in range(0, len(central_atoms_angles)):
            unit_PA_all_angles.append([])
            for j in range(0, len(central_atoms_angles[i])):
                # For the angle at central_atoms_angles[i][j,:] the
                # corresponding u_PA value
                # is found for the plane ABC and bond AB, where ABC
                # corresponds to the order
                # of the arguements. This is why the reverse order
                # was also added
                unit_PA_all_angles[i].append(
                    linear_algebra.u_PA_from_angles(
                        central_atoms_angles[i][j][0],
                        i,
                        central_atoms_angles[i][j][1],
                        coords,
                    )
                )
        # Finds the contributing factors from the other angle terms
        # scaling_factor_all_angles
        # = cell(max(max(angle_list))); %This array will contain
        # scaling factor and angle list position
        scaling_factor_all_angles = []
        for i in range(0, len(central_atoms_angles)):
            scaling_factor_all_angles.append([])
            for j in range(0, len(central_atoms_angles[i])):
                n = 1
                m = 1
                angles_around = 0
                additional_contributions = 0
                scaling_factor_all_angles[i].append([0, 0])
                # Position in angle list
                scaling_factor_all_angles[i][j][1] = central_atoms_angles[i][
                    j
                ][2]
                # Goes through the list of angles with the same central atom
                # and computes the
                # term need for the modified Seminario method
                # Forwards directions, finds the same bonds with the central atom i
                while (
                    ((j + n) < len(central_atoms_angles[i]))
                    and central_atoms_angles[i][j][0]
                    == central_atoms_angles[i][j + n][0]
                ):
                    additional_contributions = (
                        additional_contributions
                        + (
                            abs(
                                np.dot(
                                    unit_PA_all_angles[i][j][:],
                                    unit_PA_all_angles[i][j + n][:],
                                )
                            )
                        )
                        ** 2
                    )
                    n = n + 1
                    angles_around = angles_around + 1
                # Backwards direction, finds the same bonds with the central atom i
                while ((j - m) >= 0) and central_atoms_angles[i][j][
                    0
                ] == central_atoms_angles[i][j - m][0]:
                    additional_contributions = (
                        additional_contributions
                        + (
                            abs(
                                np.dot(
                                    unit_PA_all_angles[i][j][:],
                                    unit_PA_all_angles[i][j - m][:],
                                )
                            )
                        )
                        ** 2
                    )
                    m = m + 1
                    angles_around = angles_around + 1
                if n != 1 or m != 1:
                    # Finds the mean value of the additional contribution to
                    # change to normal
                    # Seminario method comment out + part
                    scaling_factor_all_angles[i][j][0] = 1 + (
                        additional_contributions / (m + n - 2)
                    )
                else:
                    scaling_factor_all_angles[i][j][0] = 1
        scaling_factors_angles_list = []
        for i in range(0, len(angle_list)):
            scaling_factors_angles_list.append([])
        # Orders the scaling factors according to the angle list
        for i in range(0, len(central_atoms_angles)):
            for j in range(0, len(central_atoms_angles[i])):
                scaling_factors_angles_list[
                    scaling_factor_all_angles[i][j][1]
                ].append(scaling_factor_all_angles[i][j][0])
        # Finds the angle force constants with the scaling factors
        # included for each angle
        for i in range(0, len(angle_list)):
            # Ensures that there is no difference when the
            # ordering is changed
            [AB_k_theta, AB_theta_0] = modified_seminario.force_angle_constant(
                angle_list[i][0],
                angle_list[i][1],
                angle_list[i][2],
                bond_lengths,
                eigenvalues,
                eigenvectors,
                coords,
                scaling_factors_angles_list[i][0],
                scaling_factors_angles_list[i][1],
            )
            [BA_k_theta, BA_theta_0] = modified_seminario.force_angle_constant(
                angle_list[i][2],
                angle_list[i][1],
                angle_list[i][0],
                bond_lengths,
                eigenvalues,
                eigenvectors,
                coords,
                scaling_factors_angles_list[i][1],
                scaling_factors_angles_list[i][0],
            )
            k_theta[i] = (AB_k_theta + BA_k_theta) / 2
            theta_0[i] = (AB_theta_0 + BA_theta_0) / 2
            # Vibrational_scaling takes into account DFT
            # deficities/ anharmonicity
            k_theta[i] = k_theta[i] * vibrational_scaling_squared
            file_angle.write(
                atom_names[angle_list[i][0]]
                + "-"
                + atom_names[angle_list[i][1]]
                + "-"
                + atom_names[angle_list[i][2]]
                + "  "
            )
            file_angle.write(
                str("%#.4g" % k_theta[i])
                + "   "
                + str("%#.4g" % theta_0[i])
                + "   "
                + str(angle_list[i][0] + 1)
                + "   "
                + str(angle_list[i][1] + 1)
                + "   "
                + str(angle_list[i][2] + 1)
            )
            file_angle.write("\n")
            unique_values_angles.append(
                [
                    atom_names[angle_list[i][0]],
                    atom_names[angle_list[i][1]],
                    atom_names[angle_list[i][2]],
                    k_theta[i],
                    theta_0[i],
                    1,
                ]
            )
        file_angle.close()

    def get_charges(self):
        """
        Saves the atomic charges in a text file obtained from
        the Gaussian log file.
        """
        log_file = self.guest_pdb[:-4] + ".log"
        with open(log_file, "r") as f:
            lines = f.readlines()
        for i in range(len(lines)):
            if "Fitting point charges to electrostatic potential" in lines[i]:
                to_begin = int(i)
            if " Sum of ESP charges =" in lines[i]:
                to_end = int(i)
        charges = lines[to_begin + 4 : to_end]
        charge_list = []
        for i in range(len(charges)):
            charge_list.append(charges[i].strip().split())
        charge_list_value = []
        atom_list = []
        for i in range(len(charge_list)):
            charge_list_value.append(charge_list[i][2])
            atom_list.append(charge_list[i][1])
        data_tuples = list(zip(atom_list, charge_list_value))
        df_charge = pd.DataFrame(data_tuples, columns=["Atom", "Charge"])
        df_charge.to_csv(
            self.charge_parameter_file, index=False, header=False, sep=" ",
        )

    def get_proper_dihedrals(self):
        """
        Saves proper dihedral angles of the ligand in a text file.
        """
        ppdb = PandasPdb()
        ppdb.read_pdb(self.guest_pdb)
        no_atoms = len(ppdb.df["ATOM"])
        atom_index_list = []
        for i in range(no_atoms):
            atom_index_list.append(i + 1)
        possible_dihedrals = []
        for dihed in itertools.permutations(atom_index_list, 4):
            possible_dihedrals.append(dihed)
        df_bonds = pd.read_csv(
            self.bond_parameter_file, header=None, delimiter=r"\s+"
        )
        df_bonds.columns = [
            "bond",
            "k_bond",
            "bond_length",
            "bond_1",
            "bond_2",
        ]
        bond1 = df_bonds["bond_1"].values.tolist()
        bond2 = df_bonds["bond_2"].values.tolist()
        bond_list_list = []
        for i in range(len(bond1)):
            args = (bond1[i], bond2[i])
            bond_list_list.append(list(args))
        reverse_bond_list_list = []
        for bonds in bond_list_list:
            reverse_bond_list_list.append(base.reverse_list(bonds))
        bond_lists = bond_list_list + reverse_bond_list_list
        proper_dihed_repeated = []
        for i in range(len(possible_dihedrals)):
            dihed_frag = (
                [possible_dihedrals[i][0], possible_dihedrals[i][1]],
                [possible_dihedrals[i][1], possible_dihedrals[i][2]],
                [possible_dihedrals[i][2], possible_dihedrals[i][3]],
            )
            a = [
                dihed_frag[0] in bond_lists,
                dihed_frag[1] in bond_lists,
                dihed_frag[2] in bond_lists,
            ]
            if a == [True, True, True]:
                proper_dihed_repeated.append(possible_dihedrals[i])

        len_repeated_dihed_list = len(proper_dihed_repeated)
        proper_dihedrals = proper_dihed_repeated
        for x in proper_dihedrals:
            z = x[::-1]
            if z in proper_dihedrals:
                proper_dihedrals.remove(z)
        len_non_repeated_dihed_list = len(proper_dihedrals)
        # print(len_repeated_dihed_list == len_non_repeated_dihed_list * 2)
        np.savetxt(self.proper_dihedral_file, proper_dihedrals, fmt="%s")
        # return(proper_dihedrals)

class ParameterizeHost:

    """
    A class used to obtain force field parameters for the QM region
    of the receptor (bond, angle and charge parameters) from QM
    calculations.

    This class contain methods to process the output files of the
    Gaussian QM output files (.chk, .fchk and .log files). Methods
    in the class extract the unprocessed hessian matrix from the
    Gaussian QM calculations, processes it and uses the Modified
    Seminario Method to ontain the bond and angle parameters. The
    class also extracts the QM charges from the log file.

    ...

    Attributes
    ----------
    xyz_file: str, optional
        XYZ file for ligand coordinates obtained from its corresponding
        formatted checkpoint file.

    coordinate_file: str, optional
        Text file containing the receptor coordinates (extracted
        from the formatted checkpoint file).

    unprocessed_hessian_file: str, optional
        Unprocessed hessian matrix of the receptor obtained from the
        formatted checkpoint file.

    bond_list_file: str, optional
        Text file containing the bond information of the receptor
        extracted from the log file.

    angle_list_file: str, optional
        Text file containing the angle information of the receptor
        extracted from the log file.

    hessian_file: str, optional
        Processed hessian matrix of the receptor.

    atom_names_file: str, optional
        Text file containing the list of atom names from the fchk file.

    bond_parameter_file: str, optional
        Text file containing the bond parameters for the receptor
        obtained using the Modified Seminario method.

    angle_parameter_file: str, optional
        Text file containing the angle parameters of the receptor.

    charge_parameter_file: str, optional
        Text file containing the QM charges of the receptor.

    host_qm_pdb: str, optional
        PDB file for the receptor's QM region.

    functional: str, optional
        Exchange/Correlation or hybrid functional to use in the Gaussian
        QM calculation.

    basis_set: str, optional
        Basis set to use for the Gaussian QM calculation.

    """

    def __init__(
        self,
        xyz_file="host_qm_coords.xyz",
        coordinate_file="host_qm_coordinates.txt",
        unprocessed_hessian_file="host_qm_unprocessed_hessian.txt",
        bond_list_file="host_qm_bond_list.txt",
        angle_list_file="host_qm_angle_list.txt",
        hessian_file="host_qm_hessian.txt",
        atom_names_file="host_qm_atom_names.txt",
        bond_parameter_file="host_qm_bonds.txt",
        angle_parameter_file="host_qm_angles.txt",
        charge_parameter_file="host_qm_surround_charges.txt",
        host_qm_pdb="host_qm.pdb",
        functional="B3LYP",
        basis_set="6-31G",
    ):

        self.xyz_file = xyz_file
        self.coordinate_file = coordinate_file
        self.unprocessed_hessian_file = unprocessed_hessian_file
        self.bond_list_file = bond_list_file
        self.angle_list_file = angle_list_file
        self.hessian_file = hessian_file
        self.atom_names_file = atom_names_file
        self.bond_parameter_file = bond_parameter_file
        self.angle_parameter_file = angle_parameter_file
        self.charge_parameter_file = charge_parameter_file
        self.host_qm_pdb = host_qm_pdb
        self.functional = functional
        self.basis_set = basis_set

    def get_xyz(self):
        """
        Saves XYZ file from the formatted checkpoint file.
        """
        fchk_file = self.host_qm_pdb[:-4] + ".fchk"
        with open(fchk_file, "r") as f:
            lines = f.readlines()
        for i in range(len(lines)):
            if "Current cartesian coordinates" in lines[i]:
                no_coordinates = re.findall(r"\d+|\d+.\d+", lines[i])
                no_coordinates = int(no_coordinates[0])
        for i in range(len(lines)):
            if "Current cartesian coordinates" in lines[i]:
                to_begin = int(i)
        cartesian_coords = lines[
            to_begin + 1 : to_begin + 1 + int(math.ceil(no_coordinates / 5))
        ]
        cartesian_list = []
        for i in range(len(cartesian_coords)):
            cartesian_list.append(cartesian_coords[i].strip().split())
        coordinates_list = [
            item for sublist in cartesian_list for item in sublist
        ]
        list_coords = [float(x) * float(0.529) for x in coordinates_list]
        for i in range(len(lines)):
            if "Atomic numbers" in lines[i]:
                to_begin = int(i)
            if "Nuclear charges" in lines[i]:
                to_end = int(i)
        atomic_numbers = lines[to_begin + 1 : to_end]
        atom_numbers = []
        for i in range(len(atomic_numbers)):
            atom_numbers.append(atomic_numbers[i].strip().split())
        numbers = [item for sublist in atom_numbers for item in sublist]
        N = int(no_coordinates / 3)
        # Opens the new xyz file
        file = open(self.xyz_file, "w")
        file.write(str(N) + "\n \n")
        coords = np.zeros((N, 3))
        n = 0
        names = []
        # Gives name for atomic number
        for x in range(0, len(numbers)):
            names.append(const.element_list[int(numbers[x]) - 1][1])
        # Print coordinates to new input_coords.xyz file
        for i in range(0, N):
            for j in range(0, 3):
                coords[i][j] = list_coords[n]
                n = n + 1
            file.write(
                names[i]
                + str(round(coords[i][0], 3))
                + " "
                + str(round(coords[i][1], 3))
                + " "
                + str(round(coords[i][2], 3))
                + "\n"
            )
        file.close()
        np.savetxt(self.coordinate_file, coords, fmt="%s")

    def get_unprocessed_hessian(self):
        """
        Saves a text file of the unprocessed hessian matrix from the
        formatted checkpoint file.
        """
        fchk_file = self.host_qm_pdb[:-4] + ".fchk"
        with open(fchk_file, "r") as f:
            lines = f.readlines()
        for i in range(len(lines)):
            if "Cartesian Force Constants" in lines[i]:
                no_hessian = re.findall(r"\d+|\d+.\d+", lines[i])
                no_hessian = int(no_hessian[0])
        for i in range(len(lines)):
            if "Cartesian Force Constants" in lines[i]:
                to_begin = int(i)
        hessian = lines[
            to_begin + 1 : to_begin + 1 + int(math.ceil(no_hessian / 5))
        ]
        hessian_list = []
        for i in range(len(hessian)):
            hessian_list.append(hessian[i].strip().split())
        unprocessed_Hessian = [
            item for sublist in hessian_list for item in sublist
        ]
        np.savetxt(
            self.unprocessed_hessian_file, unprocessed_Hessian, fmt="%s",
        )

    def get_bond_angles(self):
        """
        Saves a text file containing bonds and angles from the gaussian
        log file.
        """
        log_file = self.host_qm_pdb[:-4] + ".log"
        fid = open(log_file, "r")
        tline = fid.readline()
        bond_list = []
        angle_list = []
        n = 1
        n_bond = 1
        n_angle = 1
        tmp = "R"  # States if bond or angle
        B = []
        # Finds the bond and angles from the .log file
        while tline:
            tline = fid.readline()
            # Line starts at point when bond and angle list occurs
            if (
                len(tline) > 80
                and tline[0:81].strip()
                == "! Name  Definition              Value          Derivative Info.                !"
            ):
                tline = fid.readline()
                tline = fid.readline()
                # Stops when all bond and angles recorded
                while (tmp[0] == "R") or (tmp[0] == "A"):
                    line = tline.split()
                    tmp = line[1]
                    # Bond or angles listed as string
                    list_terms = line[2][2:-1]
                    # Bond List
                    if tmp[0] == "R":
                        x = list_terms.split(",")
                        # Subtraction due to python array indexing at 0
                        x = [(int(i) - 1) for i in x]
                        bond_list.append(x)
                        # Angle List
                    if tmp[0] == "A":
                        x = list_terms.split(",")
                        # Subtraction due to python array indexing at 0
                        x = [(int(i) - 1) for i in x]
                        angle_list.append(x)
                    tline = fid.readline()
                # Leave loop
                tline = -1
        np.savetxt(self.bond_list_file, bond_list, fmt="%s")
        np.savetxt(self.angle_list_file, angle_list, fmt="%s")

    def get_hessian(self):
        """
        Extracts hessian matrix from the unprocessed hessian matrix
        and saves into a new file.
        """
        unprocessed_Hessian = np.loadtxt(self.unprocessed_hessian_file)
        fchk_file = self.host_qm_pdb[:-4] + ".fchk"
        with open(fchk_file, "r") as f:
            lines = f.readlines()
        for i in range(len(lines)):
            if "Current cartesian coordinates" in lines[i]:
                no_coordinates = re.findall(r"\d+|\d+.\d+", lines[i])
                no_coordinates = int(no_coordinates[0])
        N = int(no_coordinates / 3)
        length_hessian = 3 * N
        hessian = np.zeros((length_hessian, length_hessian))
        m = 0
        # Write the hessian in a 2D array format
        for i in range(0, (length_hessian)):
            for j in range(0, (i + 1)):
                hessian[i][j] = unprocessed_Hessian[m]
                hessian[j][i] = unprocessed_Hessian[m]
                m = m + 1
        hessian = (hessian * (627.509391)) / (
            0.529 ** 2
        )  # Change from Hartree/bohr to kcal/mol/ang
        np.savetxt(self.hessian_file, hessian, fmt="%s")

    def get_atom_names(self):
        """
        Saves a list of atom names from the formatted checkpoint file.
        """
        fchk_file = self.host_qm_pdb[:-4] + ".fchk"
        with open(fchk_file, "r") as f:
            lines = f.readlines()
        for i in range(len(lines)):
            if "Atomic numbers" in lines[i]:
                to_begin = int(i)
            if "Nuclear charges" in lines[i]:
                to_end = int(i)
        atomic_numbers = lines[to_begin + 1 : to_end]
        atom_numbers = []
        for i in range(len(atomic_numbers)):
            atom_numbers.append(atomic_numbers[i].strip().split())
        numbers = [item for sublist in atom_numbers for item in sublist]
        names = []
        # Gives name for atomic number
        for x in range(0, len(numbers)):
            names.append(const.element_list[int(numbers[x]) - 1][1])
        atom_names = []
        for i in range(0, len(names)):
            atom_names.append(names[i].strip() + str(i + 1))
        np.savetxt(self.atom_names_file, atom_names, fmt="%s")

    def get_bond_angle_params(self):
        """
        Saves the bond and angle parameter files obtained from
        the formatted checkpoint file.
        """
        fchk_file = self.host_qm_pdb[:-4] + ".fchk"
        with open(fchk_file, "r") as f:
            lines = f.readlines()
        for i in range(len(lines)):
            if "Current cartesian coordinates" in lines[i]:
                no_coordinates = re.findall(r"\d+|\d+.\d+", lines[i])
                no_coordinates = int(no_coordinates[0])
        N = int(no_coordinates / 3)
        coords = np.loadtxt(self.coordinate_file)
        hessian = np.loadtxt(self.hessian_file)
        bond_list = np.loadtxt(self.bond_list_file, dtype=int)
        atom_names = np.loadtxt(self.atom_names_file, dtype=str)
        # Find bond lengths
        bond_lengths = np.zeros((N, N))
        for i in range(0, N):
            for j in range(0, N):
                diff_i_j = np.array(coords[i, :]) - np.array(coords[j, :])
                bond_lengths[i][j] = np.linalg.norm(diff_i_j)
        eigenvectors = np.empty((3, 3, N, N), dtype=complex)
        eigenvalues = np.empty((N, N, 3), dtype=complex)
        partial_hessian = np.zeros((3, 3))
        for i in range(0, N):
            for j in range(0, N):
                partial_hessian = hessian[
                    (i * 3) : ((i + 1) * 3), (j * 3) : ((j + 1) * 3)
                ]
                [a, b] = np.linalg.eig(partial_hessian)
                eigenvalues[i, j, :] = a
                eigenvectors[:, :, i, j] = b
        # Modified Seminario method to find the bond parameters
        # and print them to file
        file_bond = open(self.bond_parameter_file, "w")
        k_b = np.zeros(len(bond_list))
        bond_length_list = np.zeros(len(bond_list))
        unique_values_bonds = []  # Used to find average values
        for i in range(0, len(bond_list)):
            AB = modified_seminario.force_constant_bond(
                bond_list[i][0],
                bond_list[i][1],
                eigenvalues,
                eigenvectors,
                coords,
            )
            BA = modified_seminario.force_constant_bond(
                bond_list[i][1],
                bond_list[i][0],
                eigenvalues,
                eigenvectors,
                coords,
            )
            # Order of bonds sometimes causes slight differences,
            # find the mean
            k_b[i] = np.real((AB + BA) / 2)
            # Vibrational_scaling takes into account DFT deficities
            # / anharmocity
            vibrational_scaling = const.get_vibrational_scaling(
                functional=self.functional, basis_set=self.basis_set
            )
            vibrational_scaling_squared = vibrational_scaling ** 2
            k_b[i] = k_b[i] * vibrational_scaling_squared
            bond_length_list[i] = bond_lengths[bond_list[i][0]][
                bond_list[i][1]
            ]
            file_bond.write(
                atom_names[bond_list[i][0]]
                + "-"
                + atom_names[bond_list[i][1]]
                + "  "
            )
            file_bond.write(
                str("%#.5g" % k_b[i])
                + "   "
                + str("%#.4g" % bond_length_list[i])
                + "   "
                + str(bond_list[i][0] + 1)
                + "   "
                + str(bond_list[i][1] + 1)
            )
            file_bond.write("\n")
            unique_values_bonds.append(
                [
                    atom_names[bond_list[i][0]],
                    atom_names[bond_list[i][1]],
                    k_b[i],
                    bond_length_list[i],
                    1,
                ]
            )
        file_bond.close()
        angle_list = np.loadtxt(self.angle_list_file, dtype=int)
        # Modified Seminario method to find the angle parameters
        # and print them to file
        file_angle = open(self.angle_parameter_file, "w")
        k_theta = np.zeros(len(angle_list))
        theta_0 = np.zeros(len(angle_list))
        unique_values_angles = []  # Used to find average values
        # Modified Seminario part goes here ...
        # Connectivity information for Modified Seminario Method
        central_atoms_angles = []
        # A structure is created with the index giving the central
        # atom of the angle, an array then lists the angles with
        # that central atom.
        # i.e. central_atoms_angles{3} contains an array of angles
        # with central atom 3
        for i in range(0, len(coords)):
            central_atoms_angles.append([])
            for j in range(0, len(angle_list)):
                if i == angle_list[j][1]:
                    # For angle ABC, atoms A C are written to array
                    AC_array = [angle_list[j][0], angle_list[j][2], j]
                    central_atoms_angles[i].append(AC_array)
                    # For angle ABC, atoms C A are written to array
                    CA_array = [angle_list[j][2], angle_list[j][0], j]
                    central_atoms_angles[i].append(CA_array)
        # Sort rows by atom number
        for i in range(0, len(coords)):
            central_atoms_angles[i] = sorted(
                central_atoms_angles[i], key=itemgetter(0)
            )
        # Find normals u_PA for each angle
        unit_PA_all_angles = []
        for i in range(0, len(central_atoms_angles)):
            unit_PA_all_angles.append([])
            for j in range(0, len(central_atoms_angles[i])):
                # For the angle at central_atoms_angles[i][j,:] the corresponding
                # u_PA value is found for the plane ABC and bond AB,
                # where ABC corresponds to the order of the arguements.
                # This is why the reverse order was also added
                unit_PA_all_angles[i].append(
                    linear_algebra.u_PA_from_angles(
                        central_atoms_angles[i][j][0],
                        i,
                        central_atoms_angles[i][j][1],
                        coords,
                    )
                )
        # Finds the contributing factors from the other angle terms
        # scaling_factor_all_angles = cell(max(max(angle_list)));
        # This array will contain scaling factor and angle list position
        scaling_factor_all_angles = []
        for i in range(0, len(central_atoms_angles)):
            scaling_factor_all_angles.append([])
            for j in range(0, len(central_atoms_angles[i])):
                n = 1
                m = 1
                angles_around = 0
                additional_contributions = 0
                scaling_factor_all_angles[i].append([0, 0])
                # Position in angle list
                scaling_factor_all_angles[i][j][1] = central_atoms_angles[i][
                    j
                ][2]
                # Goes through the list of angles with the same central
                # atom and computes the term need for the modified Seminario method
                # Forwards directions, finds the same bonds with the central atom i
                while (
                    ((j + n) < len(central_atoms_angles[i]))
                    and central_atoms_angles[i][j][0]
                    == central_atoms_angles[i][j + n][0]
                ):
                    additional_contributions = (
                        additional_contributions
                        + (
                            abs(
                                np.dot(
                                    unit_PA_all_angles[i][j][:],
                                    unit_PA_all_angles[i][j + n][:],
                                )
                            )
                        )
                        ** 2
                    )
                    n = n + 1
                    angles_around = angles_around + 1
                # Backwards direction, finds the same bonds with the central atom i
                while ((j - m) >= 0) and central_atoms_angles[i][j][
                    0
                ] == central_atoms_angles[i][j - m][0]:
                    additional_contributions = (
                        additional_contributions
                        + (
                            abs(
                                np.dot(
                                    unit_PA_all_angles[i][j][:],
                                    unit_PA_all_angles[i][j - m][:],
                                )
                            )
                        )
                        ** 2
                    )
                    m = m + 1
                    angles_around = angles_around + 1
                if n != 1 or m != 1:
                    # Finds the mean value of the additional contribution to
                    # change to normal Seminario method comment out + part
                    scaling_factor_all_angles[i][j][0] = 1 + (
                        additional_contributions / (m + n - 2)
                    )
                else:
                    scaling_factor_all_angles[i][j][0] = 1
        scaling_factors_angles_list = []
        for i in range(0, len(angle_list)):
            scaling_factors_angles_list.append([])
        # Orders the scaling factors according to the angle list
        for i in range(0, len(central_atoms_angles)):
            for j in range(0, len(central_atoms_angles[i])):
                scaling_factors_angles_list[
                    scaling_factor_all_angles[i][j][1]
                ].append(scaling_factor_all_angles[i][j][0])
        # Finds the angle force constants with the scaling factors
        # included for each angle
        for i in range(0, len(angle_list)):
            # Ensures that there is no difference when the
            # ordering is changed
            [AB_k_theta, AB_theta_0] = modified_seminario.force_angle_constant(
                angle_list[i][0],
                angle_list[i][1],
                angle_list[i][2],
                bond_lengths,
                eigenvalues,
                eigenvectors,
                coords,
                scaling_factors_angles_list[i][0],
                scaling_factors_angles_list[i][1],
            )
            [BA_k_theta, BA_theta_0] = modified_seminario.force_angle_constant(
                angle_list[i][2],
                angle_list[i][1],
                angle_list[i][0],
                bond_lengths,
                eigenvalues,
                eigenvectors,
                coords,
                scaling_factors_angles_list[i][1],
                scaling_factors_angles_list[i][0],
            )
            k_theta[i] = (AB_k_theta + BA_k_theta) / 2
            theta_0[i] = (AB_theta_0 + BA_theta_0) / 2
            # Vibrational_scaling takes into account DFT
            # deficities / anharmonicity
            k_theta[i] = k_theta[i] * vibrational_scaling_squared
            file_angle.write(
                atom_names[angle_list[i][0]]
                + "-"
                + atom_names[angle_list[i][1]]
                + "-"
                + atom_names[angle_list[i][2]]
                + "  "
            )
            file_angle.write(
                str("%#.4g" % k_theta[i])
                + "   "
                + str("%#.4g" % theta_0[i])
                + "   "
                + str(angle_list[i][0] + 1)
                + "   "
                + str(angle_list[i][1] + 1)
                + "   "
                + str(angle_list[i][2] + 1)
            )
            file_angle.write("\n")
            unique_values_angles.append(
                [
                    atom_names[angle_list[i][0]],
                    atom_names[angle_list[i][1]],
                    atom_names[angle_list[i][2]],
                    k_theta[i],
                    theta_0[i],
                    1,
                ]
            )
        file_angle.close()

    def get_charges(self):
        """
        Saves the atomic charges in a text file obtained from
        the Gaussian log file.
        """
        log_file = self.host_qm_pdb[:-4] + ".log"
        with open(log_file, "r") as f:
            lines = f.readlines()
        for i in range(len(lines)):
            if "Fitting point charges to electrostatic potential" in lines[i]:
                to_begin = int(i)
            if " Sum of ESP charges =" in lines[i]:
                to_end = int(i)
        charges = lines[to_begin + 4 : to_end]
        charge_list = []
        for i in range(len(charges)):
            charge_list.append(charges[i].strip().split())
        charge_list_value = []
        atom_list = []
        for i in range(len(charge_list)):
            charge_list_value.append(charge_list[i][2])
            atom_list.append(charge_list[i][1])
        data_tuples = list(zip(atom_list, charge_list_value))
        df_charge = pd.DataFrame(data_tuples, columns=["Atom", "Charge"])
        df_charge.to_csv(
            self.charge_parameter_file, index=False, header=False, sep=" ",
        )



