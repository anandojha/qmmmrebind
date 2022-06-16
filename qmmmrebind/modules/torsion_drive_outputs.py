"""

"""

import re
import os

from biopandas.pdb import PandasPdb
from scipy import optimize
import numpy as np
import scipy

import modules.constants as const
import modules.base as base
import modules.file_utilities as file_utilities
import modules.file_modify as file_modify
import modules.openmm_utilities as openmm_utilities
import modules.torsion_drive_inputs as torsion_inputs

class TorsionDriveParams:

    """

    A class used to parameterize the torsional parameters
    of the ligand by fitting the torsional parameters obtained
    from torsiondrive calculations.

    Previously obtained reparameterized XML forcefield file did
    not have the torsional parameters obtained from QM calculations.
    The torsional parameters obtained from torsiondrive scans are
    fitted and a new XML forcefield file is generated.

    ...

    Attributes
    ----------
    num_charge_atoms : int, optional
        Number of charged atoms in the molecule.

    index_charge_atom_1: int, optional
        Index of the first charged atom.

    charge_atom_1 : int, optional
        Charge on the first charged atom.

    tor_dir : str, optional
        Torsiondrive directory containing separate torsiondrive folders,
        each containing files for a separate torsiondrive calculation
        for a particular dihedral angle.

    reparameterized_torsional_params_file : str, optional
        Text file containing the forcefield parameters for the
        ligand previously obtained without torsional reparameterization.

    psi_input_file : str, optional
        Input file for psi4 QM engine.

    xyz_file : str, optional
        XYZ file for ligand coordinates.

    coords_file : str, optional
        Text file containing the XYZ coordinates of the ligand.

    template_pdb: str, optional
        Ligand PDB with atoms beginning from 1 to be used as a template PDB
        to retrieve atom indices and symbols.

    system_pdb: str, optional
        PDB file for the torsiondrive torsion scans

    system_sdf : str, optional
        Maximum number of geometry optimization steps.

    system_xml : str, optional
        XML force field file for the ligand.

    qm_scan_file : str, optional
        Output scan file for the torsiondrive scans.

    load_topology : str, optional
        Argument to specify how to load the topology. Can either
        be "openmm" or "parmed".

    method : str, optional
        Minimization method for fitting of torsional
        parameters.

    dihedral_text_file : str, optional
        Dihedral information file for torsiondrive.

    system_init_sdf : str, optional
        Ligand SDF (structure-data) format file. This file will be generated
        only if the ligand is charged.

    reparameterised_system_xml_file : str, optional
        Reparameterized force field XML file obtained using
        openforcefield without torsional reparamaterization.

    reparameterised_torsional_system_xml_file : str, optional
        XML force field file for the ligand obtained with
        torsional reparamaterization.

    """
    
    def __init__(
        self,
        # TODO: some of these variables are ints, and should be initialized as ints
        num_charge_atoms="",
        index_charge_atom_1="",
        charge_atom_1="",
        tor_dir="torsion_dir",
        reparameterized_torsional_params_file="reparameterized_torsional_params.txt",
        psi_input_file="torsion_drive_input.dat",
        xyz_file="torsion_drive_input.xyz",
        coords_file="torsion_drive_input.txt",
        template_pdb="guest_init_II.pdb",
        system_pdb="torsion_drive_input.pdb",
        system_sdf="torsion_drive_input.sdf",
        system_xml="torsion_drive_input.xml",
        qm_scan_file="scan.xyz",
        load_topology="openmm",
        method="L-BFGS-B",
        dihedral_text_file="dihedrals.txt",
        system_init_sdf="torsion_drive_input_init.sdf",
        reparameterised_system_xml_file="guest_reparameterised.xml",
        reparameterised_torsional_system_xml_file="guest_torsional_reparameterized.xml",
    ):

        self.num_charge_atoms = num_charge_atoms
        self.index_charge_atom_1 = index_charge_atom_1
        self.charge_atom_1 = charge_atom_1
        self.tor_dir = tor_dir
        self.reparameterized_torsional_params_file = (
            reparameterized_torsional_params_file
        )
        self.psi_input_file = psi_input_file
        self.xyz_file = xyz_file
        self.coords_file = coords_file
        self.template_pdb = template_pdb
        self.system_pdb = system_pdb
        self.system_sdf = system_sdf
        self.system_xml = system_xml
        self.qm_scan_file = qm_scan_file
        self.method = method
        self.dihedral_text_file = dihedral_text_file
        self.system_init_sdf = system_init_sdf
        self.load_topology = load_topology
        self.reparameterised_system_xml_file = reparameterised_system_xml_file
        self.reparameterised_torsional_system_xml_file = (
            reparameterised_torsional_system_xml_file
        )

    def write_reparams_torsion_lines(self):
        """
        Saves a text file containing torsional parameters for the ligand
        obtained through openforcefield.
        """
        torsional_parameters_list = []
        parent_cwd = os.getcwd()
        # TODO: use os.path.join
        target_dir = os.path.join(parent_cwd, self.tor_dir)
        # TODO: let's use a more informative variable name than 'i'
        for i in os.listdir(target_dir):
            os.chdir(os.path.join(parent_cwd, self.tor_dir, i))
            if os.path.isfile(self.qm_scan_file):
                print("Entering directory" + " : " + os.getcwd())
                torsion_inputs.torsiondrive_input_to_xyz(
                    psi_input_file=self.psi_input_file, xyz_file=self.xyz_file,
                )
                file_modify.xyz_to_pdb(
                    xyz_file=self.xyz_file,
                    coords_file=self.coords_file,
                    template_pdb=self.template_pdb,
                    system_pdb=self.system_pdb,
                )
                file_modify.generate_xml_from_charged_pdb_sdf(
                    system_pdb=self.system_pdb,
                    system_init_sdf=self.system_init_sdf,
                    system_sdf=self.system_sdf,
                    num_charge_atoms=self.num_charge_atoms,
                    index_charge_atom_1=self.index_charge_atom_1,
                    charge_atom_1=self.charge_atom_1,
                    system_xml=self.system_xml,
                )
                torsional_lines = get_torsional_lines(
                    template_pdb=self.template_pdb,
                    system_xml=self.system_xml,
                    qm_scan_file=self.qm_scan_file,
                    load_topology=self.load_topology,
                    method=self.method,
                    dihedral_text_file=self.dihedral_text_file,
                )
                # print(torsional_lines)
                torsional_parameters_list.append(torsional_lines)
                file_utilities.remove_mm_files(qm_scan_file=self.qm_scan_file)
                os.chdir(parent_cwd)
            else:
                print("Entering directory" + " : " + os.getcwd())
                print(
                    "Torsional Scan file not found, optimization may not \
                     be complete. Existing!!"
                )
                os.chdir(parent_cwd)
        torsional_parameters = [
            item for sublist in torsional_parameters_list for item in sublist
        ]
        with open(self.reparameterized_torsional_params_file, "w") as f:
            for i in torsional_parameters:
                f.write(i + "\n")

    def write_reparams_torsion_lines_charged(self):
        """
        Saves a text file containing torsional parameters for a charged ligand
        obtained through openforcefield.
        """
        torsional_parameters_list = []
        parent_cwd = os.getcwd()
        target_dir = os.path.join(parent_cwd, self.tor_dir)
        for i in os.listdir(target_dir):
            os.chdir(os.path.join(parent_cwd, self.tor_dir, i))
            if os.path.isfile(self.qm_scan_file):
                print("Entering directory" + " : " + os.getcwd())
                torsion_inputs.torsiondrive_input_to_xyz(
                    psi_input_file=self.psi_input_file, xyz_file=self.xyz_file,
                )
                file_modify.xyz_to_pdb(
                    xyz_file=self.xyz_file,
                    coords_file=self.coords_file,
                    template_pdb=self.template_pdb,
                    system_pdb=self.system_pdb,
                )
                file_modify.generate_xml_from_charged_pdb_sdf(
                    system_pdb=self.system_pdb,
                    system_init_sdf=self.system_init_sdf,
                    system_sdf=self.system_sdf,
                    num_charge_atoms=self.num_charge_atoms,
                    index_charge_atom_1=self.index_charge_atom_1,
                    charge_atom_1=self.charge_atom_1,
                    system_xml=self.system_xml,
                )
                torsional_lines = get_torsional_lines(
                    template_pdb=self.template_pdb,
                    system_xml=self.system_xml,
                    qm_scan_file=self.qm_scan_file,
                    load_topology=self.load_topology,
                    method=self.method,
                    dihedral_text_file=self.dihedral_text_file,
                )
                # print(torsional_lines)
                torsional_parameters_list.append(torsional_lines)
                file_utilities.remove_mm_files(qm_scan_file=self.qm_scan_file)
                os.chdir(parent_cwd)
            else:
                print("Entering directory" + " : " + os.getcwd())
                print(
                    "Torsional Scan file not found, optimization may not \
                     be complete. Existing!!"
                )
                os.chdir(parent_cwd)
        torsional_parameters = [
            item for sublist in torsional_parameters_list for item in sublist
        ]
        with open(self.reparameterized_torsional_params_file, "w") as f:
            for i in torsional_parameters:
                f.write(i + "\n")

    def write_torsional_reparams(self):
        """
        Generates a XML force field file for the ligand with reparameterized
        torsional parameters.
        """
        with open(self.reparameterized_torsional_params_file, "r") as xml_tor:
            xml_tor_lines = xml_tor.readlines()
        non_zero_k_tor = []
        for i in xml_tor_lines:
            to_find = "k=" + '"' + "0.0" + '"'
            if to_find not in i:
                non_zero_k_tor.append(i)
        # print(non_zero_k_tor)
        p1 = []
        for i in range(len(non_zero_k_tor)):
            p1.append(int(re.findall("\d*\.?\d+", non_zero_k_tor[i])[2]))
        # print(p1)
        p2 = []
        for i in range(len(non_zero_k_tor)):
            p2.append(int(re.findall("\d*\.?\d+", non_zero_k_tor[i])[4]))
        # print(p2)
        p3 = []
        for i in range(len(non_zero_k_tor)):
            p3.append(int(re.findall("\d*\.?\d+", non_zero_k_tor[i])[6]))
        # print(p3)
        p4 = []
        for i in range(len(non_zero_k_tor)):
            p4.append(int(re.findall("\d*\.?\d+", non_zero_k_tor[i])[8]))
        # print(p4)
        periodicity = []
        for i in range(len(non_zero_k_tor)):
            periodicity.append(
                int(re.findall("\d*\.?\d+", non_zero_k_tor[i])[9])
            )
        # print(periodicity)
        # TODO: there may be a way to consolidate the reparametrization of 
        #  the XML file to obey the DRY principle
        xml_tor_reparams = open(self.reparameterised_system_xml_file, "r")
        xml_tor_reparams_lines = xml_tor_reparams.readlines()
        # A string template and formatting should be used here
        for j in range(len(xml_tor_reparams_lines)):
            for i in range(len(non_zero_k_tor)):
                to_find_tor = (
                    "p1="
                    + '"'
                    + str(p1[i])
                    + '"'
                    + " "
                    + "p2="
                    + '"'
                    + str(p2[i])
                    + '"'
                    + " "
                    + "p3="
                    + '"'
                    + str(p3[i])
                    + '"'
                    + " "
                    + "p4="
                    + '"'
                    + str(p4[i])
                    + '"'
                    + " "
                    + "periodicity="
                    + '"'
                    + str(periodicity[i])
                    + '"'
                )
                if to_find_tor in xml_tor_reparams_lines[j]:
                    # print(xml_tor_reparams_lines[j])
                    xml_tor_reparams_lines[j] = non_zero_k_tor[i]
        with open(self.reparameterised_torsional_system_xml_file, "w") as f:
            for i in xml_tor_reparams_lines:
                f.write(i)

def get_dihedrals(qm_scan_file):

    """
    Returns dihedrals from the torsiondrive scan file.

    Parameters
    ----------
    qm_scan_file : str
        Output scan file containing torsiondrive scans.

    Returns
    -------
    dihedrals : list
        List of all the dihedral values from the qm scan file.

    """
    with open(qm_scan_file, "r") as f:
        lines = f.readlines()
    energy_dihedral_lines = []
    for i in range(len(lines)):
        if "Dihedral" in lines[i]:
            energy_dihedral_lines.append(lines[i])
    dihedrals = []
    for i in energy_dihedral_lines:
        energy_dihedral = i
        energy_dihedral = re.findall(r"[-+]?\d+[.]?\d*", energy_dihedral)
        dihedral = float(energy_dihedral[0])
        dihedrals.append(dihedral)
    return dihedrals


def get_qm_energies(qm_scan_file):

    """
    Returns QM optimized energies from the torsiondrive
    scan file.

    Parameters
    ----------
    qm_scan_file : str
        Output scan file containing torsiondrive scans.

    Returns
    -------
    qm_energies : list
        List of all the qm optimiseed energies extracted from the torsiondrive
        scan file.
    """
    with open(qm_scan_file, "r") as f:
        lines = f.readlines()
    energy_dihedral_lines = []
    for i in range(len(lines)):
        if "Dihedral" in lines[i]:
            energy_dihedral_lines.append(lines[i])
    qm_energies = []
    for i in energy_dihedral_lines:
        energy_dihedral = i
        energy_dihedral = re.findall(r"[-+]?\d+[.]?\d*", energy_dihedral)
        energy = float(energy_dihedral[1])
        qm_energies.append(energy)
    return qm_energies


def generate_mm_pdbs(qm_scan_file, template_pdb):

    """
    Generate PDBs from the torsiondrive scan file
    based on a template PDB.

    """
    with open(qm_scan_file, "r") as f:
        lines = f.readlines()
    energy_dihedral_lines = []
    for i in range(len(lines)):
        if "Dihedral" in lines[i]:
            energy_dihedral_lines.append(lines[i])
    dihedrals = []
    for i in energy_dihedral_lines:
        energy_dihedral = i
        energy_dihedral = re.findall(r"[-+]?\d+[.]?\d*", energy_dihedral)
        dihedral = float(energy_dihedral[0])
        dihedrals.append(dihedral)
    lines_markers = []
    for i in range(len(lines)):
        if "Dihedral" in lines[i]:
            lines_markers.append(i)
    lines_markers.append(len(lines) + 1)
    for i in range(len(lines_markers) - 1):
        # pdb_file_to_write = str(dihedrals[i]) + ".pdb"
        if dihedrals[i] > 0:
            pdb_file_to_write = "plus_" + str(abs(dihedrals[i])) + ".pdb"
        if dihedrals[i] < 0:
            pdb_file_to_write = "minus_" + str(abs(dihedrals[i])) + ".pdb"
        to_begin = lines_markers[i]
        to_end = lines_markers[i + 1]
        lines_to_write = lines[to_begin + 1 : to_end - 1]
        x_coords = []
        y_coords = []
        z_coords = []
        for i in lines_to_write:
            coordinates = i
            coordinates = re.findall(r"[-+]?\d+[.]?\d*", coordinates)
            x = float(coordinates[0])
            y = float(coordinates[1])
            z = float(coordinates[2])
            x_coords.append(x)
            y_coords.append(y)
            z_coords.append(z)
        ppdb = PandasPdb()
        ppdb.read_pdb(template_pdb)
        ppdb.df["ATOM"]["x_coord"] = x_coords
        ppdb.df["ATOM"]["y_coord"] = y_coords
        ppdb.df["ATOM"]["z_coord"] = z_coords
        ppdb.to_pdb(pdb_file_to_write)

def gen_init_guess(qm_scan_file, load_topology, system_xml):

    """
    Initial guess for the torsional parameter.

    Parameters
    ----------
    qm_scan_file : str
        Output scan file containing torsiondrive scans.

    load_topology : {"openmm", "parmed"}
        Argument to speify how to load the topology.

    system_xml : str
        XML force field file for the system.

    Returns
    -------
    k_init_guess : list
        Initial guess for the torsional parameters.

    """
    x = get_dihedrals(qm_scan_file)
    y = base.scale_list(
        list_=openmm_utilities.get_mm_potential_energies(
            qm_scan_file=qm_scan_file,
            load_topology=load_topology,
            system_xml=system_xml,
        )
    )
    init_vals = [0.0, 0.0, 0.0, 0.0]
    k_init_guess, covar = scipy.optimize.curve_fit(
        torsion_inputs.dihedral_energy, x, y, p0=init_vals
    )
    for i in range(len(k_init_guess)):
        if k_init_guess[i] < 0:
            k_init_guess[i] = 0
    return k_init_guess


def objective_function(k_array, x, delta_qm):
    """
    Objective function for the torsional parameter fitting.
    """
    delta_mm = torsion_inputs.dihedral_energy(
        x, k1=k_array[0], k2=k_array[1], k3=k_array[2], k4=k_array[3]
    )
    loss_function = torsion_inputs.error_function(delta_qm, delta_mm)
    return loss_function


def fit_params(qm_scan_file, load_topology, system_xml, method):
    """
    Optimization of the objective function.
    """
    k_guess = gen_init_guess(
        qm_scan_file=qm_scan_file,
        load_topology=load_topology,
        system_xml=system_xml,
    )
    x_data = np.array(get_dihedrals(qm_scan_file))
    delta_qm = np.array(
        base.scale_list(const.list_hartree_kcal(list_=get_qm_energies(qm_scan_file)))
    )
    optimise = scipy.optimize.minimize(
        objective_function,
        k_guess,
        args=(x_data, delta_qm),
        method=method,
        bounds=[(0.00, None), (0.00, None), (0.00, None), (0.00, None),],
    )
    return optimise.x


def get_tor_params(
    qm_scan_file, template_pdb, load_topology, system_xml, method
):
    """
    Returns the fitted torsional parameters.
    """
    qm_e = get_qm_energies(qm_scan_file=qm_scan_file)
    qm_e_kcal = const.list_hartree_kcal(qm_e)
    delta_qm = base.scale_list(qm_e_kcal)
    generate_mm_pdbs(qm_scan_file=qm_scan_file, template_pdb=template_pdb)
    mm_pe_no_torsion_kcal = openmm_utilities.get_mm_potential_energies(
        qm_scan_file=qm_scan_file,
        load_topology=load_topology,
        system_xml=system_xml,
    )
    delta_mm = base.scale_list(mm_pe_no_torsion_kcal)
    opt_param = fit_params(
        qm_scan_file=qm_scan_file,
        load_topology=load_topology,
        system_xml=system_xml,
        method=method,
    )
    return opt_param


def get_torsional_lines(
    template_pdb,
    system_xml,
    qm_scan_file,
    load_topology,
    method,
    dihedral_text_file,
):
    """
    Returns the torsional lines for the XML forcefield file.
    """
    opt_param = get_tor_params(
        qm_scan_file=qm_scan_file,
        template_pdb=template_pdb,
        load_topology=load_topology,
        system_xml=system_xml,
        method=method,
    )
    dihedral_text = open(dihedral_text_file, "r")
    dihedral_text_lines = dihedral_text.readlines()
    atom_numbers = dihedral_text_lines[-1]
    atom_index_from_1 = [
        int(re.findall(r"\d+", atom_numbers)[0]),
        int(re.findall(r"\d+", atom_numbers)[1]),
        int(re.findall(r"\d+", atom_numbers)[2]),
        int(re.findall(r"\d+", atom_numbers)[3]),
    ]
    atom_index = [i - 1 for i in atom_index_from_1]
    atom_index_lines = (
        " "
        + "p1="
        + '"'
        + str(atom_index[0])
        + '"'
        + " "
        + "p2="
        + '"'
        + str(atom_index[1])
        + '"'
        + " "
        + "p3="
        + '"'
        + str(atom_index[2])
        + '"'
        + " "
        + "p4="
        + '"'
        + str(atom_index[3])
        + '"'
        + " "
    )
    tor_lines = []
    for i in range(len(opt_param)):
        line_to_append = (
            "                "
            + "<Torsion "
            + "k="
            + '"'
            + str(round(opt_param[i], 8))
            + '"'
            + atom_index_lines
            + "periodicity="
            + '"'
            + str(i + 1)
            + '"'
            + " "
            + "phase="
            + '"'
            + "0"
            + '"'
            + "/>"
        )
        # print(line_to_append)
        tor_lines.append(line_to_append)
    return tor_lines
