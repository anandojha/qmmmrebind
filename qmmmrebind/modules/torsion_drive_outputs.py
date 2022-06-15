"""

"""

import re

from biopandas.pdb import PandasPdb
import numpy as np
import scipy

import modules.constants as const
import modules.base as base
import modules.openmm_utilities as openmm_utilities
import modules.torsion_drive_inputs as torsion_inputs

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
