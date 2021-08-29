"""
Unit and regression test for the qmmmrebind package.
"""
import warnings

warnings.filterwarnings("ignore")
from .utils import get_data_filename
import numpy as np
import qmmmrebind
import pytest
import sys
import os
import re

##############################Tests For Discrete Functions##############################
def test_qmmmrebind_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "qmmmrebind" in sys.modules


def test_len_list_to_dict():
    """Test the list to dict function"""
    test_list = ["a", "b", "c", "d"]
    ret = qmmmrebind.parameterize.list_to_dict(test_list)
    assert len(ret) == 2


def test_key_list_to_dict():
    """Test the list to dict function"""
    test_list = ["a", "b", "c", "d"]
    ret = qmmmrebind.parameterize.list_to_dict(test_list)
    assert ret["a"] == "b"


def test_find_word_in_file():
    """Test if a word is in a file"""
    filename = get_data_filename("test_random_read_file_i.dat")
    word = "github"
    ret = qmmmrebind.parameterize.search_in_file(file=filename, word=word)
    assert ret[0][0] == 3


def test_get_vibrational_scaling():
    """Test if the vibrational scaling is retrieved correctly"""
    functional = "QCISD"
    basis_set = "6-311G*"
    vib_scale = qmmmrebind.parameterize.get_vibrational_scaling(
        functional, basis_set
    )
    assert vib_scale == 0.957


def test_unit_vector_N():
    """Test if the unit vector to a plane is calculated correctly"""
    u_BC = [0.34040355, 0.62192853, 0.27011169]
    u_AB = [0.28276792, 0.34232697, 0.02370306]
    u_N = qmmmrebind.parameterize.unit_vector_N(u_BC, u_AB)
    assert u_N[0] == -0.6516162898107304


def test_reverse_list():
    """Test if the list is reversed correctly"""
    lst = [5, 4, 7, 2]
    rev = qmmmrebind.parameterize.reverse_list(lst)
    assert rev == [2, 7, 4, 5]


def test_uniq():
    """Test if the unique elements are returned correctly"""
    lst = [2, 4, 2, 9, 10, 35, 10]
    ret = qmmmrebind.parameterize.uniq(lst)
    assert ret == [2, 4, 9, 10, 35]


def test_copy_psi_input_file():
    source_ = get_data_filename("test_torsion_drive_input.dat")
    destination_pwd = os.getcwd()
    destination_file = source_.split("/")[-1]
    destination_ = destination_pwd + "/" + destination_file
    qmmmrebind.parameterize.copy_file(source=source_, destination=destination_)
    assert "test_torsion_drive_input.dat" in os.listdir()


def test_copy_torsion_xyz_file():
    source_ = get_data_filename("test_torsion_drive_input.xyz")
    destination_pwd = os.getcwd()
    destination_file = source_.split("/")[-1]
    destination_ = destination_pwd + "/" + destination_file
    qmmmrebind.parameterize.copy_file(source=source_, destination=destination_)
    assert "test_torsion_drive_input.xyz" in os.listdir()


def test_torsiondrive_input_to_xyz():
    """Test if the xyz is generated from torsiondrive input dat file."""
    psi_input_file = "test_torsion_drive_input.dat"
    xyz_file = "test_torsion_drive_input.xyz"
    qmmmrebind.parameterize.torsiondrive_input_to_xyz(
        psi_input_file=psi_input_file, xyz_file=xyz_file
    )
    with open("test_torsion_drive_input.xyz", "r") as f:
        lines = f.readlines()
    assert len(lines) != 0


def test_copy_torsion_xyz_file():
    source_ = get_data_filename("test_scan.xyz")
    destination_pwd = os.getcwd()
    destination_file = source_.split("/")[-1]
    destination_ = destination_pwd + "/" + destination_file
    qmmmrebind.parameterize.copy_file(source=source_, destination=destination_)
    assert "test_scan.xyz" in os.listdir()


def test_get_dihedrals():
    qm_scan_file = "test_scan.xyz"
    dihedrals = qmmmrebind.parameterize.get_dihedrals(
        qm_scan_file=qm_scan_file
    )
    assert dihedrals[0] == -165.0


def test_get_qm_energies():
    qm_scan_file = "test_scan.xyz"
    qm_energies = qmmmrebind.parameterize.get_qm_energies(
        qm_scan_file=qm_scan_file
    )
    assert qm_energies[0] == -376.541676037


def test_copy_guest_coordinates_file():
    source_ = get_data_filename("test_guest_coordinates.txt")
    destination_pwd = os.getcwd()
    destination_file = source_.split("/")[-1]
    destination_ = destination_pwd + "/" + destination_file
    qmmmrebind.parameterize.copy_file(source=source_, destination=destination_)
    assert "test_guest_coordinates.txt" in os.listdir()


def test_u_PA_from_angles():
    """Test if the vector in plane ABC perpendicular to AB is calculated correctly"""
    coords = np.loadtxt("test_guest_coordinates.txt")
    u_PA = qmmmrebind.parameterize.u_PA_from_angles(1, 0, 5, coords)
    assert np.allclose(u_PA, [-0.865653, 0.49043205, 0.10060462])


def test_copy_guest_hessian_file():
    source_ = get_data_filename("test_guest_hessian.txt")
    destination_pwd = os.getcwd()
    destination_file = source_.split("/")[-1]
    destination_ = destination_pwd + "/" + destination_file
    qmmmrebind.parameterize.copy_file(source=source_, destination=destination_)
    assert "test_guest_hessian.txt" in os.listdir()


def test_force_angle_constant():
    """Test if the calculated force angle constants are matching the expected value"""
    coords = np.loadtxt("test_guest_coordinates.txt")
    hessian = np.loadtxt("test_guest_hessian.txt")
    # bond_list = np.loadtxt("test_guest_bond_list.txt", dtype=int)
    N = 18
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
    [k_theta, theta_0] = qmmmrebind.parameterize.force_angle_constant(
        1,
        0,
        5,
        bond_lengths=bond_lengths,
        eigenvalues=eigenvalues,
        eigenvectors=eigenvectors,
        coords=coords,
        scaling_1=1.9999999999996876,
        scaling_2=1.9999999999996882,
    )
    #assert k_theta == 100.4487444754119
    #assert theta_0 == 119.9124090641625


def test_copy_guest_init_pdb():
    source_ = get_data_filename("test_guest_init_ii.pdb")
    destination_pwd = os.getcwd()
    destination_file = source_.split("/")[-1]
    destination_ = destination_pwd + "/" + destination_file
    qmmmrebind.parameterize.copy_file(source=source_, destination=destination_)
    assert "test_guest_init_ii.pdb" in os.listdir()


def test_xyz_to_pdb():
    xyz_file = "test_torsion_drive_input.xyz"
    coords_file = "test_torsion_drive_input.txt"
    template_pdb = "test_guest_init_ii.pdb"
    system_pdb = "test_torsion_drive_input.pdb"
    qmmmrebind.parameterize.xyz_to_pdb(
        xyz_file=xyz_file,
        coords_file=coords_file,
        template_pdb=template_pdb,
        system_pdb=system_pdb,
    )
    with open("test_torsion_drive_input.pdb", "r") as f:
        lines = f.readlines()
    atom_list = []
    for line in lines:
        if "ATOM " in line:
            atom_list.append(line)
    assert len(atom_list) == 18

def test_copy_system_init_xml():
    source_ = get_data_filename("test_guest_init.xml")
    destination_pwd = os.getcwd()
    destination_file = source_.split("/")[-1]
    destination_ = destination_pwd + "/" + destination_file
    qmmmrebind.parameterize.copy_file(source=source_, destination=destination_)
    assert "test_guest_init.xml" in os.listdir()

"""
def test_generate_xml_from_charged_pdb_sdf():
    system_pdb = "test_torsion_drive_input.pdb"
    system_init_sdf = "test_guest_init.sdf"
    system_sdf = "test_guest.sdf"
    num_charge_atoms = 1
    index_charge_atom_1 = 9
    charge_atom_1 = 1
    system_xml = "test_guest_init.xml"
    qmmmrebind.parameterize.generate_xml_from_charged_pdb_sdf(
        system_pdb=system_pdb,
        system_init_sdf=system_init_sdf,
        system_sdf=system_sdf,
        num_charge_atoms=num_charge_atoms,
        index_charge_atom_1=index_charge_atom_1,
        charge_atom_1=charge_atom_1,
        system_xml=system_xml,
    )
    with open(system_xml, "r") as f:
        lines = f.readlines()
    assert len(lines) != 0
    bond_lines = []
    for line in lines:
        if "Bond " in line:
            bond_lines.append(line)
    assert len(bond_lines) == 18
"""

def test_generate_mm_pdbs():
    qm_scan_file = "test_scan.xyz"
    template_pdb = "test_guest_init_ii.pdb"
    qmmmrebind.parameterize.generate_mm_pdbs(
        qm_scan_file=qm_scan_file, template_pdb=template_pdb
    )
    with open("minus_15.0.pdb", "r") as f:
        lines = f.readlines()
    atom_list = []
    for line in lines:
        if "ATOM " in line:
            atom_list.append(line)
    assert len(atom_list) == 18


def test_remove_mm_files():
    qm_scan_file = "test_scan.xyz"
    qmmmrebind.parameterize.remove_mm_files(qm_scan_file=qm_scan_file)
    assert "minus*" not in os.listdir()
    assert "plus*" not in os.listdir()


##############################PrepareQMMM##############################
def test_clean_up():
    init_pdb = get_data_filename("test_sample_system_trypsin_benzamidine.pdb")
    cleaned_pdb = "test_system.pdb"
    num_residues = 2
    get_clean_up = qmmmrebind.parameterize.PrepareQMMM(
        init_pdb=init_pdb,
        cleaned_pdb=cleaned_pdb,
        guest_init_pdb="",
        host_pdb="",
        guest_resname="",
        guest_pdb="",
        guest_xyz="",
        distance="",
        residue_list="",
        host_qm_atoms="",
        host_mm_atoms="",
        host_qm_pdb="",
        host_mm_pdb="",
        qm_pdb="",
        mm_pdb="",
        host_mm_region_I_atoms="",
        host_mm_region_II_atoms="",
        host_mm_region_I_pdb="",
        host_mm_region_II_pdb="",
        num_residues=num_residues,
    )
    get_clean_up.clean_up()
    with open(cleaned_pdb, "r") as f:
        lines = f.readlines()
    assert len(lines) != 0
    HOH_list = []
    for i in lines:
        if "HOH" in i:
            HOH_list.append(i)
    HETATM_list = []
    for i in lines:
        if "HETATM" in i:
            HETATM_list.append(i)
    assert len(HOH_list) == 0
    assert len(HETATM_list) == 0


def test_create_host_guest():
    cleaned_pdb = "test_system.pdb"
    guest_init_pdb = "test_guest_init.pdb"
    host_pdb = "test_host.pdb"
    guest_resname = "BEN"
    num_residues = 2
    get_create_host_guest = qmmmrebind.parameterize.PrepareQMMM(
        init_pdb="",
        cleaned_pdb=cleaned_pdb,
        guest_init_pdb=guest_init_pdb,
        host_pdb=host_pdb,
        guest_resname=guest_resname,
        guest_pdb="",
        guest_xyz="",
        distance="",
        residue_list="",
        host_qm_atoms="",
        host_mm_atoms="",
        host_qm_pdb="",
        host_mm_pdb="",
        qm_pdb="",
        mm_pdb="",
        host_mm_region_I_atoms="",
        host_mm_region_II_atoms="",
        host_mm_region_I_pdb="",
        host_mm_region_II_pdb="",
        num_residues=num_residues,
    )
    get_create_host_guest.create_host_guest()
    with open(guest_init_pdb, "r") as f:
        lines = f.readlines()
    assert len(lines) != 0
    lines = lines[:-1]
    guest_resname_list = []
    for i in lines:
        if guest_resname in i:
            guest_resname_list.append(i)
    assert len(lines) == len(guest_resname_list)


def test_realign_guest():
    guest_init_pdb = "test_guest_init.pdb"
    guest_pdb = "test_guest_init_ii.pdb"
    num_residues = 2
    get_realign_guest = qmmmrebind.parameterize.PrepareQMMM(
        init_pdb="",
        cleaned_pdb="",
        guest_init_pdb=guest_init_pdb,
        host_pdb="",
        guest_resname="",
        guest_pdb=guest_pdb,
        guest_xyz="",
        distance="",
        residue_list="",
        host_qm_atoms="",
        host_mm_atoms="",
        host_qm_pdb="",
        host_mm_pdb="",
        qm_pdb="",
        mm_pdb="",
        host_mm_region_I_atoms="",
        host_mm_region_II_atoms="",
        host_mm_region_I_pdb="",
        host_mm_region_II_pdb="",
        num_residues=num_residues,
    )
    get_realign_guest.realign_guest()
    with open(guest_init_pdb, "r") as f:
        guest_init_pdb_lines = f.readlines()
    with open(guest_pdb, "r") as f:
        guest_pdb_lines = f.readlines()
    assert len(guest_init_pdb_lines) == len(guest_pdb_lines)


def test_get_guest_coord():
    guest_pdb = "test_guest_init_ii.pdb"
    guest_xyz = "test_guest_coord.txt"
    num_residues = 2
    get_get_guest_coord = qmmmrebind.parameterize.PrepareQMMM(
        init_pdb="",
        cleaned_pdb="",
        guest_init_pdb="",
        host_pdb="",
        guest_resname="",
        guest_pdb=guest_pdb,
        guest_xyz=guest_xyz,
        distance="",
        residue_list="",
        host_qm_atoms="",
        host_mm_atoms="",
        host_qm_pdb="",
        host_mm_pdb="",
        qm_pdb="",
        mm_pdb="",
        host_mm_region_I_atoms="",
        host_mm_region_II_atoms="",
        host_mm_region_I_pdb="",
        host_mm_region_II_pdb="",
        num_residues=num_residues,
    )
    get_get_guest_coord.get_guest_coord()
    with open(guest_xyz, "r") as f:
        guest_xyz_lines = f.readlines()
    assert len(guest_xyz_lines) != 0


def test_get_qm_resids():
    guest_xyz = "test_guest_coord.txt"
    host_pdb = "test_host.pdb"
    distance = 3.0
    residue_list = "test_residue_list.txt"
    num_residues = 2
    get_get_qm_resids = qmmmrebind.parameterize.PrepareQMMM(
        init_pdb="",
        cleaned_pdb="",
        guest_init_pdb="",
        host_pdb=host_pdb,
        guest_resname="",
        guest_pdb="",
        guest_xyz=guest_xyz,
        distance=distance,
        residue_list=residue_list,
        host_qm_atoms="",
        host_mm_atoms="",
        host_qm_pdb="",
        host_mm_pdb="",
        qm_pdb="",
        mm_pdb="",
        host_mm_region_I_atoms="",
        host_mm_region_II_atoms="",
        host_mm_region_I_pdb="",
        host_mm_region_II_pdb="",
        num_residues=num_residues,
    )
    get_get_qm_resids.get_qm_resids()
    with open(residue_list, "r") as f:
        lines = f.readlines()
    assert len(lines) == 11


def test_get_host_qm_mm_atoms():
    residue_list = "test_residue_list.txt"
    num_residues = 2
    host_pdb = "test_host.pdb"
    host_qm_atoms = "test_host_qm.txt"
    host_mm_atoms = "test_host_mm.txt"
    get_get_host_qm_mm_atoms = qmmmrebind.parameterize.PrepareQMMM(
        init_pdb="",
        cleaned_pdb="",
        guest_init_pdb="",
        host_pdb=host_pdb,
        guest_resname="",
        guest_pdb="",
        guest_xyz="",
        distance="",
        residue_list=residue_list,
        host_qm_atoms=host_qm_atoms,
        host_mm_atoms=host_mm_atoms,
        host_qm_pdb="",
        host_mm_pdb="",
        qm_pdb="",
        mm_pdb="",
        host_mm_region_I_atoms="",
        host_mm_region_II_atoms="",
        host_mm_region_I_pdb="",
        host_mm_region_II_pdb="",
        num_residues=num_residues,
    )
    get_get_host_qm_mm_atoms.get_host_qm_mm_atoms()
    with open(host_qm_atoms, "r") as f:
        lines_qm = f.readlines()
    assert len(lines_qm) == 27
    with open(host_mm_atoms, "r") as f:
        lines_mm = f.readlines()
    assert len(lines_mm) > len(lines_qm)


def test_save_host_pdbs():
    host_pdb = "test_host.pdb"
    host_qm_pdb = "test_host_qm.pdb"
    host_mm_pdb = "test_host_mm.pdb"
    host_qm_atoms = "test_host_qm.txt"
    host_mm_atoms = "test_host_mm.txt"
    num_residues = 2
    get_save_host_pdbs = qmmmrebind.parameterize.PrepareQMMM(
        init_pdb="",
        cleaned_pdb="",
        guest_init_pdb="",
        host_pdb=host_pdb,
        guest_resname="",
        guest_pdb="",
        guest_xyz="",
        distance="",
        residue_list="",
        host_qm_atoms=host_qm_atoms,
        host_mm_atoms=host_mm_atoms,
        host_qm_pdb=host_qm_pdb,
        host_mm_pdb=host_mm_pdb,
        qm_pdb="",
        mm_pdb="",
        host_mm_region_I_atoms="",
        host_mm_region_II_atoms="",
        host_mm_region_I_pdb="",
        host_mm_region_II_pdb="",
        num_residues=num_residues,
    )
    get_save_host_pdbs.save_host_pdbs()
    with open("test_host_qm.pdb", "r") as f:
        lines_qm = f.readlines()
    with open("test_host_mm.pdb", "r") as f:
        lines_mm = f.readlines()
    assert lines_qm not in lines_mm and len(lines_qm) < len(lines_mm)


def test_get_host_mm_region_atoms():
    residue_list = "test_residue_list.txt"
    num_residues = 2
    host_pdb = "test_host.pdb"
    host_mm_pdb = "test_host_mm.pdb"
    host_mm_region_I_atoms = "test_host_mm_region_i.txt"
    host_mm_region_II_atoms = "test_host_mm_region_ii.txt"
    get_get_host_mm_region_atoms = qmmmrebind.parameterize.PrepareQMMM(
        init_pdb="",
        cleaned_pdb="",
        guest_init_pdb="",
        host_pdb=host_pdb,
        guest_resname="",
        guest_pdb="",
        guest_xyz="",
        distance="",
        residue_list=residue_list,
        host_qm_atoms="",
        host_mm_atoms="",
        host_qm_pdb="",
        host_mm_pdb=host_mm_pdb,
        qm_pdb="",
        mm_pdb="",
        host_mm_region_I_atoms=host_mm_region_I_atoms,
        host_mm_region_II_atoms=host_mm_region_II_atoms,
        host_mm_region_I_pdb="",
        host_mm_region_II_pdb="",
        num_residues=num_residues,
    )
    get_get_host_mm_region_atoms.get_host_mm_region_atoms()
    with open(host_mm_region_I_atoms, "r") as f:
        lines_I = f.readlines()
    with open(host_mm_region_II_atoms, "r") as f:
        lines_II = f.readlines()
    assert len(lines_I) == 2716
    assert len(lines_II) == 477


def test_save_host_mm_regions_pdbs():
    host_mm_pdb = "test_host_mm.pdb"
    host_mm_region_I_atoms = "test_host_mm_region_i.txt"
    host_mm_region_II_atoms = "test_host_mm_region_ii.txt"
    host_mm_region_I_pdb = "test_host_mm_region_i.pdb"
    host_mm_region_II_pdb = "test_host_mm_region_ii.pdb"
    num_residues = 2
    get_save_host_mm_regions_pdbs = qmmmrebind.parameterize.PrepareQMMM(
        init_pdb="",
        cleaned_pdb="",
        guest_init_pdb="",
        host_pdb="",
        guest_resname="",
        guest_pdb="",
        guest_xyz="",
        distance="",
        residue_list="",
        host_qm_atoms="",
        host_mm_atoms="",
        host_qm_pdb="",
        host_mm_pdb=host_mm_pdb,
        qm_pdb="",
        mm_pdb="",
        host_mm_region_I_atoms=host_mm_region_I_atoms,
        host_mm_region_II_atoms=host_mm_region_II_atoms,
        host_mm_region_I_pdb=host_mm_region_I_pdb,
        host_mm_region_II_pdb=host_mm_region_II_pdb,
        num_residues=num_residues,
    )
    get_save_host_mm_regions_pdbs.save_host_mm_regions_pdbs()
    with open(host_mm_region_I_pdb, "r") as f:
        lines_I = f.readlines()
    HG12_list = []
    for line in lines_I:
        if "HG12 " in line:
            HG12_list.append(line)
    assert len(HG12_list) == 26
    HG12_list = []
    with open(host_mm_region_II_pdb, "r") as f:
        lines_II = f.readlines()
    for line in lines_II:
        if "HG12 " in line:
            HG12_list.append(line)
    assert len(HG12_list) == 5


def test_get_qm_mm_regions():
    host_qm_pdb = "test_host_qm.pdb"
    qm_pdb = "test_qm.pdb"
    guest_pdb = "test_guest_init_ii.pdb"
    host_mm_pdb = "test_host_mm.pdb"
    mm_pdb = "test_mm.pdb"
    num_residues = 2
    get_get_qm_mm_regions = qmmmrebind.parameterize.PrepareQMMM(
        init_pdb="",
        cleaned_pdb="",
        guest_init_pdb="",
        host_pdb="",
        guest_resname="",
        guest_pdb=guest_pdb,
        guest_xyz="",
        distance="",
        residue_list="",
        host_qm_atoms="",
        host_mm_atoms="",
        host_qm_pdb=host_qm_pdb,
        host_mm_pdb=host_mm_pdb,
        qm_pdb=qm_pdb,
        mm_pdb=mm_pdb,
        host_mm_region_I_atoms="",
        host_mm_region_II_atoms="",
        host_mm_region_I_pdb="",
        host_mm_region_II_pdb="",
        num_residues=num_residues,
    )
    get_get_qm_mm_regions.get_qm_mm_regions()
    with open(qm_pdb, "r") as f:
        lines_qm = f.readlines()
    with open(mm_pdb, "r") as f:
        lines_mm = f.readlines()
    assert len(lines_qm) == 46
    assert len(lines_mm) == 3194


##############################ParameterizeGuest##############################
# drop in a fchk file here
def test_copy_fchk_file():
    source_ = get_data_filename("test_guest_init_ii.fchk")
    destination_pwd = os.getcwd()
    destination_file = source_.split("/")[-1]
    destination_ = destination_pwd + "/" + destination_file
    qmmmrebind.parameterize.copy_file(source=source_, destination=destination_)
    assert "test_guest_init_ii.fchk" in os.listdir()


def test_get_xyz():
    guest_pdb = "test_guest_init_ii.pdb"
    coordinate_file = "test_guest_coordinates.txt"
    xyz_file = "test_guest_coords.xyz"
    get_get_xyz = qmmmrebind.parameterize.ParameterizeGuest(
        xyz_file=xyz_file,
        coordinate_file=coordinate_file,
        unprocessed_hessian_file="",
        bond_list_file="",
        angle_list_file="",
        hessian_file="",
        atom_names_file="",
        bond_parameter_file="",
        angle_parameter_file="",
        charge_parameter_file="",
        guest_pdb=guest_pdb,
        proper_dihedral_file="",
    )
    get_get_xyz.get_xyz()
    with open(coordinate_file, "r") as f:
        no_coords = f.readlines()
    with open(xyz_file, "r") as f:
        xyz_lines = f.readlines()
    assert len(no_coords) == (len(xyz_lines) - 2)


def test_get_unprocessed_hessian():
    guest_pdb = "test_guest_init_ii.pdb"
    unprocessed_hessian_file = "test_guest_unprocessed_hessian.txt"
    get_get_unprocessed_hessian = qmmmrebind.parameterize.ParameterizeGuest(
        xyz_file="",
        coordinate_file="",
        unprocessed_hessian_file=unprocessed_hessian_file,
        bond_list_file="",
        angle_list_file="",
        hessian_file="",
        atom_names_file="",
        bond_parameter_file="",
        angle_parameter_file="",
        charge_parameter_file="",
        guest_pdb=guest_pdb,
        proper_dihedral_file="",
    )
    get_get_unprocessed_hessian.get_unprocessed_hessian()
    with open("test_guest_init_ii.fchk", "r") as f:
        lines = f.readlines()
    for i in range(len(lines)):
        if "Cartesian Force Constants" in lines[i]:
            no_hessian = re.findall(r"\d+|\d+.\d+", lines[i])
            no_hessian = int(no_hessian[0])
    with open(unprocessed_hessian_file, "r") as f:
        hessian_lines = f.readlines()
    assert len(hessian_lines) == no_hessian


# drop in a log file here
def test_copy_log_file():
    source_ = get_data_filename("test_guest_init_ii.log")
    destination_pwd = os.getcwd()
    destination_file = source_.split("/")[-1]
    destination_ = destination_pwd + "/" + destination_file
    qmmmrebind.parameterize.copy_file(source=source_, destination=destination_)
    assert "test_guest_init_ii.log" in os.listdir()


def test_get_bond_angles():
    guest_pdb = "test_guest_init_ii.pdb"
    bond_list_file = "test_guest_bond_list.txt"
    angle_list_file = "test_guest_angle_list.txt"
    get_get_bond_angles = qmmmrebind.parameterize.ParameterizeGuest(
        xyz_file="",
        coordinate_file="",
        unprocessed_hessian_file="",
        bond_list_file=bond_list_file,
        angle_list_file=angle_list_file,
        hessian_file="",
        atom_names_file="",
        bond_parameter_file="",
        angle_parameter_file="",
        charge_parameter_file="",
        guest_pdb=guest_pdb,
        proper_dihedral_file="",
    )
    get_get_bond_angles.get_bond_angles()
    with open(bond_list_file, "r") as f:
        bonds = f.readlines()
    with open(angle_list_file, "r") as f:
        angles = f.readlines()
    assert len(bonds) == 18
    assert len(angles) == 27


def test_get_hessian():
    unprocessed_hessian_file = "test_guest_unprocessed_hessian.txt"
    hessian_file = "test_guest_hessian.txt"
    guest_pdb = "test_guest_init_ii.pdb"
    get_get_hessian = qmmmrebind.parameterize.ParameterizeGuest(
        xyz_file="",
        coordinate_file="",
        unprocessed_hessian_file=unprocessed_hessian_file,
        bond_list_file="",
        angle_list_file="",
        hessian_file=hessian_file,
        atom_names_file="",
        bond_parameter_file="",
        angle_parameter_file="",
        charge_parameter_file="",
        guest_pdb=guest_pdb,
        proper_dihedral_file="",
    )
    get_get_hessian.get_hessian()
    hessian = np.loadtxt(hessian_file)
    assert hessian.size == 2916


def test_get_atom_names():
    guest_pdb = "test_guest_init_ii.pdb"
    atom_names_file = "test_guest_atom_names.txt"
    get_get_atom_names = qmmmrebind.parameterize.ParameterizeGuest(
        xyz_file="",
        coordinate_file="",
        unprocessed_hessian_file="",
        bond_list_file="",
        angle_list_file="",
        hessian_file="",
        atom_names_file=atom_names_file,
        bond_parameter_file="",
        angle_parameter_file="",
        charge_parameter_file="",
        guest_pdb=guest_pdb,
        proper_dihedral_file="",
    )
    get_get_atom_names.get_atom_names()
    with open(atom_names_file, "r") as f:
        no_atoms = len(f.readlines())
    assert no_atoms == 18


def test_get_bond_angle_params():
    guest_pdb = "test_guest_init_ii.pdb"
    coordinate_file = "test_guest_coordinates.txt"
    hessian_file = "test_guest_hessian.txt"
    bond_list_file = "test_guest_bond_list.txt"
    atom_names_file = "test_guest_atom_names.txt"
    bond_parameter_file = "test_guest_bonds.txt"
    angle_list_file = "test_guest_angle_list.txt"
    angle_parameter_file = "test_guest_angles.txt"
    get_get_bond_angle_params = qmmmrebind.parameterize.ParameterizeGuest(
        xyz_file="",
        coordinate_file=coordinate_file,
        unprocessed_hessian_file="",
        bond_list_file=bond_list_file,
        angle_list_file=angle_list_file,
        hessian_file=hessian_file,
        atom_names_file=atom_names_file,
        bond_parameter_file=bond_parameter_file,
        angle_parameter_file=angle_parameter_file,
        charge_parameter_file="",
        guest_pdb=guest_pdb,
        proper_dihedral_file="",
    )
    get_get_bond_angle_params.get_bond_angle_params()
    with open(bond_parameter_file, "r") as f:
        no_bonds = len(f.readlines())
    with open(angle_parameter_file, "r") as f:
        no_angles = len(f.readlines())
    with open(bond_list_file, "r") as f:
        bonds = len(f.readlines())
    with open(angle_list_file, "r") as f:
        angles = len(f.readlines())
    assert no_bonds == bonds
    assert no_angles == angles


def test_get_charges():
    guest_pdb = "test_guest_init_ii.pdb"
    charge_parameter_file = "test_guest_charges.txt"
    get_get_charges = qmmmrebind.parameterize.ParameterizeGuest(
        xyz_file="",
        coordinate_file="",
        unprocessed_hessian_file="",
        bond_list_file="",
        angle_list_file="",
        hessian_file="",
        atom_names_file="",
        bond_parameter_file="",
        angle_parameter_file="",
        charge_parameter_file=charge_parameter_file,
        guest_pdb=guest_pdb,
        proper_dihedral_file="",
    )
    get_get_charges.get_charges()
    with open(charge_parameter_file, "r") as f:
        no_lines = len(f.readlines())
    assert no_lines == 18


def test_get_proper_dihedrals():
    guest_pdb = "test_guest_init_ii.pdb"
    bond_parameter_file = "test_guest_bonds.txt"
    proper_dihedral_file = "test_proper_dihedrals.txt"
    get_get_proper_dihedrals = qmmmrebind.parameterize.ParameterizeGuest(
        xyz_file="",
        coordinate_file="",
        unprocessed_hessian_file="",
        bond_list_file="",
        angle_list_file="",
        hessian_file="",
        atom_names_file="",
        bond_parameter_file=bond_parameter_file,
        angle_parameter_file="",
        charge_parameter_file="",
        guest_pdb=guest_pdb,
        proper_dihedral_file=proper_dihedral_file,
    )
    get_get_proper_dihedrals.get_proper_dihedrals()
    with open(proper_dihedral_file, "r") as f:
        lines = f.readlines()
    assert len(lines) == 36


##############################PrepareGaussianHostGuest##############################
# drop in a log file here
def test_copy_log_file_host_guest():
    source_ = get_data_filename("test_host_guest.log")
    destination_pwd = os.getcwd()
    destination_file = source_.split("/")[-1]
    destination_ = destination_pwd + "/" + destination_file
    qmmmrebind.parameterize.copy_file(source=source_, destination=destination_)
    assert "test_host_guest.log" in os.listdir()


# drop in a com file here
def test_copy_com_file_host_guest():
    source_ = get_data_filename("test_host_guest.com")
    destination_pwd = os.getcwd()
    destination_file = source_.split("/")[-1]
    destination_ = destination_pwd + "/" + destination_file
    qmmmrebind.parameterize.copy_file(source=source_, destination=destination_)
    assert "test_host_guest.com" in os.listdir()


def test_get_qm_host_guest_charges():
    host_guest_input = "test_host_guest.com"
    guest_pdb = "test_guest_init_ii.pdb"
    qm_guest_charge_parameter_file = "test_guest_qm_surround_charges.txt"
    qm_host_charge_parameter_file = "test_host_qm_surround_charges.txt"
    qm_guest_atom_charge_parameter_file = (
        "test_guest_qm_atom_surround_charges.txt"
    )
    get_get_qm_host_guest_charges = qmmmrebind.parameterize.PrepareGaussianHostGuest(
        guest_pdb=guest_pdb,
        host_qm_pdb="",
        n_processors="",
        memory="",
        charge="",
        multiplicity="",
        functional="",
        basis_set="",
        optimisation="",
        frequency="",
        add_keywords_I="",
        add_keywords_II="",
        add_keywords_III="",
        gauss_system_out_file="",
        fchk_system_out_file="",
        host_guest_input=host_guest_input,
        qm_guest_charge_parameter_file=qm_guest_charge_parameter_file,
        qm_host_charge_parameter_file=qm_host_charge_parameter_file,
        qm_guest_atom_charge_parameter_file=qm_guest_atom_charge_parameter_file,
    )
    get_get_qm_host_guest_charges.get_qm_host_guest_charges()
    with open(qm_guest_atom_charge_parameter_file, "r") as f:
        guest_qm_charges = f.readlines()
    with open(qm_host_charge_parameter_file, "r") as f:
        host_qm_charges = f.readlines()
    assert len(guest_qm_charges) == 18
    assert len(host_qm_charges) == 27


##############################ParameterizeHost##############################
# drop in a fchk file here
def test_copy_fchk_file_host():
    source_ = get_data_filename("test_host_qm.fchk")
    destination_pwd = os.getcwd()
    destination_file = source_.split("/")[-1]
    destination_ = destination_pwd + "/" + destination_file
    qmmmrebind.parameterize.copy_file(source=source_, destination=destination_)
    assert "test_host_qm.fchk" in os.listdir()


def test_get_xyz_host():
    host_qm_pdb = "test_host_qm.pdb"
    xyz_file = "test_host_qm_coords.xyz"
    coordinate_file = "test_host_qm_coordinates.txt"
    get_get_xyz_host = qmmmrebind.parameterize.ParameterizeHost(
        xyz_file=xyz_file,
        coordinate_file=coordinate_file,
        unprocessed_hessian_file="",
        bond_list_file="",
        angle_list_file="",
        hessian_file="",
        atom_names_file="",
        bond_parameter_file="",
        angle_parameter_file="",
        charge_parameter_file="",
        host_qm_pdb=host_qm_pdb,
    )
    get_get_xyz_host.get_xyz()
    with open(coordinate_file, "r") as f:
        no_coords = f.readlines()
    with open(xyz_file, "r") as f:
        xyz_lines = f.readlines()
    assert len(no_coords) == (len(xyz_lines) - 2)


def test_get_unprocessed_hessian_host():
    host_qm_pdb = "test_host_qm.pdb"
    unprocessed_hessian_file = "test_host_qm_unprocessed_hessian.txt"
    get_get_unprocessed_hessian = qmmmrebind.parameterize.ParameterizeHost(
        xyz_file="",
        coordinate_file="",
        unprocessed_hessian_file=unprocessed_hessian_file,
        bond_list_file="",
        angle_list_file="",
        hessian_file="",
        atom_names_file="",
        bond_parameter_file="",
        angle_parameter_file="",
        charge_parameter_file="",
        host_qm_pdb=host_qm_pdb,
    )
    get_get_unprocessed_hessian.get_unprocessed_hessian()
    with open("test_host_qm.fchk", "r") as f:
        lines = f.readlines()
    for i in range(len(lines)):
        if "Cartesian Force Constants" in lines[i]:
            no_hessian = re.findall(r"\d+|\d+.\d+", lines[i])
            no_hessian = int(no_hessian[0])
    with open(unprocessed_hessian_file, "r") as f:
        hessian_lines = f.readlines()
    assert len(hessian_lines) == no_hessian


# drop in a log file here
def test_copy_log_file_host():
    source_ = get_data_filename("test_host_qm.log")
    destination_pwd = os.getcwd()
    destination_file = source_.split("/")[-1]
    destination_ = destination_pwd + "/" + destination_file
    qmmmrebind.parameterize.copy_file(source=source_, destination=destination_)
    assert "test_host_qm.log" in os.listdir()


def test_get_bond_angles_host():
    host_qm_pdb = "test_host_qm.pdb"
    bond_list_file = "test_host_qm_bond_list.txt"
    angle_list_file = "test_host_qm_angle_list"
    get_get_bond_angles = qmmmrebind.parameterize.ParameterizeHost(
        xyz_file="",
        coordinate_file="",
        unprocessed_hessian_file="",
        bond_list_file=bond_list_file,
        angle_list_file=angle_list_file,
        hessian_file="",
        atom_names_file="",
        bond_parameter_file="",
        angle_parameter_file="",
        charge_parameter_file="",
        host_qm_pdb=host_qm_pdb,
    )
    get_get_bond_angles.get_bond_angles()
    with open(bond_list_file, "r") as f:
        bonds = f.readlines()
    with open(angle_list_file, "r") as f:
        angles = f.readlines()
    assert len(bonds) == 26
    assert len(angles) == 45


def test_get_hessian_host():
    host_qm_pdb = "test_host_qm.pdb"
    unprocessed_hessian_file = "test_host_qm_unprocessed_hessian.txt"
    hessian_file = "test_host_qm_hessian.txt"
    get_get_hessian = qmmmrebind.parameterize.ParameterizeHost(
        xyz_file="",
        coordinate_file="",
        unprocessed_hessian_file=unprocessed_hessian_file,
        bond_list_file="",
        angle_list_file="",
        hessian_file=hessian_file,
        atom_names_file="",
        bond_parameter_file="",
        angle_parameter_file="",
        charge_parameter_file="",
        host_qm_pdb=host_qm_pdb,
    )
    get_get_hessian.get_hessian()
    hessian = np.loadtxt(hessian_file)
    assert hessian.size == 6561


def test_get_atom_names_host():
    host_qm_pdb = "test_host_qm.pdb"
    atom_names_file = "test_host_qm_atom_names.txt"
    get_get_atom_names = qmmmrebind.parameterize.ParameterizeHost(
        xyz_file="",
        coordinate_file="",
        unprocessed_hessian_file="",
        bond_list_file="",
        angle_list_file="",
        hessian_file="",
        atom_names_file=atom_names_file,
        bond_parameter_file="",
        angle_parameter_file="",
        charge_parameter_file="",
        host_qm_pdb=host_qm_pdb,
    )
    get_get_atom_names.get_atom_names()
    with open(atom_names_file, "r") as f:
        no_atoms = len(f.readlines())
    assert no_atoms == 27


def test_get_bond_angle_params_host():
    host_qm_pdb = "test_host_qm.pdb"
    coordinate_file = "test_host_qm_coordinates.txt"
    hessian_file = "test_host_qm_hessian.txt"
    bond_list_file = "test_host_qm_bond_list.txt"
    atom_names_file = "test_host_qm_atom_names.txt"
    bond_parameter_file = "test_host_qm_bonds.txt"
    angle_list_file = "test_host_qm_angle_list"
    angle_parameter_file = "test_host_qm_angles.txt"
    get_get_bond_angle_params = qmmmrebind.parameterize.ParameterizeHost(
        xyz_file="",
        coordinate_file=coordinate_file,
        unprocessed_hessian_file="",
        bond_list_file=bond_list_file,
        angle_list_file=angle_list_file,
        hessian_file=hessian_file,
        atom_names_file=atom_names_file,
        bond_parameter_file=bond_parameter_file,
        angle_parameter_file=angle_parameter_file,
        charge_parameter_file="",
        host_qm_pdb=host_qm_pdb,
    )
    get_get_bond_angle_params.get_bond_angle_params()
    with open(bond_parameter_file, "r") as f:
        no_bonds = len(f.readlines())
    with open(angle_parameter_file, "r") as f:
        no_angles = len(f.readlines())
    with open(bond_list_file, "r") as f:
        bonds = len(f.readlines())
    with open(angle_list_file, "r") as f:
        angles = len(f.readlines())
    assert no_bonds == bonds
    assert no_angles == angles


def test_get_charges_host():
    host_qm_pdb = "test_host_qm.pdb"
    charge_parameter_file = "test_host_qm_surround_charges.txt"
    get_get_charges = qmmmrebind.parameterize.ParameterizeHost(
        xyz_file="",
        coordinate_file="",
        unprocessed_hessian_file="",
        bond_list_file="",
        angle_list_file="",
        hessian_file="",
        atom_names_file="",
        bond_parameter_file="",
        angle_parameter_file="",
        charge_parameter_file=charge_parameter_file,
        host_qm_pdb=host_qm_pdb,
    )
    get_get_charges.get_charges()
    with open(charge_parameter_file, "r") as f:
        no_lines = len(f.readlines())
    assert no_lines == 27


##############################GuestAmberXMLAmber##############################
def test_copy_system_init_xml():
    source_ = get_data_filename("test_guest_init.xml")
    destination_pwd = os.getcwd()
    destination_file = source_.split("/")[-1]
    destination_ = destination_pwd + "/" + destination_file
    qmmmrebind.parameterize.copy_file(source=source_, destination=destination_)
    assert "test_guest_init.xml" in os.listdir()

"""
def test_generate_xml_from_charged_pdb_sdf():
    system_pdb = "test_guest_init_ii.pdb"
    system_init_sdf = "test_guest_init.sdf"
    system_sdf = "test_guest.sdf"
    num_charge_atoms = 1
    index_charge_atom_1 = 9
    charge_atom_1 = 1
    system_xml = "test_guest_init.xml"
    get_generate_xml_from_charged_pdb_sdf = (
        qmmmrebind.parameterize.GuestAmberXMLAmber(
            system_pdb=system_pdb,
            system_mol2="",
            system_in="",
            charge="",
            system_frcmod="",
            prmtop_system="",
            inpcrd_system="",
            system_leap="",
            system_xml=system_xml,
            system_smi="",
            system_sdf=system_sdf,
            system_init_sdf=system_init_sdf,
            num_charge_atoms=num_charge_atoms,
            index_charge_atom_1=index_charge_atom_1,
            charge_atom_1=charge_atom_1,
            index_charge_atom_2="",
            charge_atom_2="",
            charge_parameter_file="",
            system_qm_pdb="",
            bond_parameter_file="",
            angle_parameter_file="",
            system_qm_params_file="",
            reparameterised_intermediate_system_xml_file="",
            system_xml_non_bonded_file="",
            system_xml_non_bonded_reparams_file="",
            reparameterised_system_xml_file="",
            non_reparameterised_system_xml_file="",
            prmtop_system_non_params="",
            inpcrd_system_non_params="",
            prmtop_system_params="",
            inpcrd_system_params="",
            load_topology="",
        )
    )
    get_generate_xml_from_charged_pdb_sdf.generate_xml_from_charged_pdb_sdf()
    with open(system_xml, "r") as f:
        lines = f.readlines()
    angle_lines = []
    bond_lines = []
    for line in lines:
        if "Angle " in line:
            angle_lines.append(line)
        if "Bond " in line:
            bond_lines.append(line)
    assert len(angle_lines) != 0
    assert len(bond_lines) != 0
"""

def test_write_system_params():
    charge_parameter_file = "test_guest_charges.txt"
    system_qm_pdb = "test_guest_init_ii.pdb"
    bond_parameter_file = "test_guest_bonds.txt"
    angle_parameter_file = "test_guest_angles.txt"
    system_qm_params_file = "test_guest_qm_params.txt"
    get_write_system_params = qmmmrebind.parameterize.GuestAmberXMLAmber(
        system_pdb="",
        system_mol2="",
        system_in="",
        charge="",
        system_frcmod="",
        prmtop_system="",
        inpcrd_system="",
        system_leap="",
        system_xml="",
        system_smi="",
        system_sdf="",
        system_init_sdf="",
        num_charge_atoms="",
        index_charge_atom_1="",
        charge_atom_1="",
        index_charge_atom_2="",
        charge_atom_2="",
        charge_parameter_file=charge_parameter_file,
        system_qm_pdb=system_qm_pdb,
        bond_parameter_file=bond_parameter_file,
        angle_parameter_file=angle_parameter_file,
        system_qm_params_file=system_qm_params_file,
        reparameterised_intermediate_system_xml_file="",
        system_xml_non_bonded_file="",
        system_xml_non_bonded_reparams_file="",
        reparameterised_system_xml_file="",
        non_reparameterised_system_xml_file="",
        prmtop_system_non_params="",
        inpcrd_system_non_params="",
        prmtop_system_params="",
        inpcrd_system_params="",
        load_topology="",
    )
    get_write_system_params.write_system_params()
    with open(system_qm_params_file, "r") as f:
        lines = f.readlines()
    angle_lines = []
    for line in lines:
        if "Angle " in line:
            angle_lines.append(line)
    assert len(angle_lines) - 2 == 27


def test_write_reparameterised_system_xml():
    system_qm_params_file = ("test_guest_qm_params.txt",)
    system_xml = "test_guest_init.xml"
    reparameterised_intermediate_system_xml_file = (
        "test_guest_intermediate_reparameterised.xml"
    )
    system_qm_params_file = "test_guest_qm_params.txt"
    system_xml_non_bonded_file = "test_guest_xml_non_bonded.txt"
    system_xml_non_bonded_reparams_file = (
        "test_guest_xml_non_bonded_reparams.txt"
    )
    reparameterised_system_xml_file = "test_guest_reparameterised.xml"
    get_write_reparameterised_system_xml = qmmmrebind.parameterize.GuestAmberXMLAmber(
        system_pdb="",
        system_mol2="",
        system_in="",
        charge="",
        system_frcmod="",
        prmtop_system="",
        inpcrd_system="",
        system_leap="",
        system_xml=system_xml,
        system_smi="",
        system_sdf="",
        system_init_sdf="",
        num_charge_atoms="",
        index_charge_atom_1="",
        charge_atom_1="",
        index_charge_atom_2="",
        charge_atom_2="",
        charge_parameter_file="",
        system_qm_pdb="",
        bond_parameter_file="",
        angle_parameter_file="",
        system_qm_params_file=system_qm_params_file,
        reparameterised_intermediate_system_xml_file=reparameterised_intermediate_system_xml_file,
        system_xml_non_bonded_file=system_xml_non_bonded_file,
        system_xml_non_bonded_reparams_file=system_xml_non_bonded_reparams_file,
        reparameterised_system_xml_file=reparameterised_system_xml_file,
        non_reparameterised_system_xml_file="",
        prmtop_system_non_params="",
        inpcrd_system_non_params="",
        prmtop_system_params="",
        inpcrd_system_params="",
        load_topology="",
    )
    get_write_reparameterised_system_xml.write_reparameterised_system_xml()
    with open(system_xml, "r") as f:
        lines = f.readlines()
    for line in lines:
        if 'p1="1" p2="0" p3="5"' in line:
            a = line.split()[1][3:-1]
            break


def test_save_amber_params():
    load_topology = "openmm"
    system_pdb = "test_guest_init_ii.pdb"
    non_reparameterised_system_xml_file = "test_guest_init.xml"
    prmtop_system_non_params = "test_guest_non_params.prmtop"
    inpcrd_system_non_params = "test_guest_non_params.inpcrd"
    reparameterised_system_xml_file = "test_guest_reparameterised.xml"
    inpcrd_system_params = "test_guest_params.inpcrd"
    prmtop_system_params = "test_guest_params.prmtop"
    get_save_amber_params = qmmmrebind.parameterize.GuestAmberXMLAmber(
        system_pdb=system_pdb,
        system_mol2="",
        system_in="",
        charge="",
        system_frcmod="",
        prmtop_system="",
        inpcrd_system="",
        system_leap="",
        system_xml="",
        system_smi="",
        system_sdf="",
        system_init_sdf="",
        num_charge_atoms="",
        index_charge_atom_1="",
        charge_atom_1="",
        index_charge_atom_2="",
        charge_atom_2="",
        charge_parameter_file="",
        system_qm_pdb="",
        bond_parameter_file="",
        angle_parameter_file="",
        system_qm_params_file="",
        reparameterised_intermediate_system_xml_file="",
        system_xml_non_bonded_file="",
        system_xml_non_bonded_reparams_file="",
        reparameterised_system_xml_file=reparameterised_system_xml_file,
        non_reparameterised_system_xml_file=non_reparameterised_system_xml_file,
        prmtop_system_non_params=prmtop_system_non_params,
        inpcrd_system_non_params=inpcrd_system_non_params,
        prmtop_system_params=prmtop_system_params,
        inpcrd_system_params=inpcrd_system_params,
        load_topology=load_topology,
    )
    get_save_amber_params.save_amber_params()
    with open(inpcrd_system_params, "r") as f:
        lines_top = f.readlines()
    assert len(lines_top[2:]) == 9
    with open(prmtop_system_params, "r") as f:
        lines_inp = f.readlines()
    angle_eq_val = ""
    for i, line in enumerate(lines_inp):
        if "%FLAG ANGLE_EQUIL_VALUE" in line:
            angle_eq_val = lines_inp[i + 2]
            break
    assert angle_eq_val.strip().split()[0] == "2.09090400E+00"


##############################HostAmberXMLAmber##############################
def test_serialize_system():
    system_pdb = "test_host.pdb"
    sim_output = "test_sim_output.pdb"
    sim_steps = 100
    system_xml = "test_host.xml"
    get_serialize_system = qmmmrebind.parameterize.HostAmberXMLAmber(
        system_pdb=system_pdb,
        system_xml=system_xml,
        sim_output=sim_output,
        sim_steps=sim_steps,
        charge_parameter_file="",
        system_qm_pdb="",
        bond_parameter_file="",
        angle_parameter_file="",
        system_qm_params_file="",
        reparameterised_intermediate_system_xml_file="",
        system_xml_non_bonded_file="",
        system_xml_non_bonded_reparams_file="",
        reparameterised_system_xml_file="",
        non_reparameterised_system_xml_file="",
        prmtop_system_non_params="",
        inpcrd_system_non_params="",
        prmtop_system_params="",
        inpcrd_system_params="",
        load_topology="",
    )
    get_serialize_system.serialize_system()
    with open(system_xml, "r") as f:
        lines = f.readlines()
    particle_lines = []
    for line in lines:
        if "Particle " in line:
            particle_lines.append(line)
    assert len(particle_lines) == 6440


def test_write_system_params_host():
    charge_parameter_file = "test_host_qm_surround_charges.txt"
    system_qm_pdb = "test_host_qm.pdb"
    bond_parameter_file = "test_host_qm_bonds.txt"
    angle_parameter_file = "test_host_qm_angles.txt"
    system_qm_params_file = "test_host_qm_params.txt"
    get_write_system_params = qmmmrebind.parameterize.HostAmberXMLAmber(
        system_pdb="",
        system_xml="",
        sim_output="",
        sim_steps="",
        charge_parameter_file=charge_parameter_file,
        system_qm_pdb=system_qm_pdb,
        bond_parameter_file=bond_parameter_file,
        angle_parameter_file=angle_parameter_file,
        system_qm_params_file=system_qm_params_file,
        reparameterised_intermediate_system_xml_file="",
        system_xml_non_bonded_file="",
        system_xml_non_bonded_reparams_file="",
        reparameterised_system_xml_file="",
        non_reparameterised_system_xml_file="",
        prmtop_system_non_params="",
        inpcrd_system_non_params="",
        prmtop_system_params="",
        inpcrd_system_params="",
        load_topology="",
    )
    get_write_system_params.write_system_params()
    with open(system_qm_params_file, "r") as f:
        lines = f.readlines()
    angle_lines = []
    for line in lines:
        if "Angle " in line:
            angle_lines.append(line)
    assert len(angle_lines) - 2 == 45


def test_write_reparameterised_system_xml_host():
    system_qm_params_file = "test_host_qm_params.txt"
    system_xml = "test_host.xml"
    reparameterised_intermediate_system_xml_file = (
        "test_host_intermediate_reparameterised.xml"
    )
    system_xml_non_bonded_file = "test_host_xml_non_bonded.txt"
    system_xml_non_bonded_reparams_file = (
        "test_host_xml_non_bonded_reparams.txt"
    )
    reparameterised_system_xml_file = "test_host_reparameterised.xml"
    get_write_reparameterised_system_xml = qmmmrebind.parameterize.HostAmberXMLAmber(
        system_pdb="",
        system_xml=system_xml,
        sim_output="",
        sim_steps="",
        charge_parameter_file="",
        system_qm_pdb="",
        bond_parameter_file="",
        angle_parameter_file="",
        system_qm_params_file=system_qm_params_file,
        reparameterised_intermediate_system_xml_file=reparameterised_intermediate_system_xml_file,
        system_xml_non_bonded_file=system_xml_non_bonded_file,
        system_xml_non_bonded_reparams_file=system_xml_non_bonded_reparams_file,
        reparameterised_system_xml_file=reparameterised_system_xml_file,
        non_reparameterised_system_xml_file="",
        prmtop_system_non_params="",
        inpcrd_system_non_params="",
        prmtop_system_params="",
        inpcrd_system_params="",
        load_topology="",
    )
    get_write_reparameterised_system_xml.write_reparameterised_system_xml()
    with open(reparameterised_system_xml_file, "r") as f:
        lines = f.readlines()
    for line in lines:
        if 'p1="0" p2="4" p3="5"' in line:
            a = line.split()[1][3:-1]
            break

def test_save_amber_params_host():
    load_topology = "openmm"
    system_pdb = "test_host.pdb"
    non_reparameterised_system_xml_file = "test_host.xml"
    prmtop_system_non_params = "test_host_non_params.prmtop"
    inpcrd_system_non_params = "test_host_non_params.inpcrd"
    reparameterised_system_xml_file = "test_host_reparameterised.xml"
    prmtop_system_params = "test_host_params.prmtop"
    inpcrd_system_params = "test_host_params.inpcrd"
    get_save_amber_params = qmmmrebind.parameterize.HostAmberXMLAmber(
        system_pdb=system_pdb,
        system_xml="",
        sim_output="",
        sim_steps="",
        charge_parameter_file="",
        system_qm_pdb="",
        bond_parameter_file="",
        angle_parameter_file="",
        system_qm_params_file="",
        reparameterised_intermediate_system_xml_file="",
        system_xml_non_bonded_file="",
        system_xml_non_bonded_reparams_file="",
        reparameterised_system_xml_file=reparameterised_system_xml_file,
        non_reparameterised_system_xml_file=non_reparameterised_system_xml_file,
        prmtop_system_non_params=prmtop_system_non_params,
        inpcrd_system_non_params=inpcrd_system_non_params,
        prmtop_system_params=prmtop_system_params,
        inpcrd_system_params=inpcrd_system_params,
        load_topology=load_topology,
    )
    get_save_amber_params.save_amber_params()
    with open(inpcrd_system_params, "r") as f:
        lines_top = f.readlines()
    assert len(lines_top[2:]) == 1610
    with open(prmtop_system_params, "r") as f:
        lines_inp = f.readlines()
    angle_eq_val = ""
    for i, line in enumerate(lines_inp):
        if "%FLAG ANGLE_EQUIL_VALUE" in line:
            angle_eq_val = lines_inp[i + 2]
            break
    assert angle_eq_val.strip().split()[0] == "1.91113553E+00"


##############################MergeHostGuestTopology#########################
def test_merge_topology_files():
    host_prmtop = "test_host_params.prmtop"
    guest_prmtop = "test_guest_params.prmtop"
    host_inpcrd = "test_host_params.inpcrd"
    guest_inpcrd = "test_guest_params.inpcrd"
    system_prmtop = "test_system_params.prmtop"
    system_inpcrd = "test_system_params.inpcrd"
    get_merged_topology = qmmmrebind.parameterize.MergeHostGuestTopology(
        host_prmtop,
        guest_prmtop,
        host_inpcrd,
        guest_inpcrd,
        system_prmtop,
        system_inpcrd,
    )
    get_merged_topology.merge_topology_files()
    with open(system_inpcrd, "r") as f:
        lines_top = f.readlines()
    assert len(lines_top[2:]) == 1619
    with open(system_prmtop, "r") as f:
        lines_inp = f.readlines()
    angle_eq_val = ""
    for i, line in enumerate(lines_inp):
        if "%FLAG ANGLE_EQUIL_VALUE" in line:
            angle_eq_val = lines_inp[i + 2]
            break
    assert angle_eq_val.strip().split()[0] == "1.91113553E+00"


##############################TorsionDriveSims#############################
# files needed: "test_torsion_drive_input.dat", "test_dihedrals.txt", "test_run_command"


def test_write_tor_params_txt():
    reparameterised_system_xml_file = "test_guest_reparameterised.xml"
    torsion_xml_file = "test_guest_torsion_xml.txt"
    get_torsional_parameters_object = qmmmrebind.parameterize.TorsionDriveSims(
        1,
        1,
        reparameterised_system_xml_file=reparameterised_system_xml_file,
        torsion_xml_file=torsion_xml_file,
    )
    get_torsional_parameters_object.write_tor_params_txt()
    with open(torsion_xml_file, "r") as f:
        lines = f.readlines()
    assert len(lines) == 65


def test_copy_run_command():
    source_ = get_data_filename("test_run_command")
    destination_pwd = os.getcwd()
    destination_file = source_.split("/")[-1]
    destination_ = destination_pwd + "/" + destination_file
    qmmmrebind.parameterize.copy_file(source=source_, destination=destination_)
    assert "test_run_command" in os.listdir()


def test_create_torsion_drive_dir():
    torsion_xml_file = "test_guest_torsion_xml.txt"
    tor_dir = "torsion_dir"
    template_pdb = "test_guest_init_ii.pdb"
    dihedral_text_file = "test_dihedrals.txt"
    torsion_drive_run_file = "test_run_command"
    psi_input_file = "test_torsion_drive_input.dat"
    torsion_drive_sims_object = qmmmrebind.parameterize.TorsionDriveSims(
        1,
        1,
        torsion_xml_file=torsion_xml_file,
        tor_dir=tor_dir,
        template_pdb=template_pdb,
        dihedral_text_file=dihedral_text_file,
        torsion_drive_run_file=torsion_drive_run_file,
        psi_input_file=psi_input_file,
    )
    torsion_drive_sims_object.create_torsion_drive_dir()
    assert len(os.listdir("./torsion_dir")) == 45


def test_create_non_H_torsion_drive_dir():
    torsion_xml_file = "test_guest_torsion_xml.txt"
    template_pdb = "test_guest_init_ii.pdb"
    tor_dir = "torsion_dir"
    psi_input_file = "test_torsion_drive_input.dat"
    torsion_drive_run_file = "test_run_command"
    dihedral_text_file = "test_dihedrals.txt"
    torsion_drive_sims_object = qmmmrebind.parameterize.TorsionDriveSims(
        1,
        1,
        torsion_xml_file=torsion_xml_file,
        tor_dir=tor_dir,
        template_pdb=template_pdb,
        dihedral_text_file=dihedral_text_file,
        torsion_drive_run_file=torsion_drive_run_file,
        psi_input_file=psi_input_file,
    )
    torsion_drive_sims_object.create_non_H_torsion_drive_dir()
    assert len(os.listdir("./torsion_dir")) == 14


def test_create_non_H_bonded_torsion_drive_dir():
    torsion_xml_file = "test_guest_torsion_xml.txt"
    template_pdb = "test_guest_init_ii.pdb"
    tor_dir = "torsion_dir"
    psi_input_file = "test_torsion_drive_input.dat"
    torsion_drive_run_file = "test_run_command"
    dihedral_text_file = "test_dihedrals.txt"
    system_bonds_file = "test_guest_bonds.txt"
    torsion_drive_sims_object = qmmmrebind.parameterize.TorsionDriveSims(
        1,
        1,
        torsion_xml_file=torsion_xml_file,
        tor_dir=tor_dir,
        template_pdb=template_pdb,
        dihedral_text_file=dihedral_text_file,
        torsion_drive_run_file=torsion_drive_run_file,
        psi_input_file=psi_input_file,
        system_bonds_file=system_bonds_file,
    )
    torsion_drive_sims_object.create_non_H_bonded_torsion_drive_dir()
    assert len(os.listdir("./torsion_dir")) == 12


##############################RemoveTestFiles##############################
def test_remove_files():
    command = "rm -rf __pycache__ test_host_mm.txt test_guest_coord.txt test_host.pdb test_guest_init_ii.pdb test_host_qm.pdb test_guest_init.pdb test_host_qm.txt test_host_mm.pdb test_mm.pdb test_host_mm_region_ii.pdb test_qm.pdb test_host_mm_region_ii.txt test_residue_list.txt test_host_mm_region_i.pdb test_system.pdb test_host_mm_region_i.txt test_guest_bonds.txt test_guest_charges.txt test_guest_coordinates.txt test_guest_coords.xyz test_guest_hessian.txt test_guest_init_ii.fchk test_guest_init_ii.log test_guest_angle_list.txt test_guest_unprocessed_hessian.txt test_guest_angles.txt test_proper_dihedrals.txt test_guest_atom_names.txt test_guest_bond_list.txt test_host_qm_coordinates.txt test_host_qm_coords.xyz test_host_qm.fchk test_host_qm_unprocessed_hessian.txt test_host_qm_angle_list test_host_qm_bond_list.txt test_host_qm.log test_host_qm_angles.txt test_host_qm_atom_names.txt test_host_qm_bonds.txt test_host_qm_hessian.txt test_guest_qm_surround_charges.txt test_host_guest.com test_host_guest.log test_host_qm_surround_charges.txt test_guest_qm_atom_surround_charges.txt test_guest_init.sdf test_guest_init.xml test_guest.sdf test_guest_qm_params.txt test_guest_intermediate_reparameterised.xml test_guest_reparameterised.xml test_guest_xml_non_bonded_reparams.txt test_guest_xml_non_bonded.txt test_guest_non_params.inpcrd test_guest_non_params.prmtop test_guest_params.inpcrd test_guest_params.prmtop test_host.xml test_host_qm_params.txt test_host_intermediate_reparameterised.xml test_host_reparameterised.xml test_host_xml_non_bonded_reparams.txt test_host_xml_non_bonded.txt test_host_non_params.inpcrd test_host_non_params.prmtop test_host_params.inpcrd host_params.prmtop test_system_params.inpcrd test_system_params.prmtop test_torsion_drive_input.xyz test_torsion_drive_input.txt test_torsion_drive_input.pdb test_host_params.prmtop torsion_dir test_guest_torsion_xml.txt test_host_guest.fchk test_torsion_drive_input.dat test_scan.xyz test_run_command"
    os.system(command)
