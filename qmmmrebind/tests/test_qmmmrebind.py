"""
Unit and regression test for the qmmmrebind package.
"""
import warnings
warnings.filterwarnings("ignore")
from .utils import get_data_filename
import qmmmrebind
import pytest
import sys
import os
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
    ret = qmmmrebind.parameterize.search_in_file(file = filename, word = word)
    assert ret[0][0] == 3

def test_find_word_in_file():
    """Test if a word is in a file"""
    filename = get_data_filename("test_random_read_file_i.dat") 
    word = "github"
    ret = qmmmrebind.parameterize.search_in_file(file = filename, word = word)
    assert ret[0][0] == 3
##############################PrepareQMMM############################## 
def test_clean_up():
    init_pdb = get_data_filename("test_sample_system_trypsin_benzamidine.pdb")
    cleaned_pdb = "test_system.pdb"
    num_residues = 2
    get_clean_up = qmmmrebind.parameterize.PrepareQMMM(init_pdb = init_pdb, cleaned_pdb = cleaned_pdb, guest_init_pdb = "", host_pdb = "", guest_resname = "", guest_pdb = "", guest_xyz = "", distance = "", residue_list = "", host_qm_atoms = "", host_mm_atoms = "", host_qm_pdb = "", host_mm_pdb = "", qm_pdb = "", mm_pdb = "", host_mm_region_I_atoms = "", host_mm_region_II_atoms = "", host_mm_region_I_pdb = "", host_mm_region_II_pdb = "", num_residues = num_residues)
    get_clean_up.clean_up()
    with open (cleaned_pdb, "r") as f:
        lines = f.readlines()
    assert len(lines) != 0
    HOH_list = []
    for i in lines:
        if "HOH"in i:
            HOH_list.append(i)
    HETATM_list = []
    for i in lines:
        if "HETATM"in i:
            HETATM_list.append(i)
    assert len(HOH_list) == 0
    assert len(HETATM_list) == 0

def test_create_host_guest():  
    cleaned_pdb = "test_system.pdb"
    guest_init_pdb = "test_guest_init.pdb"
    host_pdb = "test_host.pdb"    
    guest_resname  = "BEN"  
    num_residues = 2 
    get_create_host_guest = qmmmrebind.parameterize.PrepareQMMM(init_pdb = "", cleaned_pdb = cleaned_pdb, guest_init_pdb = guest_init_pdb, host_pdb = host_pdb, guest_resname = guest_resname, guest_pdb = "", guest_xyz = "", distance = "", residue_list = "", host_qm_atoms = "", host_mm_atoms = "", host_qm_pdb = "", host_mm_pdb = "", qm_pdb = "", mm_pdb = "", host_mm_region_I_atoms = "", host_mm_region_II_atoms = "", host_mm_region_I_pdb = "", host_mm_region_II_pdb = "", num_residues = num_residues)  
    get_create_host_guest.create_host_guest()  
    with open (guest_init_pdb, "r") as f:
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
    get_realign_guest = qmmmrebind.parameterize.PrepareQMMM(init_pdb = "", cleaned_pdb = "", guest_init_pdb = guest_init_pdb, host_pdb = "", guest_resname = "", guest_pdb = guest_pdb, guest_xyz = "", distance = "", residue_list = "", host_qm_atoms = "", host_mm_atoms = "", host_qm_pdb = "", host_mm_pdb = "", qm_pdb = "", mm_pdb = "", host_mm_region_I_atoms = "", host_mm_region_II_atoms = "", host_mm_region_I_pdb = "", host_mm_region_II_pdb = "", num_residues = num_residues)
    get_realign_guest.realign_guest()  
    with open (guest_init_pdb, "r") as f:
        guest_init_pdb_lines = f.readlines()
    with open (guest_pdb, "r") as f:
        guest_pdb_lines = f.readlines()
    assert len(guest_init_pdb_lines) == len(guest_pdb_lines)

def test_get_guest_coord():   
    guest_pdb = "test_guest_init_ii.pdb"
    guest_xyz = "test_guest_coord.txt"
    num_residues = 2
    get_get_guest_coord = qmmmrebind.parameterize.PrepareQMMM(init_pdb = "", cleaned_pdb = "", guest_init_pdb = "", host_pdb = "", guest_resname = "", guest_pdb = guest_pdb, guest_xyz = guest_xyz, distance = "", residue_list = "", host_qm_atoms = "", host_mm_atoms = "", host_qm_pdb = "", host_mm_pdb = "", qm_pdb = "", mm_pdb = "", host_mm_region_I_atoms = "", host_mm_region_II_atoms = "", host_mm_region_I_pdb = "", host_mm_region_II_pdb = "", num_residues = num_residues)
    get_get_guest_coord.get_guest_coord()     

def test_get_qm_resids():   
    guest_xyz = "test_guest_coord.txt"
    host_pdb = "test_host.pdb"
    distance = 3.0
    residue_list = "test_residue_list.txt"
    num_residues = 2
    get_get_qm_resids = qmmmrebind.parameterize.PrepareQMMM(init_pdb = "", cleaned_pdb = "", guest_init_pdb = "", host_pdb = host_pdb, guest_resname = "", guest_pdb = "", guest_xyz = guest_xyz, distance = distance, residue_list = residue_list, host_qm_atoms = "", host_mm_atoms = "", host_qm_pdb = "", host_mm_pdb = "", qm_pdb = "", mm_pdb = "", host_mm_region_I_atoms = "", host_mm_region_II_atoms = "", host_mm_region_I_pdb = "", host_mm_region_II_pdb = "", num_residues = num_residues)
    get_get_qm_resids.get_qm_resids()     

def test_get_host_qm_mm_atoms():   
    residue_list = "test_residue_list.txt"
    num_residues = 2
    host_pdb = "test_host.pdb"
    host_qm_atoms = "test_host_qm.txt" 
    host_mm_atoms = "test_host_mm.txt"
    get_get_host_qm_mm_atoms = qmmmrebind.parameterize.PrepareQMMM(init_pdb = "", cleaned_pdb = "", guest_init_pdb = "", host_pdb = host_pdb, guest_resname = "", guest_pdb = "", guest_xyz = "", distance = "", residue_list = residue_list, host_qm_atoms = host_qm_atoms, host_mm_atoms = host_mm_atoms, host_qm_pdb = "", host_mm_pdb = "", qm_pdb = "", mm_pdb = "", host_mm_region_I_atoms = "", host_mm_region_II_atoms = "", host_mm_region_I_pdb = "", host_mm_region_II_pdb = "", num_residues = num_residues)
    get_get_host_qm_mm_atoms.get_host_qm_mm_atoms()   

def test_save_host_pdbs():   
    host_pdb = "test_host.pdb"
    host_qm_pdb = "test_host_qm.pdb" 
    host_mm_pdb = "test_host_mm.pdb"
    host_qm_atoms = "test_host_qm.txt" 
    host_mm_atoms = "test_host_mm.txt"
    num_residues = 2
    get_save_host_pdbs = qmmmrebind.parameterize.PrepareQMMM(init_pdb = "", cleaned_pdb = "", guest_init_pdb = "", host_pdb = host_pdb, guest_resname = "", guest_pdb = "", guest_xyz = "", distance = "", residue_list = "", host_qm_atoms = host_qm_atoms, host_mm_atoms = host_mm_atoms, host_qm_pdb = host_qm_pdb, host_mm_pdb = host_mm_pdb, qm_pdb = "", mm_pdb = "", host_mm_region_I_atoms = "", host_mm_region_II_atoms = "", host_mm_region_I_pdb = "", host_mm_region_II_pdb = "", num_residues = num_residues)
    get_save_host_pdbs.save_host_pdbs()  

def test_get_host_mm_region_atoms():   
    residue_list = "test_residue_list.txt"
    num_residues = 2
    host_pdb = "test_host.pdb"
    host_mm_pdb = "test_host_mm.pdb"
    host_mm_region_I_atoms = "test_host_mm_region_i.txt"
    host_mm_region_II_atoms = "test_host_mm_region_ii.txt"
    get_get_host_mm_region_atoms = qmmmrebind.parameterize.PrepareQMMM(init_pdb = "", cleaned_pdb = "", guest_init_pdb = "", host_pdb = host_pdb, guest_resname = "", guest_pdb = "", guest_xyz = "", distance = "", residue_list = residue_list, host_qm_atoms = "", host_mm_atoms = "", host_qm_pdb = "", host_mm_pdb = host_mm_pdb, qm_pdb = "", mm_pdb = "", host_mm_region_I_atoms = host_mm_region_I_atoms, host_mm_region_II_atoms = host_mm_region_II_atoms, host_mm_region_I_pdb = "", host_mm_region_II_pdb = "", num_residues = num_residues)
    get_get_host_mm_region_atoms.get_host_mm_region_atoms()  

def test_save_host_mm_regions_pdbs():   
    host_mm_pdb = "test_host_mm.pdb"
    host_mm_region_I_atoms = "test_host_mm_region_i.txt"
    host_mm_region_II_atoms = "test_host_mm_region_ii.txt"
    host_mm_region_I_pdb = "test_host_mm_region_i.pdb"
    host_mm_region_II_pdb = "test_host_mm_region_ii.pdb"
    num_residues = 2
    get_save_host_mm_regions_pdbs = qmmmrebind.parameterize.PrepareQMMM(init_pdb = "", cleaned_pdb = "", guest_init_pdb = "", host_pdb = "", guest_resname = "", guest_pdb = "", guest_xyz = "", distance = "", residue_list = "", host_qm_atoms = "", host_mm_atoms = "", host_qm_pdb = "", host_mm_pdb = host_mm_pdb, qm_pdb = "", mm_pdb = "", host_mm_region_I_atoms = host_mm_region_I_atoms, host_mm_region_II_atoms = host_mm_region_II_atoms, host_mm_region_I_pdb = host_mm_region_I_pdb, host_mm_region_II_pdb = host_mm_region_II_pdb, num_residues = num_residues)
    get_save_host_mm_regions_pdbs.save_host_mm_regions_pdbs()  

def test_get_qm_mm_regions(): 
    host_qm_pdb = "test_host_qm.pdb" 
    qm_pdb = "test_qm.pdb"
    guest_pdb = "test_guest_init_ii.pdb"
    host_mm_pdb = "test_host_mm.pdb"
    mm_pdb = "test_mm.pdb"
    num_residues = 2
    get_get_qm_mm_regions = qmmmrebind.parameterize.PrepareQMMM(init_pdb = "", cleaned_pdb = "", guest_init_pdb = "", host_pdb = "", guest_resname = "", guest_pdb = guest_pdb, guest_xyz = "", distance = "", residue_list = "", host_qm_atoms = "", host_mm_atoms = "", host_qm_pdb = host_qm_pdb, host_mm_pdb = host_mm_pdb, qm_pdb = qm_pdb, mm_pdb = mm_pdb, host_mm_region_I_atoms = "", host_mm_region_II_atoms = "", host_mm_region_I_pdb = "", host_mm_region_II_pdb = "", num_residues = num_residues)
    get_get_qm_mm_regions.get_qm_mm_regions() 

##############################ParameterizeGuest############################## 
#drop in a fchk file here
def test_copy_fchk_file():
    source_ = get_data_filename("test_guest_init_ii.fchk") 
    destination_pwd = os.getcwd()
    destination_file = source_.split('/')[-1]
    destination_ = destination_pwd + "/" + destination_file
    qmmmrebind.parameterize.copy_file(source = source_, destination = destination_)

def test_get_xyz(): 
    guest_pdb = "test_guest_init_ii.pdb"
    coordinate_file = "test_guest_coordinates.txt"
    xyz_file = "test_guest_coords.xyz"
    get_get_xyz = qmmmrebind.parameterize.ParameterizeGuest(xyz_file = xyz_file, coordinate_file = coordinate_file, unprocessed_hessian_file = "", bond_list_file = "", angle_list_file = "", hessian_file = "", atom_names_file = "", bond_parameter_file = "", vibrational_scaling = "", angle_parameter_file = "", charge_parameter_file = "", guest_pdb = guest_pdb, proper_dihedral_file = "")
    get_get_xyz.get_xyz() 

def test_get_unprocessed_hessian(): 
    guest_pdb = "test_guest_init_ii.pdb"
    unprocessed_hessian_file = "test_guest_unprocessed_hessian.txt"
    get_get_unprocessed_hessian = qmmmrebind.parameterize.ParameterizeGuest(xyz_file = "", coordinate_file = "", unprocessed_hessian_file = unprocessed_hessian_file, bond_list_file = "", angle_list_file = "", hessian_file = "", atom_names_file = "", bond_parameter_file = "", vibrational_scaling = "", angle_parameter_file = "", charge_parameter_file = "", guest_pdb = guest_pdb, proper_dihedral_file = "")
    get_get_unprocessed_hessian.get_unprocessed_hessian() 

#drop in a log file here
def test_copy_log_file():
    source_ = get_data_filename("test_guest_init_ii.log") 
    destination_pwd = os.getcwd()
    destination_file = source_.split('/')[-1]
    destination_ = destination_pwd + "/" + destination_file
    qmmmrebind.parameterize.copy_file(source = source_, destination = destination_)

def test_get_bond_angles(): 
    guest_pdb = "test_guest_init_ii.pdb"
    bond_list_file = "test_guest_bond_list.txt"
    angle_list_file = "test_guest_angle_list.txt"
    get_get_bond_angles = qmmmrebind.parameterize.ParameterizeGuest(xyz_file = "", coordinate_file = "", unprocessed_hessian_file = "", bond_list_file = bond_list_file, angle_list_file = angle_list_file, hessian_file = "", atom_names_file = "", bond_parameter_file = "", vibrational_scaling = "", angle_parameter_file = "", charge_parameter_file = "", guest_pdb = guest_pdb, proper_dihedral_file = "")
    get_get_bond_angles.get_bond_angles() 

def test_get_hessian(): 
    unprocessed_hessian_file = "test_guest_unprocessed_hessian.txt"
    hessian_file = "test_guest_hessian.txt"
    guest_pdb = "test_guest_init_ii.pdb"
    get_get_hessian = qmmmrebind.parameterize.ParameterizeGuest(xyz_file = "", coordinate_file = "", unprocessed_hessian_file = unprocessed_hessian_file, bond_list_file = "", angle_list_file = "", hessian_file = hessian_file, atom_names_file = "", bond_parameter_file = "", vibrational_scaling = "", angle_parameter_file = "", charge_parameter_file = "", guest_pdb = guest_pdb, proper_dihedral_file = "")
    get_get_hessian.get_hessian() 

def test_get_atom_names(): 
    guest_pdb = "test_guest_init_ii.pdb"
    atom_names_file = "test_guest_atom_names.txt"
    get_get_atom_names = qmmmrebind.parameterize.ParameterizeGuest(xyz_file = "", coordinate_file = "", unprocessed_hessian_file = "", bond_list_file = "", angle_list_file = "", hessian_file = "", atom_names_file = atom_names_file, bond_parameter_file = "", vibrational_scaling = "", angle_parameter_file = "", charge_parameter_file = "", guest_pdb = guest_pdb, proper_dihedral_file = "")
    get_get_atom_names.get_atom_names() 

def test_get_bond_angle_params(): 
    guest_pdb = "test_guest_init_ii.pdb"
    coordinate_file = "test_guest_coordinates.txt"
    hessian_file = "test_guest_hessian.txt"
    bond_list_file = "test_guest_bond_list.txt"
    atom_names_file = "test_guest_atom_names.txt"
    bond_parameter_file = "test_guest_bonds.txt"
    vibrational_scaling = 1.00
    angle_list_file = "test_guest_angle_list.txt"
    angle_parameter_file = "test_guest_angles.txt"
    get_get_bond_angle_params = qmmmrebind.parameterize.ParameterizeGuest(xyz_file = "", coordinate_file = coordinate_file, unprocessed_hessian_file = "", bond_list_file = bond_list_file, angle_list_file = angle_list_file, hessian_file = hessian_file, atom_names_file = atom_names_file, bond_parameter_file = bond_parameter_file, vibrational_scaling = vibrational_scaling, angle_parameter_file = angle_parameter_file, charge_parameter_file = "", guest_pdb = guest_pdb, proper_dihedral_file = "")
    get_get_bond_angle_params.get_bond_angle_params() 

def test_get_charges(): 
    guest_pdb = "test_guest_init_ii.pdb"
    charge_parameter_file = "test_guest_charges.txt"
    get_get_charges = qmmmrebind.parameterize.ParameterizeGuest(xyz_file = "", coordinate_file = "", unprocessed_hessian_file = "", bond_list_file = "", angle_list_file = "", hessian_file = "", atom_names_file = "", bond_parameter_file = "", vibrational_scaling = "", angle_parameter_file = "", charge_parameter_file = charge_parameter_file, guest_pdb = guest_pdb, proper_dihedral_file = "")
    get_get_charges.get_charges() 

def test_get_proper_dihedrals(): 
    guest_pdb = "test_guest_init_ii.pdb"
    bond_parameter_file = "test_guest_bonds.txt"
    proper_dihedral_file = "test_proper_dihedrals.txt"
    get_get_proper_dihedrals = qmmmrebind.parameterize.ParameterizeGuest(xyz_file = "", coordinate_file = "", unprocessed_hessian_file = "", bond_list_file = "", angle_list_file = "", hessian_file = "", atom_names_file = "", bond_parameter_file = bond_parameter_file, vibrational_scaling = "", angle_parameter_file = "", charge_parameter_file = "", guest_pdb = guest_pdb, proper_dihedral_file = proper_dihedral_file)
    get_get_proper_dihedrals.get_proper_dihedrals() 

##############################PrepareGaussianHostGuest############################## 
#drop in a log file here
def test_copy_log_file_host_guest():
    source_ = get_data_filename("test_host_guest.log") 
    destination_pwd = os.getcwd()
    destination_file = source_.split('/')[-1]
    destination_ = destination_pwd + "/" + destination_file
    qmmmrebind.parameterize.copy_file(source = source_, destination = destination_)

#drop in a com file here
def test_copy_com_file_host_guest():
    source_ = get_data_filename("test_host_guest.com") 
    destination_pwd = os.getcwd()
    destination_file = source_.split('/')[-1]
    destination_ = destination_pwd + "/" + destination_file
    qmmmrebind.parameterize.copy_file(source = source_, destination = destination_)

def test_get_qm_host_guest_charges(): 
    host_guest_input = "test_host_guest.com"
    guest_pdb = "test_guest_init_ii.pdb"
    qm_guest_charge_parameter_file = "test_guest_qm_surround_charges.txt"
    qm_host_charge_parameter_file = "test_host_qm_surround_charges.txt"
    qm_guest_atom_charge_parameter_file = "test_guest_qm_atom_surround_charges.txt"
   
    get_get_qm_host_guest_charges = qmmmrebind.parameterize.PrepareGaussianHostGuest(guest_pdb = guest_pdb, host_qm_pdb = "", n_processors = "", memory = "", charge = "", multiplicity = "", functional = "", basis_set = "", optimisation = "", frequency = "", add_keywords_I = "", add_keywords_II = "", add_keywords_III = "", gauss_system_out_file = "", fchk_system_out_file = "", host_guest_input = host_guest_input, qm_guest_charge_parameter_file = qm_guest_charge_parameter_file, qm_host_charge_parameter_file = qm_host_charge_parameter_file, qm_guest_atom_charge_parameter_file = qm_guest_atom_charge_parameter_file)
    get_get_qm_host_guest_charges.get_qm_host_guest_charges() 
##############################ParameterizeHost############################## 
#drop in a fchk file here
def test_copy_fchk_file_host():
    source_ = get_data_filename("test_host_qm.fchk") 
    destination_pwd = os.getcwd()
    destination_file = source_.split('/')[-1]
    destination_ = destination_pwd + "/" + destination_file
    qmmmrebind.parameterize.copy_file(source = source_, destination = destination_)

def test_get_xyz_host(): 
    host_qm_pdb = "test_host_qm.pdb"
    xyz_file = "test_host_qm_coords.xyz"
    coordinate_file = "test_host_qm_coordinates.txt"
    get_get_xyz_host = qmmmrebind.parameterize.ParameterizeHost(xyz_file = xyz_file, coordinate_file = coordinate_file, unprocessed_hessian_file = "", bond_list_file = "", angle_list_file = "", hessian_file= "", atom_names_file = "", bond_parameter_file = "", vibrational_scaling = "", angle_parameter_file = "", charge_parameter_file = "", host_qm_pdb = host_qm_pdb)
    get_get_xyz_host.get_xyz() 

def test_get_unprocessed_hessian_host(): 
    host_qm_pdb = "test_host_qm.pdb"
    unprocessed_hessian_file = "test_host_qm_unprocessed_hessian.txt"
    get_get_unprocessed_hessian = qmmmrebind.parameterize.ParameterizeHost(xyz_file = "", coordinate_file = "", unprocessed_hessian_file = unprocessed_hessian_file, bond_list_file = "", angle_list_file = "", hessian_file= "", atom_names_file = "", bond_parameter_file = "", vibrational_scaling = "", angle_parameter_file = "", charge_parameter_file = "", host_qm_pdb = host_qm_pdb)
    get_get_unprocessed_hessian.get_unprocessed_hessian()

#drop in a log file here
def test_copy_log_file_host():
    source_ = get_data_filename("test_host_qm.log") 
    destination_pwd = os.getcwd()
    destination_file = source_.split('/')[-1]
    destination_ = destination_pwd + "/" + destination_file
    qmmmrebind.parameterize.copy_file(source = source_, destination = destination_)
 
def test_get_bond_angles_host(): 
    host_qm_pdb = "test_host_qm.pdb"
    bond_list_file = "test_host_qm_bond_list.txt"
    angle_list_file = "test_host_qm_angle_list"
    get_get_bond_angles = qmmmrebind.parameterize.ParameterizeHost(xyz_file = "", coordinate_file = "", unprocessed_hessian_file = "", bond_list_file = bond_list_file, angle_list_file = angle_list_file, hessian_file= "", atom_names_file = "", bond_parameter_file = "", vibrational_scaling = "", angle_parameter_file = "", charge_parameter_file = "", host_qm_pdb = host_qm_pdb)
    get_get_bond_angles.get_bond_angles() 

def test_get_hessian_host(): 
    host_qm_pdb = "test_host_qm.pdb"
    unprocessed_hessian_file = "test_host_qm_unprocessed_hessian.txt"
    hessian_file = "test_host_qm_hessian.txt"
    get_get_hessian = qmmmrebind.parameterize.ParameterizeHost(xyz_file = "", coordinate_file = "", unprocessed_hessian_file = unprocessed_hessian_file, bond_list_file = "", angle_list_file = "", hessian_file = hessian_file, atom_names_file = "", bond_parameter_file = "", vibrational_scaling = "", angle_parameter_file = "", charge_parameter_file = "", host_qm_pdb = host_qm_pdb)
    get_get_hessian.get_hessian() 

def test_get_atom_names_host(): 
    host_qm_pdb = "test_host_qm.pdb"
    atom_names_file = "test_host_qm_atom_names.txt"
    get_get_atom_names = qmmmrebind.parameterize.ParameterizeHost(xyz_file = "", coordinate_file = "", unprocessed_hessian_file = "", bond_list_file = "", angle_list_file = "", hessian_file = "", atom_names_file = atom_names_file, bond_parameter_file = "", vibrational_scaling = "", angle_parameter_file = "", charge_parameter_file = "", host_qm_pdb = host_qm_pdb)
    get_get_atom_names.get_atom_names() 

def test_get_bond_angle_params_host(): 
    host_qm_pdb = "test_host_qm.pdb"
    coordinate_file = "test_host_qm_coordinates.txt"
    hessian_file = "test_host_qm_hessian.txt"
    bond_list_file = "test_host_qm_bond_list.txt"
    atom_names_file = "test_host_qm_atom_names.txt"
    bond_parameter_file = "test_host_qm_bonds.txt"
    vibrational_scaling = 1.00
    angle_list_file = "test_host_qm_angle_list"
    angle_parameter_file = "test_host_qm_angles.txt"
    get_get_bond_angle_params = qmmmrebind.parameterize.ParameterizeHost(xyz_file = "", coordinate_file = coordinate_file, unprocessed_hessian_file = "", bond_list_file = bond_list_file, angle_list_file = angle_list_file, hessian_file = hessian_file, atom_names_file = atom_names_file, bond_parameter_file = bond_parameter_file, vibrational_scaling = vibrational_scaling, angle_parameter_file = angle_parameter_file, charge_parameter_file = "", host_qm_pdb = host_qm_pdb)
    get_get_bond_angle_params.get_bond_angle_params() 

def test_get_charges_host(): 
    host_qm_pdb = "test_host_qm.pdb"
    charge_parameter_file = "test_host_qm_surround_charges.txt"
    get_get_charges = qmmmrebind.parameterize.ParameterizeHost(xyz_file = "", coordinate_file = "", unprocessed_hessian_file = "", bond_list_file = "", angle_list_file = "", hessian_file= "", atom_names_file = "", bond_parameter_file = "", vibrational_scaling = "", angle_parameter_file = "", charge_parameter_file = charge_parameter_file, host_qm_pdb = host_qm_pdb)
    get_get_charges.get_charges() 

##############################GuestAmberXMLAmber############################## 
def test_generate_xml_from_charged_pdb_sdf(): 
    system_pdb = "test_guest_init_ii.pdb"
    system_init_sdf = "test_guest_init.sdf"
    system_sdf = "test_guest.sdf"
    num_charge_atoms = 1
    index_charge_atom_1 = 9
    charge_atom_1 = 1
    system_xml = "test_guest_init.xml"
    get_generate_xml_from_charged_pdb_sdf = qmmmrebind.parameterize.GuestAmberXMLAmber(system_pdb = system_pdb, system_mol2 = "", system_in = "", charge = "", system_frcmod = "", prmtop_system = "", inpcrd_system = "", system_leap = "", system_xml = system_xml, system_smi = "", system_sdf = system_sdf, system_init_sdf = system_init_sdf, num_charge_atoms = num_charge_atoms, index_charge_atom_1 = index_charge_atom_1, charge_atom_1 = charge_atom_1, index_charge_atom_2 = "", charge_atom_2 = "", charge_parameter_file = "", system_qm_pdb = "", bond_parameter_file = "", angle_parameter_file = "", system_qm_params_file = "", reparameterised_intermediate_system_xml_file = "", system_xml_non_bonded_file = "", system_xml_non_bonded_reparams_file = "", reparameterised_system_xml_file = "", non_reparameterised_system_xml_file = "", prmtop_system_non_params = "", inpcrd_system_non_params = "", prmtop_system_params = "", inpcrd_system_params = "", load_topology = "")
    get_generate_xml_from_charged_pdb_sdf.generate_xml_from_charged_pdb_sdf() 

def test_write_system_params(): 
    charge_parameter_file = "test_guest_charges.txt"
    system_qm_pdb = "test_guest_init_ii.pdb"
    bond_parameter_file = "test_guest_bonds.txt"
    angle_parameter_file = "test_guest_angles.txt"
    system_qm_params_file = "test_guest_qm_params.txt"
    get_write_system_params = qmmmrebind.parameterize.GuestAmberXMLAmber(system_pdb = "", system_mol2 = "", system_in = "", charge = "", system_frcmod = "", prmtop_system = "", inpcrd_system = "", system_leap = "", system_xml = "", system_smi = "", system_sdf = "", system_init_sdf = "", num_charge_atoms = "", index_charge_atom_1 = "", charge_atom_1 = "", index_charge_atom_2 = "", charge_atom_2 = "", charge_parameter_file = charge_parameter_file, system_qm_pdb = system_qm_pdb, bond_parameter_file = bond_parameter_file, angle_parameter_file = angle_parameter_file, system_qm_params_file = system_qm_params_file, reparameterised_intermediate_system_xml_file = "", system_xml_non_bonded_file = "", system_xml_non_bonded_reparams_file = "", reparameterised_system_xml_file = "", non_reparameterised_system_xml_file = "", prmtop_system_non_params = "", inpcrd_system_non_params = "", prmtop_system_params = "", inpcrd_system_params = "", load_topology = "")
    get_write_system_params.write_system_params() 

def test_write_reparameterised_system_xml(): 
    system_qm_params_file = "test_guest_qm_params.txt",
    system_xml = "test_guest_init.xml" 
    reparameterised_intermediate_system_xml_file =  "test_guest_intermediate_reparameterised.xml"
    system_qm_params_file = "test_guest_qm_params.txt" 
    system_xml_non_bonded_file = "test_guest_xml_non_bonded.txt"
    system_xml_non_bonded_reparams_file = "test_guest_xml_non_bonded_reparams.txt"
    reparameterised_system_xml_file = "test_guest_reparameterised.xml"
    get_write_reparameterised_system_xml = qmmmrebind.parameterize.GuestAmberXMLAmber(system_pdb = "", system_mol2 = "", system_in = "", charge = "", system_frcmod = "", prmtop_system = "", inpcrd_system = "", system_leap = "", system_xml = system_xml, system_smi = "", system_sdf = "", system_init_sdf = "", num_charge_atoms = "", index_charge_atom_1 = "", charge_atom_1 = "", index_charge_atom_2 = "", charge_atom_2 = "", charge_parameter_file = "", system_qm_pdb = "", bond_parameter_file = "", angle_parameter_file = "", system_qm_params_file = system_qm_params_file, reparameterised_intermediate_system_xml_file = reparameterised_intermediate_system_xml_file, system_xml_non_bonded_file = system_xml_non_bonded_file, system_xml_non_bonded_reparams_file = system_xml_non_bonded_reparams_file  , reparameterised_system_xml_file = reparameterised_system_xml_file, non_reparameterised_system_xml_file = "", prmtop_system_non_params = "", inpcrd_system_non_params = "", prmtop_system_params = "", inpcrd_system_params = "", load_topology = "")
    get_write_reparameterised_system_xml.write_reparameterised_system_xml() 

def test_save_amber_params(): 
    load_topology = "openmm"
    system_pdb = "test_guest_init_ii.pdb"
    non_reparameterised_system_xml_file = "test_guest_init.xml"
    prmtop_system_non_params = "test_guest_non_params.prmtop"
    inpcrd_system_non_params =  "test_guest_non_params.inpcrd"
    reparameterised_system_xml_file = "test_guest_reparameterised.xml"
    inpcrd_system_params = "test_guest_params.inpcrd"
    prmtop_system_params = "test_guest_params.prmtop"
    get_save_amber_params = qmmmrebind.parameterize.GuestAmberXMLAmber(system_pdb = system_pdb, system_mol2 = "", system_in = "", charge = "", system_frcmod = "", prmtop_system = "", inpcrd_system = "", system_leap = "", system_xml = "", system_smi = "", system_sdf = "", system_init_sdf = "", num_charge_atoms = "", index_charge_atom_1 = "", charge_atom_1 = "", index_charge_atom_2 = "", charge_atom_2 = "", charge_parameter_file = "", system_qm_pdb = "", bond_parameter_file = "", angle_parameter_file = "", system_qm_params_file = "", reparameterised_intermediate_system_xml_file = "", system_xml_non_bonded_file = "", system_xml_non_bonded_reparams_file = "", reparameterised_system_xml_file = reparameterised_system_xml_file, non_reparameterised_system_xml_file = non_reparameterised_system_xml_file, prmtop_system_non_params = prmtop_system_non_params, inpcrd_system_non_params = inpcrd_system_non_params, prmtop_system_params = prmtop_system_params, inpcrd_system_params = inpcrd_system_params, load_topology = load_topology)
    get_save_amber_params.save_amber_params() 

##############################HostAmberXMLAmber############################## 
def test_serialize_system(): 
    system_pdb = "test_host.pdb"
    sim_output = "test_sim_output.pdb"
    sim_steps = 100
    system_xml = "test_host.xml"
    get_serialize_system = qmmmrebind.parameterize.HostAmberXMLAmber(system_pdb = system_pdb, system_xml = system_xml, sim_output = sim_output, sim_steps = sim_steps, charge_parameter_file = "", system_qm_pdb = "", bond_parameter_file = "", angle_parameter_file = "", system_qm_params_file = "", reparameterised_intermediate_system_xml_file = "", system_xml_non_bonded_file = "", system_xml_non_bonded_reparams_file = "", reparameterised_system_xml_file = "", non_reparameterised_system_xml_file = "", prmtop_system_non_params = "", inpcrd_system_non_params = "", prmtop_system_params = "", inpcrd_system_params = "", load_topology = "")
    get_serialize_system.serialize_system() 

def test_write_system_params_host(): 
    charge_parameter_file = "test_host_qm_surround_charges.txt"
    system_qm_pdb = "test_host_qm.pdb"
    bond_parameter_file = "test_host_qm_bonds.txt"
    angle_parameter_file = "test_host_qm_angles.txt"
    system_qm_params_file = "test_host_qm_params.txt"
    get_write_system_params = qmmmrebind.parameterize.HostAmberXMLAmber(system_pdb = "", system_xml = "", sim_output = "", sim_steps = "", charge_parameter_file = charge_parameter_file, system_qm_pdb = system_qm_pdb, bond_parameter_file = bond_parameter_file, angle_parameter_file = angle_parameter_file, system_qm_params_file = system_qm_params_file, reparameterised_intermediate_system_xml_file = "", system_xml_non_bonded_file = "", system_xml_non_bonded_reparams_file = "", reparameterised_system_xml_file = "", non_reparameterised_system_xml_file = "", prmtop_system_non_params = "", inpcrd_system_non_params = "", prmtop_system_params = "", inpcrd_system_params = "", load_topology = "")
    get_write_system_params.write_system_params() 

def test_write_reparameterised_system_xml_host(): 
    system_qm_params_file = "test_host_qm_params.txt"
    system_xml = "test_host.xml"
    reparameterised_intermediate_system_xml_file = "test_host_intermediate_reparameterised.xml"
    system_xml_non_bonded_file = "test_host_xml_non_bonded.txt"
    system_xml_non_bonded_reparams_file = "test_host_xml_non_bonded_reparams.txt"
    reparameterised_system_xml_file = "test_host_reparameterised.xml"
    get_write_reparameterised_system_xml = qmmmrebind.parameterize.HostAmberXMLAmber(system_pdb = "", system_xml = system_xml, sim_output = "", sim_steps = "", charge_parameter_file = "", system_qm_pdb = "", bond_parameter_file = "", angle_parameter_file = "", system_qm_params_file = system_qm_params_file, reparameterised_intermediate_system_xml_file = reparameterised_intermediate_system_xml_file, system_xml_non_bonded_file = system_xml_non_bonded_file, system_xml_non_bonded_reparams_file = system_xml_non_bonded_reparams_file, reparameterised_system_xml_file = reparameterised_system_xml_file, non_reparameterised_system_xml_file = "", prmtop_system_non_params = "", inpcrd_system_non_params = "", prmtop_system_params = "", inpcrd_system_params = "", load_topology = "")
    get_write_reparameterised_system_xml.write_reparameterised_system_xml() 

def test_save_amber_params_host(): 
    load_topology = "openmm"
    system_pdb = "test_host.pdb"
    non_reparameterised_system_xml_file = "test_host.xml"
    prmtop_system_non_params = "test_host_non_params.prmtop"
    inpcrd_system_non_params = "test_host_non_params.inpcrd"
    reparameterised_system_xml_file = "test_host_reparameterised.xml"
    prmtop_system_params = "host_params.prmtop"
    inpcrd_system_params = "test_host_params.inpcrd"
    get_save_amber_params = qmmmrebind.parameterize.HostAmberXMLAmber(system_pdb = system_pdb, system_xml = "", sim_output = "", sim_steps = "", charge_parameter_file = "", system_qm_pdb = "", bond_parameter_file = "", angle_parameter_file = "", system_qm_params_file = "", reparameterised_intermediate_system_xml_file = "", system_xml_non_bonded_file = "", system_xml_non_bonded_reparams_file = "", reparameterised_system_xml_file = reparameterised_system_xml_file, non_reparameterised_system_xml_file = non_reparameterised_system_xml_file, prmtop_system_non_params = prmtop_system_non_params, inpcrd_system_non_params = inpcrd_system_non_params, prmtop_system_params = prmtop_system_params, inpcrd_system_params = inpcrd_system_params, load_topology = load_topology)
    get_save_amber_params.save_amber_params() 

##############################RemoveTestFiles############################## 
def test_remove_files():
    command = "rm -rf __pycache__ test_host_mm.txt test_guest_coord.txt test_host.pdb test_guest_init_ii.pdb test_host_qm.pdb test_guest_init.pdb test_host_qm.txt test_host_mm.pdb test_mm.pdb test_host_mm_region_ii.pdb test_qm.pdb test_host_mm_region_ii.txt test_residue_list.txt test_host_mm_region_i.pdb test_system.pdb test_host_mm_region_i.txt test_guest_bonds.txt test_guest_charges.txt test_guest_coordinates.txt test_guest_coords.xyz test_guest_hessian.txt test_guest_init_ii.fchk test_guest_init_ii.log test_guest_angle_list.txt test_guest_unprocessed_hessian.txt test_guest_angles.txt test_proper_dihedrals.txt test_guest_atom_names.txt test_guest_bond_list.txt test_host_qm_coordinates.txt test_host_qm_coords.xyz test_host_qm.fchk test_host_qm_unprocessed_hessian.txt test_host_qm_angle_list test_host_qm_bond_list.txt test_host_qm.log test_host_qm_angles.txt test_host_qm_atom_names.txt test_host_qm_bonds.txt test_host_qm_hessian.txt test_guest_qm_surround_charges.txt test_host_guest.com test_host_guest.log test_host_qm_surround_charges.txt test_guest_qm_atom_surround_charges.txt test_guest_init.sdf test_guest_init.xml test_guest.sdf test_guest_qm_params.txt test_guest_intermediate_reparameterised.xml test_guest_reparameterised.xml test_guest_xml_non_bonded_reparams.txt test_guest_xml_non_bonded.txt test_guest_non_params.inpcrd test_guest_non_params.prmtop test_guest_params.inpcrd test_guest_params.prmtop test_host.xml test_host_qm_params.txt test_host_intermediate_reparameterised.xml test_host_reparameterised.xml test_host_xml_non_bonded_reparams.txt test_host_xml_non_bonded.txt test_host_non_params.inpcrd test_host_non_params.prmtop test_host_params.inpcrd host_params.prmtop"
    os.system(command) 









