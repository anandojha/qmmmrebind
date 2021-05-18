"""
Unit and regression test for the qmmmrebind package.
"""
# Import package, test suite, and other packages as needed
import qmmmrebind
import pytest
import sys
import os
from .utils import get_data_filename

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
    filename = get_data_filename("data_input_file_first.dat") 
    word = "github"
    ret = qmmmrebind.parameterize.search_in_file(file = filename, word = word)
    assert ret[0][0] == 3

def test_find_word_in_file():
    """Test if a word is in a file"""
    filename = get_data_filename("data_input_file_first.dat") 
    word = "github"
    ret = qmmmrebind.parameterize.search_in_file(file = filename, word = word)
    assert ret[0][0] == 3

def test_method_a():
    ret = qmmmrebind.parameterize.test_class(a = 1, b = 2, c = " ")
    c = ret.method_a()
    assert c == 3

def test_clean_up():
    init_pdb = get_data_filename("data_input_file_second.pdb")
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


def test_realign_guest():   
    guest_init_pdb = "test_guest_init.pdb"
    guest_pdb = "test_guest_init_II.pdb"
    num_residues = 2
    get_realign_guest = qmmmrebind.parameterize.PrepareQMMM(init_pdb = "", cleaned_pdb = "", guest_init_pdb = guest_init_pdb, host_pdb = "", guest_resname = "", guest_pdb = guest_pdb, guest_xyz = "", distance = "", residue_list = "", host_qm_atoms = "", host_mm_atoms = "", host_qm_pdb = "", host_mm_pdb = "", qm_pdb = "", mm_pdb = "", host_mm_region_I_atoms = "", host_mm_region_II_atoms = "", host_mm_region_I_pdb = "", host_mm_region_II_pdb = "", num_residues = num_residues)
    get_realign_guest.realign_guest()    















