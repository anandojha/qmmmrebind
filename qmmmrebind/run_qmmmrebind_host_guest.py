# Host-Guest Systems Force Field Parameterization
# Note 1 : Functions to perform gaussian and torsiondrive calculations have been hashed since they were computed elsewhere.
# Note 2 : Place all the host and guest gaussian calculations in the directory called, "qm_data" and copy them to the parent directory while executing this script.
# Note 3 : All the torsion drive calculations are performed prior to running this script and placed in the torsion_dir.
# Note 4 : Topology files  (PDB and parm7/prmtop files are placed in the directory called "pdb_top_files".
##################################################################################################################################################################################################################
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import argparse
import shutil
import sys
import os
import re

import modules.base as base
import modules.constants as const
import modules.file_modify as file_modify
import modules.file_utilities as file_utilities
import modules.forcefield_utils as ff_utils
import modules.gaussian_inputs as gaussian_inputs
import modules.gaussian_outputs as gaussian_outputs
import modules.linear_algebra as linear_algebra
import modules.modified_seminario as modified_seminario
import modules.openmm_utilities as openmm_utilities
import modules.stage_inputs as stage_inputs
import modules.template_ff as template_ff
import modules.torsion_drive_inputs as torsion_inputs
import modules.torsion_drive_outputs as torsion_outputs
import modules.torsion_drive_run as torsion_run

pwd_qmmmrebind = "/home/aaojha/qmmmrebind/"  # PWD of the directory where qmmmrebind is installed
path_join = pwd_qmmmrebind + "qmmmrebind/"
module_path = os.path.abspath(os.path.join(path_join))
if module_path not in sys.path:
    sys.path.append(module_path)
#from parameterize import *
parser = argparse.ArgumentParser()
parser.add_argument("--pdbfile", type=str, help="System's PDB file to be parameterized")
parser.add_argument("--prmtopfile", type=str, help="System's initial topology file")
args = parser.parse_args()
##################################################################################################################################################################################################################
## Step 0 : PDB standardization and/or relaxation
qmmmrebindpdb = "qmmmrebind_init.pdb"
file_modify.singular_resid(pdbfile=args.pdbfile, qmmmrebind_init_file=qmmmrebindpdb)
#relax_init_structure(pdbfile=args.pdbfile, prmtopfile=args.prmtopfile, sim_output="output.pdb", sim_steps=10000, qmmmrebindpdb=qmmmrebindpdb)
##################################################################################################################################################################################################################
## Step I : Defining QM and MM regions
system = stage_inputs.PrepareQMMM(init_pdb=qmmmrebindpdb, distance=6.0, num_residues=8, guest_resname="APN")
## Step II : Ligand QM calculation
qm_guest = gaussian_inputs.PrepareGaussianGuest(functional="B3LYP", basis_set="cc-pVDZ")
## Step III : QM calculation of the receptor-ligand region for QM-derived charges
qm_system = gaussian_inputs.PrepareGaussianHostGuest(basis_set="6-31G(d,p)", functional="B3LYP", add_keywords_II="POP(CHELPG, REGULAR)", add_keywords_III="IOP(6/33=2) SCRF=PCM")
## Step IV : Ligand QM reparameterization
params_guest = gaussian_outputs.ParameterizeGuest(functional="B3LYP", basis_set="cc-pVDZ")
## Step V : Generation of reparameterized prmtop & inpcrd files for the ligand
system_guest = template_ff.GuestAmberXMLAmber(charge_parameter_file="guest_qm_surround_charges.txt") # Use charge_parameter_file="guest_qm_atom_surround_charges.txt" for bound state polarised charges and charge_parameter_file="guest_qm_surround_charges.txt" for non-polarised charges
## Step VI : Reparameterization of the torsion angles for the ligand
guest_torsion_params = torsion_run.TorsionDriveSims(method_torsion_drive="geometric")
## Step VII : Generation of torsional parameters for the ligand
guest_write_params = torsion_outputs.TorsionDriveParams()
## Step VIII : Receptor (Residues Surrounding the ligand) QM calculation
qm_host = gaussian_inputs.PrepareGaussianHost()
## Step IX : Receptor (Residues surrounding the ligand) reparameterization
params_host = gaussian_outputs.ParameterizeHost()
## Step X : System reparameterization
#hostguest_system = SystemAmberSystem(system_pdb = qmmmrebindpdb)
hostguest_system = ff_utils.SystemGuestAmberSystem(system_pdb = qmmmrebindpdb)
##################################################################################################################################################################################################################
## Step I : Defining QM and MM regions
system.clean_up()
system.create_host_guest()
system.realign_guest()
system.get_guest_coord()
system.get_qm_resids()
system.get_host_qm_mm_atoms()
system.save_host_pdbs()
system.get_host_mm_region_atoms()
system.save_host_mm_regions_pdbs()
system.get_qm_mm_regions()
##################################################################################################################################################################################################################
## Step II : Ligand QM calculation
qm_guest.write_input()
#qm_guest.run_gaussian()
#qm_guest.get_fchk()
##################################################################################################################################################################################################################
## Step III : QM calculation of the receptor-ligand region for QM-derived charges
qm_system.write_input()
#qm_system.run_gaussian()
#qm_system.get_fchk()
qm_system.get_qm_host_guest_charges()
##################################################################################################################################################################################################################
## Step IV : Ligand QM reparameterization
params_guest.get_xyz()
params_guest.get_unprocessed_hessian()
params_guest.get_bond_angles()
params_guest.get_hessian()
params_guest.get_atom_names()
params_guest.get_bond_angle_params()
params_guest.get_charges()
params_guest.get_proper_dihedrals()
##################################################################################################################################################################################################################
## Step V : Generation of reparameterized prmtop & inpcrd files for the ligand
system_guest.generate_xml_from_pdb_sdf()
system_guest.write_system_params()
system_guest.write_reparameterised_system_xml()
system_guest.save_amber_params()
system_guest.analyze_diff_energies()
##################################################################################################################################################################################################################
## Step VI : Reparameterization of the torsion angles for the ligand
#guest_torsion_params.write_torsion_drive_run_file()
#guest_torsion_params.write_tor_params_txt()
#guest_torsion_params.write_psi4_input()
#guest_torsion_params.create_non_H_bonded_torsion_drive_dir()
#guest_torsion_params.run_torsion_sim()
##################################################################################################################################################################################################################
## Step VII : Generation of torsional parameters for the ligand
guest_write_params.write_reparams_torsion_lines()
guest_write_params.write_torsional_reparams()
##################################################################################################################################################################################################################
## Step VIII : Receptor (Residues Surrounding the ligand) QM calculation
qm_host.write_input()
#qm_host.run_gaussian()
#qm_host.get_fchk()
##################################################################################################################################################################################################################
## Step IX : Receptor (Residues surrounding the ligand) reparameterization
params_host.get_xyz()
params_host.get_unprocessed_hessian()
params_host.get_bond_angles()
params_host.get_hessian()
params_host.get_atom_names()
params_host.get_bond_angle_params()
params_host.get_charges()
##################################################################################################################################################################################################################
## Step X : System reparameterization

"""
# Reparameterize the bound state residues and ligand (angle, bond, torsion and charge)
hostguest_system.generate_xml_from_prmtop()
hostguest_system.write_guest_params_non_zero()
hostguest_system.write_host_params()
hostguest_system.merge_qm_params()
remove_bad_angle_params(angle=1.00, k_angle=500)
hostguest_system.write_reparameterised_system_xml()
hostguest_system.write_torsional_reparams()
hostguest_system.save_amber_params()
"""

"""
# Reparameterize the bound state residues and ligand (angle, bond, and torsion)
hostguest_system.generate_xml_from_prmtop()
hostguest_system.write_guest_params_non_zero()
hostguest_system.write_host_params()
hostguest_system.merge_qm_params()
remove_bad_angle_params(angle=1.00, k_angle=500)
hostguest_system.write_intermediate_reparameterised_system_xml()
hostguest_system.write_torsional_reparams_intermediate()
hostguest_system.save_amber_params_non_qm_charges()
"""

# Reparameterize the ligand (angle, bond, torsion and charge)
hostguest_system = ff_utils.SystemGuestAmberSystem(system_pdb = qmmmrebindpdb)
hostguest_system.generate_xml_from_prmtop()
hostguest_system.write_guest_params_non_zero()
file_modify.remove_bad_angle_params(angle=1.00, k_angle=500)
hostguest_system.write_reparameterised_system_xml()
hostguest_system.write_torsional_reparams()
hostguest_system.save_amber_params()

"""
# Reparameterize the ligand (angle, bond, and torsion)
hostguest_system = SystemGuestAmberSystem(system_pdb = qmmmrebindpdb)
hostguest_system.generate_xml_from_prmtop()
hostguest_system.write_guest_params_non_zero()
remove_bad_angle_params(angle=1.00, k_angle=500)
hostguest_system.write_intermediate_reparameterised_system_xml()
hostguest_system.write_torsional_reparams_intermediate()
hostguest_system.save_amber_params_non_qm_charges()
"""
##################################################################################################################################################################################################################
## Step XI : Standardisation and OpenMM Run
file_utilities.change_names(inpcrd_file = "hostguest_params.inpcrd", prmtop_file = "hostguest_params.prmtop", pdb_file = qmmmrebindpdb)
file_modify.add_dim_prmtop(pdbfile = "system_qmmmrebind.pdb", prmtopfile = "system_qmmmrebind.prmtop")
#prmtop_calibration(prmtopfile = "system_qmmmrebind.prmtop", inpcrdfile = "system_qmmmrebind.inpcrd")
openmm_utilities.run_openmm_prmtop_pdb()
openmm_utilities.run_openmm_prmtop_inpcrd()
##################################################################################################################################################################################################################
## Step XII: Assignment of designated directories
#move_qmmmmrebind_files()
#move_qm_files()
#move_qmmmrebind_files()
print("Done")
##################################################################################################################################################################################################################
