---
description: |
    API documentation for modules: parameterize.

lang: en

classoption: oneside
geometry: margin=1in
papersize: a4

linkcolor: blue
links-as-notes: true
...


    
# Module `parameterize` {#parameterize}

QMMMReBind : Quantum Mechanics – Molecular Mechanics ( *QMMM* ) forcefield *Re*paramaterisation of the *Bind*ing site for the receptor-ligand complexes




    
## Functions


    
### Function `OPLS_LJ` {#parameterize.OPLS_LJ}




>     def OPLS_LJ(
>         system
>     )




    
### Function `copy_file` {#parameterize.copy_file}




>     def copy_file(
>         source,
>         destination
>     )




    
### Function `dihedral_energy` {#parameterize.dihedral_energy}




>     def dihedral_energy(
>         x,
>         k1,
>         k2,
>         k3,
>         k4=0
>     )




    
### Function `dot_product` {#parameterize.dot_product}




>     def dot_product(
>         u_PA,
>         eig_AB
>     )




    
### Function `error_function` {#parameterize.error_function}




>     def error_function(
>         delta_qm,
>         delta_mm
>     )




    
### Function `error_function_boltzmann` {#parameterize.error_function_boltzmann}




>     def error_function_boltzmann(
>         delta_qm,
>         delta_mm,
>         T
>     )




    
### Function `fit_params` {#parameterize.fit_params}




>     def fit_params(
>         qm_scan_file,
>         load_topology,
>         system_xml,
>         method
>     )




    
### Function `force_angle_constant` {#parameterize.force_angle_constant}




>     def force_angle_constant(
>         atom_A,
>         atom_B,
>         atom_C,
>         bond_lengths,
>         eigenvalues,
>         eigenvectors,
>         coords,
>         scaling_1,
>         scaling_2
>     )




    
### Function `force_angle_constant_special_case` {#parameterize.force_angle_constant_special_case}




>     def force_angle_constant_special_case(
>         atom_A,
>         atom_B,
>         atom_C,
>         bond_lengths,
>         eigenvalues,
>         eigenvectors,
>         coords,
>         scaling_1,
>         scaling_2
>     )




    
### Function `force_constant_bond` {#parameterize.force_constant_bond}




>     def force_constant_bond(
>         atom_A,
>         atom_B,
>         eigenvalues,
>         eigenvectors,
>         coords
>     )




    
### Function `gen_init_guess` {#parameterize.gen_init_guess}




>     def gen_init_guess(
>         qm_scan_file,
>         load_topology,
>         system_xml
>     )




    
### Function `generate_mm_pdbs` {#parameterize.generate_mm_pdbs}




>     def generate_mm_pdbs(
>         qm_scan_file,
>         template_pdb
>     )




    
### Function `generate_xml_from_charged_pdb_sdf` {#parameterize.generate_xml_from_charged_pdb_sdf}




>     def generate_xml_from_charged_pdb_sdf(
>         system_pdb,
>         system_init_sdf,
>         system_sdf,
>         num_charge_atoms,
>         index_charge_atom_1,
>         charge_atom_1,
>         system_xml
>     )


This function generates an openforcefield xml file from the pdb file via SDF file and openforcefield.

    
### Function `generate_xml_from_pdb_sdf` {#parameterize.generate_xml_from_pdb_sdf}




>     def generate_xml_from_pdb_sdf(
>         system_pdb,
>         system_sdf,
>         system_xml
>     )


This function generates an openforcefield xml file from the pdb file

    
### Function `get_dihedrals` {#parameterize.get_dihedrals}




>     def get_dihedrals(
>         qm_scan_file
>     )




    
### Function `get_mm_potential_energies` {#parameterize.get_mm_potential_energies}




>     def get_mm_potential_energies(
>         qm_scan_file,
>         load_topology,
>         system_xml
>     )




    
### Function `get_non_torsion_mm_energy` {#parameterize.get_non_torsion_mm_energy}




>     def get_non_torsion_mm_energy(
>         system_pdb,
>         load_topology,
>         system_xml
>     )




    
### Function `get_qm_energies` {#parameterize.get_qm_energies}




>     def get_qm_energies(
>         qm_scan_file
>     )




    
### Function `get_tor_params` {#parameterize.get_tor_params}




>     def get_tor_params(
>         qm_scan_file,
>         template_pdb,
>         load_topology,
>         system_xml,
>         method
>     )




    
### Function `get_torsional_lines` {#parameterize.get_torsional_lines}




>     def get_torsional_lines(
>         template_pdb,
>         system_xml,
>         qm_scan_file,
>         load_topology,
>         method,
>         dihedral_text_file
>     )




    
### Function `list_diff` {#parameterize.list_diff}




>     def list_diff(
>         list_1,
>         list_2
>     )




    
### Function `list_hartree_kcal` {#parameterize.list_hartree_kcal}




>     def list_hartree_kcal(
>         list_
>     )




    
### Function `list_kJ_kcal` {#parameterize.list_kJ_kcal}




>     def list_kJ_kcal(
>         list_
>     )




    
### Function `list_to_dict` {#parameterize.list_to_dict}




>     def list_to_dict(
>         lst
>     )




    
### Function `objective_function` {#parameterize.objective_function}




>     def objective_function(
>         k_array,
>         x,
>         delta_qm
>     )




    
### Function `remove_mm_files` {#parameterize.remove_mm_files}




>     def remove_mm_files(
>         qm_scan_file
>     )




    
### Function `reverse_list` {#parameterize.reverse_list}




>     def reverse_list(
>         lst
>     )




    
### Function `scale_list` {#parameterize.scale_list}




>     def scale_list(
>         list_
>     )




    
### Function `search_in_file` {#parameterize.search_in_file}




>     def search_in_file(
>         file: str,
>         word: str
>     ) ‑> list


Search for the given string in file and return lines containing that string along with line numbers

    
### Function `torsiondrive_input_to_xyz` {#parameterize.torsiondrive_input_to_xyz}




>     def torsiondrive_input_to_xyz(
>         psi_input_file,
>         xyz_file
>     )




    
### Function `u_PA_from_angles` {#parameterize.u_PA_from_angles}




>     def u_PA_from_angles(
>         atom_A,
>         atom_B,
>         atom_C,
>         coords
>     )




    
### Function `uniq` {#parameterize.uniq}




>     def uniq(
>         input_
>     )




    
### Function `unit_vector_N` {#parameterize.unit_vector_N}




>     def unit_vector_N(
>         u_BC,
>         u_AB
>     )




    
### Function `xyz_to_pdb` {#parameterize.xyz_to_pdb}




>     def xyz_to_pdb(
>         xyz_file,
>         coords_file,
>         template_pdb,
>         system_pdb
>     )





    
## Classes


    
### Class `GuestAmberXMLAmber` {#parameterize.GuestAmberXMLAmber}




>     class GuestAmberXMLAmber(
>         charge,
>         num_charge_atoms,
>         charge_atom_1,
>         index_charge_atom_1,
>         system_pdb='guest_init_II.pdb',
>         system_mol2='guest.mol2',
>         system_in='guest.in',
>         system_frcmod='guest.frcmod',
>         prmtop_system='guest.prmtop',
>         inpcrd_system='guest.inpcrd',
>         system_leap='guest.leap',
>         system_xml='guest_init.xml',
>         system_smi='guest.smi',
>         system_sdf='guest.sdf',
>         system_init_sdf='guest_init.sdf',
>         index_charge_atom_2=' ',
>         charge_atom_2=' ',
>         charge_parameter_file='guest_charges.txt',
>         system_qm_pdb='guest_init_II.pdb',
>         bond_parameter_file='guest_bonds.txt',
>         angle_parameter_file='guest_angles.txt',
>         system_qm_params_file='guest_qm_params.txt',
>         reparameterised_intermediate_system_xml_file='guest_intermediate_reparameterised.xml',
>         system_xml_non_bonded_file='guest_xml_non_bonded.txt',
>         system_xml_non_bonded_reparams_file='guest_xml_non_bonded_reparams.txt',
>         reparameterised_system_xml_file='guest_reparameterised.xml',
>         non_reparameterised_system_xml_file='guest_init.xml',
>         prmtop_system_non_params='guest_non_params.prmtop',
>         inpcrd_system_non_params='guest_non_params.inpcrd',
>         prmtop_system_params='guest_params.prmtop',
>         inpcrd_system_params='guest_params.inpcrd',
>         load_topology='openmm'
>     )










    
#### Methods


    
##### Method `analyze_diff_energies` {#parameterize.GuestAmberXMLAmber.analyze_diff_energies}




>     def analyze_diff_energies(
>         self
>     )




    
##### Method `generate_xml_antechamber` {#parameterize.GuestAmberXMLAmber.generate_xml_antechamber}




>     def generate_xml_antechamber(
>         self
>     )


This function generates an xml file from the pdb file through antechamber

    
##### Method `generate_xml_from_charged_pdb_sdf` {#parameterize.GuestAmberXMLAmber.generate_xml_from_charged_pdb_sdf}




>     def generate_xml_from_charged_pdb_sdf(
>         self
>     )


This function generates an openforcefield xml file from the pdb file via SDF file and openforcefield.

    
##### Method `generate_xml_from_doubly_charged_pdb_sdf` {#parameterize.GuestAmberXMLAmber.generate_xml_from_doubly_charged_pdb_sdf}




>     def generate_xml_from_doubly_charged_pdb_sdf(
>         self
>     )


This function generates an openforcefield xml file from the pdb file via SDF file and openforcefield.

    
##### Method `generate_xml_from_pdb_sdf` {#parameterize.GuestAmberXMLAmber.generate_xml_from_pdb_sdf}




>     def generate_xml_from_pdb_sdf(
>         self
>     )


This function generates an openforcefield xml file from the pdb file

    
##### Method `generate_xml_from_pdb_smi` {#parameterize.GuestAmberXMLAmber.generate_xml_from_pdb_smi}




>     def generate_xml_from_pdb_smi(
>         self
>     )


This function generates an openforcefield xml file from the pdb file

    
##### Method `save_amber_params` {#parameterize.GuestAmberXMLAmber.save_amber_params}




>     def save_amber_params(
>         self
>     )




    
##### Method `write_reparameterised_system_xml` {#parameterize.GuestAmberXMLAmber.write_reparameterised_system_xml}




>     def write_reparameterised_system_xml(
>         self
>     )




    
##### Method `write_system_params` {#parameterize.GuestAmberXMLAmber.write_system_params}




>     def write_system_params(
>         self
>     )


This function saves the parameters obtained from the QM log files in a text file.

    
### Class `HostAmberXMLAmber` {#parameterize.HostAmberXMLAmber}




>     class HostAmberXMLAmber(
>         system_pdb='host.pdb',
>         system_xml='host.xml',
>         sim_output='sim_output.pdb',
>         sim_steps=1000,
>         charge_parameter_file='host_qm_surround_charges.txt',
>         system_qm_pdb='host_qm.pdb',
>         bond_parameter_file='host_qm_bonds.txt',
>         angle_parameter_file='host_qm_angles.txt',
>         system_qm_params_file='host_qm_params.txt',
>         reparameterised_intermediate_system_xml_file='host_intermediate_reparameterised.xml',
>         system_xml_non_bonded_file='host_xml_non_bonded.txt',
>         system_xml_non_bonded_reparams_file='host_xml_non_bonded_reparams.txt',
>         reparameterised_system_xml_file='host_reparameterised.xml',
>         non_reparameterised_system_xml_file='host.xml',
>         prmtop_system_non_params='host_non_params.prmtop',
>         inpcrd_system_non_params='host_non_params.inpcrd',
>         prmtop_system_params='host_params.prmtop',
>         inpcrd_system_params='host_params.inpcrd',
>         load_topology='openmm'
>     )










    
#### Methods


    
##### Method `analyze_diff_energies` {#parameterize.HostAmberXMLAmber.analyze_diff_energies}




>     def analyze_diff_energies(
>         self
>     )




    
##### Method `save_amber_params` {#parameterize.HostAmberXMLAmber.save_amber_params}




>     def save_amber_params(
>         self
>     )




    
##### Method `serialize_system` {#parameterize.HostAmberXMLAmber.serialize_system}




>     def serialize_system(
>         self
>     )




    
##### Method `write_reparameterised_system_xml` {#parameterize.HostAmberXMLAmber.write_reparameterised_system_xml}




>     def write_reparameterised_system_xml(
>         self
>     )




    
##### Method `write_system_params` {#parameterize.HostAmberXMLAmber.write_system_params}




>     def write_system_params(
>         self
>     )


This function saves the parameters obtained from the QM log files in a text file.

    
### Class `MergeHostGuestTopology` {#parameterize.MergeHostGuestTopology}




>     class MergeHostGuestTopology(
>         host_prmtop,
>         guest_prmtop,
>         host_inpcrd,
>         guest_inpcrd,
>         system_prmtop,
>         system_inpcrd
>     )










    
#### Methods


    
##### Method `merge_topology_files` {#parameterize.MergeHostGuestTopology.merge_topology_files}




>     def merge_topology_files(
>         self
>     )




    
### Class `ParameterizeGuest` {#parameterize.ParameterizeGuest}




>     class ParameterizeGuest(
>         vibrational_scaling,
>         xyz_file='guest_coords.xyz',
>         coordinate_file='guest_coordinates.txt',
>         unprocessed_hessian_file='guest_unprocessed_hessian.txt',
>         bond_list_file='guest_bond_list.txt',
>         angle_list_file='guest_angle_list.txt',
>         hessian_file='guest_hessian.txt',
>         atom_names_file='guest_atom_names.txt',
>         bond_parameter_file='guest_bonds.txt',
>         angle_parameter_file='guest_angles.txt',
>         charge_parameter_file='guest_charges.txt',
>         guest_pdb='guest_init_II.pdb',
>         proper_dihedral_file='proper_dihedrals.txt'
>     )










    
#### Methods


    
##### Method `get_atom_names` {#parameterize.ParameterizeGuest.get_atom_names}




>     def get_atom_names(
>         self
>     )


This function saves a list of atom names from the formatted checkpoint file.

    
##### Method `get_bond_angle_params` {#parameterize.ParameterizeGuest.get_bond_angle_params}




>     def get_bond_angle_params(
>         self
>     )


This function saves the bond and angle parameter files obtained from the formatted checkpoint file.

    
##### Method `get_bond_angles` {#parameterize.ParameterizeGuest.get_bond_angles}




>     def get_bond_angles(
>         self
>     )


This function saves a text file of the bonds and angles from the gaussian log file.

    
##### Method `get_charges` {#parameterize.ParameterizeGuest.get_charges}




>     def get_charges(
>         self
>     )


This function saves the charges in a text file obtained from the Gaussian log file.

    
##### Method `get_hessian` {#parameterize.ParameterizeGuest.get_hessian}




>     def get_hessian(
>         self
>     )


This function extracts hessian from the unprocessed hessian and saves into a new file.

    
##### Method `get_proper_dihedrals` {#parameterize.ParameterizeGuest.get_proper_dihedrals}




>     def get_proper_dihedrals(
>         self
>     )


This function saves proper dihedral angles of the guest molecule in a text file.

    
##### Method `get_unprocessed_hessian` {#parameterize.ParameterizeGuest.get_unprocessed_hessian}




>     def get_unprocessed_hessian(
>         self
>     )


This function saves a text file of the unprocessed hessian from the formatted checkpoint file.

    
##### Method `get_xyz` {#parameterize.ParameterizeGuest.get_xyz}




>     def get_xyz(
>         self
>     )


This function saves a xyz file from the formatted checkpoint file.

    
### Class `ParameterizeHost` {#parameterize.ParameterizeHost}




>     class ParameterizeHost(
>         vibrational_scaling,
>         xyz_file='host_qm_coords.xyz',
>         coordinate_file='host_qm_coordinates.txt',
>         unprocessed_hessian_file='host_qm_unprocessed_hessian.txt',
>         bond_list_file='host_qm_bond_list.txt',
>         angle_list_file='host_qm_angle_list.txt',
>         hessian_file='host_qm_hessian.txt',
>         atom_names_file='host_qm_atom_names.txt',
>         bond_parameter_file='host_qm_bonds.txt',
>         angle_parameter_file='host_qm_angles.txt',
>         charge_parameter_file='host_qm_surround_charges.txt',
>         host_qm_pdb='host_qm.pdb'
>     )










    
#### Methods


    
##### Method `get_atom_names` {#parameterize.ParameterizeHost.get_atom_names}




>     def get_atom_names(
>         self
>     )


This function saves a list of atom names from the formatted checkpoint file.

    
##### Method `get_bond_angle_params` {#parameterize.ParameterizeHost.get_bond_angle_params}




>     def get_bond_angle_params(
>         self
>     )


This function saves the bond and angle parameter files obtained from the formatted checkpoint file.

    
##### Method `get_bond_angles` {#parameterize.ParameterizeHost.get_bond_angles}




>     def get_bond_angles(
>         self
>     )


This function saves a text file of the bonds and angles from the gaussian log file.

    
##### Method `get_charges` {#parameterize.ParameterizeHost.get_charges}




>     def get_charges(
>         self
>     )


This function saves the charges in a text file obtained from the Gaussian log file.

    
##### Method `get_hessian` {#parameterize.ParameterizeHost.get_hessian}




>     def get_hessian(
>         self
>     )


This function extracts hessian from the unprocessed hessian and saves into a new file.

    
##### Method `get_unprocessed_hessian` {#parameterize.ParameterizeHost.get_unprocessed_hessian}




>     def get_unprocessed_hessian(
>         self
>     )


This function saves a text file of the unprocessed hessian from the formatted checkpoint file.

    
##### Method `get_xyz` {#parameterize.ParameterizeHost.get_xyz}




>     def get_xyz(
>         self
>     )


This function saves a xyz file from the formatted checkpoint file.

    
### Class `PrepareGaussianGuest` {#parameterize.PrepareGaussianGuest}




>     class PrepareGaussianGuest(
>         charge,
>         multiplicity,
>         guest_pdb='guest_init_II.pdb',
>         n_processors=12,
>         memory=50,
>         functional='B3LYP',
>         basis_set='6-31G',
>         optimisation='OPT',
>         frequency='FREQ',
>         add_keywords_I='Integral=(Grid=UltraFine)',
>         add_keywords_II='Pop(MK,ReadRadii)',
>         add_keywords_III='IOp(6/33=2,6/42=6)',
>         gauss_out_file='guest.out',
>         fchk_out_file='guest_fchk.out'
>     )










    
#### Methods


    
##### Method `get_fchk` {#parameterize.PrepareGaussianGuest.get_fchk}




>     def get_fchk(
>         self
>     )


This function converts the checkpoint file file into the formatted chechkpoint file.

    
##### Method `run_gaussian` {#parameterize.PrepareGaussianGuest.run_gaussian}




>     def run_gaussian(
>         self
>     )


This function runs the gaussian QM calculation.

    
##### Method `write_input` {#parameterize.PrepareGaussianGuest.write_input}




>     def write_input(
>         self
>     )


This function prints out the commands section of the gaussian input file.

    
### Class `PrepareGaussianHost` {#parameterize.PrepareGaussianHost}




>     class PrepareGaussianHost(
>         charge,
>         multiplicity,
>         host_qm_pdb='host_qm.pdb',
>         n_processors=12,
>         memory=50,
>         functional='B3LYP',
>         basis_set='6-31G',
>         optimisation='OPT',
>         frequency='FREQ',
>         add_keywords_I='Integral=(Grid=UltraFine)',
>         add_keywords_II='Pop(MK,ReadRadii)',
>         add_keywords_III='IOp(6/33=2,6/42=6)',
>         gauss_out_file='host_qm.out',
>         fchk_out_file='host_qm_fchk.out'
>     )










    
#### Methods


    
##### Method `get_fchk` {#parameterize.PrepareGaussianHost.get_fchk}




>     def get_fchk(
>         self
>     )


This function converts the checkpoint file file into the formatted chechkpoint file.

    
##### Method `run_gaussian` {#parameterize.PrepareGaussianHost.run_gaussian}




>     def run_gaussian(
>         self
>     )


This function runs the gaussian QM calculation.

    
##### Method `write_input` {#parameterize.PrepareGaussianHost.write_input}




>     def write_input(
>         self
>     )


This function prints out the commands section of the gaussian input file.

    
### Class `PrepareGaussianHostGuest` {#parameterize.PrepareGaussianHostGuest}




>     class PrepareGaussianHostGuest(
>         charge,
>         multiplicity,
>         guest_pdb='guest_init_II.pdb',
>         host_qm_pdb='host_qm.pdb',
>         n_processors=12,
>         memory=50,
>         functional='B3LYP',
>         basis_set='6-31G',
>         optimisation='',
>         frequency='',
>         add_keywords_I='Integral=(Grid=UltraFine)',
>         add_keywords_II='Pop(MK,ReadRadii)',
>         add_keywords_III='IOp(6/33=2,6/42=6)',
>         gauss_system_out_file='system_qm.out',
>         fchk_system_out_file='system_qm_fchk.out',
>         host_guest_input='host_guest.com',
>         qm_guest_charge_parameter_file='guest_qm_surround_charges.txt',
>         qm_host_charge_parameter_file='host_qm_surround_charges.txt',
>         qm_guest_atom_charge_parameter_file='guest_qm_atom_surround_charges.txt'
>     )










    
#### Methods


    
##### Method `get_fchk` {#parameterize.PrepareGaussianHostGuest.get_fchk}




>     def get_fchk(
>         self
>     )


This function converts the checkpoint file file into the formatted chechkpoint file.

    
##### Method `get_qm_host_guest_charges` {#parameterize.PrepareGaussianHostGuest.get_qm_host_guest_charges}




>     def get_qm_host_guest_charges(
>         self
>     )


This function extracts charges and saves them separately for the host and guest

    
##### Method `run_gaussian` {#parameterize.PrepareGaussianHostGuest.run_gaussian}




>     def run_gaussian(
>         self
>     )


This function runs the gaussian QM calculation.

    
##### Method `write_input` {#parameterize.PrepareGaussianHostGuest.write_input}




>     def write_input(
>         self
>     )


This function prints out the commands section of the gaussian input file.

    
### Class `PrepareQMMM` {#parameterize.PrepareQMMM}




>     class PrepareQMMM(
>         init_pdb,
>         distance,
>         num_residues,
>         guest_resname,
>         cleaned_pdb='system.pdb',
>         guest_init_pdb='guest_init.pdb',
>         host_pdb='host.pdb',
>         guest_pdb='guest_init_II.pdb',
>         guest_xyz='guest_coord.txt',
>         residue_list='residue_list.txt',
>         host_qm_atoms='host_qm.txt',
>         host_mm_atoms='host_mm.txt',
>         host_qm_pdb='host_qm.pdb',
>         host_mm_pdb='host_mm.pdb',
>         qm_pdb='qm.pdb',
>         mm_pdb='mm.pdb',
>         host_mm_region_I_atoms='host_mm_region_I.txt',
>         host_mm_region_II_atoms='host_mm_region_II.txt',
>         host_mm_region_I_pdb='host_mm_region_I.pdb',
>         host_mm_region_II_pdb='host_mm_region_II.pdb'
>     )


A class used to segregate the QM and MM regions.

This class contains methods to remove the solvent, ions and all
entities that are exclusive of receptor and the ligand. It also
defines the Quantum Mechanical (QM) region and the Molecular
Mechanical (MM) region based upon the distance of the ligand
from the receptor and the chosen number of receptor residues. It
is also assumed that the initial PDB file will have the receptor
followed by the ligand.

...

#### Attributes

**```init_pdb```** :&ensp;<code>str</code>
:   Initial PDB file containing the receptor-ligand complex with
    solvent, ions, etc.


**```cleaned_pdb```** :&ensp;<code>str</code>
:   Formatted PDB file containing only the receptor and the ligand.
    (This file will be saved in the current working directory)


**```guest_init_pdb```** :&ensp;<code>str</code>
:   A separate ligand PDB file with atom numbers not beginning from 1.
    (This file will be saved in the current working directory)


**```host_pdb```** :&ensp;<code>str</code>
:   A separate receptor PDB file with atom numbers beginning from 1.


**```guest_resname```** :&ensp;<code>str</code>
:   Three letter residue ID for the ligand


**```guest_pdb```** :&ensp;<code>str</code>
:   Ligand PDB file with atom numbers beginning from 1.
    (This file will be saved in the current working directory)


**```guest_xyz```** :&ensp;<code>str</code>
:   A text file of the XYZ corordinates of the ligand.
    (This file will be saved in the current working directory)


**```distance```** :&ensp;<code>float</code>
:   The distance required to define the QM region of the receptor.
    This is the distance between the atoms of the ligand and the
    atoms of the receptor.


**```residue_list```** :&ensp;<code>str</code>
:   A text file of the residue numbers of the receptor within the
    proximity (as defined by the distance) from the ligand.
    (This file will be saved in the current working directory)


**```host_qm_atoms```** :&ensp;<code>str</code>
:   A text file of the atom numbers of the receptors in the QM
    region.
    (This file will be saved in the current working directory)


**```host_mm_atoms```** :&ensp;<code>str</code>
:   A text file of the atom numbers of the receptors in the MM
    region (all atoms except atoms in the QM region)
    (This file will be saved in the current working directory)


**```host_qm_pdb```** :&ensp;<code>str</code>
:   PDB file for the receptor's QM region.
    (This file will be saved in the current working directory)


**```host_mm_pdb```** :&ensp;<code>str</code>
:   PDB file for the receptor's MM region.
    (This file will be saved in the current working directory)


**```qm_pdb```** :&ensp;<code>str</code>
:   PDB file for the QM region (receptor's QM region and the
    ligand).
    (This file will be saved in the current working directory)


**```mm_pdb```** :&ensp;<code>str</code>
:   PDB file for the MM region.
    (This file will be saved in the current working directory)


**```host_mm_region_I_atoms```** :&ensp;<code>str</code>
:   A text file of the atom numbers of the receptors in the MM
    region preceeding the QM region.
    (This file will be saved in the current working directory)


**```host_mm_region_II_atoms```** :&ensp;<code>str</code>
:   A text file of the atom numbers of the receptors in the MM
    region following the QM region.
    (This file will be saved in the current working directory)


**```host_mm_region_I_pdb```** :&ensp;<code>str</code>
:   PDB file of the receptor in the MM region preceeding the
    QM region.
    (This file will be saved in the current working directory)


**```host_mm_region_II_pdb```** :&ensp;<code>str</code>
:   PDB file of the receptor in the MM region following the
    QM region.
    (This file will be saved in the current working directory)


**```num_residues```** :&ensp;<code>int</code>
:   Number of residues required in the QM region of the receptor.

#### Parameters

**```init_pdb```** :&ensp;<code>str</code>
:   Initial PDB file containing the receptor-ligand complex with
    solvent, ions, etc.


**```cleaned_pdb```** :&ensp;<code>str</code>
:   Formatted PDB file containing only the receptor and the ligand.
    (This file will be saved in the current working directory)


**```guest_init_pdb```** :&ensp;<code>str</code>
:   A separate ligand PDB file with atom numbers not beginning from 1.
    (This file will be saved in the current working directory)


**```host_pdb```** :&ensp;<code>str</code>
:   A separate receptor PDB file with atom numbers beginning from 1.


**```guest_resname```** :&ensp;<code>str</code>
:   Three letter residue ID for the ligand.


**```guest_pdb```** :&ensp;<code>str</code>
:   Ligand PDB file with atom numbers beginning from 1.
    (This file will be saved in the current working directory)


**```guest_xyz```** :&ensp;<code>str</code>
:   A text file of the XYZ corordinates of the ligand.
    (This file will be saved in the current working directory)


**```distance```** :&ensp;<code>float</code>
:   The distance required to define the QM region of the receptor.
    This is the distance between the atoms of the ligand and the
    atoms of the receptor.


**```residue_list```** :&ensp;<code>str</code>
:   A text file of the residue numbers of the receptor within the
    proximity (as defined by the distance) from the ligand.
    (This file will be saved in the current working directory)


**```host_qm_atoms```** :&ensp;<code>str</code>
:   A text file of the atom numbers of the receptors in the QM
    region.
    (This file will be saved in the current working directory)


**```host_mm_atoms```** :&ensp;<code>str</code>
:   A text file of the atom numbers of the receptors in the MM
    region (all atoms except atoms in the QM region)
    (This file will be saved in the current working directory)


**```host_qm_pdb```** :&ensp;<code>str</code>
:   PDB file for the receptor's QM region.
    (This file will be saved in the current working directory)


**```host_mm_pdb```** :&ensp;<code>str</code>
:   PDB file for the receptor's MM region.
    (This file will be saved in the current working directory)


**```qm_pdb```** :&ensp;<code>str</code>
:   PDB file for the QM region (receptor's QM region and the
    ligand).
    (This file will be saved in the current working directory)


**```mm_pdb```** :&ensp;<code>str</code>
:   PDB file for the MM region.
    (This file will be saved in the current working directory)


**```host_mm_region_I_atoms```** :&ensp;<code>str</code>
:   A text file of the atom numbers of the receptors in the MM
    region preceeding the QM region.
    (This file will be saved in the current working directory)


**```host_mm_region_II_atoms```** :&ensp;<code>str</code>
:   A text file of the atom numbers of the receptors in the MM
    region following the QM region.
    (This file will be saved in the current working directory)


**```host_mm_region_I_pdb```** :&ensp;<code>str</code>
:   PDB file of the receptor in the MM region preceeding the
    QM region.
    (This file will be saved in the current working directory)


**```host_mm_region_II_pdb```** :&ensp;<code>str</code>
:   PDB file of the receptor in the MM region following the
    QM region.
    (This file will be saved in the current working directory)


**```num_residues```** :&ensp;<code>int</code>
:   Number of residues required in the QM region of the receptor.









    
#### Methods


    
##### Method `clean_up` {#parameterize.PrepareQMMM.clean_up}




>     def clean_up(
>         self
>     )


Reads the given PDB file, removes all entities except the
receptor and ligand and saves a new pdb file.

    
##### Method `create_host_guest` {#parameterize.PrepareQMMM.create_host_guest}




>     def create_host_guest(
>         self
>     )


Saves separate receptor and ligand PDB files.

    
##### Method `get_guest_coord` {#parameterize.PrepareQMMM.get_guest_coord}




>     def get_guest_coord(
>         self
>     )


Saves a text file of the XYZ corordinates of the ligand.

    
##### Method `get_host_mm_region_atoms` {#parameterize.PrepareQMMM.get_host_mm_region_atoms}




>     def get_host_mm_region_atoms(
>         self
>     )


Saves a text file for the atoms of the receptor's MM region
preceding the QM region and saves another text file for the
atoms of the receptor's MM region folllowing the QM region.

    
##### Method `get_host_qm_mm_atoms` {#parameterize.PrepareQMMM.get_host_qm_mm_atoms}




>     def get_host_qm_mm_atoms(
>         self
>     )


Saves a text file of the atom numbers of the receptors in the QM
region and MM region separately.

    
##### Method `get_qm_mm_regions` {#parameterize.PrepareQMMM.get_qm_mm_regions}




>     def get_qm_mm_regions(
>         self
>     )


Saves separate PDB files for the QM and MM regions.
QM regions comprise the QM region of the receptor
and the entire ligand where the MM region comprise
the non-selected QM regions of the receptor.

    
##### Method `get_qm_resids` {#parameterize.PrepareQMMM.get_qm_resids}




>     def get_qm_resids(
>         self
>     )


Saves a text file of the residue numbers of the receptor within the
proximity (as defined by the distance) from the ligand.

    
##### Method `realign_guest` {#parameterize.PrepareQMMM.realign_guest}




>     def realign_guest(
>         self
>     )


Saves a ligand PDB file with atom numbers beginning from 1.

    
##### Method `save_host_mm_regions_pdbs` {#parameterize.PrepareQMMM.save_host_mm_regions_pdbs}




>     def save_host_mm_regions_pdbs(
>         self
>     )


Saves a PDB file for the receptor's MM region preceding
the QM region and saves another PDB file for the receptor's
MM region folllowing the QM region.

    
##### Method `save_host_pdbs` {#parameterize.PrepareQMMM.save_host_pdbs}




>     def save_host_pdbs(
>         self
>     )


Saves a PDB file for the receptor's QM region and MM
region separately.

    
### Class `RunOpenMMSims` {#parameterize.RunOpenMMSims}




>     class RunOpenMMSims(
>         system_prmtop,
>         system_inpcrd,
>         system_pdb,
>         system_output='sim_output.pdb',
>         sim_steps=1000
>     )










    
#### Methods


    
##### Method `run_openmm_prmtop_inpcrd` {#parameterize.RunOpenMMSims.run_openmm_prmtop_inpcrd}




>     def run_openmm_prmtop_inpcrd(
>         self
>     )




    
##### Method `run_openmm_prmtop_pdb` {#parameterize.RunOpenMMSims.run_openmm_prmtop_pdb}




>     def run_openmm_prmtop_pdb(
>         self
>     )




    
### Class `TorsionDriveParams` {#parameterize.TorsionDriveParams}




>     class TorsionDriveParams(
>         num_charge_atoms,
>         index_charge_atom_1,
>         charge_atom_1,
>         tor_dir='torsion_dir',
>         reparameterized_torsional_params_file='reparameterized_torsional_params.txt',
>         psi_input_file='torsion_drive_input.dat',
>         xyz_file='torsion_drive_input.xyz',
>         coords_file='torsion_drive_input.txt',
>         template_pdb='guest_init_II.pdb',
>         system_pdb='torsion_drive_input.pdb',
>         system_sdf='torsion_drive_input.sdf',
>         system_xml='torsion_drive_input.xml',
>         qm_scan_file='scan.xyz',
>         load_topology='openmm',
>         method='L-BFGS-B',
>         dihedral_text_file='dihedrals.txt',
>         system_init_sdf='torsion_drive_input_init.sdf',
>         reparameterised_system_xml_file='guest_reparameterised.xml',
>         reparameterised_torsional_system_xml_file='guest_torsional_reparameterized.xml'
>     )










    
#### Methods


    
##### Method `write_reparams_torsion_lines` {#parameterize.TorsionDriveParams.write_reparams_torsion_lines}




>     def write_reparams_torsion_lines(
>         self
>     )




    
##### Method `write_reparams_torsion_lines_charged` {#parameterize.TorsionDriveParams.write_reparams_torsion_lines_charged}




>     def write_reparams_torsion_lines_charged(
>         self
>     )




    
##### Method `write_torsional_reparams` {#parameterize.TorsionDriveParams.write_torsional_reparams}




>     def write_torsional_reparams(
>         self
>     )




    
### Class `TorsionDriveSims` {#parameterize.TorsionDriveSims}




>     class TorsionDriveSims(
>         charge,
>         multiplicity,
>         reparameterised_system_xml_file='guest_reparameterised.xml',
>         torsion_xml_file='guest_torsion_xml.txt',
>         xyz_file='guest_coords.xyz',
>         psi_input_file='torsion_drive_input.dat',
>         memory=50,
>         basis_set='STO-3G',
>         functional='BLYP',
>         iterations=2000,
>         method_torsion_drive='native_opt',
>         system_bonds_file='guest_bonds.txt',
>         tor_dir='torsion_dir',
>         dihedral_text_file='dihedrals.txt',
>         template_pdb='guest_init_II.pdb',
>         torsion_drive_run_file='run_command',
>         dihedral_interval=15,
>         engine='psi4',
>         energy_threshold=0.001
>     )










    
#### Methods


    
##### Method `create_non_H_bonded_torsion_drive_dir` {#parameterize.TorsionDriveSims.create_non_H_bonded_torsion_drive_dir}




>     def create_non_H_bonded_torsion_drive_dir(
>         self
>     )




    
##### Method `create_non_H_torsion_drive_dir` {#parameterize.TorsionDriveSims.create_non_H_torsion_drive_dir}




>     def create_non_H_torsion_drive_dir(
>         self
>     )




    
##### Method `create_torsion_drive_dir` {#parameterize.TorsionDriveSims.create_torsion_drive_dir}




>     def create_torsion_drive_dir(
>         self
>     )




    
##### Method `run_torsion_sim` {#parameterize.TorsionDriveSims.run_torsion_sim}




>     def run_torsion_sim(
>         self
>     )




    
##### Method `write_psi4_input` {#parameterize.TorsionDriveSims.write_psi4_input}




>     def write_psi4_input(
>         self
>     )




    
##### Method `write_tor_params_txt` {#parameterize.TorsionDriveSims.write_tor_params_txt}




>     def write_tor_params_txt(
>         self
>     )




    
##### Method `write_torsion_drive_run_file` {#parameterize.TorsionDriveSims.write_torsion_drive_run_file}




>     def write_torsion_drive_run_file(
>         self
>     )





-----
Generated by *pdoc* 0.9.2 (<https://pdoc3.github.io>).
