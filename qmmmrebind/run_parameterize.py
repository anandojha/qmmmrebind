from parameterize import *
########################################################################################################################################################################################################
system = PrepareQMMM(init_pdb = "test_system.pdb", cleaned_pdb = "system.pdb", guest_init_pdb = "guest_init.pdb", host_pdb = "host.pdb", guest_resname = "BEN", guest_pdb = "guest_init_II.pdb", guest_xyz = "guest_coord.txt", distance = 6.0, residue_list = "residue_list.txt", host_qm_atoms = "host_qm.txt", host_mm_atoms = "host_mm.txt", host_qm_pdb = "host_qm.pdb", host_mm_pdb = "host_mm.pdb", qm_pdb = "qm.pdb", mm_pdb = "mm.pdb", host_mm_region_I_atoms = "host_mm_region_I.txt", host_mm_region_II_atoms = "host_mm_region_II.txt", host_mm_region_I_pdb = "host_mm_region_I.pdb", host_mm_region_II_pdb = "host_mm_region_II.pdb", num_residues = 2)
########################################################################################################################################################################################################
qm_guest = PrepareGaussianGuest(guest_pdb = "guest_init_II.pdb", n_processors = 12, memory = 50, charge = 1, multiplicity = 1, functional = "B3LYP", basis_set = "6-31G", optimisation = "OPT", frequency = "FREQ", add_keywords_I = "Integral=(Grid=UltraFine)", add_keywords_II = "Pop(MK,ReadRadii)", add_keywords_III = "IOp(6/33=2,6/42=6)", gauss_out_file = "guest.out", fchk_out_file = "guest_fchk.out")
########################################################################################################################################################################################################
qm_system = PrepareGaussianHostGuest(guest_pdb = "guest_init_II.pdb", host_qm_pdb = "host_qm.pdb", n_processors = 12, memory = 50, charge = 1, multiplicity = 1, functional = "B3LYP", basis_set = "6-31G", optimisation = "", frequency = "", add_keywords_I = "Integral=(Grid=UltraFine)", add_keywords_II = "Pop(MK,ReadRadii)", add_keywords_III = "IOp(6/33=2,6/42=6)", gauss_system_out_file = "system_qm.out", fchk_system_out_file = "system_qm_fchk.out", host_guest_input = "host_guest.com", qm_guest_charge_parameter_file = "guest_qm_surround_charges.txt", qm_host_charge_parameter_file = "host_qm_surround_charges.txt", qm_guest_atom_charge_parameter_file = "guest_qm_atom_surround_charges.txt")
########################################################################################################################################################################################################
params_guest = ParameterizeGuest(xyz_file = "guest_coords.xyz", coordinate_file = "guest_coordinates.txt", unprocessed_hessian_file = "guest_unprocessed_hessian.txt", bond_list_file = "guest_bond_list.txt", angle_list_file = "guest_angle_list.txt", hessian_file = "guest_hessian.txt", atom_names_file = "guest_atom_names.txt", bond_parameter_file = "guest_bonds.txt", vibrational_scaling = 1.00, angle_parameter_file = "guest_angles.txt", charge_parameter_file = "guest_charges.txt", guest_pdb = "guest_init_II.pdb", proper_dihedral_file = "proper_dihedrals.txt")
########################################################################################################################################################################################################
off_xml_guest = OffXMLGuest(guest_pdb = "guest_init_II.pdb", sdf_file = "guest.sdf", xml_off_file = "guest_off.xml", atom_names_file = "guest_atom_names.txt", torsion_off_file = "guest_torsions_off.txt", lj_file = "guest_lj_off.txt", charge = 1)
########################################################################################################################################################################################################
xml_guest = XMLGuest(charge_parameter_file = "guest_qm_surround_charges.txt", atom_charge_file = "guest_qm_atom_surround_charges.txt", xml_file = "guest.xml", guest_pdb = "guest_init_II.pdb", bond_parameter_file = "guest_bonds.txt", angle_parameter_file = "guest_angles.txt", torsion_off_file = "guest_torsions_off.txt", lj_file = "guest_lj_off.txt", coulomb14scale = "0.833333", lj14scale = "0.5", residue_name = "BEN")
########################################################################################################################################################################################################
connect_guest = ConnectPdbGuest(xyz_file = "guest_coords.xyz", bond_parameter_file = "guest_bonds.txt", conect_file = "guest_conect.txt", guest_pdb = "guest_init_II.pdb", conect_guest_pdb = "guest.pdb", xml_file = "guest.xml", sim_output = "guest_openmm_sim.out", sim_steps = 10000)
########################################################################################################################################################################################################
torsion_drive_guest = TorsionDriveSims(xyz_file = "guest_coords.xyz", memory = 50, charge = 1, multiplicity = 1, functional = "B3LYP", basis_set = "6-31G", psi_input_file = "input.dat", dihedral_text_file = "dihedrals.txt", proper_dihedral_file = "proper_dihedrals.txt", tor_dir = "torsion_dir", dihedral_interval = 15, engine = "openmm", xml_file = "guest.xml", conect_guest_pdb = "guest.pdb", atom_names_file = "guest_atom_names.txt" ) 
########################################################################################################################################################################################################
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
########################################################################################################################################################################################################
qm_guest.write_input()
qm_guest.run_gaussian()
qm_guest.get_fchk()
########################################################################################################################################################################################################
qm_system.write_input()
qm_system.run_gaussian()
qm_system.get_fchk()
qm_system.get_qm_host_guest_charges()
########################################################################################################################################################################################################
params_guest.get_xyz()
params_guest.get_unprocessed_hessian()
params_guest.get_bond_angles()
params_guest.get_hessian()
params_guest.get_atom_names()
params_guest.get_bond_angle_params()
params_guest.get_charges()
params_guest.get_proper_dihedrals()
########################################################################################################################################################################################################
off_xml_guest.generate_off_xml_amber()
off_xml_guest.get_torsion_off()
off_xml_guest.get_lj_off()
########################################################################################################################################################################################################
xml_guest.get_qm_charges()
xml_guest.write_xml()
########################################################################################################################################################################################################
connect_guest.create_conect()
connect_guest.run_openmm()
########################################################################################################################################################################################################
torsion_drive_guest.create_tor_dir()
#torsion_drive_guest.run_torsion_drive()
########################################################################################################################################################################################################
qm_host = PrepareGaussianHost(host_qm_pdb = "host_qm.pdb", n_processors = 12, memory = 50, charge = 0, multiplicity = 1, functional = "B3LYP", basis_set = "6-31G", optimisation = "OPT", frequency = "FREQ", add_keywords_I = "Integral=(Grid=UltraFine)", add_keywords_II = "Pop(MK,ReadRadii)", add_keywords_III = "IOp(6/33=2,6/42=6)", gauss_out_file = "host_qm.out", fchk_out_file = "host_qm_fchk.out")
########################################################################################################################################################################################################
params_host = ParameterizeHost(xyz_file = "host_qm_coords.xyz", coordinate_file = "host_qm_coordinates.txt", unprocessed_hessian_file = "host_qm_unprocessed_hessian.txt", bond_list_file = "host_qm_bond_list.txt", angle_list_file = "host_qm_angle_list.txt", hessian_file = "host_qm_hessian.txt", atom_names_file = "host_qm_atom_names.txt", bond_parameter_file = "host_qm_bonds.txt", vibrational_scaling = 1.00, angle_parameter_file = "host_qm_angles.txt", charge_parameter_file = "host_qm_surround_charges.txt", host_qm_pdb = "host_qm.pdb", host_qm_params_file = "host_qm_params.txt", host_pdb = "host.pdb", ffxml = "amber14-all.xml", host_xml = "host_off.xml", sim_output = "sim_output.pdb", sim_steps = 10000, host_residue_name = "HST", host_singular_file = "host_HST.pdb", conect_pdb_txt = "host_connect.txt", conect_pdb_file = "host_connect.pdb", host_singular_xml_file = "host_HST.xml", coulomb14scale = "0.833333", lj14scale = "0.5", reparameterised_host_xml_file = "host_reparameterized.xml", reparams_host_file = "host_final.xml")
########################################################################################################################################################################################################
qm_host.write_input()
qm_host.run_gaussian()
qm_host.get_fchk()
########################################################################################################################################################################################################
params_host.get_xyz()
params_host.get_unprocessed_hessian()
params_host.get_bond_angles()
params_host.get_hessian()
params_host.get_atom_names()
params_host.get_bond_angle_params()
params_host.get_charges()
params_host.write_host_params()
params_host.serialise_host()
params_host.create_host_xml()
params_host.run_host_mm()
params_host.write_reparameterised_host_xml()
params_host.run_host_mm_qm()
########################################################################################################################################################################################################
os.system("rm -rf __pycache__")
########################################################################################################################################################################################################

