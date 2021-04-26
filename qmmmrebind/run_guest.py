from parameterize_guest import *
########################################################################################################################################################################################################
system = PrepareQMMM(init_pdb = "test_system.pdb", cleaned_pdb = "system.pdb", guest_init_pdb = "guest_init.pdb", host_pdb = "host.pdb", guest_resname = "BEN", guest_pdb = "guest_init_II.pdb", guest_xyz = "guest_coord.txt", distance = 6.0, residue_list = "residue_list.txt", host_qm_atoms = "host_qm.txt", host_mm_atoms = "host_mm.txt", host_qm_pdb = "host_qm.pdb", host_mm_pdb = "host_mm.pdb", qm_pdb = "qm.pdb", mm_pdb = "mm.pdb", host_mm_region_I_atoms = "host_mm_region_I.txt", host_mm_region_II_atoms = "host_mm_region_II.txt", host_mm_region_I_pdb = "host_mm_region_I.pdb", host_mm_region_II_pdb = "host_mm_region_II.pdb", num_residues = 2)
########################################################################################################################################################################################################
qm_system = PrepareGaussian(guest_pdb = "guest_init_II.pdb", n_processors = 12, memory = 50, charge = 1, multiplicity = 1, functional = "B3LYP", basis_set = "6-31G", optimisation = "OPT", frequency = "FREQ", add_keywords_I = "Integral=(Grid=UltraFine)", add_keywords_II = "Pop(MK,ReadRadii)", add_keywords_III = "IOp(6/33=2,6/42=6)", gauss_out_file = "guest.out", fchk_out_file = "guest_fchk.out")
########################################################################################################################################################################################################
params = Parameterize(xyz_file = "guest_coords.xyz", coordinate_file = "guest_coordinates.txt", unprocessed_hessian_file = "guest_unprocessed_hessian.txt", bond_list_file = "guest_bond_list.txt", angle_list_file = "guest_angle_list.txt", hessian_file = "guest_hessian.txt", atom_names_file = "guest_atom_names.txt", bond_parameter_file = "guest_bonds.txt", vibrational_scaling = 1.00, angle_parameter_file = "guest_angles.txt", charge_parameter_file = "guest_charges.txt", guest_pdb = "guest_init_II.pdb", proper_dihedral_file = "proper_dihedrals.txt")
########################################################################################################################################################################################################
off_xml = OffXML(guest_pdb = "guest_init_II.pdb", sdf_file = "guest.sdf", xml_off_file = "guest_off.xml", atom_names_file = "guest_atom_names.txt", torsion_off_file = "guest_torsions_off.txt", lj_file = "guest_lj_off.txt", charge = 1)
########################################################################################################################################################################################################
xml = XML(charge_parameter_file = "guest_charges.txt", atom_charge_file = "guest_atom_charges.txt", xml_file = "guest.xml", guest_pdb = "guest_init_II.pdb", bond_parameter_file = "guest_bonds.txt", angle_parameter_file = "guest_angles.txt", torsion_off_file = "guest_torsions_off.txt", lj_file = "guest_lj_off.txt", coulomb14scale = "0.833333", lj14scale = "0.5", residue_name = "BEN")
########################################################################################################################################################################################################
connect = ConnectPdb(xyz_file = "guest_coords.xyz", bond_parameter_file = "guest_bonds.txt", conect_file = "guest_conect.txt", guest_pdb = "guest_init_II.pdb", conect_guest_pdb = "guest.pdb", xml_file = "guest.xml", sim_output = "guest_openmm_sim.out", sim_steps = 10000)
########################################################################################################################################################################################################
torsion_drive = TorsionDriveSims(xyz_file = "guest_coords.xyz", memory = 50, charge = 1, multiplicity = 1, functional = "B3LYP", basis_set = "6-31G", psi_input_file = "input.dat", dihedral_text_file = "dihedrals.txt", proper_dihedral_file = "proper_dihedrals.txt", tor_dir = "torsion_dir", dihedral_interval = 15, engine = "openmm", xml_file = "guest.xml", conect_guest_pdb = "guest.pdb", atom_names_file = "guest_atom_names.txt" ) 
########################################################################################################################################################################################################
system_qm = PrepareGaussianHostGuest(guest_pdb = "guest_init_II.pdb", host_qm_pdb = "host_qm.pdb", n_processors = 12, memory = 50, charge = 1, multiplicity = 1, functional = "B3LYP", basis_set = "6-31G", optimisation = "OPT", frequency = "", add_keywords_I = "Integral=(Grid=UltraFine)", add_keywords_II = "Pop(MK,ReadRadii)", add_keywords_III = "IOp(6/33=2,6/42=6)", gauss_system_out_file = "system_qm.out", fchk_system_out_file = "system_qm_fchk.out", host_guest_input = "host_guest.com")
########################################################################################################################################################################################################
#system.clean_up()
#system.create_host_guest()
#system.realign_guest()
#system.get_guest_coord()
#system.get_qm_resids()
#system.get_host_qm_mm_atoms()
#system.save_host_pdbs()
#system.get_host_mm_region_atoms()
#system.save_host_mm_regions_pdbs()
#system.get_qm_mm_regions()
########################################################################################################################################################################################################
#qm_system.write_input()
#qm_system.run_gaussian()
#qm_system.get_fchk()
########################################################################################################################################################################################################
#params.get_xyz()
#params.get_unprocessed_hessian()
#params.get_bond_angles()
#params.get_hessian()
#params.get_atom_names()
#params.get_bond_angle_params()
#params.get_charges()
#params.get_proper_dihedrals()
########################################################################################################################################################################################################
#off_xml.generate_off_xml_amber()
#off_xml.get_torsion_off()
#off_xml.get_lj_off()
########################################################################################################################################################################################################
#xml.get_qm_charges()
#xml.write_xml()
########################################################################################################################################################################################################
#connect.create_conect()
#connect.run_openmm()
########################################################################################################################################################################################################
#torsion_drive.create_tor_dir()
#torsion_drive.run_torsion_drive()
########################################################################################################################################################################################################
system_qm.write_input()
system_qm.run_gaussian()
system_qm.get_fchk()
########################################################################################################################################################################################################
os.system("rm -rf __pycache__")
########################################################################################################################################################################################################
