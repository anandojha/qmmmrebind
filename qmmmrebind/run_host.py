from parameterize_host import *
########################################################################################################################################################################################################
qm_system_host = PrepareGaussian(host_qm_pdb = "host_qm.pdb", n_processors = 12, memory = 50, charge = 0, multiplicity = 1, functional = "B3LYP", basis_set = "6-31G", optimisation = "OPT", frequency = "FREQ", add_keywords_I = "Integral=(Grid=UltraFine)", add_keywords_II = "Pop(MK,ReadRadii)", add_keywords_III = "IOp(6/33=2,6/42=6)", gauss_out_file = "host_qm.out", fchk_out_file = "host_qm_fchk.out")
########################################################################################################################################################################################################
params_host = Parameterize(xyz_file = "host_qm_coords.xyz", coordinate_file = "host_qm_coordinates.txt", unprocessed_hessian_file = "host_qm_unprocessed_hessian.txt", bond_list_file = "host_qm_bond_list.txt", angle_list_file = "host_qm_angle_list.txt", hessian_file = "host_qm_hessian.txt", atom_names_file = "host_qm_atom_names.txt", bond_parameter_file = "host_qm_bonds.txt", vibrational_scaling = 1.00, angle_parameter_file = "host_qm_angles.txt", charge_parameter_file = "host_qm_charges.txt", host_qm_pdb = "host_qm.pdb", host_qm_params_file = "host_qm_params.txt", host_pdb = "host.pdb", ffxml = "amber14-all.xml", host_xml = "host_off.xml", sim_output = "sim_output.pdb", sim_steps = 10000, host_residue_name = "HST", host_singular_file = "host_HST.pdb", conect_pdb_txt = "host_connect.txt", conect_pdb_file = "host_connect.pdb", host_singular_xml_file = "host_HST.xml", coulomb14scale = "0.833333", lj14scale = "0.5", reparameterised_host_xml_file = "host_reparameterized.xml", reparams_host_file = "host_final.xml")
########################################################################################################################################################################################################
#qm_system_host.write_input()
#qm_system_host.run_gaussian()
#qm_system_host.get_fchk()
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
