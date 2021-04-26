from biopandas.pdb import PandasPdb
import matplotlib.pyplot as plt
from simtk.openmm.app import *
from scipy import optimize
#from simtk.openmm import *
#from simtk.unit import *
import pandas as pd
import numpy as np
import simtk
import scipy
import math
import re
import os
####################################################################################################################################################################################
def get_dihedrals(qm_scan_file):
    with open(qm_scan_file, "r") as f:
        lines = f.readlines()
    energy_dihedral_lines = []
    for i in range(len(lines)):
        if "Dihedral" in lines[i]:
            energy_dihedral_lines.append(lines[i])
    dihedrals = []
    for i in energy_dihedral_lines:
        energy_dihedral = i
        energy_dihedral = re.findall(r'[-+]?\d+[.]?\d*', energy_dihedral)
        dihedral = float(energy_dihedral[0])
        dihedrals.append(dihedral)
    return(dihedrals)
####################################################################################################################################################################################
def get_qm_energies(qm_scan_file):
    with open(qm_scan_file, "r") as f:
        lines = f.readlines()
    energy_dihedral_lines = []
    for i in range(len(lines)):
        if "Dihedral" in lines[i]:
            energy_dihedral_lines.append(lines[i])
    qm_energies = []
    for i in energy_dihedral_lines:
        energy_dihedral = i 
        energy_dihedral = re.findall(r'[-+]?\d+[.]?\d*', energy_dihedral)
        energy = float(energy_dihedral[1])
        qm_energies.append(energy)
    return(qm_energies)
####################################################################################################################################################################################
def generate_mm_pdbs(qm_scan_file, pdb_file):
    with open (qm_scan_file, "r") as f:
        lines = f.readlines()
    energy_dihedral_lines = []
    for i in range(len(lines)):
        if "Dihedral" in lines[i]:
            energy_dihedral_lines.append(lines[i])
    dihedrals = []
    for i in energy_dihedral_lines:
        energy_dihedral = i 
        energy_dihedral = re.findall(r'[-+]?\d+[.]?\d*', energy_dihedral)
        dihedral = float(energy_dihedral[0])
        dihedrals.append(dihedral)
    lines_markers = []
    for i in range(len(lines)):
        if "Dihedral" in lines[i]:
            lines_markers.append(i)
    lines_markers.append(len(lines) + 1)
    for i in range(len(lines_markers) - 1):
        #pdb_file_to_write = str(dihedrals[i]) + ".pdb"
        if dihedrals[i] > 0 :
            pdb_file_to_write = "plus_" + str(abs(dihedrals[i])) + ".pdb"
        if dihedrals[i] < 0 :
            pdb_file_to_write = "minus_" + str(abs(dihedrals[i])) + ".pdb"        
        to_begin = lines_markers[i]
        to_end = lines_markers[i+1]
        lines_to_write = lines[to_begin + 1:to_end - 1]
        x_coords = []
        y_coords = []
        z_coords = []
        for i in lines_to_write:
            coordinates = i 
            coordinates = re.findall(r'[-+]?\d+[.]?\d*', coordinates)
            x = float(coordinates[0])
            y = float(coordinates[1])
            z = float(coordinates[2])
            x_coords.append(x)
            y_coords.append(y)
            z_coords.append(z)
        ppdb = PandasPdb()
        ppdb.read_pdb(pdb_file)
        ppdb.df["ATOM"]["x_coord"] = x_coords
        ppdb.df["ATOM"]["y_coord"] = y_coords
        ppdb.df["ATOM"]["z_coord"] = z_coords    
        ppdb.to_pdb(pdb_file_to_write)  
####################################################################################################################################################################################       
def get_mm_energy(pdb_file, forcefield_file):
    pdb = simtk.openmm.app.PDBFile(pdb_file)
    forcefield = simtk.openmm.app.ForceField(forcefield_file)
    system = forcefield.createSystem(pdb.topology)
    integrator = simtk.openmm.LangevinIntegrator(300 * simtk.unit.kelvin, 1 / simtk.unit.picosecond, 0.002 * simtk.unit.picoseconds)
    simulation = simtk.openmm.app.Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    #simulation.minimizeEnergy()
    state = simulation.context.getState(getEnergy = True, getParameters = True, getForces = True)
    potential_energy = state.getPotentialEnergy()
    potential_energy = str(potential_energy)
    potential_energy = re.findall(r'[-+]?\d+[.]?\d*', potential_energy)
    return(float(potential_energy[0])) 
####################################################################################################################################################################################
def create_xml_no_bond(forcefield_file, exclude_bond_forces_file):
    with open (forcefield_file, "r") as f:
        lines = f.readlines()
    for i in range(len(lines)):
        if "<HarmonicBondForce>" in lines[i]:
            to_begin = int(i) - 1
        if "</HarmonicBondForce>" in lines[i]:
            to_end = int(i) + 1 
    bond_forces = lines[to_begin + 1:to_end]  
    lines_to_write_1 = lines[ : to_begin +1 ]
    lines_to_write_2 = lines[to_end : ]
    with open (exclude_bond_forces_file, "w") as f:
        f.writelines(lines_to_write_1)
        f.writelines(lines_to_write_2)   
####################################################################################################################################################################################
def create_xml_no_angle(forcefield_file, exclude_angle_forces_file):
    with open (forcefield_file, "r") as f:
        lines = f.readlines()
    for i in range(len(lines)):
        if "<HarmonicAngleForce>" in lines[i]:
            to_begin = int(i) - 1
        if "</HarmonicAngleForce>" in lines[i]:
            to_end = int(i) + 1 
    bond_forces = lines[to_begin + 1:to_end]  
    lines_to_write_1 = lines[ : to_begin +1 ]
    lines_to_write_2 = lines[to_end : ]
    with open (exclude_angle_forces_file, "w") as f:
        f.writelines(lines_to_write_1)
        f.writelines(lines_to_write_2)
####################################################################################################################################################################################       
def create_xml_no_torsion(forcefield_file, exclude_torsion_forces_file):
    with open (forcefield_file, "r") as f:
        lines = f.readlines()
    for i in range(len(lines)):
        if "<PeriodicTorsionForce>" in lines[i]:
            to_begin = int(i) - 1
        if "</PeriodicTorsionForce>" in lines[i]:
            to_end = int(i) + 1 
    bond_forces = lines[to_begin + 1:to_end]  
    lines_to_write_1 = lines[ : to_begin +1 ]
    lines_to_write_2 = lines[to_end : ]
    with open (exclude_torsion_forces_file, "w") as f:
        f.writelines(lines_to_write_1)
        f.writelines(lines_to_write_2)
####################################################################################################################################################################################
def create_xml_no_non_bonded(forcefield_file, exclude_non_bonded_forces_file):
    with open (forcefield_file, "r") as f:
        lines = f.readlines()
    for i in range(len(lines)):
        if "<NonbondedForce" in lines[i]:
            to_begin = int(i) - 1
        if "</NonbondedForce" in lines[i]:
            to_end = int(i) + 1 
    bond_forces = lines[to_begin + 1:to_end]  
    lines_to_write_1 = lines[ : to_begin +1 ]
    lines_to_write_2 = lines[to_end : ]
    with open (exclude_non_bonded_forces_file, "w") as f:
        f.writelines(lines_to_write_1)
        f.writelines(lines_to_write_2)
####################################################################################################################################################################################
def scale_list(list_):
    scaled_list = [i - min(list_) for i in list_]
    return (scaled_list)
####################################################################################################################################################################################
def list_hartree_kcal(list_):
    converted_list = [i * 627.5094 for i in list_]
    return (converted_list)
####################################################################################################################################################################################
def list_kJ_kcal(list_):
    converted_list = [i / 4.184 for i in list_]
    return (converted_list)
####################################################################################################################################################################################
def get_mm_bond_energies(qm_scan_file, forcefield_file, exclude_bond_forces_file):
    mm_pdb_list = []
    for i in get_dihedrals(qm_scan_file):
        if i > 0 :
            pdb_file = "plus_" + str(abs(i)) + ".pdb"
        if i < 0 :
            pdb_file = "minus_" + str(abs(i)) + ".pdb"  
        mm_pdb_list.append(pdb_file)
    for i in mm_pdb_list:
        mm_pdb_file = i
    mm_bond_energies = []
    for i in mm_pdb_list:
        mm_pdb_file = i
        mm_energy = get_mm_energy(pdb_file = mm_pdb_file, forcefield_file = forcefield_file)
        create_xml_no_bond(forcefield_file, exclude_bond_forces_file)
        angle_torsion_nonbonded_energy = get_mm_energy(pdb_file = mm_pdb_file, forcefield_file = exclude_bond_forces_file)
        bond_energy = mm_energy - angle_torsion_nonbonded_energy
        mm_bond_energies.append(bond_energy)
        remove_command = "rm -rf " + exclude_bond_forces_file
        os.system(remove_command)   
    return(mm_bond_energies)
####################################################################################################################################################################################
def get_mm_angle_energies(qm_scan_file, forcefield_file, exclude_angle_forces_file):
    mm_pdb_list = []
    for i in get_dihedrals(qm_scan_file):
        if i > 0 :
            pdb_file = "plus_" + str(abs(i)) + ".pdb"
        if i < 0 :
            pdb_file = "minus_" + str(abs(i)) + ".pdb"  
        mm_pdb_list.append(pdb_file)
    for i in mm_pdb_list:
        mm_pdb_file = i
    mm_angle_energies = []
    for i in mm_pdb_list:
        mm_pdb_file = i
        mm_energy = get_mm_energy(pdb_file = mm_pdb_file, forcefield_file = forcefield_file)
        create_xml_no_angle(forcefield_file, exclude_angle_forces_file)
        bond_torsion_nonbonded_energy = get_mm_energy(pdb_file = mm_pdb_file, forcefield_file = exclude_angle_forces_file)
        angle_energy = mm_energy - bond_torsion_nonbonded_energy
        mm_angle_energies.append(angle_energy)
        remove_command = "rm -rf " + exclude_angle_forces_file
        os.system(remove_command)   
    return(mm_angle_energies)
####################################################################################################################################################################################
def get_mm_torsion_energies(qm_scan_file, forcefield_file, exclude_torsion_forces_file):
    mm_pdb_list = []
    for i in get_dihedrals(qm_scan_file):
        if i > 0 :
            pdb_file = "plus_" + str(abs(i)) + ".pdb"
        if i < 0 :
            pdb_file = "minus_" + str(abs(i)) + ".pdb"  
        mm_pdb_list.append(pdb_file)
    for i in mm_pdb_list:
        mm_pdb_file = i
    mm_torsion_energies = []
    for i in mm_pdb_list:
        mm_pdb_file = i
        mm_energy = get_mm_energy(pdb_file = mm_pdb_file, forcefield_file = forcefield_file)
        create_xml_no_torsion(forcefield_file, exclude_torsion_forces_file)
        bond_angle_nonbonded_energy = get_mm_energy(pdb_file = mm_pdb_file, forcefield_file = exclude_torsion_forces_file)
        torsion_energy = mm_energy - bond_angle_nonbonded_energy
        mm_torsion_energies.append(torsion_energy)
        remove_command = "rm -rf " + exclude_torsion_forces_file
        os.system(remove_command)   
    return(mm_torsion_energies)
####################################################################################################################################################################################
def get_mm_non_bonded_energies(qm_scan_file, forcefield_file, exclude_non_bonded_forces_file):
    mm_pdb_list = []
    for i in get_dihedrals(qm_scan_file):
        if i > 0 :
            pdb_file = "plus_" + str(abs(i)) + ".pdb"
        if i < 0 :
            pdb_file = "minus_" + str(abs(i)) + ".pdb"  
        mm_pdb_list.append(pdb_file)
    for i in mm_pdb_list:
        mm_pdb_file = i
    mm_non_bonded_energies = []
    for i in mm_pdb_list:
        mm_pdb_file = i
        mm_energy = get_mm_energy(pdb_file = mm_pdb_file, forcefield_file = forcefield_file)
        create_xml_no_non_bonded(forcefield_file, exclude_non_bonded_forces_file)
        bond_angle_torsion_energy = get_mm_energy(pdb_file = mm_pdb_file, forcefield_file = exclude_non_bonded_forces_file)
        nonbonded_energy = mm_energy - bond_angle_torsion_energy
        mm_non_bonded_energies.append(nonbonded_energy)
        remove_command = "rm -rf " + exclude_non_bonded_forces_file
        os.system(remove_command)   
    return(mm_non_bonded_energies)
####################################################################################################################################################################################
def get_mm_potential_energies(qm_scan_file, forcefield_file):
    mm_pdb_list = []
    for i in get_dihedrals(qm_scan_file):
        if i > 0 :
            pdb_file = "plus_" + str(abs(i)) + ".pdb"
        if i < 0 :
            pdb_file = "minus_" + str(abs(i)) + ".pdb"  
        mm_pdb_list.append(pdb_file)
    for i in mm_pdb_list:
        mm_pdb_file = i
    mm_potential_energies = []
    for i in mm_pdb_list:
        mm_pdb_file = i
        mm_energy = get_mm_energy(pdb_file = mm_pdb_file, forcefield_file = forcefield_file)
        mm_potential_energies.append(mm_energy)
    return(mm_potential_energies)
####################################################################################################################################################################################
def remove_mm_pdbs(qm_scan_file):
    mm_pdb_list = []
    for i in get_dihedrals(qm_scan_file):
        if i > 0 :
            pdb_file = "plus_" + str(abs(i)) + ".pdb"
        if i < 0 :
            pdb_file = "minus_" + str(abs(i)) + ".pdb"  
        mm_pdb_list.append(pdb_file)
    for i in mm_pdb_list:
        command = "rm -rf  " + i
        os.system(command) 
####################################################################################################################################################################################
def list_diff(list_1, list_2):
    diff_list = []
    zipped_list = zip(list_1, list_2)
    for list1_i, list2_i in zipped_list:
        diff_list.append(list1_i-list2_i)
    return(diff_list)
####################################################################################################################################################################################
def dihedral_energy (x, k1, k2, k3, k4 = 0):
    energy_1 = k1*(1 + np.cos (1*x*0.01745))
    energy_2 = k2*(1 - np.cos (2*x*0.01745))
    energy_3 = k3*(1 + np.cos (3*x*0.01745))
    energy_4 = k4*(1 - np.cos (4*x*0.01745))
    dihedral_energy = energy_1 + energy_2 + energy_3 + energy_4
    return (dihedral_energy)
####################################################################################################################################################################################
def error_function(delta_qm, delta_mm):
    squared_error = np.square(np.subtract(delta_qm, delta_mm))
    mean_squared_error = squared_error.mean()
    root_mean_squared_error = math.sqrt(mean_squared_error)
    return (root_mean_squared_error)
####################################################################################################################################################################################
def error_function_boltzmann(delta_qm, delta_mm, T):
    kb = 3.297623483* 10**(-24) #in cal/K
    delta_qm_boltzmann_weighted = [np.exp(-i / (kb * T)) for i in delta_qm]
    squared_error = np.square(np.subtract(delta_qm, delta_mm)) * delta_qm_boltzmann_weighted
    mean_squared_error = squared_error.mean()
    root_mean_squared_error = math.sqrt(mean_squared_error)
    return (root_mean_squared_error)
####################################################################################################################################################################################
def gen_init_guess(qm_scan_file, forcefield_file):
    x = get_dihedrals(qm_scan_file)
    y  = scale_list(list_kJ_kcal(list_ = get_mm_potential_energies(qm_scan_file,forcefield_file)))  
    init_vals = [0.0, 0.0, 0.0, 0.0]
    k_init_guess, covar = scipy.optimize.curve_fit(dihedral_energy, x, y, p0 = init_vals)
    for i in range(len(k_init_guess)):
        if k_init_guess[i] < 0:
            k_init_guess[i] = 0
    return(k_init_guess)  
####################################################################################################################################################################################
def objective_function(k_array, x, delta_qm):
    delta_mm = dihedral_energy(x, k1 = k_array[0], k2 = k_array[1], k3 = k_array[2],  k4 = k_array[3])    
    loss_function = error_function(delta_qm, delta_mm)
    return(loss_function)
####################################################################################################################################################################################
def fit_params(qm_scan_file, forcefield_file, method):
    k_guess = gen_init_guess(qm_scan_file, forcefield_file)
    x_data = np.array(get_dihedrals(qm_scan_file))
    delta_qm = np.array(scale_list(list_hartree_kcal(list_ = get_qm_energies(qm_scan_file)))) 
    optimise = scipy.optimize.minimize(objective_function, k_guess, args=(x_data,delta_qm), method = method, bounds = [(0.00, None), (0.00, None), (0.00, None), (0.00, None)])
    return(optimise.x)
####################################################################################################################################################################################
def get_tor_params(qm_scan_file, pdb_file, forcefield_file, exclude_torsion_forces_file, method):
    remove_mm_pdbs(qm_scan_file)
    qm_e = get_qm_energies(qm_scan_file)
    qm_e_kcal = list_hartree_kcal(list_ = qm_e)
    delta_qm = scale_list(qm_e_kcal)
    generate_mm_pdbs(qm_scan_file, pdb_file)
    mm_torsion_e = get_mm_torsion_energies(qm_scan_file, forcefield_file, exclude_torsion_forces_file)
    mm_pe = get_mm_potential_energies(qm_scan_file, forcefield_file)
    mm_pe_no_torsion = list_diff(mm_pe, mm_torsion_e)
    mm_pe_no_torsion_kcal = list_kJ_kcal(list_ = mm_pe_no_torsion)
    delta_mm = scale_list(mm_pe_no_torsion_kcal)  
    opt_param = fit_params(qm_scan_file, forcefield_file, method)
    remove_mm_pdbs(qm_scan_file)
    return(opt_param)
####################################################################################################################################################################################
def get_tor_param_lines(qm_scan_file, pdb_file, forcefield_file, exclude_torsion_forces_file, method, 
                        dihedral_text_file, atom_names_file):
    opt_param = get_tor_params(qm_scan_file, pdb_file,forcefield_file, exclude_torsion_forces_file, method)
    df_atoms = pd.read_csv( atom_names_file, header = None, delimiter = r"\s+")
    df_atoms['atom_numbers'] = df_atoms.index
    df_atoms['atom_numbers'] = df_atoms['atom_numbers'] + 1
    df_atoms.columns = ["atom", "atom_numbers"]
    df_atoms.head()
    class_list = []
    for i in range(len(opt_param)):
        class_list.append("class" + str(i+1) + "=" + '"' + df_atoms['atom'].iloc[int(np.loadtxt(dihedral_text_file)[i]) - 1] + '"')
    tor_param_lines = []
    for i in range(len(opt_param)):
        tor_param_lines.append("<Proper" + " " + "k" + "=" + '"' + str(round(opt_param[i], 8)) +  '"' + " " + class_list[0] + " " + class_list[1] + " " + class_list[2] + " "  + class_list[3] + " " +" periodicity" + "=" + '"' + str(i+1) + '"' + " " + "phase" + "=" + "0.0" + '"' + "/>")
    return(tor_param_lines)
####################################################################################################################################################################################
#print(get_tor_params(qm_scan_file = "scan.xyz", pdb_file = "guest.pdb", forcefield_file = "guest.xml", exclude_torsion_forces_file = "guest_no_torsion_forces.xml", method = "L-BFGS-B"))
#print(get_tor_param_lines(qm_scan_file = "scan.xyz", pdb_file = "guest.pdb", forcefield_file = "guest.xml", exclude_torsion_forces_file = "guest_no_torsion_forces.xml", method = "L-BFGS-B", dihedral_text_file = "dihedrals.txt", atom_names_file = "guest_atom_names.txt"))
