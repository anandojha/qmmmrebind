from operator import itemgetter
import pandas as pd
import numpy as np
import time
import math
import os
import re
########################################################################################################################################################################################################
element_list = [['1 ', 'H ', 'Hydrogen'], ['2 ', 'He', 'Helium'], ['3 ', 'Li', 'Lithium'], ['4 ', 'Be', 'Beryllium'], ['5 ', 'B ', 'Boron'], ['6 ', 'C ', 'Carbon'], ['7 ', 'N ', 'Nitrogen'], ['8 ', 'O ', 'Oxygen'], ['9 ', 'F ', 'Fluorine'], ['10', 'Ne', 'Neon'], ['11', 'Na', 'Sodium'], ['12', 'Mg', 'Magnesium'], ['13', 'Al', 'Aluminum'], ['14', 'Si', 'Silicon'], ['15', 'P ', 'Phosphorus'], ['16', 'S ', 'Sulfur'], ['17', 'Cl', 'Chlorine'], ['18', 'Ar', 'Argon'], ['19', 'K ', 'Potassium'], ['20', 'Ca', 'Calcium'], ['21', 'Sc', 'Scandium'], ['22', 'Ti', 'Titanium'], ['23', 'V ', 'Vanadium'], ['24', 'Cr', 'Chromium'], ['25', 'Mn', 'Manganese'], ['26', 'Fe', 'Iron'], ['27', 'Co', 'Cobalt'], ['28', 'Ni', 'Nickel'], ['29', 'Cu', 'Copper'], ['30', 'Zn', 'Zinc'], ['31', 'Ga', 'Gallium'], ['32', 'Ge', 'Germanium'], ['33', 'As', 'Arsenic'], ['34', 'Se', 'Selenium'], ['35', 'Br', 'Bromine'], ['36', 'Kr', 'Krypton'], ['37', 'Rb', 'Rubidium'], ['38', 'Sr', 'Strontium'], ['39', 'Y ', 'Yttrium'], ['40', 'Zr', 'Zirconium'], ['41', 'Nb', 'Niobium'], ['42', 'Mo', 'Molybdenum'], ['43', 'Tc', 'Technetium'], ['44', 'Ru', 'Ruthenium'], ['45', 'Rh', 'Rhodium'], ['46', 'Pd', 'Palladium'], ['47', 'Ag', 'Silver'], ['48', 'Cd', 'Cadmium'], ['49', 'In', 'Indium'], ['50', 'Sn', 'Tin'], ['51', 'Sb', 'Antimony'], ['52', 'Te', 'Tellurium'], ['53', 'I ', 'Iodine'], ['54', 'Xe', 'Xenon'], ['55', 'Cs', 'Cesium'], ['56', 'Ba', 'Barium'], ['57', 'La', 'Lanthanum'], ['58', 'Ce', 'Cerium'], ['59', 'Pr', 'Praseodymium'], ['60', 'Nd', 'Neodymium'], ['61', 'Pm', 'Promethium'], ['62', 'Sm', 'Samarium'], ['63', 'Eu', 'Europium'], ['64', 'Gd', 'Gadolinium'], ['65', 'Tb', 'Terbium'], ['66', 'Dy', 'Dysprosium'], ['67', 'Ho', 'Holmium'], ['68', 'Er', 'Erbium'], ['69', 'Tm', 'Thulium'], ['70', 'Yb', 'Ytterbium'], ['71', 'Lu', 'Lutetium'], ['72', 'Hf', 'Hafnium'], ['73', 'Ta', 'Tantalum'], ['74', 'W ', 'Tungsten'], ['75', 'Re', 'Rhenium'], ['76', 'Os', 'Osmium'], ['77', 'Ir', 'Iridium'], ['78', 'Pt', 'Platinum'], ['79', 'Au', 'Gold'], ['80', 'Hg', 'Mercury'], ['81', 'Tl', 'Thallium'], ['82', 'Pb', 'Lead'], ['83', 'Bi', 'Bismuth'], ['84', 'Po', 'Polonium'], ['85', 'At', 'Astatine'], ['86', 'Rn', 'Radon'], ['87', 'Fr', 'Francium'], ['88', 'Ra', 'Radium'], ['89', 'Ac', 'Actinium'], ['90', 'Th', 'Thorium'], ['91', 'Pa', 'Protactinium'], ['92', 'U ', 'Uranium'], ['93', 'Np', 'Neptunium'], ['94', 'Pu', 'Plutonium'], ['95', 'Am', 'Americium'], ['96', 'Cm', 'Curium'], ['97', 'Bk', 'Berkelium'], ['98', 'Cf', 'Californium'], ['99', 'Es', 'Einsteinium']]
########################################################################################################################################################################################################
def unit_vector_N(u_BC, u_AB):
    # Calculates unit normal vector which is perpendicular to plane ABC
    cross_product = np.cross(u_BC, u_AB)
    norm_u_N = np.linalg.norm(cross_product)
    u_N = cross_product / norm_u_N
    return u_N
########################################################################################################################################################################################################
def  u_PA_from_angles(atom_A, atom_B, atom_C, coords):
    # Gives the vector in the plane A,B,C and perpendicular to A to B
    diff_AB = coords[atom_B,:] - coords[atom_A,:]
    norm_diff_AB = np.linalg.norm(diff_AB)
    u_AB = diff_AB / norm_diff_AB
    diff_CB = coords[atom_B,:] - coords[atom_C,:]
    norm_diff_CB = np.linalg.norm(diff_CB)
    u_CB = diff_CB / norm_diff_CB
    u_N = unit_vector_N( u_CB, u_AB )
    u_PA = np.cross(u_N,  u_AB)
    norm_PA = np.linalg.norm(u_PA)
    u_PA = u_PA /norm_PA;
    return u_PA
########################################################################################################################################################################################################
def force_angle_constant( atom_A, atom_B, atom_C, bond_lengths, eigenvalues, eigenvectors, coords, scaling_1, scaling_2 ):
    # Force Constant- Equation 14 of seminario calculation paper - gives force constant for angle (in kcal/mol/rad^2) and equilibrium angle in degrees
    # Vectors along bonds calculated
    diff_AB = coords[atom_B,:] - coords[atom_A,:]
    norm_diff_AB = np.linalg.norm(diff_AB)
    u_AB = diff_AB / norm_diff_AB
    diff_CB = coords[atom_B,:] - coords[atom_C,:]
    norm_diff_CB = np.linalg.norm(diff_CB)
    u_CB = diff_CB / norm_diff_CB
    # Bond lengths and eigenvalues found
    bond_length_AB = bond_lengths[atom_A,atom_B]
    eigenvalues_AB = eigenvalues[atom_A, atom_B, :]
    eigenvectors_AB = eigenvectors[0:3, 0:3, atom_A, atom_B]
    bond_length_BC = bond_lengths[atom_B,atom_C]
    eigenvalues_CB = eigenvalues[atom_C, atom_B,  :]
    eigenvectors_CB = eigenvectors[0:3, 0:3, atom_C, atom_B]
    # Normal vector to angle plane found
    u_N = unit_vector_N( u_CB, u_AB )
    u_PA = np.cross(u_N,  u_AB)
    norm_u_PA = np.linalg.norm(u_PA)
    u_PA = u_PA / norm_u_PA
    u_PC = np.cross(u_CB, u_N)
    norm_u_PC = np.linalg.norm(u_PC)
    u_PC = u_PC / norm_u_PC
    sum_first = 0 
    sum_second = 0
    # Projections of eigenvalues
    for i in range(0,3):
        eig_AB_i = eigenvectors_AB[:,i]
        eig_BC_i = eigenvectors_CB[:,i]
        sum_first = sum_first + ( eigenvalues_AB[i] * abs(dot_product(u_PA , eig_AB_i) ) ) 
        sum_second = sum_second +  ( eigenvalues_CB[i] * abs(dot_product(u_PC, eig_BC_i) ) ) 
    # Scaling due to additional angles - Modified Seminario Part
    sum_first = sum_first/scaling_1 
    sum_second = sum_second/scaling_2
    # Added as two springs in series
    k_theta = ( 1 / ( (bond_length_AB**2) * sum_first) ) + ( 1 / ( (bond_length_BC**2) * sum_second) ) 
    k_theta = 1/k_theta
    k_theta = - k_theta #Change to OPLS form
    k_theta = abs(k_theta * 0.5) #Change to OPLS form
    # Equilibrium Angle
    theta_0 = math.degrees(math.acos(np.dot(u_AB, u_CB)))
    # If the vectors u_CB and u_AB are linearly dependent u_N cannot be defined
    # This case is dealt with here
    if abs(sum((u_CB) - (u_AB))) < 0.01 or ( abs(sum((u_CB) - (u_AB))) > 1.99 and abs(sum((u_CB) - (u_AB))) < 2.01):
        scaling_1 = 1 
        scaling_2 = 1
        [ k_theta, theta_0 ] = force_angle_constant_special_case( atom_A, atom_B, atom_C, bond_lengths, eigenvalues, eigenvectors, coords, scaling_1, scaling_2 )
    return k_theta, theta_0
########################################################################################################################################################################################################
def dot_product(u_PA , eig_AB):
    x = 0     
    for i in range(0,3):
        x = x + u_PA[i] * eig_AB[i].conjugate()
    return x 
########################################################################################################################################################################################################
def force_angle_constant_special_case( atom_A, atom_B, atom_C, bond_lengths, eigenvalues, eigenvectors, coords, scaling_1, scaling_2 ):
    # Force Constant- Equation 14 of seminario calculation paper - gives force constant for angle (in kcal/mol/rad^2) and equilibrium angle in degrees
    # Deals with cases when u_N cannot be defined and instead takes samples of u_N across a unit sphere. 
    # Vectors along bonds calculated
    diff_AB = coords[atom_B,:] - coords[atom_A,:]
    norm_diff_AB = np.linalg.norm(diff_AB)
    u_AB = diff_AB / norm_diff_AB
    diff_CB = coords[atom_B,:] - coords[atom_C,:]
    norm_diff_CB = np.linalg.norm(diff_CB)
    u_CB = diff_CB / norm_diff_CB
    # Bond lengths and eigenvalues found
    bond_length_AB = bond_lengths[atom_A,atom_B]
    eigenvalues_AB = eigenvalues[atom_A, atom_B, :]
    eigenvectors_AB = eigenvectors[0:3, 0:3, atom_A, atom_B]
    bond_length_BC = bond_lengths[atom_B,atom_C]
    eigenvalues_CB = eigenvalues[atom_C, atom_B,  :]
    eigenvectors_CB = eigenvectors[0:3, 0:3, atom_C, atom_B]
    k_theta_array = np.zeros((180, 360))
    # Find force constant with varying u_N (with vector uniformly sampled across a sphere)
    for theta in range(0,180):
        for phi in range(0,360):   
            r = 1
            u_N = [r * math.sin(math.radians(theta))*math.cos(math.radians(theta)),r * math.sin(math.radians(theta))*math.sin(math.radians(theta)), r*math.cos(math.radians(theta)) ]
            u_PA = np.cross(u_N,  u_AB)
            u_PA = u_PA / np.linalg.norm(u_PA)
            u_PC = np.cross(u_CB, u_N)
            u_PC = u_PC / np.linalg.norm(u_PC)
            sum_first = 0
            sum_second = 0
            # Projections of eigenvalues
            for i in range(0,3):
                eig_AB_i = eigenvectors_AB[:,i]
                eig_BC_i = eigenvectors_CB[:,i]
                sum_first = sum_first + ( eigenvalues_AB[i] * abs(dot_product(u_PA , eig_AB_i) ) ) 
                sum_second = sum_second +  ( eigenvalues_CB[i] * abs(dot_product(u_PC, eig_BC_i) ) ) 
            # Added as two springs in series
            k_theta_ij = ( 1 / ( (bond_length_AB**2) * sum_first) ) + ( 1 / ( (bond_length_BC**2) * sum_second) ) 
            k_theta_ij = 1/k_theta_ij
            k_theta_ij = - k_theta_ij #Change to OPLS form
            k_theta_ij = abs(k_theta_ij * 0.5) #Change to OPLS form
            k_theta_array[theta, phi] = k_theta_ij
    # Removes cases where u_N was linearly dependent of u_CB or u_AB
    # Force constant used is taken as the mean
    k_theta = np.mean(np.mean(k_theta_array))
    # Equilibrium Angle independent of u_N
    theta_0 = math.degrees(math.cos(np.dot(u_AB, u_CB)))
    return k_theta, theta_0
########################################################################################################################################################################################################
def force_constant_bond(atom_A, atom_B, eigenvalues, eigenvectors, coords):
    # Force Constant - Equation 10 of Seminario paper - gives force constant for bond
    # Eigenvalues and eigenvectors calculated 
    eigenvalues_AB = eigenvalues[atom_A,atom_B,:]
    eigenvectors_AB = eigenvectors[:,:,atom_A,atom_B]
    # Vector along bond 
    diff_AB = np.array(coords[atom_B,:]) - np.array(coords[atom_A,:])
    norm_diff_AB = np.linalg.norm(diff_AB)
    unit_vectors_AB = diff_AB / norm_diff_AB
    k_AB = 0 
    # Projections of eigenvalues 
    for i in range(0,3):
        dot_product = abs(np.dot(unit_vectors_AB, eigenvectors_AB[:,i]))
        k_AB = k_AB + ( eigenvalues_AB[i] * dot_product )
    k_AB = -k_AB * 0.5 # Convert to OPLS form
    return k_AB
########################################################################################################################################################################################################
def  u_PA_from_angles(atom_A, atom_B, atom_C, coords):
    # Gives the vector in the plane A,B,C and perpendicular to A to B
    diff_AB = coords[atom_B,:] - coords[atom_A,:]
    norm_diff_AB = np.linalg.norm(diff_AB)
    u_AB = diff_AB / norm_diff_AB
    diff_CB = coords[atom_B,:] - coords[atom_C,:]
    norm_diff_CB = np.linalg.norm(diff_CB)
    u_CB = diff_CB / norm_diff_CB
    u_N = unit_vector_N( u_CB, u_AB )
    u_PA = np.cross(u_N,  u_AB)
    norm_PA = np.linalg.norm(u_PA)
    u_PA = u_PA /norm_PA;
    return u_PA
########################################################################################################################################################################################################
def reverse_list(lst): 
    reversed_list = lst[::-1] 
    return (reversed_list)
########################################################################################################################################################################################################
def uniq(input_):
    output = []
    for x in input_:
        if x not in output:
            output.append(x)
    return output
########################################################################################################################################################################################################
def search_in_file(file: str, word: str) -> list:
    """Search for the given string in file and return lines containing that string along with line numbers"""
    line_number = 0
    list_of_results = []
    with open(file, 'r') as f:
        for line in f:
            line_number += 1
            if word in line:
                list_of_results.append((line_number, line.rstrip()))
    return(list_of_results)
########################################################################################################################################################################################################
def OPLS_LJ(system):
    forces = {system.getForce(index).__class__.__name__: system.getForce(index) for index in range(system.getNumForces())}
    nonbonded_force = forces['NonbondedForce']
    lorentz = simtk.openmm.openmm.CustomNonbondedForce('4*epsilon*((sigma/r)^12-(sigma/r)^6); sigma=sqrt(sigma1*sigma2); epsilon=sqrt(epsilon1*epsilon2)')
    lorentz.setNonbondedMethod(nonbonded_force.getNonbondedMethod())
    lorentz.addPerParticleParameter('sigma')
    lorentz.addPerParticleParameter('epsilon')
    lorentz.setCutoffDistance(nonbonded_force.getCutoffDistance())
    system.addForce(lorentz)
    LJset = {}
    for index in range(nonbonded_force.getNumParticles()):
        charge, sigma, epsilon = nonbonded_force.getParticleParameters(index)
        LJset[index] = (sigma, epsilon)
        lorentz.addParticle([sigma, epsilon])
        nonbonded_force.setParticleParameters(index, charge, sigma, epsilon * 0)
    for i in range(nonbonded_force.getNumExceptions()):
        (p1, p2, q, sig, eps) = nonbonded_force.getExceptionParameters(i)
        lorentz.addExclusion(p1, p2)
        if eps._value != 0.0:
            #print (p1,p2,sig,eps)
            sig14 = math.sqrt(LJset[p1][0].value_in_unit(simtk.unit.nanometer) * LJset[p2][0].value_in_unit(simtk.unit.nanometer))
            eps14 = math.sqrt(LJset[p1][1].value_in_unit(simtk.unit.kilojoule/simtk.unit.mole) * LJset[p2][1].value_in_unit(simtk.unit.kilojoule/simtk.unit.mole))
            nonbonded_force.setExceptionParameters(i, p1, p2, q, sig14, eps)
    return (system)
########################################################################################################################################################################################################
def list_to_dict(lst):
    res_dct = {lst[i]: lst[i + 1] for i in range(0, len(lst), 2)}
    return (res_dct)
########################################################################################################################################################################################################
