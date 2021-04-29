from biopandas.pdb import PandasPdb
from operator import itemgetter
from mendeleev import element
from simtk.openmm import app
import subprocess as sp
from functions import * 
from sys import stdout
import pandas as pd
import numpy as np
import statistics
import itertools
import parmed
import pickle
import shutil
import simtk
import time
import math
import sys
import re
import os
####################################################################################################################################################################################
sys.setrecursionlimit(10000)
####################################################################################################################################################################################
class PrepareQMMM:
 
    def __init__(self, init_pdb, cleaned_pdb, guest_init_pdb, host_pdb, guest_resname, guest_pdb, guest_xyz, distance, residue_list, host_qm_atoms, host_mm_atoms, host_qm_pdb, host_mm_pdb, qm_pdb, mm_pdb, host_mm_region_I_atoms, host_mm_region_II_atoms, host_mm_region_I_pdb, host_mm_region_II_pdb, num_residues):
        self.init_pdb = init_pdb
        self.cleaned_pdb = cleaned_pdb
        self.guest_init_pdb = guest_init_pdb
        self.host_pdb = host_pdb
        self.guest_resname = guest_resname
        self.guest_pdb = guest_pdb
        self.guest_xyz  = guest_xyz
        self.distance = distance
        self.residue_list = residue_list
        self.host_qm_atoms = host_qm_atoms
        self.host_mm_atoms = host_mm_atoms
        self.host_qm_pdb = host_qm_pdb
        self.host_mm_pdb = host_mm_pdb
        self.qm_pdb = qm_pdb
        self.mm_pdb = mm_pdb
        self.host_mm_region_I_atoms = host_mm_region_I_atoms
        self.host_mm_region_II_atoms = host_mm_region_II_atoms
        self.host_mm_region_I_pdb = host_mm_region_I_pdb
        self.host_mm_region_II_pdb = host_mm_region_II_pdb
        self.num_residues = num_residues

    def clean_up(self):   
        """
        This function removes everything from the initial pdb file except the host and guest molecules from the system and saves a pdb file.
        """
        ions = ["Na+", "Cs+", "K+","Li+", "Rb+", "Cl-", "Br-", "F-", "I-", "Ca2"]
        intermediate_file_1 = self.cleaned_pdb[:-4] + "_intermediate_1.pdb"
        intermediate_file_2 = self.cleaned_pdb[:-4] + "_intermediate_2.pdb"
        command = "pdb4amber -i " + self.init_pdb + " -o " + intermediate_file_1 + " --noter --dry"
        os.system(command)
        to_delete = intermediate_file_1[:-4] + "_nonprot.pdb" , intermediate_file_1[:-4] + "_renum.txt", intermediate_file_1[:-4] + "_sslink", intermediate_file_1[:-4] + "_water.pdb"   
        os.system("rm -rf " + to_delete[0] + " " + to_delete[1] + " " + to_delete[2] + " " + to_delete[3])  
        with open(intermediate_file_1) as f1, open(intermediate_file_2, 'w') as f2:
            for line in f1:
                if not any(ion in line for ion in ions):
                    f2.write(line)                
        with open(intermediate_file_2, 'r') as f1 :
            filedata = f1.read()
        filedata = filedata.replace('HETATM', 'ATOM  ')
        with open(self.cleaned_pdb, 'w') as f2:
            f2.write(filedata)
        command = "rm -rf " + intermediate_file_1 + " " + intermediate_file_2
        os.system(command)
    
    def create_host_guest(self):   
        """
        This function creates separate host and guest pdb files.
        """
        with open(self.cleaned_pdb) as f1, open(self.host_pdb, 'w') as f2:
            for line in f1:
                if not self.guest_resname in line and not "CRYST1" in line:
                    f2.write(line)           
        with open(self.cleaned_pdb) as f1, open(self.guest_init_pdb, 'w') as f2:
            for line in f1:
                if self.guest_resname in line or "END" in line:
                    f2.write(line)
                    
    def realign_guest(self):   
        """
        This function realigns the atom numbers in the initial guest pdb file given the assumption that the original pdb file had the host molecule followed by the guest molecule.
        """
        ppdb = PandasPdb()
        ppdb.read_pdb(self.guest_init_pdb)
        to_substract = min(ppdb.df["ATOM"]["atom_number"]) - 1 
        ppdb.df["ATOM"]["atom_number"] = ppdb.df["ATOM"]["atom_number"] - to_substract
        intermediate_file_1 = self.guest_pdb[:-4] + "_intermediate_1.pdb"
        intermediate_file_2 = self.guest_pdb[:-4] + "_intermediate_2.pdb"
        ppdb.to_pdb(path = intermediate_file_1)
        command = "pdb4amber -i " + intermediate_file_1 +  " -o " + intermediate_file_2
        os.system(command)
        to_delete = intermediate_file_2[:-4] + "_nonprot.pdb" , intermediate_file_2[:-4] + "_renum.txt", intermediate_file_2[:-4] + "_sslink"
        os.system("rm -rf " + to_delete[0] + " " + to_delete[1] + " " + to_delete[2])  
        with open(intermediate_file_2, 'r') as f1 :
            filedata = f1.read()
        filedata = filedata.replace('HETATM', 'ATOM  ')
        with open(self.guest_pdb, 'w') as f2:
            f2.write(filedata)
        command = "rm -rf " + intermediate_file_1 + " " + intermediate_file_2
        os.system(command)
        
    def get_guest_coord(self):   
        """
        This function saves a list of xyz coordinates from the pdb file.
        """
        ppdb = PandasPdb()
        ppdb.read_pdb(self.guest_pdb)
        xyz = ppdb.df['ATOM'][["x_coord","y_coord","z_coord"]]
        xyz_to_list = xyz.values.tolist()
        np.savetxt(self.guest_xyz, xyz_to_list)
        
    def get_qm_resids(self):  
        """
        This function saves a list of residue numbers of the residues in the host molecule surrounding the guest within a mentioned distance from any of the atom in the guest molecule.
        """
        guest_coord_list = np.loadtxt(self.guest_xyz)
        host_atom_list = []
        for i in range(len(guest_coord_list)):
            reference_point = guest_coord_list[i]
            ppdb = PandasPdb()
            ppdb.read_pdb(self.host_pdb)
            distances = ppdb.distance(xyz=reference_point, records=('ATOM'))
            all_within_distance = ppdb.df['ATOM'][distances < float(self.distance)]
            host_df = all_within_distance["atom_number"]
            host_list = host_df.values.tolist()
            host_atom_list.append(host_list)
        host_atom_list = list(itertools.chain(*host_atom_list))
        host_atom_list = set(host_atom_list)
        host_atom_list = list(host_atom_list)
        host_atom_list.sort() 
        ppdb = PandasPdb()
        ppdb.read_pdb(self.host_pdb)
        df = ppdb.df['ATOM'][['atom_number','residue_number','residue_name']]
        index_list = []
        for i in host_atom_list:
            indices = np.where(df["atom_number"] == i)
            indices = list(indices)[0]
            indices = list(indices)
            index_list.append(indices)
        index_list = list(itertools.chain.from_iterable(index_list))
        df1 = df.iloc[index_list,] 
        resid_num = list(df1.residue_number.unique())
        np.savetxt(self.residue_list, resid_num)
        
    def get_host_qm_mm_atoms(self):  
        """
        This function saves a list of atom numbers of the residues in the host molecule surrounding the guest within a certain mentioned distance from any of the atom in the guest molecule.
        """
        resid_num = np.loadtxt(self.residue_list)
        #approximated_res_list = [int(i) for i in resid_num]
        approximated_res_list = []
        for i in range(int(statistics.median(resid_num)) - int(int(self.num_residues)/2), int(statistics.median(resid_num)) + int(int(self.num_residues)/2)):
            approximated_res_list.append(i)
        ppdb = PandasPdb()
        ppdb.read_pdb(self.host_pdb)
        df = ppdb.df['ATOM'][['atom_number','residue_number','residue_name']]
        host_index_nested_list = []
        for i in approximated_res_list:
            indices = np.where(df["residue_number"] == i)
            indices = list(indices)[0]
            indices = list(indices)
            host_index_nested_list.append(indices)  
        host_index_list = list(itertools.chain.from_iterable(host_index_nested_list))
        df_atom = df.iloc[host_index_list] 
        df_atom_number = df_atom['atom_number']
        host_atom_list = df_atom_number.values.tolist()
        selected_atoms = []
        selected_atoms.extend(host_atom_list)
        ppdb = PandasPdb()
        ppdb.read_pdb(self.host_pdb)
        len_atoms = []
        for i in range(len(ppdb.df['ATOM'])):
            len_atoms.append(i+1)
        non_selected_atoms = list(set(len_atoms).difference(selected_atoms))
        if len(non_selected_atoms) + len(selected_atoms) == len(len_atoms):
            print("Sum of the selected and non-selected region equals the length of list of total atoms")
        else:
            print("Error")
        np.savetxt(self.host_qm_atoms, selected_atoms)  
        np.savetxt(self.host_mm_atoms, non_selected_atoms) 
        
    def save_host_pdbs(self):  
        """
        This function saves a host pdb file of selected region and another host pdb file non-selected region. 
        """
        selected_atoms = np.loadtxt(self.host_qm_atoms)
        selected_atoms = [int(i) for i in selected_atoms]
        ppdb = PandasPdb()
        ppdb.read_pdb(self.host_pdb)
        for i in selected_atoms:
            ppdb.df['ATOM'] = ppdb.df['ATOM'][ppdb.df['ATOM']['atom_number'] != i]
        ppdb.to_pdb(path = self.host_mm_pdb, records = None, gz = False, append_newline = True)
        non_selected_atoms = np.loadtxt(self.host_mm_atoms)
        non_selected_atoms = [int(i) for i in non_selected_atoms]
        ppdb = PandasPdb()
        ppdb.read_pdb(self.host_pdb)
        for i in non_selected_atoms:
            ppdb.df['ATOM'] = ppdb.df['ATOM'][ppdb.df['ATOM']['atom_number'] != i]
        ppdb.to_pdb(path = self.host_qm_pdb, records = None, gz = False, append_newline = True)

    def get_host_mm_region_atoms(self):  
        """
        This function divides the host MM region into two sections, i.e. one preceding the QM region and one successing the QM region.
        """
        resid_num = np.loadtxt(self.residue_list)
        approximated_res_list = []
        for i in range(int(statistics.median(resid_num)) - int(int(self.num_residues)/2), int(statistics.median(resid_num)) + int(int(self.num_residues)/2)):
            approximated_res_list.append(i)
        #print(approximated_res_list)
        ppdb = PandasPdb()
        ppdb.read_pdb(self.host_pdb)
        df = ppdb.df['ATOM'][['residue_number']]
        res_list = list(set(df['residue_number'].to_list()))
        res_mm_list = list(set(res_list).difference(approximated_res_list))
        #print(res_mm_list)
        res_mm_region_I_list = []
        for i in res_mm_list:
            for j in approximated_res_list:
                if i < j:
                    res_mm_region_I_list.append(i)
        res_mm_region_I_list = list(set(res_mm_region_I_list))
        #print(res_mm_region_I_list)
        res_mm_region_II_list = list(set(res_mm_list).difference(res_mm_region_I_list))
        #print(res_mm_region_II_list)
        ppdb.read_pdb(self.host_mm_pdb)
        df = ppdb.df['ATOM'][['atom_number','residue_number','residue_name']]
        mm_region_I_index_nested_list = []
        for i in res_mm_region_I_list:
            indices = np.where(df["residue_number"] == i)
            indices = list(indices)[0]
            indices = list(indices)
            mm_region_I_index_nested_list.append(indices)  
        mm_region_I_index_list = list(itertools.chain.from_iterable(mm_region_I_index_nested_list))
        df_atom = df.iloc[mm_region_I_index_list] 
        df_atom_number = df_atom['atom_number']
        mm_region_I_atom_list = df_atom_number.values.tolist()
        mm_region_I_atoms = []
        mm_region_I_atoms.extend(mm_region_I_atom_list)
        mm_region_II_index_nested_list = []
        for i in res_mm_region_II_list:
            indices = np.where(df["residue_number"] == i)
            indices = list(indices)[0]
            indices = list(indices)
            mm_region_II_index_nested_list.append(indices)  
        mm_region_II_index_list = list(itertools.chain.from_iterable(mm_region_II_index_nested_list))
        df_atom = df.iloc[mm_region_II_index_list] 
        df_atom_number = df_atom['atom_number']
        mm_region_II_atom_list = df_atom_number.values.tolist()
        mm_region_II_atoms = []
        mm_region_II_atoms.extend(mm_region_II_atom_list)
        ppdb.read_pdb(self.host_mm_pdb)
        len_atoms = []
        for i in range(len(ppdb.df['ATOM'])):
            len_atoms.append(i+1)
        if len(mm_region_I_atoms) + len(mm_region_II_atoms) == len(len_atoms):
            print("Sum of the MM region I atoms and  MM region II atoms equals the length of list of total MM atoms")
        else:
            print("Error")
        np.savetxt(self.host_mm_region_I_atoms, mm_region_I_atoms)  
        np.savetxt(self.host_mm_region_II_atoms, mm_region_II_atoms) 

    def save_host_mm_regions_pdbs(self):  
        """
        This function saves two pdb files, one pdb preceding the QM region and one pdb file succeeding the QM region.
        """
        mm_region_I_atoms = np.loadtxt(self.host_mm_region_I_atoms)
        mm_region_I_atoms = [int(i) for i in mm_region_I_atoms]
        mm_region_II_atoms = np.loadtxt(self.host_mm_region_II_atoms)
        mm_region_II_atoms = [int(i) for i in mm_region_II_atoms]
        ppdb = PandasPdb()
        ppdb.read_pdb(self.host_mm_pdb)
        for i in mm_region_II_atoms:
            ppdb.df['ATOM'] = ppdb.df['ATOM'][ppdb.df['ATOM']['atom_number'] != i]
        ppdb.to_pdb(path = self.host_mm_region_I_pdb, records = None, gz = False, append_newline = True)
        ppdb = PandasPdb()
        ppdb.read_pdb(self.host_mm_pdb)
        for i in mm_region_I_atoms:
            ppdb.df['ATOM'] = ppdb.df['ATOM'][ppdb.df['ATOM']['atom_number'] != i]
        ppdb.to_pdb(path = self.host_mm_region_II_pdb, records = None, gz = False, append_newline = True)        
        
    def get_qm_mm_regions(self):  
        """
        This function saves a qm region comprising of the guest molecule and the selected atoms of the host molecule and a mm region comprising of the unselected atoms of the host molecule. When merging two pdb files, all lines beginning from "ATOM" will be selected and it would end with the "END"
        """
        with open(self.host_qm_pdb) as f1, open(self.qm_pdb, 'w') as f2:
            for line in f1:
                if "ATOM" in line:
                    f2.write(line)  
        with open(self.guest_pdb) as f1, open(self.qm_pdb, 'a') as f2:
            for line in f1:
                if "ATOM" in line:
                    f2.write(line) 
        with open(self.qm_pdb, 'a') as f:
            f.write("END")
        with open(self.host_mm_pdb) as f1, open(self.mm_pdb, 'w') as f2:
            for line in f1:
                if "ATOM" in line:
                    f2.write(line)  
        with open(self.mm_pdb, 'a') as f:
            f.write("END")  

####################################################################################################################################################################################
class PrepareGaussianGuest:
    
    def __init__(self, guest_pdb, n_processors, memory, charge, multiplicity, functional, basis_set, optimisation, frequency, add_keywords_I, add_keywords_II, add_keywords_III, gauss_out_file, fchk_out_file):
        self.guest_pdb = guest_pdb
        self.n_processors = n_processors
        self.memory = memory
        self.charge = charge
        self.multiplicity = multiplicity
        self.functional = functional
        self.basis_set  = basis_set
        self.optimisation = optimisation
        self.frequency = frequency
        self.gauss_out_file = gauss_out_file
        self.fchk_out_file = fchk_out_file
        self.add_keywords_I = add_keywords_I
        self.add_keywords_II = add_keywords_II
        self.add_keywords_III = add_keywords_III
        
    def write_input(self):   
        """
        This function prints out the commands section of the gaussian input file.
        """
        command_line_1 = "%Chk = " + self.guest_pdb[:-4] + ".chk" 
        command_line_2 = "%Mem = " + str(self.memory) + "GB"
        command_line_3 = "%NProcShared = " + str(self.n_processors)
        command_line_4 = "# " + self.functional + " " + self.basis_set + " " + self.optimisation + " " + self.frequency + " " + self.add_keywords_I + " " + self.add_keywords_II + " " + self.add_keywords_III
        command_line_5 = " "
        command_line_6 = self.guest_pdb[:-4] + " "  + "gaussian input file"
        command_line_7 = " "
        command_line_8 = str(self.charge) + " "  + str(self.multiplicity)
        ppdb = PandasPdb()
        ppdb.read_pdb(self.guest_pdb)
        df = ppdb.df['ATOM']
        df_1 = ppdb.df["ATOM"]["element_symbol"]
        df_1.columns = ['atom']
        df_2 = df[['x_coord', 'y_coord', 'z_coord']]
        df_merged = pd.concat([df_1, df_2], axis=1)
        command_line_9 = df_merged.to_string(header = False, index = False)
        command_line_10 = " "
        command = [command_line_1, command_line_2, command_line_3, command_line_4, 
                   command_line_5, command_line_6, command_line_7, command_line_8, 
                   command_line_9, command_line_10]
        commands = '\n'.join(command)
        with open(self.guest_pdb[:-4] + ".com", "w") as f:
            f.write(commands)
            
    def run_gaussian(self):   
        """
        This function runs the gaussian QM calculation.
        """
        execute_command = "g16" + " < " + self.guest_pdb[:-4] + ".com"  + " > "  + self.guest_pdb[:-4] + ".log" 
        with open(self.gauss_out_file, "w+") as f:
            sp.run(execute_command, shell = True, stdout = f, stderr = sp.STDOUT)
        
    def get_fchk(self):   
        """
        This function converts the checkpoint file file into the formatted chechkpoint file.
        """
        execute_command = "formchk"+ " " + self.guest_pdb[:-4] + ".chk" + " " + self.guest_pdb[:-4] + ".fchk"
        with open(self.fchk_out_file, "w+") as f:
            sp.run(execute_command, shell = True, stdout = f, stderr = sp.STDOUT)
####################################################################################################################################################################################
class PrepareGaussianHostGuest:
    
    def __init__(self, guest_pdb, host_qm_pdb, n_processors, memory, charge, multiplicity, functional, basis_set, optimisation, frequency, add_keywords_I, add_keywords_II, add_keywords_III, gauss_system_out_file, fchk_system_out_file, host_guest_input, qm_guest_charge_parameter_file, qm_host_charge_parameter_file, qm_guest_atom_charge_parameter_file):
        self.guest_pdb = guest_pdb
        self.host_qm_pdb = host_qm_pdb
        self.n_processors = n_processors
        self.memory = memory
        self.charge = charge
        self.multiplicity = multiplicity
        self.functional = functional
        self.basis_set  = basis_set
        self.optimisation = optimisation
        self.frequency = frequency
        self.add_keywords_I = add_keywords_I
        self.add_keywords_II = add_keywords_II
        self.add_keywords_III = add_keywords_III
        self.gauss_system_out_file = gauss_system_out_file
        self.fchk_system_out_file = fchk_system_out_file
        self.host_guest_input = host_guest_input
        self.qm_guest_charge_parameter_file = qm_guest_charge_parameter_file
        self.qm_host_charge_parameter_file = qm_host_charge_parameter_file
        self.qm_guest_atom_charge_parameter_file = qm_guest_atom_charge_parameter_file
        
    def write_input(self):   
        """
        This function prints out the commands section of the gaussian input file.
        """
        command_line_1 = "%Chk = " + self.host_guest_input[:-4] + ".chk" 
        command_line_2 = "%Mem = " + str(self.memory) + "GB"
        command_line_3 = "%NProcShared = " + str(self.n_processors)
        command_line_4 = "# " + self.functional + " " + self.basis_set + " " + self.optimisation + " " + self.frequency + " " + self.add_keywords_I + " " + self.add_keywords_II + " " + self.add_keywords_III
        command_line_5 = " "
        command_line_6 = "Gaussian Input File"
        command_line_7 = " "
        command_line_8 = str(self.charge) + " "  + str(self.multiplicity)
        ppdb = PandasPdb()
        ppdb.read_pdb(self.guest_pdb)
        df = ppdb.df['ATOM']
        df_1 = ppdb.df["ATOM"]["element_symbol"]
        df_1.columns = ['atom']
        df_3 = df[['x_coord', 'y_coord', 'z_coord']]
        df_2 = pd.Series(["0"] * len(df), name='decide_freeze') 
        df_merged_1 = pd.concat([df_1, df_2, df_3], axis=1)
        ppdb = PandasPdb()
        ppdb.read_pdb(self.host_qm_pdb)
        df = ppdb.df['ATOM']
        df_1 = ppdb.df["ATOM"]["element_symbol"]
        df_1.columns = ['atom']
        df_3 = df[['x_coord', 'y_coord', 'z_coord']]
        df_2 = pd.Series(["0"] * len(df), name='decide_freeze') 
        df_merged_2 = pd.concat([df_1, df_2, df_3], axis=1)    
        df_merged = pd.concat([df_merged_1, df_merged_2], axis=0)
        command_line_9 = df_merged.to_string(header = False, index = False)
        command_line_10 = " "
        command = [command_line_1, command_line_2, command_line_3, command_line_4, 
                   command_line_5, command_line_6, command_line_7, command_line_8, 
                   command_line_9, command_line_10]
        commands = '\n'.join(command)
       
        with open(self.host_guest_input, "w") as f:
            f.write(commands)
           
    def run_gaussian(self):   
        """
        This function runs the gaussian QM calculation.
        """
        execute_command = "g16" + " < " + self.host_guest_input + " > "  + self.host_guest_input[:-4] + ".log" 
        with open(self.gauss_system_out_file, "w+") as f:
            sp.run(execute_command, shell = True, stdout = f, stderr = sp.STDOUT)
            
    def get_fchk(self):   
        """
        This function converts the checkpoint file file into the formatted chechkpoint file.
        """
        execute_command = "formchk"+ " " + self.host_guest_input[:-4] + ".chk" + " " + self.host_guest_input[:-4] + ".fchk"
        with open(self.fchk_system_out_file, "w+") as f:
            sp.run(execute_command, shell = True, stdout = f, stderr = sp.STDOUT)

    def get_qm_host_guest_charges(self):
        """This function extracts charges and saves them separately for the host and guest"""
        log_file = self.host_guest_input[:-4] + ".log"
        with open (log_file, "r") as f:
            lines = f.readlines()
        for i in range(len(lines)):
            if "Fitting point charges to electrostatic potential" in lines[i]:
                to_begin = int(i)   
            if " Sum of ESP charges =" in lines[i]:
                to_end = int(i)  
        charges = lines[to_begin + 4: to_end]
        charge_list = []
        for i in range(len(charges)):
            charge_list.append(charges[i].strip().split())
        charge_list_value = []
        atom_list = []
        for i in range(len(charge_list)):
            charge_list_value.append(charge_list[i][2])
            atom_list.append(charge_list[i][1])
        ppdb = PandasPdb()
        ppdb.read_pdb(self.guest_pdb)
        df_guest = ppdb.df['ATOM']
        number_guest_atoms = df_guest.shape[0]
        data_tuples = list(zip(atom_list,charge_list_value))
        df_charge = pd.DataFrame(data_tuples, columns=['Atom','Charge'])
        number_host_atoms = df_charge.shape[0] - number_guest_atoms
        df_charge_guest = df_charge.head(number_guest_atoms)
        df_charge_host = df_charge.tail(number_host_atoms)
        df_charge_only_guest = df_charge_guest["Charge"]
        df_charge_guest.to_csv (self.qm_guest_charge_parameter_file, index = False, header = False, sep = ' ')
        df_charge_host.to_csv (self.qm_host_charge_parameter_file, index = False, header = False, sep = ' ')
        df_charge_only_guest.to_csv (self.qm_guest_atom_charge_parameter_file, index = False, header = False, sep = ' ')
####################################################################################################################################################################################
class ParameterizeGuest:
    
    def __init__(self, xyz_file, coordinate_file, unprocessed_hessian_file, bond_list_file, angle_list_file, hessian_file, atom_names_file, bond_parameter_file, vibrational_scaling, angle_parameter_file, charge_parameter_file, guest_pdb, proper_dihedral_file):
        self.xyz_file = xyz_file
        self.coordinate_file = coordinate_file
        self.unprocessed_hessian_file = unprocessed_hessian_file
        self.bond_list_file = bond_list_file
        self.angle_list_file = angle_list_file
        self.hessian_file = hessian_file
        self.atom_names_file = atom_names_file
        self.bond_parameter_file = bond_parameter_file
        self.vibrational_scaling = vibrational_scaling
        self.angle_parameter_file = angle_parameter_file
        self.charge_parameter_file = charge_parameter_file
        self.guest_pdb = guest_pdb
        self.proper_dihedral_file = proper_dihedral_file
        
    def get_xyz(self):   
        """
        This function saves a xyz file from the formatted checkpoint file.
        """
        fchk_file = self.guest_pdb[:-4] + ".fchk"
        with open (fchk_file, "r") as f:
            lines = f.readlines()
        for i in range(len(lines)):
            if "Current cartesian coordinates" in lines[i]:
                no_coordinates = re.findall(r'\d+|\d+.\d+', lines[i])
                no_coordinates = int(no_coordinates[0])
        for i in range(len(lines)):
            if "Current cartesian coordinates" in lines[i]:
                to_begin = int(i)     
        cartesian_coords = lines[to_begin + 1:to_begin + 1 + int(math.ceil(no_coordinates/5))]
        cartesian_list = []
        for i in range(len(cartesian_coords)):
            cartesian_list.append(cartesian_coords[i].strip().split())
        coordinates_list = [item for sublist in cartesian_list for item in sublist]
        list_coords = [float(x)*float(0.529) for x in coordinates_list]
        for i in range(len(lines)):
            if "Atomic numbers" in lines[i]:
                to_begin = int(i) 
            if "Nuclear charges" in lines[i]:
                to_end = int(i)    
        atomic_numbers = lines[to_begin + 1:to_end]
        atom_numbers = []
        for i in range(len(atomic_numbers)):
            atom_numbers.append(atomic_numbers[i].strip().split())
        numbers = [item for sublist in atom_numbers for item in sublist]
        N = int(no_coordinates/3)
        # Opens the new xyz file 
        file = open(self.xyz_file, "w")
        file.write(str(N) + '\n \n')
        coords = np.zeros((N,3))
        n = 0
        names = []   
        # Gives name for atomic number
        for x in range(0,len(numbers)):
            names.append(element_list[int(numbers[x]) - 1][1]) 
        # Print coordinates to new input_coords.xyz file
        for i in range(0, N):
            for j in range(0,3):
                coords[i][j] = list_coords[n]
                n = n + 1
            file.write(names[i] + str(round(coords[i][0],3)) + ' ' + str(round(coords[i][1],3)) + ' ' + str(round(coords[i][2], 3)) + '\n')
        file.close()
        np.savetxt(self.coordinate_file, coords, fmt='%s')  
        
    def get_unprocessed_hessian(self):   
        """
        This function saves a text file of the unprocessed hessian from the formatted checkpoint file.
        """    
        fchk_file = self.guest_pdb[:-4] + ".fchk"
        with open (fchk_file, "r") as f:
            lines = f.readlines()
        for i in range(len(lines)):
            if "Cartesian Force Constants" in lines[i]:
                no_hessian = re.findall(r'\d+|\d+.\d+', lines[i])
                no_hessian = int(no_hessian[0])
        for i in range(len(lines)):
            if "Cartesian Force Constants" in lines[i]:
                to_begin = int(i)     
        hessian = lines[to_begin + 1:to_begin + 1 + int(math.ceil(no_hessian/5))]
        hessian_list = []
        for i in range(len(hessian)):
            hessian_list.append(hessian[i].strip().split())
        unprocessed_Hessian = [item for sublist in hessian_list for item in sublist]
        np.savetxt(self.unprocessed_hessian_file, unprocessed_Hessian, fmt='%s')
        
    def get_bond_angles(self):   
        """
        This function saves a text file of the bonds and angles from the gaussian log file.
        """ 
        log_file = self.guest_pdb[:-4] + ".log"
        fid = open(log_file, "r")
        tline = fid.readline()
        bond_list = []
        angle_list = []
        n = 1
        n_bond = 1
        n_angle = 1
        tmp = 'R' # States if bond or angle
        B = []
        # Finds the bond and angles from the .log file
        while tline:
            tline = fid.readline()
            # Line starts at point when bond and angle list occurs
            if len(tline) > 80 and tline[0:81].strip() == '! Name  Definition              Value          Derivative Info.                !':
                tline = fid.readline()
                tline = fid.readline()
                # Stops when all bond and angles recorded 
                while ( ( tmp[0] == 'R' ) or (tmp[0] == 'A') ):
                    line = tline.split()
                    tmp = line[1]
                    # Bond or angles listed as string
                    list_terms = line[2][2:-1]
                    # Bond List 
                    if ( tmp[0] == 'R' ): 
                        x = list_terms.split(',')
                        # Subtraction due to python array indexing at 0
                        x = [(int(i) - 1 ) for i in x]
                        bond_list.append(x)
                        # Angle List 
                    if (tmp[0] == 'A' ): 
                        x = list_terms.split(',')
                        # Subtraction due to python array indexing at 0
                        x = [(int(i) - 1 ) for i in x]
                        angle_list.append(x)
                    tline = fid.readline()
                # Leave loop
                tline = -1
        np.savetxt(self.bond_list_file, bond_list, fmt='%s')
        np.savetxt(self.angle_list_file, angle_list, fmt='%s')
        
    def get_hessian(self):
        """
        This function extracts hessian from the unprocessed hessian and saves into a new file.
        """       
        unprocessed_Hessian = np.loadtxt(self.unprocessed_hessian_file) 
        fchk_file = self.guest_pdb[:-4] + ".fchk"     
        with open (fchk_file, "r") as f:
            lines = f.readlines()
        for i in range(len(lines)):
            if "Current cartesian coordinates" in lines[i]:
                no_coordinates = re.findall(r'\d+|\d+.\d+', lines[i])
                no_coordinates = int(no_coordinates[0])
        N = int(no_coordinates/3)    
        length_hessian = 3 * N
        hessian = np.zeros((length_hessian, length_hessian))
        m = 0
        # Write the hessian in a 2D array format 
        for i in range (0,(length_hessian)):
            for j in range (0,(i + 1)):
                hessian[i][j] = unprocessed_Hessian[m]
                hessian[j][i] = unprocessed_Hessian[m]
                m = m + 1
        hessian = (hessian * (627.509391))/ (0.529**2)  # Change from Hartree/bohr to kcal/mol/ang
        np.savetxt(self.hessian_file, hessian, fmt='%s')

    def get_atom_names(self):
        """
        This function saves a list of atom names from the formatted checkpoint file.
        """
        fchk_file = self.guest_pdb[:-4] + ".fchk"
        with open (fchk_file, "r") as f:
            lines = f.readlines()    
        for i in range(len(lines)):
            if "Atomic numbers" in lines[i]:
                to_begin = int(i) 
            if "Nuclear charges" in lines[i]:
                to_end = int(i) 
        atomic_numbers = lines[to_begin + 1:to_end]
        atom_numbers = []
        for i in range(len(atomic_numbers)):
            atom_numbers.append(atomic_numbers[i].strip().split())
        numbers = [item for sublist in atom_numbers for item in sublist]
        names = []   
        # Gives name for atomic number
        for x in range(0,len(numbers)):
            names.append(element_list[int(numbers[x]) - 1][1]) 
        atom_names = []
        for i in range(0,len(names)):
            atom_names.append(names[i].strip() + str(i + 1))
        np.savetxt(self.atom_names_file, atom_names, fmt ='%s')
        
    def get_bond_angle_params(self):  
        """
        This function saves the bond and angle parameter files obtained from the formatted checkpoint file. 
        """
        fchk_file = self.guest_pdb[:-4] + ".fchk"
        with open (fchk_file, "r") as f:
            lines = f.readlines()
        for i in range(len(lines)):
            if "Current cartesian coordinates" in lines[i]:
                no_coordinates = re.findall(r'\d+|\d+.\d+', lines[i])
                no_coordinates = int(no_coordinates[0])    
        N = int(no_coordinates/3)
        coords = np.loadtxt(self.coordinate_file)
        hessian = np.loadtxt(self.hessian_file)
        bond_list = np.loadtxt(self.bond_list_file, dtype = int)
        atom_names = np.loadtxt(self.atom_names_file, dtype = str) 
        # Find bond lengths
        bond_lengths = np.zeros((N, N))
        for i in range (0,N):
            for j in range(0,N):
                diff_i_j = np.array(coords[i,:]) - np.array(coords[j,:])
                bond_lengths[i][j] =  np.linalg.norm(diff_i_j)
        eigenvectors = np.empty((3, 3, N, N), dtype=complex)
        eigenvalues = np.empty((N, N, 3), dtype=complex)
        partial_hessian = np.zeros((3, 3))
        for i in range(0,N):
            for j in range(0,N):
                partial_hessian = hessian[(i * 3):((i + 1)*3),(j * 3):((j + 1)*3)]
                [a, b] = np.linalg.eig(partial_hessian)
                eigenvalues[i,j,:] = (a)
                eigenvectors[:,:,i,j] = (b)
        # Modified Seminario method to find the bond parameters and print them to file
        file_bond = open(self.bond_parameter_file, 'w')
        k_b = np.zeros(len(bond_list))
        bond_length_list = np.zeros(len(bond_list))
        unique_values_bonds = [] # Used to find average values 
        for i in range(0, len(bond_list)):
            AB = force_constant_bond(bond_list[i][0], bond_list[i][1],eigenvalues, eigenvectors, coords)
            BA = force_constant_bond(bond_list[i][1], bond_list[i][0],eigenvalues, eigenvectors, coords)
            # Order of bonds sometimes causes slight differences, find the mean
            k_b[i] = np.real(( AB + BA ) /2); 
            # Vibrational_scaling takes into account DFT deficities/ anharmocity   
            vibrational_scaling_squared = self.vibrational_scaling**2        
            k_b[i] = k_b[i] * vibrational_scaling_squared
            bond_length_list[i] =  bond_lengths[bond_list[i][0]][bond_list[i][1]]
            file_bond.write(atom_names[bond_list[i][0]] + '-' + atom_names[bond_list[i][1]] + '  ')
            file_bond.write(str("%#.5g" % k_b[i])+ '   ' + str("%#.4g" % bond_length_list[i]) +  '   ' + str(bond_list[i][0] + 1) +  '   ' + str(bond_list[i][1] + 1))
            file_bond.write('\n')
            unique_values_bonds.append([atom_names[bond_list[i][0]], atom_names[bond_list[i][1]], k_b[i], bond_length_list[i], 1 ])
        file_bond.close()
        angle_list = np.loadtxt(self.angle_list_file, dtype = int)
        # Modified Seminario method to find the angle parameters and print them to file
        file_angle = open(self.angle_parameter_file, 'w')
        k_theta = np.zeros(len(angle_list))
        theta_0 = np.zeros(len(angle_list))
        unique_values_angles = [] # Used to find average values
        # Modified Seminario part goes here ...
        # Connectivity information for Modified Seminario Method
        central_atoms_angles = []
        # A structure is created with the index giving the central atom of the angle, an array then lists the angles with that central atom. 
        # i.e. central_atoms_angles{3} contains an array of angles with central atom 3
        for i in range(0, len(coords)):
            central_atoms_angles.append([])
            for j in range(0, len(angle_list)):
                if i == angle_list[j][1]:
                    # For angle ABC, atoms A C are written to array
                    AC_array = [angle_list[j][0],  angle_list[j][2], j]
                    central_atoms_angles[i].append(AC_array)
                    # For angle ABC, atoms C A are written to array
                    CA_array = [angle_list[j][2],  angle_list[j][0], j]
                    central_atoms_angles[i].append(CA_array)
        # Sort rows by atom number
        for i in range(0, len(coords)):
            central_atoms_angles[i] = sorted(central_atoms_angles[i], key=itemgetter(0))
        # Find normals u_PA for each angle
        unit_PA_all_angles = []
        for i in range(0,len(central_atoms_angles)):
            unit_PA_all_angles.append([])
            for j in range(0, len(central_atoms_angles[i])):
                # For the angle at central_atoms_angles[i][j,:] the corresponding u_PA value is found for the plane ABC and bond AB, where ABC corresponds to the order of the arguements. This is why the reverse order was also added
                unit_PA_all_angles[i].append(u_PA_from_angles(central_atoms_angles[i][j][0], i, central_atoms_angles[i][j][1], coords))
        # Finds the contributing factors from the other angle terms scaling_factor_all_angles = cell(max(max(angle_list))); %This array will contain scaling factor and angle list position
        scaling_factor_all_angles = []
        for i in range(0,len(central_atoms_angles)):
            scaling_factor_all_angles.append([])
            for j in range(0,len(central_atoms_angles[i])):
                n = 1
                m = 1
                angles_around = 0 
                additional_contributions = 0 
                scaling_factor_all_angles[i].append([0,0]) 
                # Position in angle list
                scaling_factor_all_angles[i][j][1] =  central_atoms_angles[i][j][2] 
                # Goes through the list of angles with the same central atom and computes the term need for the modified Seminario method    
                # Forwards directions, finds the same bonds with the central atom i  
                while( ( (j + n ) < len(central_atoms_angles[i]) ) and central_atoms_angles[i][j][0] == central_atoms_angles[i][j+n][0] ):
                    additional_contributions = additional_contributions + (abs(np.dot(unit_PA_all_angles[i][j][:], unit_PA_all_angles[i][j + n][:])))**2 
                    n = n + 1
                    angles_around = angles_around + 1
                 # Backwards direction, finds the same bonds with the central atom i   
                while( ( (j - m ) >= 0 ) and central_atoms_angles[i][j][0] == central_atoms_angles[i][j-m][0] ):
                    additional_contributions = additional_contributions + (abs(np.dot(unit_PA_all_angles[i][j][:], unit_PA_all_angles[i][j - m][:] ) ) )**2
                    m = m + 1
                    angles_around =  angles_around + 1
                if (n != 1 or m != 1):
                    # Finds the mean value of the additional contribution to change to normal Seminario method comment out + part 
                    scaling_factor_all_angles[i][j][0] = 1 + ( additional_contributions / (m  + n - 2) )  
                else:
                    scaling_factor_all_angles[i][j][0] = 1
        scaling_factors_angles_list = []
        for i in range(0,len(angle_list) ):
            scaling_factors_angles_list.append([]) 
        # Orders the scaling factors according to the angle list
        for i in range(0,len(central_atoms_angles)):
            for j in range(0,len(central_atoms_angles[i]) ):
                scaling_factors_angles_list[scaling_factor_all_angles[i][j][1]].append(scaling_factor_all_angles[i][j][0]) 
        # Finds the angle force constants with the scaling factors included for each angle
        for i in range(0,len(angle_list) ):
            # Ensures that there is no difference when the ordering is changed 
            [AB_k_theta, AB_theta_0] = force_angle_constant( angle_list[i][0], angle_list[i][1], angle_list[i][2], bond_lengths, eigenvalues, eigenvectors, coords, scaling_factors_angles_list[i][0], scaling_factors_angles_list[i][1] ) 
            [BA_k_theta, BA_theta_0] = force_angle_constant( angle_list[i][2], angle_list[i][1], angle_list[i][0], bond_lengths, eigenvalues, eigenvectors, coords, scaling_factors_angles_list[i][1], scaling_factors_angles_list[i][0] ) 
            k_theta[i] = (AB_k_theta + BA_k_theta ) / 2
            theta_0[i] = (AB_theta_0 +  BA_theta_0 ) / 2
            # Vibrational_scaling takes into account DFT deficities/ anharmonicity 
            k_theta[i] =  k_theta[i] * vibrational_scaling_squared
            file_angle.write(atom_names[angle_list[i][0]] + '-' + atom_names[angle_list[i][1] ] + '-' + atom_names[angle_list[i][2]] + '  ' )
            file_angle.write(str("%#.4g" % k_theta[i]) + '   ' + str("%#.4g" % theta_0[i]) + '   ' + str(angle_list[i][0] + 1)  + '   ' + str(angle_list[i][1] + 1) + '   ' + str(angle_list[i][2] + 1))
            file_angle.write('\n')
            unique_values_angles.append([atom_names[angle_list[i][0]], atom_names[angle_list[i][1]], atom_names[angle_list[i][2]], k_theta[i], theta_0[i], 1 ])
        file_angle.close()
        
    def get_charges(self):
        """
        This function saves the charges in a text file obtained from the Gaussian log file. 
        """
        log_file = self.guest_pdb[:-4] + ".log"
        with open (log_file, "r") as f:
            lines = f.readlines()
        for i in range(len(lines)):
            if "Fitting point charges to electrostatic potential" in lines[i]:
                to_begin = int(i)   
            if " Sum of ESP charges =" in lines[i]:
                to_end = int(i)  
        charges = lines[to_begin + 4: to_end]
        charge_list = []
        for i in range(len(charges)):
            charge_list.append(charges[i].strip().split())
        charge_list_value = []
        atom_list = []
        for i in range(len(charge_list)):
            charge_list_value.append(charge_list[i][2])
            atom_list.append(charge_list[i][1])
        data_tuples = list(zip(atom_list,charge_list_value))
        df_charge = pd.DataFrame(data_tuples, columns=['Atom','Charge'])
        df_charge.to_csv (self.charge_parameter_file, index = False, header = False, sep = ' ')
        
    def get_proper_dihedrals(self):
        """
        This function saves proper dihedral angles of the guest molecule in a text file.
        """
        ppdb = PandasPdb()
        ppdb.read_pdb(self.guest_pdb)
        no_atoms = len(ppdb.df["ATOM"])
        atom_index_list = []
        for i in range(no_atoms):
            atom_index_list.append(i+1)
        possible_dihedrals = []
        for dihed in itertools.permutations(atom_index_list, 4):
            possible_dihedrals.append(dihed)
        df_bonds = pd.read_csv(self.bond_parameter_file, header = None, delimiter = r"\s+")
        df_bonds.columns = ["bond", "k_bond", "bond_length", "bond_1", "bond_2"]
        bond1 = df_bonds["bond_1"].values.tolist()
        bond2 = df_bonds["bond_2"].values.tolist()
        bond_list_list = []
        for i in range(len(bond1)):
            args = (bond1[i],bond2[i])
            bond_list_list.append(list(args))   
        reverse_bond_list_list = []
        for bonds in bond_list_list:
            reverse_bond_list_list.append(reverse_list(bonds))
        bond_lists = bond_list_list + reverse_bond_list_list
        proper_dihed_repeated  = []
        for i in range(len(possible_dihedrals)):
            dihed_frag = [possible_dihedrals[i][0],possible_dihedrals[i][1]], [possible_dihedrals[i][1],possible_dihedrals[i][2]], [possible_dihedrals[i][2],possible_dihedrals[i][3]]
            a = [dihed_frag[0] in bond_lists,dihed_frag[1] in bond_lists,dihed_frag[2] in bond_lists]
            if a == [True, True, True] :
                proper_dihed_repeated.append(possible_dihedrals[i])
        
        len_repeated_dihed_list = len(proper_dihed_repeated)
        proper_dihedrals = proper_dihed_repeated
        for x in proper_dihedrals:
            z = x[::-1]
            if z in proper_dihedrals:
                proper_dihedrals.remove(z)
        len_non_repeated_dihed_list = len(proper_dihedrals)
        #print(len_repeated_dihed_list == len_non_repeated_dihed_list * 2)  
        np.savetxt(self.proper_dihedral_file, proper_dihedrals, fmt='%s')
        #return(proper_dihedrals)
####################################################################################################################################################################################
class OffXMLGuest:
    
    def __init__(self, guest_pdb, sdf_file, xml_off_file, atom_names_file, torsion_off_file, lj_file, charge):
        self.guest_pdb = guest_pdb
        self.sdf_file = sdf_file
        self.xml_off_file = xml_off_file
        self.atom_names_file = atom_names_file
        self.torsion_off_file = torsion_off_file
        self.lj_file = lj_file
        self.charge = charge

    def generate_off_xml_amber(self):
        """
        This function generates an openforcefield xml file from the pdb file through antechamber.
        """
        mol2_file = self.guest_pdb[:-4] + ".mol2"
        in_file = self.guest_pdb[:-4] + ".in"
        frcmod_file = self.guest_pdb[:-4] + ".frcmod"
        leap_file = self.guest_pdb[:-4] + ".leap"
        prmtop_file = self.guest_pdb[:-4] + ".prmtop"
        inpcrd_file = self.guest_pdb[:-4] + ".inpcrd"
        babel_command = "babel -ipdb " + self.guest_pdb + " -omol2 " + mol2_file
        os.system(babel_command)
        antechamber_command = "antechamber -i " + mol2_file + " -fi mol2 -o " + in_file + " -fo prepi -c bcc -nc " + str(self.charge)
        os.system(antechamber_command)
        os.system("rm -rf *ANTECHAMBER* *sqm* PREP.INF NEWPDB.PDB ATOMTYPE.INF")
        parmchk2_command = "parmchk2 -i " + in_file + " -o " + frcmod_file + " -f prepi -a Y"
        os.system(parmchk2_command)
        os.system("rm -rf *ANTECHAMBER*")
        line_0 = " "
        line_1 = "loadamberprep " + in_file
        line_2 = "loadamberparams " + frcmod_file
        line_3 = "source leaprc.gaff"
        line_4 = "pdb = loadpdb " + self.guest_pdb
        line_5 = "saveamberparm pdb " + prmtop_file + " " + inpcrd_file
        line_6 = "quit"
        with open(leap_file, 'w') as f:
            f.write(line_0 + "\n")
            f.write(line_1 + "\n")
            f.write(line_2 + "\n")
            f.write(line_3 + "\n")
            f.write(line_4 + "\n")
            f.write(line_5 + "\n")
            f.write(line_6 + "\n")
        leap_command = "tleap -f " + leap_file
        os.system(leap_command)
        os.system("rm -rf leap.log")
        prmtop = simtk.openmm.app.AmberPrmtopFile(prmtop_file)
        system = prmtop.createSystem(nonbondedMethod = simtk.openmm.app.NoCutoff, constraints = None)
        with open(self.xml_off_file, "w+") as out:
            out.write(simtk.openmm.XmlSerializer.serializeSystem(system))
            
    def get_torsion_off(self):
        """
        This function saves the torsional parameters from the openforcefield file into a text file
        """
        xml_off = open(self.xml_off_file, 'r') 
        xml_off_lines = xml_off.readlines() 
        for i in range(len(xml_off_lines)):
            if "<Torsions>" in xml_off_lines[i]:
                to_begin = int(i)
            if "</Torsions>" in xml_off_lines[i]:
                to_end = int(i)        
        torsion_params = xml_off_lines[to_begin + 1:to_end]
        for i in range(len(torsion_params)):
            torsion_params[i] = torsion_params[i].strip()
            torsion_params[i] = torsion_params[i].replace("p1","class1")
            torsion_params[i] = torsion_params[i].replace("p2","class2")
            torsion_params[i] = torsion_params[i].replace("p3","class3")
            torsion_params[i] = torsion_params[i].replace("p4","class4")
            torsion_params[i] = torsion_params[i].replace("Torsion","Proper")
        k_list_off = []
        for i in range(len(torsion_params)):
            k_list_off.append(float(re.findall('\d*\.?\d+',torsion_params[i])[0]))
        k_list_off = [round(num, 10) for num in k_list_off]
        #print(k_list_off)
        class1 = []
        for i in range(len(torsion_params)):
            class1.append(int(re.findall('\d*\.?\d+',torsion_params[i])[2]))
        #print(class1)
        class2 = []
        for i in range(len(torsion_params)):
            class2.append(int(re.findall('\d*\.?\d+',torsion_params[i])[4]))
        #print(class2)
        class3 = []
        for i in range(len(torsion_params)):
            class3.append(int(re.findall('\d*\.?\d+',torsion_params[i])[6]))
        #print(class3)
        class4 = []
        for i in range(len(torsion_params)):
            class4.append(int(re.findall('\d*\.?\d+',torsion_params[i])[8]))
        #print(class4)
        periodicity = []
        for i in range(len(torsion_params)):
            periodicity.append(int(re.findall('\d*\.?\d+',torsion_params[i])[9]))
        #print(periodicity)
        phase = []
        for i in range(len(torsion_params)):
            phase.append(float(re.findall('\d*\.?\d+',torsion_params[i])[10]))
        phase = [round(num, 8) for num in phase]
        #print(phase)
        data_tuples = list(zip(k_list_off, class1, class2, class3, class4, periodicity, phase))
        df_tor = pd.DataFrame(data_tuples, columns=['k','class1','class2','class3','class4','periodicity','phase'])
        #print(df_tor.head())
        class_1_list = df_tor["class1"].to_list()
        class_2_list = df_tor["class2"].to_list()
        class_3_list = df_tor["class3"].to_list()
        class_4_list = df_tor["class4"].to_list()
        df_atoms = pd.read_csv(self.atom_names_file, header = None, delimiter = r"\s+")
        df_atoms.columns = ["atom"]
        class_1_atoms = []
        for i in class_1_list:
            class_1_atoms.append(df_atoms['atom'].iloc[i])
        class_2_atoms = []
        for i in class_2_list:
            class_2_atoms.append(df_atoms['atom'].iloc[i])
        class_3_atoms = []
        for i in class_3_list:
            class_3_atoms.append(df_atoms['atom'].iloc[i])
        class_4_atoms = []
        for i in class_4_list:
            class_4_atoms.append(df_atoms['atom'].iloc[i])
        data_tuples = list(zip(class_1_atoms, class_2_atoms, class_3_atoms, class_4_atoms))
        df_tor_atoms = pd.DataFrame(data_tuples, columns=['class1_atoms','class2_atoms','class3_atoms','class4_atoms'])
        #print(df_tor_atoms.head())
        frames = [df_tor, df_tor_atoms]  
        df_concat = pd.concat(frames, axis = 1)
        #print(df_concat.head())
        df_concat.to_csv (self.torsion_off_file, index = False, header = False,sep = ' ')
        #print(df_tor)  
        
    def get_lj_off(self):
        """
        This function saves the torsional parameters from the openforcefield file into a text file
        """       
        # Non-Bonded Parameters from OFF file
        with open (self.xml_off_file, "r") as f:
            lines = f.readlines()
        for i in range(len(lines)):
            if "<Particles>" in lines[i]:
                to_begin = int(i)
            if "</Particles>" in lines[i]:
                to_end = int(i)        
        lj_params = lines[to_begin + 1:to_end]    
        for i in range(len(lj_params)):
            lj_params[i] = lj_params[i].strip()
        epsilon_list_off = []
        for i in range(len(lj_params)):
            epsilon_list_off.append(float(re.findall('\d*\.?\d+',lj_params[i])[0]))
        epsilon_list_off = [round(num, 6) for num in epsilon_list_off]
        #print(epsilon_list_off)
        charge_list_off = []
        for i in range(len(lj_params)):
            charge_list_off.append(float(re.findall('\d*\.?\d+',lj_params[i])[1]))
        charge_list_off = [round(num, 6) for num in charge_list_off]
        #print(charge_list_off)
        sigma_list_off = []
        for i in range(len(lj_params)):
            sigma_list_off.append(float(re.findall('\d*\.?\d+',lj_params[i])[2]))
        sigma_list_off = [round(num, 6) for num in sigma_list_off]
        #print(sigma_list_off)    
        data_tuples = list(zip(charge_list_off,epsilon_list_off,sigma_list_off))
        df_lj = pd.DataFrame(data_tuples, columns=['Charge','Epsilon','Sigma'])
        df_lj.to_csv (self.lj_file, index = False, header=False,sep=' ')
        #print(df_lj)  
####################################################################################################################################################################################
class XMLGuest:
    
    def __init__(self, charge_parameter_file, atom_charge_file, xml_file, guest_pdb, bond_parameter_file, angle_parameter_file, torsion_off_file, lj_file, coulomb14scale, lj14scale, residue_name):
        self.charge_parameter_file = charge_parameter_file
        self.atom_charge_file = atom_charge_file
        self.guest_pdb = guest_pdb
        self.bond_parameter_file = bond_parameter_file
        self.angle_parameter_file = angle_parameter_file
        self.torsion_off_file = torsion_off_file
        self.lj_file = lj_file
        self.coulomb14scale = coulomb14scale
        self.lj14scale = lj14scale
        self.residue_name = residue_name
        self.xml_file = xml_file
        
    def get_qm_charges(self):
        """
        This function saves the charges in the list format obtained from the gaussian log file.
        """
        df_charges = pd.read_csv(self.charge_parameter_file, header = None, delimiter = r"\s+")
        df_charges.columns = ["atom", "charges"]
        qm_charges = df_charges["charges"].values.tolist()
        qm_charges = [round(num, 6) for num in qm_charges]
        #print(qm_charges)
        np.savetxt(self.atom_charge_file, qm_charges, fmt='%s')   
        
    def write_xml(self):
        """ 
        This function writes a xml forcefield file.
        """
        qm_charges = np.loadtxt(self.atom_charge_file)  
        #Getting the name_list
        ppdb = PandasPdb()
        ppdb.read_pdb(self.guest_pdb)
        df = ppdb.df["ATOM"]["element_symbol"]
        atom_list = df.values.tolist()
        #print(atom_list)
        list_atom_range = []
        for i in range(df.shape[0]):
            list_atom_range.append(i)
        list_atom_name_number = [atom_list[i] + str(list_atom_range[i] + 1) for i in range(len(list_atom_range))] 
        #print(list_atom_name_number)
        #print(atom_list)
        #print(list_atom_range)
        mass_number_list = []
        for i in atom_list:
            mass_number_list.append(element(i).mass)
        #print(mass_number_list) 
        name_list = []
        for i in list_atom_range:
            name_list.append("ffxml_" + str(i + 1))
        #print(name_list) 
        # Non-Bonded Parameters OFF file
        df = pd.read_csv(self.lj_file, header = None, delimiter = r"\s+")
        df.columns = ['Charge','Epsilon','Sigma']
        #print(df.head())
        charge_list_off = df["Charge"].values.tolist()
        charge_list_off = [round(num, 6) for num in charge_list_off]
        #print(charge_list_off)
        epsilon_list_off = df["Epsilon"].values.tolist()
        epsilon_list_off = [round(num, 6) for num in epsilon_list_off]
        #print(epsilon_list_off)
        sigma_list_off = df["Sigma"].values.tolist()
        sigma_list_off = [round(num, 6) for num in sigma_list_off]
        #print(sigma_list_off)
        # Bond Parameters from QM files
        df = pd.read_csv(self.bond_parameter_file, header = None, delimiter = r"\s+")
        df.columns = ["bond", "k_bond", "bond_length", "bond_1", "bond_2"]
        #print(df.head())
        bond_1_list = df["bond_1"].values.tolist()
        bond_1_list = [x - 1 for x in bond_1_list]
        bond_2_list = df["bond_2"].values.tolist()
        bond_2_list = [x - 1 for x in bond_2_list]
        #print(bond_1_list)
        #print(bond_2_list)
        k_bond_list = df["k_bond"].values.tolist()
        k_bond_list = [i* 1000.00 for i in k_bond_list]
        k_bond_list = [round(num, 10) for num in k_bond_list]
        #print(k_bond_list)
        bond_length_list = df["bond_length"].values.tolist()
        bond_length_list = [i/10.00 for i in bond_length_list]
        bond_length_list = [round(num, 6) for num in bond_length_list]
        #print(bond_length_list)
        bond_names_list = df["bond"].values.tolist()
        #print(bond_names_list)
        bond_1_names_list = []
        for i in range(len(bond_names_list)):
            bond_1_names_list.append(bond_names_list[i].partition('-')[0])
        #print(bond_1_names_list)
        bond_2_names_list = []
        for i in range(len(bond_names_list)):
            bond_2_names_list.append(bond_names_list[i].partition('-')[2])     
        #print(bond_2_names_list)
        # Angle Parameters from QM files
        df = pd.read_csv(self.angle_parameter_file, header = None, delimiter = r"\s+")
        df.columns = ["angle", "k_angle", "angle_degrees", "angle_1", "angle_2", "angle_3"]
        #print(df.head())
        angle_names_list = df["angle"].values.tolist()
        #print(angle_names_list)
        angle_1_names_list = []
        for i in range(len(angle_names_list)):
            angle_1_names_list.append(angle_names_list[i].partition('-')[0])
        #print(angle_1_names_list)
        angle_int_names_list = []
        for i in range(len(angle_names_list)):
            angle_int_names_list.append(angle_names_list[i].partition('-')[2])
        #print(angle_int_names_list)
        angle_2_names_list = []
        for i in range(len(angle_int_names_list)):
            angle_2_names_list.append(angle_int_names_list[i].partition('-')[0])           
        #print(angle_2_names_list)
        angle_3_names_list = []
        for i in range(len(angle_int_names_list)):
            angle_3_names_list.append(angle_int_names_list[i].partition('-')[2])         
        #print(angle_3_names_list)
        k_angle_list = df["k_angle"].values.tolist()
        k_angle_list = [round(num, 6) for num in k_angle_list]
        #print(k_angle_list)
        angle_list = df["angle_degrees"].values.tolist()
        angle_list = [(i * math.pi)/180.00 for i in angle_list]
        angle_list = [round(num, 6) for num in angle_list]
        #print(angle_list)
        # Torsion Parameters OFF file
        df = pd.read_csv(self.torsion_off_file, header = None, delimiter = r"\s+")
        df.columns = ['k','class1','class2','class3','class4','periodicity','phase','class1_atoms','class2_atoms','class3_atoms','class4_atoms']
        #print(df.head(20))
        k_list_off = df["k"].values.tolist()
        k_list_off = [round(num, 6) for num in k_list_off]
        #print(k_list_off)
        class1_list_off = df["class1_atoms"].values.tolist()
        #print(class1_list_off)
        class2_list_off = df["class2_atoms"].values.tolist()
        #print(class2_list_off)
        class3_list_off = df["class3_atoms"].values.tolist()
        #print(class3_list_off)
        class4_list_off = df["class4_atoms"].values.tolist()
        #print(class4_list_off)
        periodicity_list_off = df["periodicity"].values.tolist()
        #print(periodicity_list_off)
        phase_list_off = df["phase"].values.tolist()
        phase_list_off = [round(num, 6) for num in phase_list_off]
        #print(phase_list_off)
        # Writing the xml forcefield file
        xml = open(self.xml_file, "w")
        xml.write("<?xml version=" + '"' + "1.0" + '"' + " " + "?>" "\n")
        xml.write("<ForceField>" + "\n")
        xml.write("<AtomTypes>" + "\n")
        for i in range(len(list_atom_range)):
            xml.write("<Type" + " " 
                            + "class=" + '"' + list_atom_name_number[i] + '"' + " " 
                            + "element=" + '"' + atom_list[i] + '"' + " " 
                            + "mass=" + '"' + str(mass_number_list[i]) + '"' + " " 
                            + "name=" + '"' + name_list[i] + '"' 
                            + "/>"  + "\n")
        xml.write("</AtomTypes>" + "\n")
        xml.write("<Residues>" + "\n")
        xml.write("<Residue name=" + '"' + self.residue_name + '"' + ">" + "\n")
        for i in range(len(list_atom_range)):
            xml.write("<Atom" + " " 
                            + "name=" + '"' + list_atom_name_number[i] + '"' + " " 
                            + "type=" + '"' + name_list[i] + '"' 
                            + "/>"  + "\n")
        for i in range(len(bond_1_list)):
            xml.write("<Bond" + " " 
                           + "from=" + '"' + str(bond_1_list[i]) + '"' + " " 
                           + "to=" + '"' + str(bond_2_list[i]) + '"' 
                           + "/>"  + "\n")
        xml.write("</Residue>" + "\n")
        xml.write("</Residues>" + "\n")
        xml.write("<HarmonicBondForce>" + "\n")
        for i in range(len(bond_names_list)):
            xml.write("<Bond" + " " 
                           + "class1=" + '"' + bond_1_names_list[i] + '"' + " " 
                           + "class2=" + '"' + bond_2_names_list[i] + '"' + " " 
                           + "k=" + '"' + str(k_bond_list[i]) + '"' + " " 
                           + "length=" + '"' + str(bond_length_list[i]) + '"' 
                           + "/>"  + "\n")
        xml.write("</HarmonicBondForce>" + "\n")
        xml.write("<HarmonicAngleForce>" + "\n")
        for i in range(len(angle_names_list)):
            xml.write("<Angle" + " " 
                           + "angle=" + '"' + str(angle_list[i]) + '"' + " " 
                           + "class1=" + '"' + angle_1_names_list[i] + '"' + " " 
                           + "class2=" + '"' + angle_2_names_list[i] + '"' + " " 
                           + "class3=" + '"' + angle_3_names_list[i] + '"' + " " 
                           + "k=" + '"' + str(k_angle_list[i]) + '"' 
                           + "/>"  + "\n")
        xml.write("</HarmonicAngleForce>" + "\n")
        xml.write("<PeriodicTorsionForce>" + "\n")
        for i in range(len(k_list_off)):
            xml.write("<Proper" + " " 
                           + "k=" + '"' + str(k_list_off[i]) + '"' + " " 
                           + "class1=" + '"' + str(class1_list_off[i]) + '"' + " " 
                           + "class2=" + '"' + str(class2_list_off[i]) + '"' + " " 
                           + "class3=" + '"' + str(class3_list_off[i]) + '"' + " " 
                           + "class4=" + '"' + str(class4_list_off[i]) + '"' + " " 
                           + "periodicity=" + '"' + str(periodicity_list_off[i]) + '"' + " "                                              
                           + "phase=" + '"' + str(phase_list_off[i]) + '"' 
                           + "/>"  + "\n")                   
        xml.write("</PeriodicTorsionForce>" + "\n")
        xml.write("<NonbondedForce" + " " 
                       + "coulomb14scale=" + '"' + self.coulomb14scale + '"' + " " 
                       + "lj14scale=" + '"' + self.lj14scale + '"' 
                       + ">" + "\n")
        for i in range(len(charge_list_off)):
            xml.write("<Atom" + " " 
                           + "charge=" + '"' + str(qm_charges[i]) + '"' + " " 
                           + "epsilon=" + '"' + str(epsilon_list_off[i]) + '"' + " " 
                           + "sigma=" + '"' + str(sigma_list_off[i]) + '"' + " "   
                           + "type=" + '"' + name_list[i] + '"' 
                           + "/>"  + "\n")
        xml.write("</NonbondedForce>" + "\n")
        xml.write("</ForceField>" )
        xml.close()
####################################################################################################################################################################################
class ConnectPdbGuest:

    def __init__(self, xyz_file, bond_parameter_file, conect_file, guest_pdb, conect_guest_pdb, xml_file, sim_output, sim_steps):
        self.xyz_file = xyz_file
        self.bond_parameter_file = bond_parameter_file
        self.conect_file = conect_file
        self.guest_pdb = guest_pdb
        self.conect_guest_pdb = conect_guest_pdb
        self.xml_file = xml_file
        self.sim_output = sim_output
        self.sim_steps = sim_steps
        
    def create_conect(self):
        """
        This function saves a pdb file of the system with connect information.
        """    
        df_xyz = pd.read_csv(self.xyz_file, header = None, skiprows = 2, delimiter = r"\s+")
        df_xyz.columns = ["element", "x_coord", "y_coord", "z_coord"]
        no_atoms = df_xyz.shape[0]
        #print(no_atoms)
        df_bonds = pd.read_csv(self.bond_parameter_file, header = None, delimiter = r"\s+")
        df_bonds.columns = ["bond", "k_bond", "bond_length", "bond_1", "bond_2"]
        bond1 = df_bonds["bond_1"].values.tolist()
        bond2 = df_bonds["bond_2"].values.tolist()
        bond_list_list = []
        for i in range(len(bond1)):
            args = (bond1[i],bond2[i])
            bond_list_list.append(args)
        #print(bond_list_list)
        list_list_connect = []
        for i in range(no_atoms):
            atom_index = i + 1
            list_connect = []
            for j in bond_list_list:
                if j[0] == atom_index or j[1] == atom_index:
                    list_connect.append(j)
            list_list_connect.append(list_connect)
        #print(list_list_connect)
        connect_list = []
        for i in list_list_connect:
            if len(i) > 1:
                connect_list.append(i)        
        #print(connect_list)
        connect_list_merged = []
        for i in range(len(connect_list)):
            merged = [k for z in connect_list[i] for k in z]
            connect_list_merged.append(merged)
        #print(connect_list_merged)
        connect_list_merged_sorted = []
        for i in range(len(connect_list_merged)):
            sorted_list = sorted(connect_list_merged[i],key = connect_list_merged[i].count,reverse = True)
            connect_list_merged_sorted.append(sorted_list)
        #print(connect_list_merged_sorted)
        connect_list_merged_sorted_unique = []
        for i in range(len(connect_list_merged_sorted)):
            unique_list = uniq(connect_list_merged_sorted[i])
            connect_list_merged_sorted_unique.append(unique_list)
        #print(connect_list_merged_sorted_unique)
        len_list = []
        for i in range(len(connect_list_merged_sorted_unique)):
            len_list.append(len(connect_list_merged_sorted_unique[i]))
        max_len = max(len_list)    
        for i in range(len(connect_list_merged_sorted_unique)):
            if len(connect_list_merged_sorted_unique[i]) < max_len:
                diff = max_len - len(connect_list_merged_sorted_unique[i])
                for j in range(diff):
                    connect_list_merged_sorted_unique[i].append("X")            
        #print(connect_list_merged_sorted_unique)
        f = open(self.conect_file, 'w')
        for i in range(len(connect_list_merged_sorted_unique)):
            string_line = '    '.join([str(elem) for elem in connect_list_merged_sorted_unique[i]]) 
            f.write("CONECT" + "    " + string_line + "\n")
        f.close()
        conect = open(self.conect_file, 'r')
        conect_lines = conect.readlines()
        for i in range(len(conect_lines)): 
            words = conect_lines[i].split()
            if len(words) == 6 :
                conect_lines[i] = '{:>0} {:>4} {:>4} {:>4} {:>4} {:>4}'.format(*words)  
            if len(words) == 5 :
                conect_lines[i] = '{:>0} {:>4} {:>4} {:>4} {:>4}'.format(*words)  
            if len(words) == 4 :
                conect_lines[i] = '{:>0} {:>4} {:>4} {:>4}'.format(*words)  
            if len(words) == 3 :
                conect_lines[i] = '{:>0} {:>4} {:>4}'.format(*words)  
            if len(words) == 2:
                conect_lines[i] = '{:>0} {:>4}'.format(*words)  
        f = open(self.conect_file, 'w')
        for i in range(len(conect_lines)):
            f.write(conect_lines[i] + "\n")
        f.close()
        conect = open(self.conect_file, 'r')
        conect_lines = conect.readlines()
        for i in range(len(conect_lines)): 
            conect_lines[i] = conect_lines[i].replace("X", "")
        f = open(self.conect_file, 'w')
        for i in range(len(conect_lines)):
            f.write(conect_lines[i])
        f.close()
        readFile = open(self.guest_pdb)
        pdblines = readFile.readlines()
        readFile.close()
        w = open(self.conect_guest_pdb,'w')
        w.writelines([item for item in pdblines[:-1]])
        w.close()
        conect = open(self.conect_file, 'r')
        conect_lines = conect.readlines()
        with open (self.conect_guest_pdb, "a") as f:
            for i in conect_lines:
                f.write(i)
        with open (self.conect_guest_pdb, "a") as f:
            f.write("END")
            
    def run_openmm(self):
        """
        This function runs a test openmm MD simulation for the given guest pdb file 
        with the given guest xml file.
        """            
        pdb = simtk.openmm.app.PDBFile(self.conect_guest_pdb)
        forcefield = simtk.openmm.app.ForceField(self.xml_file)
        system = forcefield.createSystem(pdb.topology)
        integrator = simtk.openmm.LangevinIntegrator(300 * simtk.unit.kelvin, 1 / simtk.unit.picosecond, 0.002 * simtk.unit.picoseconds)
        simulation = simtk.openmm.app.Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy()
        simulation.reporters.append(simtk.openmm.app.PDBReporter(self.sim_output, self.sim_steps/10))
        simulation.reporters.append(simtk.openmm.app.StateDataReporter(stdout, reportInterval = int(self.sim_steps/10), step = True, potentialEnergy = True, temperature = True))
        simulation.step(self.sim_steps)
        command = "rm -rf " + str(self.sim_output)
        os.system(command)
####################################################################################################################################################################################
class TorsionDriveSims:
    """
    This class creates a directory where all dihedrals of the molecule undergo TorsionDrive simulations
    """
    def __init__(self, xyz_file, memory, charge, multiplicity, functional, basis_set, psi_input_file, dihedral_text_file, proper_dihedral_file, tor_dir, dihedral_interval, engine, xml_file, conect_guest_pdb, atom_names_file):    
        self.xyz_file = xyz_file
        self.memory = memory
        self.charge = charge
        self.multiplicity = multiplicity
        self.functional = functional
        self.basis_set = basis_set
        self.psi_input_file = psi_input_file
        self.dihedral_text_file = dihedral_text_file
        self.proper_dihedral_file = proper_dihedral_file
        self.tor_dir = tor_dir
        self.dihedral_interval = dihedral_interval
        self.engine = engine
        self.xml_file = xml_file
        self.conect_guest_pdb = conect_guest_pdb
        self.atom_names_file = atom_names_file
        
    def create_tor_dir(self):
        xyz_lines = open(self.xyz_file, 'r').readlines()[2:]
        with open(self.psi_input_file, "w") as f: 
            f.write("memory" + " " + str(self.memory)  + " " + "GB" + "\n") 
            f.write("molecule" + " " + "{" + "\n")
            f.write(str(self.charge) + " " + str(self.multiplicity) + "\n")
            for line in xyz_lines:
                f.write(line)
            f.write("}" + "\n")
            f.write("set"  + " " + "{" + "basis" + " " + self.basis_set + "}" + "\n")
            f.write("gradient" + "(" + "'" + self.functional + "'"")""\n")
        df_dihedrals = pd.read_csv(self.proper_dihedral_file, header = None, delimiter = r"\s+")
        df_dihedrals.columns = ["atom_1", "atom_2", "atom_3", "atom_4"]
        dihedrals_list_list = []
        for i in range(len(df_dihedrals)):
            dihedrals_list_list.append(df_dihedrals.iloc[i].values.tolist())
        os.system("rm -rf " + self.tor_dir)
        os.system("mkdir " + self.tor_dir)
        parent_cwd = os.getcwd()
        shutil.move(parent_cwd + "/" + self.psi_input_file, parent_cwd + "/" + self.tor_dir + "/" + self.psi_input_file)
        os.chdir(parent_cwd + "/" + self.tor_dir)
        torsion_drive_dir = os.getcwd()
        for i in range(len(dihedrals_list_list)):
            dir_name = "torsion_drive" + "_" + str(i)
            os.system("rm -rf " + dir_name)
            os.system("mkdir "  + dir_name)
            os.chdir (torsion_drive_dir + "/" + dir_name)
            with open(self.dihedral_text_file, "w") as f: 
                f.write("# dihedral definition by atom indices starting from 1" + "\n")
                f.write("# i     j     k     l" + "\n")
                i = dihedrals_list_list[i][0]
                j = dihedrals_list_list[i][1]
                k = dihedrals_list_list[i][2]
                l = dihedrals_list_list[i][3]
                f.write(" " +  "{:< 6d}".format(i) + "{:< 6d}".format(j) + "{:< 6d}".format(k) + "{:< 6d}".format(l) +  "\n")
                shutil.copy (torsion_drive_dir + "/" + self.psi_input_file, torsion_drive_dir + "/" + dir_name + "/" + self.psi_input_file)
                shutil.copy (parent_cwd + "/" + self.xml_file, torsion_drive_dir + "/" + dir_name + "/" + self.xml_file)
                shutil.copy (parent_cwd + "/" + self.conect_guest_pdb, torsion_drive_dir + "/" + dir_name + "/" + self.conect_guest_pdb)
                shutil.copy (parent_cwd + "/" + self.atom_names_file, torsion_drive_dir + "/" + dir_name + "/" + self.atom_names_file)
                os.chdir(torsion_drive_dir)   
        os.system("rm -rf "  + self.psi_input_file)
        os.chdir(parent_cwd)  

    def run_torsion_drive(self):
        df_dihedrals = pd.read_csv(self.proper_dihedral_file, header = None, delimiter = r"\s+")
        num_folders = len(df_dihedrals)
        parent_cwd = os.getcwd()
        for i in range(num_folders):
            dir_ = "torsion_drive" + "_" + str(i)
            os.chdir(parent_cwd + "/" + self.tor_dir + "/" + dir_)
            if self.engine == "psi4":
                torsion_command = "torsiondrive-launch" + " " + self.psi_input_file + " " + self.dihedral_text_file + " " + "-g" + " " + str(self.dihedral_interval) + " " + "-e" + " " + self.engine + " " + "-v"
            if self.engine == "openmm":
                torsion_command = "torsiondrive-launch" + " " + self.conect_guest_pdb + " " + self.dihedral_text_file + " " + "-g" + " " + str(self.dihedral_interval) + " " + "-e" + " " + self.engine + " " + "-v"
            os.system(torsion_command)
            print(torsion_command)
            os.chdir(parent_cwd)  
####################################################################################################################################################################################
class PrepareGaussianHost:
    
    def __init__(self, host_qm_pdb, n_processors, memory, charge, multiplicity, functional, basis_set, optimisation, frequency, add_keywords_I, add_keywords_II, add_keywords_III, gauss_out_file, fchk_out_file):
        self.host_qm_pdb = host_qm_pdb
        self.n_processors = n_processors
        self.memory = memory
        self.charge = charge
        self.multiplicity = multiplicity
        self.functional = functional
        self.basis_set  = basis_set
        self.optimisation = optimisation
        self.frequency = frequency
        self.gauss_out_file = gauss_out_file
        self.fchk_out_file = fchk_out_file
        self.add_keywords_I = add_keywords_I
        self.add_keywords_II = add_keywords_II
        self.add_keywords_III = add_keywords_III
        
    def write_input(self):   
        """
        This function prints out the commands section of the gaussian input file.
        """
        command_line_1 = "%Chk = " + self.host_qm_pdb[:-4] + ".chk" 
        command_line_2 = "%Mem = " + str(self.memory) + "GB"
        command_line_3 = "%NProcShared = " + str(self.n_processors)
        command_line_4 = "# " + self.functional + " " + self.basis_set + " " + self.optimisation + " " + self.frequency + " " + self.add_keywords_I + " " + self.add_keywords_II + " " + self.add_keywords_III
        command_line_5 = " "
        command_line_6 = self.host_qm_pdb[:-4] + " "  + "gaussian input file"
        command_line_7 = " "
        command_line_8 = str(self.charge) + " "  + str(self.multiplicity)
        ppdb = PandasPdb()
        ppdb.read_pdb(self.host_qm_pdb)
        df = ppdb.df['ATOM']
        df_1 = ppdb.df["ATOM"]["element_symbol"]
        df_1.columns = ['atom']
        df_2 = df[['x_coord', 'y_coord', 'z_coord']]
        df_merged = pd.concat([df_1, df_2], axis=1)
        command_line_9 = df_merged.to_string(header = False, index = False)
        command_line_10 = " "
        command = [command_line_1, command_line_2, command_line_3, command_line_4, 
                   command_line_5, command_line_6, command_line_7, command_line_8, 
                   command_line_9, command_line_10]
        commands = '\n'.join(command)
        with open(self.host_qm_pdb[:-4] + ".com", "w") as f:
            f.write(commands)
            
    def run_gaussian(self):   
        """
        This function runs the gaussian QM calculation.
        """
        execute_command = "g16" + " < " + self.host_qm_pdb[:-4] + ".com"  + " > "  + self.host_qm_pdb[:-4] + ".log" 
        with open(self.gauss_out_file, "w+") as f:
            sp.run(execute_command, shell = True, stdout = f, stderr = sp.STDOUT)
        
    def get_fchk(self):   
        """
        This function converts the checkpoint file file into the formatted chechkpoint file.
        """
        execute_command = "formchk"+ " " + self.host_qm_pdb[:-4] + ".chk" + " " + self.host_qm_pdb[:-4] + ".fchk"
        with open(self.fchk_out_file, "w+") as f:
            sp.run(execute_command, shell = True, stdout = f, stderr = sp.STDOUT)
####################################################################################################################################################################################
class ParameterizeHost:
    
    def __init__(self, xyz_file, coordinate_file, unprocessed_hessian_file, bond_list_file, angle_list_file, hessian_file, atom_names_file, bond_parameter_file, vibrational_scaling, angle_parameter_file, charge_parameter_file, host_qm_pdb, host_qm_params_file, host_pdb, ffxml, host_xml, sim_output, sim_steps, host_residue_name, host_singular_file, conect_pdb_txt, conect_pdb_file, host_singular_xml_file, coulomb14scale, lj14scale, reparameterised_host_xml_file, reparams_host_file):
        self.xyz_file = xyz_file
        self.coordinate_file = coordinate_file
        self.unprocessed_hessian_file = unprocessed_hessian_file
        self.bond_list_file = bond_list_file
        self.angle_list_file = angle_list_file
        self.hessian_file = hessian_file
        self.atom_names_file = atom_names_file
        self.bond_parameter_file = bond_parameter_file
        self.vibrational_scaling = vibrational_scaling
        self.angle_parameter_file = angle_parameter_file
        self.charge_parameter_file = charge_parameter_file
        self.host_qm_pdb = host_qm_pdb
        self.host_qm_params_file = host_qm_params_file
        self.host_pdb = host_pdb
        self.ffxml = ffxml
        self.host_xml = host_xml
        self.sim_output = sim_output
        self.sim_steps = sim_steps
        self.host_residue_name = host_residue_name
        self.host_singular_file = host_singular_file
        self.conect_pdb_txt = conect_pdb_txt
        self.conect_pdb_file = conect_pdb_file
        self.host_singular_xml_file = host_singular_xml_file
        self.coulomb14scale = coulomb14scale
        self.lj14scale = lj14scale
        self.reparameterised_host_xml_file = reparameterised_host_xml_file
        self.reparams_host_file = reparams_host_file

    def get_xyz(self):   
        """
        This function saves a xyz file from the formatted checkpoint file.
        """
        fchk_file = self.host_qm_pdb[:-4] + ".fchk"
        with open (fchk_file, "r") as f:
            lines = f.readlines()
        for i in range(len(lines)):
            if "Current cartesian coordinates" in lines[i]:
                no_coordinates = re.findall(r'\d+|\d+.\d+', lines[i])
                no_coordinates = int(no_coordinates[0])
        for i in range(len(lines)):
            if "Current cartesian coordinates" in lines[i]:
                to_begin = int(i)     
        cartesian_coords = lines[to_begin + 1:to_begin + 1 + int(math.ceil(no_coordinates/5))]
        cartesian_list = []
        for i in range(len(cartesian_coords)):
            cartesian_list.append(cartesian_coords[i].strip().split())
        coordinates_list = [item for sublist in cartesian_list for item in sublist]
        list_coords = [float(x)*float(0.529) for x in coordinates_list]
        for i in range(len(lines)):
            if "Atomic numbers" in lines[i]:
                to_begin = int(i) 
            if "Nuclear charges" in lines[i]:
                to_end = int(i)    
        atomic_numbers = lines[to_begin + 1:to_end]
        atom_numbers = []
        for i in range(len(atomic_numbers)):
            atom_numbers.append(atomic_numbers[i].strip().split())
        numbers = [item for sublist in atom_numbers for item in sublist]
        N = int(no_coordinates/3)
        # Opens the new xyz file 
        file = open(self.xyz_file, "w")
        file.write(str(N) + '\n \n')
        coords = np.zeros((N,3))
        n = 0
        names = []   
        # Gives name for atomic number
        for x in range(0,len(numbers)):
            names.append(element_list[int(numbers[x]) - 1][1]) 
        # Print coordinates to new input_coords.xyz file
        for i in range(0, N):
            for j in range(0,3):
                coords[i][j] = list_coords[n]
                n = n + 1
            file.write(names[i] + str(round(coords[i][0],3)) + ' ' + str(round(coords[i][1],3)) + ' ' + str(round(coords[i][2], 3)) + '\n')
        file.close()
        np.savetxt(self.coordinate_file, coords, fmt='%s')  
        
    def get_unprocessed_hessian(self):   
        """
        This function saves a text file of the unprocessed hessian from the formatted checkpoint file.
        """    
        fchk_file = self.host_qm_pdb[:-4] + ".fchk"
        with open (fchk_file, "r") as f:
            lines = f.readlines()
        for i in range(len(lines)):
            if "Cartesian Force Constants" in lines[i]:
                no_hessian = re.findall(r'\d+|\d+.\d+', lines[i])
                no_hessian = int(no_hessian[0])
        for i in range(len(lines)):
            if "Cartesian Force Constants" in lines[i]:
                to_begin = int(i)     
        hessian = lines[to_begin + 1:to_begin + 1 + int(math.ceil(no_hessian/5))]
        hessian_list = []
        for i in range(len(hessian)):
            hessian_list.append(hessian[i].strip().split())
        unprocessed_Hessian = [item for sublist in hessian_list for item in sublist]
        np.savetxt(self.unprocessed_hessian_file, unprocessed_Hessian, fmt='%s')
        
    def get_bond_angles(self):   
        """
        This function saves a text file of the bonds and angles from the gaussian log file.
        """ 
        log_file = self.host_qm_pdb[:-4] + ".log"
        fid = open(log_file, "r")
        tline = fid.readline()
        bond_list = []
        angle_list = []
        n = 1
        n_bond = 1
        n_angle = 1
        tmp = 'R' # States if bond or angle
        B = []
        # Finds the bond and angles from the .log file
        while tline:
            tline = fid.readline()
            # Line starts at point when bond and angle list occurs
            if len(tline) > 80 and tline[0:81].strip() == '! Name  Definition              Value          Derivative Info.                !':
                tline = fid.readline()
                tline = fid.readline()
                # Stops when all bond and angles recorded 
                while ( ( tmp[0] == 'R' ) or (tmp[0] == 'A') ):
                    line = tline.split()
                    tmp = line[1]
                    # Bond or angles listed as string
                    list_terms = line[2][2:-1]
                    # Bond List 
                    if ( tmp[0] == 'R' ): 
                        x = list_terms.split(',')
                        # Subtraction due to python array indexing at 0
                        x = [(int(i) - 1 ) for i in x]
                        bond_list.append(x)
                        # Angle List 
                    if (tmp[0] == 'A' ): 
                        x = list_terms.split(',')
                        # Subtraction due to python array indexing at 0
                        x = [(int(i) - 1 ) for i in x]
                        angle_list.append(x)
                    tline = fid.readline()
                # Leave loop
                tline = -1
        np.savetxt(self.bond_list_file, bond_list, fmt='%s')
        np.savetxt(self.angle_list_file, angle_list, fmt='%s')
        
    def get_hessian(self):
        """
        This function extracts hessian from the unprocessed hessian and saves into a new file.
        """       
        unprocessed_Hessian = np.loadtxt(self.unprocessed_hessian_file) 
        fchk_file = self.host_qm_pdb[:-4] + ".fchk"     
        with open (fchk_file, "r") as f:
            lines = f.readlines()
        for i in range(len(lines)):
            if "Current cartesian coordinates" in lines[i]:
                no_coordinates = re.findall(r'\d+|\d+.\d+', lines[i])
                no_coordinates = int(no_coordinates[0])
        N = int(no_coordinates/3)    
        length_hessian = 3 * N
        hessian = np.zeros((length_hessian, length_hessian))
        m = 0
        # Write the hessian in a 2D array format 
        for i in range (0,(length_hessian)):
            for j in range (0,(i + 1)):
                hessian[i][j] = unprocessed_Hessian[m]
                hessian[j][i] = unprocessed_Hessian[m]
                m = m + 1
        hessian = (hessian * (627.509391))/ (0.529**2)  # Change from Hartree/bohr to kcal/mol/ang
        np.savetxt(self.hessian_file, hessian, fmt='%s')

    def get_atom_names(self):
        """
        This function saves a list of atom names from the formatted checkpoint file.
        """
        fchk_file = self.host_qm_pdb[:-4] + ".fchk"
        with open (fchk_file, "r") as f:
            lines = f.readlines()    
        for i in range(len(lines)):
            if "Atomic numbers" in lines[i]:
                to_begin = int(i) 
            if "Nuclear charges" in lines[i]:
                to_end = int(i) 
        atomic_numbers = lines[to_begin + 1:to_end]
        atom_numbers = []
        for i in range(len(atomic_numbers)):
            atom_numbers.append(atomic_numbers[i].strip().split())
        numbers = [item for sublist in atom_numbers for item in sublist]
        names = []   
        # Gives name for atomic number
        for x in range(0,len(numbers)):
            names.append(element_list[int(numbers[x]) - 1][1]) 
        atom_names = []
        for i in range(0,len(names)):
            atom_names.append(names[i].strip() + str(i + 1))
        np.savetxt(self.atom_names_file, atom_names, fmt ='%s')
        
    def get_bond_angle_params(self):  
        """
        This function saves the bond and angle parameter files obtained from the formatted checkpoint file. 
        """
        fchk_file = self.host_qm_pdb[:-4] + ".fchk"
        with open (fchk_file, "r") as f:
            lines = f.readlines()
        for i in range(len(lines)):
            if "Current cartesian coordinates" in lines[i]:
                no_coordinates = re.findall(r'\d+|\d+.\d+', lines[i])
                no_coordinates = int(no_coordinates[0])    
        N = int(no_coordinates/3)
        coords = np.loadtxt(self.coordinate_file)
        hessian = np.loadtxt(self.hessian_file)
        bond_list = np.loadtxt(self.bond_list_file, dtype = int)
        atom_names = np.loadtxt(self.atom_names_file, dtype = str) 
        # Find bond lengths
        bond_lengths = np.zeros((N, N))
        for i in range (0,N):
            for j in range(0,N):
                diff_i_j = np.array(coords[i,:]) - np.array(coords[j,:])
                bond_lengths[i][j] =  np.linalg.norm(diff_i_j)
        eigenvectors = np.empty((3, 3, N, N), dtype=complex)
        eigenvalues = np.empty((N, N, 3), dtype=complex)
        partial_hessian = np.zeros((3, 3))
        for i in range(0,N):
            for j in range(0,N):
                partial_hessian = hessian[(i * 3):((i + 1)*3),(j * 3):((j + 1)*3)]
                [a, b] = np.linalg.eig(partial_hessian)
                eigenvalues[i,j,:] = (a)
                eigenvectors[:,:,i,j] = (b)
        # Modified Seminario method to find the bond parameters and print them to file
        file_bond = open(self.bond_parameter_file, 'w')
        k_b = np.zeros(len(bond_list))
        bond_length_list = np.zeros(len(bond_list))
        unique_values_bonds = [] # Used to find average values 
        for i in range(0, len(bond_list)):
            AB = force_constant_bond(bond_list[i][0], bond_list[i][1],eigenvalues, eigenvectors, coords)
            BA = force_constant_bond(bond_list[i][1], bond_list[i][0],eigenvalues, eigenvectors, coords)
            # Order of bonds sometimes causes slight differences, find the mean
            k_b[i] = np.real(( AB + BA ) /2); 
            # Vibrational_scaling takes into account DFT deficities/ anharmocity   
            vibrational_scaling_squared = self.vibrational_scaling**2        
            k_b[i] = k_b[i] * vibrational_scaling_squared
            bond_length_list[i] =  bond_lengths[bond_list[i][0]][bond_list[i][1]]
            file_bond.write(atom_names[bond_list[i][0]] + '-' + atom_names[bond_list[i][1]] + '  ')
            file_bond.write(str("%#.5g" % k_b[i])+ '   ' + str("%#.4g" % bond_length_list[i]) +  '   ' + str(bond_list[i][0] + 1) +  '   ' + str(bond_list[i][1] + 1))
            file_bond.write('\n')
            unique_values_bonds.append([atom_names[bond_list[i][0]], atom_names[bond_list[i][1]], k_b[i], bond_length_list[i], 1 ])
        file_bond.close()
        angle_list = np.loadtxt(self.angle_list_file, dtype = int)
        # Modified Seminario method to find the angle parameters and print them to file
        file_angle = open(self.angle_parameter_file, 'w')
        k_theta = np.zeros(len(angle_list))
        theta_0 = np.zeros(len(angle_list))
        unique_values_angles = [] # Used to find average values
        # Modified Seminario part goes here ...
        # Connectivity information for Modified Seminario Method
        central_atoms_angles = []
        # A structure is created with the index giving the central atom of the angle, an array then lists the angles with that central atom. 
        # i.e. central_atoms_angles{3} contains an array of angles with central atom 3
        for i in range(0, len(coords)):
            central_atoms_angles.append([])
            for j in range(0, len(angle_list)):
                if i == angle_list[j][1]:
                    # For angle ABC, atoms A C are written to array
                    AC_array = [angle_list[j][0],  angle_list[j][2], j]
                    central_atoms_angles[i].append(AC_array)
                    # For angle ABC, atoms C A are written to array
                    CA_array = [angle_list[j][2],  angle_list[j][0], j]
                    central_atoms_angles[i].append(CA_array)
        # Sort rows by atom number
        for i in range(0, len(coords)):
            central_atoms_angles[i] = sorted(central_atoms_angles[i], key=itemgetter(0))
        # Find normals u_PA for each angle
        unit_PA_all_angles = []
        for i in range(0,len(central_atoms_angles)):
            unit_PA_all_angles.append([])
            for j in range(0, len(central_atoms_angles[i])):
                # For the angle at central_atoms_angles[i][j,:] the corresponding u_PA value is found for the plane ABC and bond AB, where ABC corresponds to the order of the arguements. This is why the reverse order was also added
                unit_PA_all_angles[i].append(u_PA_from_angles(central_atoms_angles[i][j][0], i, central_atoms_angles[i][j][1], coords))
        # Finds the contributing factors from the other angle terms scaling_factor_all_angles = cell(max(max(angle_list))); %This array will contain scaling factor and angle list position
        scaling_factor_all_angles = []
        for i in range(0,len(central_atoms_angles)):
            scaling_factor_all_angles.append([])
            for j in range(0,len(central_atoms_angles[i])):
                n = 1
                m = 1
                angles_around = 0 
                additional_contributions = 0 
                scaling_factor_all_angles[i].append([0,0]) 
                # Position in angle list
                scaling_factor_all_angles[i][j][1] =  central_atoms_angles[i][j][2] 
                # Goes through the list of angles with the same central atom and computes the term need for the modified Seminario method    
                # Forwards directions, finds the same bonds with the central atom i  
                while( ( (j + n ) < len(central_atoms_angles[i]) ) and central_atoms_angles[i][j][0] == central_atoms_angles[i][j+n][0] ):
                    additional_contributions = additional_contributions + (abs(np.dot(unit_PA_all_angles[i][j][:], unit_PA_all_angles[i][j + n][:])))**2 
                    n = n + 1
                    angles_around = angles_around + 1
                 # Backwards direction, finds the same bonds with the central atom i   
                while( ( (j - m ) >= 0 ) and central_atoms_angles[i][j][0] == central_atoms_angles[i][j-m][0] ):
                    additional_contributions = additional_contributions + (abs(np.dot(unit_PA_all_angles[i][j][:], unit_PA_all_angles[i][j - m][:] ) ) )**2
                    m = m + 1
                    angles_around =  angles_around + 1
                if (n != 1 or m != 1):
                    # Finds the mean value of the additional contribution to change to normal Seminario method comment out + part 
                    scaling_factor_all_angles[i][j][0] = 1 + ( additional_contributions / (m  + n - 2) )  
                else:
                    scaling_factor_all_angles[i][j][0] = 1
        scaling_factors_angles_list = []
        for i in range(0,len(angle_list) ):
            scaling_factors_angles_list.append([]) 
        # Orders the scaling factors according to the angle list
        for i in range(0,len(central_atoms_angles)):
            for j in range(0,len(central_atoms_angles[i]) ):
                scaling_factors_angles_list[scaling_factor_all_angles[i][j][1]].append(scaling_factor_all_angles[i][j][0]) 
        # Finds the angle force constants with the scaling factors included for each angle
        for i in range(0,len(angle_list) ):
            # Ensures that there is no difference when the ordering is changed 
            [AB_k_theta, AB_theta_0] = force_angle_constant( angle_list[i][0], angle_list[i][1], angle_list[i][2], bond_lengths, eigenvalues, eigenvectors, coords, scaling_factors_angles_list[i][0], scaling_factors_angles_list[i][1] ) 
            [BA_k_theta, BA_theta_0] = force_angle_constant( angle_list[i][2], angle_list[i][1], angle_list[i][0], bond_lengths, eigenvalues, eigenvectors, coords, scaling_factors_angles_list[i][1], scaling_factors_angles_list[i][0] ) 
            k_theta[i] = (AB_k_theta + BA_k_theta ) / 2
            theta_0[i] = (AB_theta_0 +  BA_theta_0 ) / 2
            # Vibrational_scaling takes into account DFT deficities/ anharmonicity 
            k_theta[i] =  k_theta[i] * vibrational_scaling_squared
            file_angle.write(atom_names[angle_list[i][0]] + '-' + atom_names[angle_list[i][1] ] + '-' + atom_names[angle_list[i][2]] + '  ' )
            file_angle.write(str("%#.4g" % k_theta[i]) + '   ' + str("%#.4g" % theta_0[i]) + '   ' + str(angle_list[i][0] + 1)  + '   ' + str(angle_list[i][1] + 1) + '   ' + str(angle_list[i][2] + 1))
            file_angle.write('\n')
            unique_values_angles.append([atom_names[angle_list[i][0]], atom_names[angle_list[i][1]], atom_names[angle_list[i][2]], k_theta[i], theta_0[i], 1 ])
        file_angle.close()
        
    def get_charges(self):
        """
        This function saves the charges in a text file obtained from the Gaussian log file. 
        """
        log_file = self.host_qm_pdb[:-4] + ".log"
        with open (log_file, "r") as f:
            lines = f.readlines()
        for i in range(len(lines)):
            if "Fitting point charges to electrostatic potential" in lines[i]:
                to_begin = int(i)   
            if " Sum of ESP charges =" in lines[i]:
                to_end = int(i)  
        charges = lines[to_begin + 4: to_end]
        charge_list = []
        for i in range(len(charges)):
            charge_list.append(charges[i].strip().split())
        charge_list_value = []
        atom_list = []
        for i in range(len(charge_list)):
            charge_list_value.append(charge_list[i][2])
            atom_list.append(charge_list[i][1])
        data_tuples = list(zip(atom_list,charge_list_value))
        df_charge = pd.DataFrame(data_tuples, columns=['Atom','Charge'])
        df_charge.to_csv (self.charge_parameter_file, index = False, header = False, sep = ' ')

    def write_host_params(self):
        """
        This function saves the parameters obtained from the QM log files in a text file. 
        """
        # Charges from QM files
        df_charges = pd.read_csv(self.charge_parameter_file, header = None, delimiter = r"\s+")
        df_charges.columns = ["atom", "charges"]
        qm_charges = df_charges["charges"].values.tolist()
        qm_charges = [round(num, 6) for num in qm_charges]
        #print(qm_charges)
        # Bond Parameters from QM files
        ppdb = PandasPdb()
        ppdb.read_pdb(self.host_qm_pdb)
        atom_name_list = ppdb.df["ATOM"]["atom_number"].values.tolist()
        #print(atom_name_list)
        df = pd.read_csv(self.bond_parameter_file, header = None, delimiter = r"\s+")
        df.columns = ["bond", "k_bond", "bond_length", "bond_1", "bond_2"]
        #print(df.head())
        bond_1_list = df["bond_1"].values.tolist()
        bond_1_list = [x - 1 + min(atom_name_list) for x in bond_1_list]
        bond_2_list = df["bond_2"].values.tolist()
        bond_2_list = [x - 1 + min(atom_name_list) for x in bond_2_list]
        #print(bond_1_list)
        #print(bond_2_list)
        k_bond_list = df["k_bond"].values.tolist()
        k_bond_list = [i* 1000.00 for i in k_bond_list]
        k_bond_list = [round(num, 10) for num in k_bond_list]
        #print(k_bond_list)
        bond_length_list = df["bond_length"].values.tolist()
        bond_length_list = [i/10.00 for i in bond_length_list]
        bond_length_list = [round(num, 6) for num in bond_length_list]
        #print(bond_length_list)
        # Angle Parameters from QM files
        ppdb = PandasPdb()
        ppdb.read_pdb(self.host_qm_pdb)
        atom_name_list = ppdb.df["ATOM"]["atom_number"].values.tolist()
        #print(atom_name_list)
        df = pd.read_csv(self.angle_parameter_file, header = None, delimiter = r"\s+")
        df.columns = ["angle", "k_angle", "angle_degrees", "angle_1", "angle_2", "angle_3"]
        #print(df.head())
        angle_1_list = df["angle_1"].values.tolist()
        angle_1_list = [x - 1 + min(atom_name_list) for x in angle_1_list]
        #print(angle_1_list)
        angle_2_list = df["angle_2"].values.tolist()
        angle_2_list = [x - 1 + min(atom_name_list) for x in angle_2_list]
        #print(angle_2_list)
        angle_3_list = df["angle_3"].values.tolist()
        angle_3_list = [x - 1 + min(atom_name_list) for x in angle_3_list]
        #print(angle_3_list)
        k_angle_list = df["k_angle"].values.tolist()
        k_angle_list = [round(num, 6) for num in k_angle_list]
        #print(k_angle_list)
        angle_list = df["angle_degrees"].values.tolist()
        angle_list = [(i * math.pi)/180.00 for i in angle_list]
        angle_list = [round(num, 6) for num in angle_list]
        #print(angle_list)
        xml = open(self.host_qm_params_file, "w")
        xml.write("Begin writing the Bond Parameters" + "\n")
        for i in range(len(k_bond_list)):
            xml.write("<Bond"  + " " 
                               + "class1=" + '"' + str(bond_1_list[i]) + '"' + " " 
                               + "class2=" + '"' + str(bond_2_list[i]) + '"' + " " 
                               + "k=" + '"' + str(k_bond_list[i]) + '"' + " " 
                               + "length=" + '"' + str(bond_length_list[i]) + '"' 
                               + "/>"  + "\n")  
        xml.write("Finish writing the Bond Parameters" + "\n")
        xml.write("Begin writing the Angle Parameters" + "\n")
        for i in range(len(k_angle_list)):
            xml.write("<Angle" + " " 
                               + "angle=" + '"' + str(angle_list[i]) + '"' + " " 
                               + "class1=" + '"' + str(angle_1_list[i]) + '"' + " " 
                               + "class2=" + '"' + str(angle_2_list[i]) + '"' + " " 
                               + "class3=" + '"' + str(angle_3_list[i]) + '"' + " " 
                               + "k=" + '"' + str(k_angle_list[i]) + '"' 
                               + "/>"  + "\n")
        xml.write("Finish writing the Angle Parameters" + "\n")    
        xml.write("Begin writing the Charge Parameters" + "\n")
        for i in range(len(qm_charges)):
            xml.write("<Atom" + " " 
                              + "charge=" + '"' + str(qm_charges[i]) + '"' + " " 
                              + "epsilon=" + '"' + str(0.00) + '"' + " " 
                              + "sigma=" + '"' + str(0.00) + '"' + " "   
                              + "type=" + '"' + "host_" + str(atom_name_list[i]) + '"' 
                              + "/>"  + "\n")
        xml.write("Finish writing the Charge Parameters" + "\n")    
        xml.close()

    def serialise_host(self):
        pdb = simtk.openmm.app.PDBFile(self.host_pdb)
        forcefield = simtk.openmm.app.ForceField(self.ffxml)
        system = forcefield.createSystem(pdb.topology)
        integrator = simtk.openmm.LangevinIntegrator(300 * simtk.unit.kelvin, 1 / simtk.unit.picosecond, 0.002 * simtk.unit.picoseconds)
        simulation = simtk.openmm.app.Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy()
        state = simulation.context.getState(getEnergy=True)
        energy = state.getPotentialEnergy()
        print("The potential energy of the system is ", energy)
        simulation.reporters.append(simtk.openmm.app.PDBReporter(self.sim_output, self.sim_steps/10))
        simulation.reporters.append(simtk.openmm.app.StateDataReporter(stdout, reportInterval = int(self.sim_steps/10), step = True, potentialEnergy = True, temperature = True))
        simulation.step(self.sim_steps)
        command = "rm -rf " + str(self.sim_output)
        os.system(command)
        with open(self.host_xml, 'w') as f:
            f.write(simtk.openmm.XmlSerializer.serialize(system))
            
    def create_host_xml(self):
        """ 
        This function takes in the GAFF file and generates an xml forcefield file.
        """
        ppdb = PandasPdb()
        ppdb.read_pdb(self.host_pdb)
        atom_numbers  = ppdb.df["ATOM"]["atom_number"].values.tolist()
        atom_numbers = [str(i) for i in atom_numbers]
        df_atom_names = pd.DataFrame(atom_numbers,columns=['atom_name'])
        res_name_list = [self.host_residue_name] * df_atom_names.shape[0]
        df_res_name = pd.DataFrame(res_name_list,columns=['residue_name'])
        res_number_list = ["1"] * df_atom_names.shape[0]
        df_res_number = pd.DataFrame(res_number_list,columns=['residue_number'])
        ppdb.df["ATOM"]["atom_name"] = df_atom_names["atom_name"]
        ppdb.df["ATOM"]["residue_name"] = df_res_name["residue_name"]
        ppdb.df["ATOM"]["residue_number"] = df_res_number["residue_number"]
        #print(ppdb.df["ATOM"].head())
        ppdb.to_pdb(path = self.host_singular_file)
        intermediate_pdb_file = self.host_singular_file[:-3] + "_intermediate.pdb"
        with open(self.host_singular_file, 'r') as f1, open(intermediate_pdb_file, 'w') as f2:
            for line in f1.readlines():
                if not (line.startswith('CONECT')):
                    f2.write(line)
        command = "rm -rf " + self.host_singular_file
        os.system(command)
        command = "mv " + intermediate_pdb_file + " " + self.host_singular_file
        os.system(command)
        ppdb = PandasPdb()
        ppdb.read_pdb(self.host_singular_file)
        #print(ppdb.df["ATOM"].tail())
        xml_off = open(self.host_xml, 'r') 
        xml_off_lines = xml_off.readlines() 
        for i in range(len(xml_off_lines)):
            if "</PeriodicBoxVectors>" in xml_off_lines[i]:
                to_begin = int(i)
            if "<Constraints/>" in xml_off_lines[i]:
                to_end = int(i)   
        atom_masses = xml_off_lines[to_begin + 2 : to_end - 1]
        list_atom_masses = []
        for i in range(len(atom_masses)):
            list_atom_masses.append(float(re.findall('\d*\.?\d+',atom_masses[i])[0]))
        #print(list_atom_masses)
        #print(len(list_atom_masses))
        ppdb = PandasPdb()
        ppdb.read_pdb(self.host_singular_file)
        element_symbol = ppdb.df["ATOM"]["element_symbol"].values.tolist()
        #print(element_symbol)
        #print(len(element_symbol))
        class_list = ppdb.df["ATOM"]["atom_name"].values.tolist()
        #print(class_list)
        #print(len(class_list))
        name_list = ["host_" + i for i in class_list]
        #print(name_list)
        #print(len(name_list))
        #Bond Parameters
        xml_off = open(self.host_xml, 'r') 
        xml_off_lines = xml_off.readlines() 
        for i in range(len(xml_off_lines)):
            if "<Bonds>" in xml_off_lines[i]:
                to_begin = int(i)
            if "</Bonds>" in xml_off_lines[i]:
                to_end = int(i)  
        bond_params = xml_off_lines[to_begin + 1:to_end]
        for i in range(len(bond_params)):
            bond_params[i] = bond_params[i].strip()
            bond_params[i] = bond_params[i].replace("p1","class1")
            bond_params[i] = bond_params[i].replace("p2","class2")
            bond_params[i] = bond_params[i].replace("d=","length=") 
        k_list_bonds_off = []
        for i in range(len(bond_params)):
            k_list_bonds_off.append(float(re.findall('\d*\.?\d+',bond_params[i])[1]))
        k_list_bonds_off = [round(num, 10) for num in k_list_bonds_off]
        #print(k_list_bonds_off)
        #print(len(k_list_bonds_off))
        class1_bonds = []
        for i in range(len(bond_params)):
            class1_bonds.append(int(re.findall('\d*\.?\d+',bond_params[i])[3]))
        #print(class1_bonds)
        #print(len(class1_bonds))
        class2_bonds = []
        for i in range(len(bond_params)):
            class2_bonds.append(int(re.findall('\d*\.?\d+',bond_params[i])[5]))
        #print(class2_bonds)
        #print(len(class2_bonds))
        bond_list_off = []
        for i in range(len(bond_params)):
            bond_list_off.append(float(re.findall('\d*\.?\d+',bond_params[i])[0]))
        bond_list_off = [round(num, 10) for num in bond_list_off]
        #print(bond_list_off)
        #print(len(bond_list_off))
        data_tuples = list(zip(bond_list_off, class1_bonds, class2_bonds, k_list_bonds_off))
        df_bond = pd.DataFrame(data_tuples, columns=['bond','class1','class2','k'])
        #print(df_bond.head())
        class_1_list_bonds = df_bond["class1"].to_list()
        class_2_list_bonds = df_bond["class2"].to_list()
        df_atoms = pd.DataFrame(class_list,columns=['atom'])
        class_1_atoms_bonds = []
        for i in class_1_list_bonds:
            class_1_atoms_bonds.append(df_atoms['atom'].iloc[i])
        class_2_atoms_bonds = []
        for i in class_2_list_bonds:
            class_2_atoms_bonds.append(df_atoms['atom'].iloc[i])
        data_tuples = list(zip(class_1_atoms_bonds, class_2_atoms_bonds))
        df_bond_atoms = pd.DataFrame(data_tuples, columns=['class1_atoms','class2_atoms'])
        #print(df_bond_atoms.head()
        frames = [df_bond, df_bond_atoms]  
        df_concat_bonds = pd.concat(frames, axis = 1)   
        #print(df_concat_bonds.tail())
        k_list_bonds_off = df_concat_bonds["k"].values.tolist()
        k_list_bonds_off = [round(num, 6) for num in k_list_bonds_off]
        #print(k_list_bonds_off)
        #print(len(k_list_bonds_off))
        class1_bond_list_off = df_concat_bonds["class1_atoms"].values.tolist()
        #print(class1_bond_list_off)
        #print(len(class1_bond_list_off))
        class2_bond_list_off = df_concat_bonds["class2_atoms"].values.tolist()
        #print(class2_bond_list_off)
        #print(len(class2_bond_list_off))
        bond_list_off = df_concat_bonds["bond"].values.tolist()
        bond_list_off = [round(num, 6) for num in bond_list_off]
        #print(bond_list_off)
        #print(len(bond_list_off))
        #Angle Parameters
        xml_off = open(self.host_xml, 'r') 
        xml_off_lines = xml_off.readlines() 
        for i in range(len(xml_off_lines)):
            if "<Angles>" in xml_off_lines[i]:
                to_begin = int(i)
            if "</Angles>" in xml_off_lines[i]:
                to_end = int(i)  
        angle_params = xml_off_lines[to_begin + 1:to_end]
        for i in range(len(angle_params)):
            angle_params[i] = angle_params[i].strip()
            angle_params[i] = angle_params[i].replace("p1","class1")
            angle_params[i] = angle_params[i].replace("p2","class2")
            angle_params[i] = angle_params[i].replace("p3","class3")
            angle_params[i] = angle_params[i].replace("a=","angle=")  
        k_list_angles_off = []
        for i in range(len(angle_params)):
            k_list_angles_off.append(float(re.findall('\d*\.?\d+',angle_params[i])[1]))
        k_list_angles_off = [round(num, 10) for num in k_list_angles_off]
        #print(k_list_angles_off)
        #print(len(k_list_angles_off))
        class1_angles = []
        for i in range(len(angle_params)):
            class1_angles.append(int(re.findall('\d*\.?\d+',angle_params[i])[3]))
        #print(class1_angles)
        class2_angles = []
        for i in range(len(angle_params)):
            class2_angles.append(int(re.findall('\d*\.?\d+',angle_params[i])[5]))
        #print(class2_angles)
        class3_angles = []
        for i in range(len(angle_params)):
            class3_angles.append(int(re.findall('\d*\.?\d+',angle_params[i])[7]))
        #print(class3_angles)
        angle_list_off = []
        for i in range(len(angle_params)):
            angle_list_off.append(float(re.findall('\d*\.?\d+',angle_params[i])[0]))
        angle_list_off = [round(num, 10) for num in angle_list_off]
        #print(angle_list_off)
        data_tuples = list(zip(angle_list_off, class1_angles, class2_angles, class3_angles, k_list_angles_off))
        df_angle = pd.DataFrame(data_tuples, columns=['angle','class1','class2','class3','k'])
        #print(df_angle.head())
        class_1_list_angles = df_angle["class1"].to_list()
        class_2_list_angles = df_angle["class2"].to_list()
        class_3_list_angles = df_angle["class3"].to_list()
        df_atoms = pd.DataFrame(class_list,columns=['atom'])
        class_1_atoms_angles = []
        for i in class_1_list_angles:
            class_1_atoms_angles.append(df_atoms['atom'].iloc[i])
        class_2_atoms_angles = []
        for i in class_2_list_angles:
            class_2_atoms_angles.append(df_atoms['atom'].iloc[i])
        class_3_atoms_angles = []
        for i in class_3_list_angles:
            class_3_atoms_angles.append(df_atoms['atom'].iloc[i])
        data_tuples = list(zip(class_1_atoms_angles, class_2_atoms_angles, class_3_atoms_angles))
        df_angle_atoms = pd.DataFrame(data_tuples, columns=['class1_atoms','class2_atoms','class3_atoms'])
        #print(df_angle_atoms.head())
        frames = [df_angle, df_angle_atoms]  
        df_concat_angles = pd.concat(frames, axis = 1)   
        #print(df_concat_angles.tail())
        k_list_angles_off = df_concat_angles["k"].values.tolist()
        k_list_angles_off = [round(num, 6) for num in k_list_angles_off]
        #print(k_list_angles_off)
        #print(len(k_list_angles_off))
        class1_angle_list_off = df_concat_angles["class1_atoms"].values.tolist()
        #print(class1_angle_list_off)
        #print(len(class1_angle_list_off))
        class2_angle_list_off = df_concat_angles["class2_atoms"].values.tolist()
        #print(class2_angle_list_off)
        #print(len(class2_angle_list_off))
        class3_angle_list_off = df_concat_angles["class3_atoms"].values.tolist()
        #print(class3_angle_list_off)
        #print(len(class3_angle_list_off))
        angle_list_off = df_concat_angles["angle"].values.tolist()
        angle_list_off = [round(num, 6) for num in angle_list_off]
        #print(angle_list_off)
        #print(len(angle_list_off))
        #Torsion Parameters
        xml_off = open(self.host_xml, 'r') 
        xml_off_lines = xml_off.readlines() 
        for i in range(len(xml_off_lines)):
            if "<Torsions>" in xml_off_lines[i]:
                to_begin = int(i)
            if "</Torsions>" in xml_off_lines[i]:
                to_end = int(i)  
        torsion_params = xml_off_lines[to_begin + 1:to_end]
        for i in range(len(torsion_params)):
            torsion_params[i] = torsion_params[i].strip()
            torsion_params[i] = torsion_params[i].replace("p1","class1")
            torsion_params[i] = torsion_params[i].replace("p2","class2")
            torsion_params[i] = torsion_params[i].replace("p3","class3")
            torsion_params[i] = torsion_params[i].replace("p4","class4")
            torsion_params[i] = torsion_params[i].replace("Torsion","Proper")  
        k_list_torsions_off = []
        for i in range(len(torsion_params)):
            k_list_torsions_off.append(float(re.findall('\d*\.?\d+',torsion_params[i])[0]))
        k_list_torsions_off = [round(num, 10) for num in k_list_torsions_off]
        #print(k_list_torsions_off)
        class1_torsions = []
        for i in range(len(torsion_params)):
            class1_torsions.append(int(re.findall('\d*\.?\d+',torsion_params[i])[2]))
        #print(class1_torsions)
        class2_torsions = []
        for i in range(len(torsion_params)):
            class2_torsions.append(int(re.findall('\d*\.?\d+',torsion_params[i])[4]))
        #print(class2_torsions)
        class3_torsions = []
        for i in range(len(torsion_params)):
            class3_torsions.append(int(re.findall('\d*\.?\d+',torsion_params[i])[6]))
        #print(class3_torsions)
        class4_torsions = []
        for i in range(len(torsion_params)):
            class4_torsions.append(int(re.findall('\d*\.?\d+',torsion_params[i])[8]))
        #print(class4_torsions)
        periodicity_torsions = []
        for i in range(len(torsion_params)):
            periodicity_torsions.append(int(re.findall('\d*\.?\d+',torsion_params[i])[9]))
        #print(periodicity_torsions)
        phase_torsions = []
        for i in range(len(torsion_params)):
            phase_torsions.append(float(re.findall('\d*\.?\d+',torsion_params[i])[10]))
        phase_torsions = [round(num, 8) for num in phase_torsions]
        #print(phase_torsions)
        data_tuples = list(zip(k_list_torsions_off, class1_torsions, class2_torsions, class3_torsions, class4_torsions, periodicity_torsions, phase_torsions))
        df_tor = pd.DataFrame(data_tuples, columns=['k','class1','class2','class3','class4','periodicity','phase'])
        #print(df_tor.head())
        class_1_list_torsions = df_tor["class1"].to_list()
        class_2_list_torsions = df_tor["class2"].to_list()
        class_3_list_torsions = df_tor["class3"].to_list()
        class_4_list_torsions = df_tor["class4"].to_list()
        df_atoms = pd.DataFrame(class_list, columns=['atom'])
        class_1_atoms_torsions = []
        for i in class_1_list_torsions:
            class_1_atoms_torsions.append(df_atoms['atom'].iloc[i])
        class_2_atoms_torsions = []
        for i in class_2_list_torsions:
            class_2_atoms_torsions.append(df_atoms['atom'].iloc[i])
        class_3_atoms_torsions = []
        for i in class_3_list_torsions:
            class_3_atoms_torsions.append(df_atoms['atom'].iloc[i])
        class_4_atoms_torsions = []
        for i in class_4_list_torsions:
            class_4_atoms_torsions.append(df_atoms['atom'].iloc[i])
        data_tuples = list(zip(class_1_atoms_torsions, class_2_atoms_torsions, class_3_atoms_torsions, class_4_atoms_torsions))
        df_tor_atoms = pd.DataFrame(data_tuples, columns=['class1_atoms','class2_atoms','class3_atoms','class4_atoms'])
        #print(df_tor_atoms.head())
        frames = [df_tor, df_tor_atoms]  
        df_concat = pd.concat(frames, axis = 1)   
        #print(df_concat.tail())
        k_list_torsions_off = df_concat["k"].values.tolist()
        k_list_torsions_off = [round(num, 6) for num in k_list_torsions_off]
        #print(k_list_torsions_off)
        #print(len(k_list_torsions_off))
        class1_torsion_list_off = df_concat["class1_atoms"].values.tolist()
        #print(class1_torsion_list_off)
        #print(len(class1_torsion_list_off))
        class2_torsion_list_off = df_concat["class2_atoms"].values.tolist()
        #print(class2_torsion_list_off)
        #print(len(class2_torsion_list_off))
        class3_torsion_list_off = df_concat["class3_atoms"].values.tolist()
        #print(class3_torsion_list_off)
        #print(len(class3_torsion_list_off))
        class4_torsion_list_off = df_concat["class4_atoms"].values.tolist()
        #print(class4_torsion_list_off)
        #print(len(class4_torsion_list_off))
        periodicity_torsion_list_off = df_concat["periodicity"].values.tolist()
        #print(periodicity_torsion_list_off)
        #print(len(periodicity_torsion_list_off))
        phase_torsion_list_off = df_concat["phase"].values.tolist()
        phase_torsion_list_off = [round(num, 6) for num in phase_torsion_list_off]
        #print(phase_torsion_list_off)
        #print(len(phase_torsion_list_off))
        #NonBond Parameters
        xml_off = open(self.host_xml, 'r') 
        xml_off_lines = xml_off.readlines() 
        for i in range(len(xml_off_lines)):
            if "<GlobalParameters/>" in xml_off_lines[i]:
                to_begin = int(i)
            if "<Exceptions>" in xml_off_lines[i]:
                to_end = int(i)  
        nonbond_params = xml_off_lines[to_begin + 4 : to_end - 1]
        for i in range(len(nonbond_params)):
            nonbond_params[i] = nonbond_params[i].strip()
            nonbond_params[i] = nonbond_params[i].replace("sig","sigma")
            nonbond_params[i] = nonbond_params[i].replace("q=","charge=")
            nonbond_params[i] = nonbond_params[i].replace("eps","epsilon") 
        epsilon = []
        for i in range(len(nonbond_params)):
            epsilon.append(float(re.findall('[-+]?\d*\.\d+|\d+',nonbond_params[i])[0]))
        epsilon = [round(num, 10) for num in epsilon]
        #print(epsilon)
        #print(len(epsilon))
        charge = []
        for i in range(len(nonbond_params)):
            charge.append(float(re.findall('[-+]?\d*\.\d+|\d+',nonbond_params[i])[1]))
        charge = [round(num, 10) for num in charge]
        #print(charge)
        #print(len(charge))
        sigma = []
        for i in range(len(nonbond_params)):
            sigma.append(float(re.findall('[-+]?\d*\.\d+|\d+',nonbond_params[i])[2]))
        sigma = [round(num, 10) for num in sigma]
        #print(sigma)
        #print(len(sigma))
        #Create connect information   
        ppdb = PandasPdb()
        ppdb.read_pdb(self.host_singular_file)
        no_atoms = ppdb.df["ATOM"].shape[0]
        class1_bond_list_connect_off = [int(i) for i in class1_bond_list_off]
        class2_bond_list_connect_off = [int(i) for i in class2_bond_list_off]
        bond_list_list = []
        for i in range(len(class1_bond_list_connect_off)):
            args = (class1_bond_list_connect_off[i],class2_bond_list_connect_off[i])
            bond_list_list.append(args)
        #print(bond_list_list)
        list_list_connect = []
        for i in range(no_atoms):
            atom_index = i + 1
            list_connect = []
            for j in bond_list_list:
                if j[0] == atom_index or j[1] == atom_index:
                    list_connect.append(j)
            list_list_connect.append(list_connect)
        #print(list_list_connect)
        connect_list = []
        for i in list_list_connect:
            if len(i) > 1:
                connect_list.append(i)        
        #print(connect_list)
        connect_list_merged = []
        for i in range(len(connect_list)):
            merged = [k for z in connect_list[i] for k in z]
            connect_list_merged.append(merged)
        #print(connect_list_merged)
        connect_list_merged_sorted = []
        for i in range(len(connect_list_merged)):
            sorted_list = sorted(connect_list_merged[i],key = connect_list_merged[i].count,reverse = True)
            connect_list_merged_sorted.append(sorted_list)
        #print(connect_list_merged_sorted)
        connect_list_merged_sorted_unique = []
        for i in range(len(connect_list_merged_sorted)):
            unique_list = uniq(connect_list_merged_sorted[i])
            connect_list_merged_sorted_unique.append(unique_list)
        #print(connect_list_merged_sorted_unique)
        len_list = []
        for i in range(len(connect_list_merged_sorted_unique)):
            len_list.append(len(connect_list_merged_sorted_unique[i]))
        max_len = max(len_list)    
        for i in range(len(connect_list_merged_sorted_unique)):
            if len(connect_list_merged_sorted_unique[i]) < max_len:
                diff = max_len - len(connect_list_merged_sorted_unique[i])
                for j in range(diff):
                    connect_list_merged_sorted_unique[i].append("X")            
        #print(connect_list_merged_sorted_unique)
        f = open(self.conect_pdb_txt, 'w')
        for i in range(len(connect_list_merged_sorted_unique)):
            string_line = '    '.join([str(elem) for elem in connect_list_merged_sorted_unique[i]]) 
            f.write("CONECT" + "    " + string_line + "\n")
        f.close()
        conect = open(self.conect_pdb_txt, 'r')
        conect_lines = conect.readlines()
        for i in range(len(conect_lines)): 
            words = conect_lines[i].split()
            if len(words) == 6 :
                conect_lines[i] = '{:>0} {:>4} {:>4} {:>4} {:>4} {:>4}'.format(*words)  
            if len(words) == 5 :
                conect_lines[i] = '{:>0} {:>4} {:>4} {:>4} {:>4}'.format(*words)  
            if len(words) == 4 :
                conect_lines[i] = '{:>0} {:>4} {:>4} {:>4}'.format(*words)  
            if len(words) == 3 :
                conect_lines[i] = '{:>0} {:>4} {:>4}'.format(*words)  
            if len(words) == 2:
                conect_lines[i] = '{:>0} {:>4}'.format(*words)  
        f = open(self.conect_pdb_txt, 'w')
        for i in range(len(conect_lines)):
            f.write(conect_lines[i] + "\n")
        f.close()
        conect = open(self.conect_pdb_txt, 'r')
        conect_lines = conect.readlines()
        for i in range(len(conect_lines)): 
            conect_lines[i] = conect_lines[i].replace("X", "")
        f = open(self.conect_pdb_txt, 'w')
        for i in range(len(conect_lines)):
            f.write(conect_lines[i])
        f.close()
        readFile = open(self.host_singular_file)
        pdblines = readFile.readlines()
        readFile.close()
        w = open(self.conect_pdb_file,'w')
        w.writelines([item for item in pdblines[:-1]])
        w.close()
        conect = open(self.conect_pdb_txt, 'r')
        conect_lines = conect.readlines()
        with open (self.conect_pdb_file, "a") as f:
            for i in conect_lines:
                f.write(i)
        with open (self.conect_pdb_file, "a") as f:
            f.write("END")    
        #Exceptions
        xml_off = open(self.host_xml, 'r') 
        xml_off_lines = xml_off.readlines() 
        for i in range(len(xml_off_lines)):
            if "<Exceptions>" in xml_off_lines[i]:
                to_begin = int(i)
            if "</Exceptions>" in xml_off_lines[i]:
                to_end = int(i)  
        exception_params = xml_off_lines[to_begin + 1 : to_end]
        for i in range(len(exception_params)):
            exception_params[i] = exception_params[i].strip()
            exception_params[i] = exception_params[i].replace("p1","class1") 
            exception_params[i] = exception_params[i].replace("p2","class2") 
            exception_params[i] = exception_params[i].replace("sig","sigma") 
            exception_params[i] = exception_params[i].replace("eps","epsilon")
            exception_params[i] = exception_params[i].replace("q=","charge=")
        epsilon_exceptions = []
        for i in range(len(exception_params)):
            epsilon_exceptions.append(float(re.findall('[-+]?\d*\.\d+|\d+',exception_params[i])[0]))
        epsilon_exceptions = [round(num, 10) for num in epsilon_exceptions]
        #print(epsilon_exceptions)
        #print(len(epsilon_exceptions))
        charge_exceptions = []
        for i in range(len(exception_params)):
            charge_exceptions.append(float(re.findall('[-+]?\d*\.\d+|\d+',exception_params[i])[5]))
        charge_exceptions = [round(num, 10) for num in charge_exceptions]
        #print(charge_exceptions)
        #print(len(charge_exceptions))
        sigma_exceptions = []
        for i in range(len(exception_params)):
            sigma_exceptions.append(float(re.findall('[-+]?\d*\.\d+|\d+',exception_params[i])[6]))
        sigma_exceptions = [round(num, 10) for num in sigma_exceptions]
        #print(sigma_exceptions)
        #print(len(sigma_exceptions))
        class1_exceptions = []
        for i in range(len(exception_params)):
            class1_exceptions.append(int(re.findall('\d*\.?\d+',exception_params[i])[2]))
        #print(class1_exceptions)
        class2_exceptions = []
        for i in range(len(exception_params)):
            class2_exceptions.append(int(re.findall('\d*\.?\d+',exception_params[i])[4]))
        #print(class2_exceptions)
        data_tuples = list(zip(epsilon_exceptions, charge_exceptions, sigma_exceptions, class1_exceptions, class2_exceptions))
        df_exception = pd.DataFrame(data_tuples, columns=['epsilon','charge','sigma','class1','class2'])
        #print(df_exception.head())
        class_1_list_exceptions = df_exception["class1"].to_list()
        class_2_list_exceptions = df_exception["class2"].to_list()
        df_atoms = pd.DataFrame(class_list, columns=['atom'])
        class_1_atoms_exceptions = []
        for i in class_1_list_exceptions:
            class_1_atoms_exceptions.append(df_atoms['atom'].iloc[i])
        class_2_atoms_exceptions = []
        for i in class_2_list_exceptions:
            class_2_atoms_exceptions.append(df_atoms['atom'].iloc[i])
        data_tuples = list(zip(class_1_atoms_exceptions, class_2_atoms_exceptions))
        df_exceptions_atoms = pd.DataFrame(data_tuples, columns=['class1_atoms','class2_atoms'])
        #print(df_exceptions_atoms.head())
        frames = [df_exception, df_exceptions_atoms]  
        df_concat = pd.concat(frames, axis = 1)   
        #print(df_concat.tail())
        class1_exception_list_off = df_concat["class1_atoms"].values.tolist()
        #print(class1_exception_list_off)
        #print(len(class1_exception_list_off))
        class2_exception_list_off = df_concat["class2_atoms"].values.tolist()
        #print(class2_exception_list_off)
        #print(len(class2_exception_list_off))
        epsilon_exception_list_off = df_concat["epsilon"].values.tolist()
        epsilon_exception_list_off = [round(num, 6) for num in epsilon_exception_list_off]
        #print(epsilon_exception_list_off)
        #print(len(epsilon_exception_list_off))
        charge_exception_list_off = df_concat["charge"].values.tolist()
        charge_exception_list_off = [round(num, 6) for num in charge_exception_list_off]
        #print(charge_exception_list_off)
        #print(len(charge_exception_list_off))
        sigma_exception_list_off = df_concat["sigma"].values.tolist()
        sigma_exception_list_off = [round(num, 6) for num in sigma_exception_list_off]
        #print(sigma_exception_list_off)
        #print(len(sigma_exception_list_off))
        xml = open(self.host_singular_xml_file, "w")
        xml.write("<?xml version=" + '"' + "1.0" + '"' + " " + "?>" "\n")
        xml.write("<ForceField>" + "\n")
        xml.write("<AtomTypes>" + "\n")
        for i in range(len(list_atom_masses)):
            xml.write("<Type" + " " 
                              + "class=" + '"' + class_list[i] + '"' + " " 
                              + "element=" + '"' + element_symbol[i] + '"' + " " 
                              + "mass=" + '"' + str(list_atom_masses[i]) + '"' + " " 
                              + "name=" + '"' + name_list[i] + '"' 
                              + "/>"  + "\n")
        xml.write("</AtomTypes>" + "\n")
        xml.write("<Residues>" + "\n")
        xml.write("<Residue name=" + '"' + self.host_residue_name + '"' + ">" + "\n")
        for i in range(len(list_atom_masses)):
            xml.write("<Atom" + " " 
                              + "name=" + '"' + class_list[i] + '"' + " " 
                              + "type=" + '"' + name_list[i] + '"' 
                              + "/>"  + "\n")
        for i in range(len(class1_bonds)):
            xml.write("<Bond" + " " 
                              + "from=" + '"' + str(class1_bonds[i]) + '"' + " " 
                              + "to=" + '"' + str(class2_bonds[i]) + '"' 
                              + "/>"  + "\n")
        xml.write("</Residue>" + "\n")
        xml.write("</Residues>" + "\n")
        xml.write("<HarmonicBondForce>" + "\n")
        for i in range(len(bond_list_off)):
            xml.write("<Bond" + " " 
                              + "class1=" + '"' + class1_bond_list_off[i] + '"' + " " 
                              + "class2=" + '"' + class2_bond_list_off[i] + '"' + " " 
                              + "k=" + '"' + str(k_list_bonds_off[i]) + '"' + " " 
                              + "length=" + '"' + str(bond_list_off[i]) + '"' 
                              + "/>"  + "\n")
        xml.write("</HarmonicBondForce>" + "\n")
        xml.write("<HarmonicAngleForce>" + "\n")
        for i in range(len(angle_list_off)):
            xml.write("<Angle" + " " 
                               + "angle=" + '"' + str(angle_list_off[i]) + '"' + " " 
                               + "class1=" + '"' + class1_angle_list_off[i] + '"' + " " 
                               + "class2=" + '"' + class2_angle_list_off[i] + '"' + " " 
                               + "class3=" + '"' + class3_angle_list_off[i] + '"' + " " 
                               + "k=" + '"' + str(k_list_angles_off[i]) + '"' 
                               + "/>"  + "\n")
        xml.write("</HarmonicAngleForce>" + "\n")
        xml.write("<PeriodicTorsionForce>" + "\n")
        for i in range(len(k_list_torsions_off)):
            xml.write("<Proper" + " " 
                                + "k=" + '"' + str(k_list_torsions_off[i]) + '"' + " " 
                                + "class1=" + '"' + str(class1_torsion_list_off[i]) + '"' + " " 
                                + "class2=" + '"' + str(class2_torsion_list_off[i]) + '"' + " " 
                                + "class3=" + '"' + str(class3_torsion_list_off[i]) + '"' + " " 
                                + "class4=" + '"' + str(class4_torsion_list_off[i]) + '"' + " " 
                                + "periodicity=" + '"' + str(periodicity_torsion_list_off[i]) + '"' + " "                                              
                                + "phase=" + '"' + str(phase_torsion_list_off[i]) + '"' 
                                + "/>"  + "\n") 
        xml.write("</PeriodicTorsionForce>" + "\n")
        xml.write("<NonbondedForce" + " " 
                                + "coulomb14scale=" + '"' + self.coulomb14scale + '"' + " " 
                                + "lj14scale=" + '"' + self.lj14scale + '"' 
                                + ">" + "\n")
        for i in range(len(name_list)):
            xml.write("<Atom" + " " 
                                + "charge=" + '"' + str(charge[i]) + '"' + " " 
                                + "epsilon=" + '"' + str(epsilon[i]) + '"' + " " 
                                + "sigma=" + '"' + str(sigma[i]) + '"' + " "   
                                + "type=" + '"' + name_list[i] + '"' 
                                + "/>"  + "\n")   
        xml.write("</NonbondedForce>" + "\n")
        xml.write("</ForceField>" )
        xml.close()  
        
    def run_host_mm(self):
        pdb = simtk.openmm.app.PDBFile(self.conect_pdb_file)
        forcefield = simtk.openmm.app.ForceField(self.host_singular_xml_file)
        system = forcefield.createSystem(pdb.topology)
        integrator = simtk.openmm.LangevinIntegrator(300 * simtk.unit.kelvin, 1 / simtk.unit.picosecond, 0.002 * simtk.unit.picoseconds)
        simulation = simtk.openmm.app.Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy()
        state = simulation.context.getState(getEnergy=True)
        energy = state.getPotentialEnergy()
        print(energy)
        simulation.reporters.append(simtk.openmm.app.PDBReporter(self.sim_output, self.sim_steps/10))
        simulation.reporters.append(simtk.openmm.app.StateDataReporter(stdout, reportInterval = int(self.sim_steps/10), step = True, potentialEnergy = True, temperature = True))
        simulation.step(self.sim_steps)
        command = "rm -rf " + str(self.sim_output)
        os.system(command)
        
    def write_reparameterised_host_xml(self):
        """
        This function replaces the angle, bond and charge parameters in the host forcefield file 
        with the QM generated features.
        """
        f = open(self.host_qm_params_file, 'r')
        lines_params = f.readlines()
        # Bond Parameters
        for i in range(len(lines_params)):
            if "Begin writing the Bond Parameters" in lines_params[i]:
                to_begin = int(i)
            if "Finish writing the Bond Parameters" in lines_params[i]:
                to_end = int(i)  
        bond_params = lines_params[to_begin + 1:to_end]
        index_search_replace_bond = []
        for i in bond_params:
            bond_line_to_replace = i
            #print(bond_line_to_replace)
            atom_number_list = [re.findall('\d*\.?\d+', i)[1], re.findall('\d*\.?\d+', i)[3]]
            #print(atom_number_list)
            comb_1 = "<Bond"  + " " + "class1=" + '"' + atom_number_list[0] + '"' + " " + "class2=" + '"' + atom_number_list[1] + '"'
            comb_2 = "<Bond"  + " " + "class1=" + '"' + atom_number_list[1] + '"' + " " + "class2=" + '"' + atom_number_list[0] + '"'
            comb_list_bond = [comb_1, comb_2]
            #print(comb_list_bond)
            list_search_bond = [search_in_file(file = self.host_singular_xml_file, word = comb_1),search_in_file(file = self.host_singular_xml_file, word = comb_2)]
            #print(list_search_bond)
            for j in range(len(list_search_bond)):
                if list_search_bond[j] != []:
                    to_add = (list_search_bond[j],i)
                    #print(to_add)
                    index_search_replace_bond.append(to_add)
        # Angle Parameters
        for i in range(len(lines_params)):
            if "Begin writing the Angle Parameters" in lines_params[i]:
                to_begin = int(i)
            if "Finish writing the Angle Parameters" in lines_params[i]:
                to_end = int(i)  
        angle_params = lines_params[to_begin + 1:to_end]
        index_search_replace_angle = []
        for i in angle_params:
            angle_line_to_replace = i
            #print(angle_line_to_replace)
            atom_number_list = [re.findall('\d*\.?\d+', i)[2], re.findall('\d*\.?\d+', i)[4], re.findall('\d*\.?\d+', i)[6]]
            #print(atom_number_list)
            comb_1 = "class1=" + '"' + atom_number_list[0] + '"' + " " + "class2=" + '"' + atom_number_list[1] + '"' + " " + "class3=" + '"' + atom_number_list[2] + '"' + " " + "k=" 
            comb_2 = "class1=" + '"' + atom_number_list[0] + '"' + " " + "class2=" + '"' + atom_number_list[2] + '"' + " " + "class3=" + '"' + atom_number_list[1] + '"' + " " + "k="
            comb_3 = "class1=" + '"' + atom_number_list[1] + '"' + " " + "class2=" + '"' + atom_number_list[0] + '"' + " " + "class3=" + '"' + atom_number_list[2] + '"' + " " + "k="
            comb_4 = "class1=" + '"' + atom_number_list[1] + '"' + " " + "class2=" + '"' + atom_number_list[2] + '"' + " " + "class3=" + '"' + atom_number_list[0] + '"' + " " + "k=" 
            comb_5 = "class1=" + '"' + atom_number_list[2] + '"' + " " + "class2=" + '"' + atom_number_list[0] + '"' + " " + "class3=" + '"' + atom_number_list[1] + '"' + " " + "k="
            comb_6 = "class1=" + '"' + atom_number_list[2] + '"' + " " + "class2=" + '"' + atom_number_list[1] + '"' + " " + "class3=" + '"' + atom_number_list[0] + '"' + " " + "k="
            comb_list_angle = [comb_1, comb_2, comb_3, comb_4, comb_5, comb_6]
            #print(comb_list_angle)    
            list_search_angle = [search_in_file(file = self.host_singular_xml_file, word = comb_1), 
                                 search_in_file(file = self.host_singular_xml_file, word = comb_2),
                                 search_in_file(file = self.host_singular_xml_file, word = comb_3),
                                 search_in_file(file = self.host_singular_xml_file, word = comb_4),
                                 search_in_file(file = self.host_singular_xml_file, word = comb_5),
                                 search_in_file(file = self.host_singular_xml_file, word = comb_6)]
            #print(list_search_angle)
            for j in range(len(list_search_angle)):
                if list_search_angle[j] != []:
                    to_add = (list_search_angle[j],i)
                    #print(to_add)
                    index_search_replace_angle.append(to_add)       
        f_org = open(self.host_singular_xml_file)
        lines = f_org.readlines()
        for i in range(len(index_search_replace_bond)):
            line_number = index_search_replace_bond[i][0][0][0] - 1
            line_to_replace = index_search_replace_bond[i][0][0][1]
            line_to_replace_with = index_search_replace_bond[i][1]
            lines[line_number] = line_to_replace_with 
        for i in range(len(index_search_replace_angle)):
            line_number = index_search_replace_angle[i][0][0][0] - 1
            line_to_replace = index_search_replace_angle[i][0][0][1]
            line_to_replace_with = index_search_replace_angle[i][1]
            lines[line_number] = line_to_replace_with 
        f_cop = open(self.reparameterised_host_xml_file, "w")
        for i in lines:
            f_cop.write(i)
        f_cop.close()
        f = open(self.host_qm_params_file, 'r')
        lines_params = f.readlines()
        # Charge Parameters
        for i in range(len(lines_params)):
            if "Begin writing the Charge Parameters" in lines_params[i]:
                to_begin = int(i)
            if "Finish writing the Charge Parameters" in lines_params[i]:
                to_end = int(i)  
        charge_params = lines_params[to_begin + 1:to_end]
        charge_atom_list = []
        for i in charge_params:
            charge_atom_list.append([float(re.findall('[-+]?\d*\.\d+|\d+',i)[0]), int(re.findall('\d*\.?\d+', i)[3])])
        f = open(self.reparameterised_host_xml_file, 'r')
        lines_params = f.readlines()
        to_replace = []
        to_replace_with = []
        for j in lines_params:
            for i in range(len(charge_atom_list)):
                if "<Atom charge" in j and "host_" + str(charge_atom_list[i][1]) in j :
                    #print(j)
                    to_replace.append(j)
                    to_replace_ = re.findall('[-+]?\d*\.\d+|\d+',j)[0]
                    to_replace_with_ = j.replace(to_replace_, str(charge_atom_list[i][0]))
                    to_replace_with.append(to_replace_with_)
                    #print(to_replace_with_)
        line_numbers = []
        for i in to_replace:
            line_number = (search_in_file(file = self.reparameterised_host_xml_file, word = i)[0][0]) - 1
            line_numbers.append(line_number)
        f_in = open(self.reparameterised_host_xml_file)
        lines = f_in.readlines()
        for i in range(len(line_numbers)):
            lines[line_numbers[i]] = to_replace_with[i]
        f_out = open(self.reparams_host_file, "w")
        for i in lines:
            f_out.write(i)
        f_out.close()
        
    def run_host_mm_qm(self):
        pdb = simtk.openmm.app.PDBFile(self.conect_pdb_file)
        forcefield = simtk.openmm.app.ForceField(self.reparams_host_file)
        system = forcefield.createSystem(pdb.topology)
        integrator = simtk.openmm.LangevinIntegrator(300 * simtk.unit.kelvin, 1 / simtk.unit.picosecond, 0.002 * simtk.unit.picoseconds)
        simulation = simtk.openmm.app.Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy()
        state = simulation.context.getState(getEnergy=True)
        energy = state.getPotentialEnergy()
        print(energy)
        simulation.reporters.append(simtk.openmm.app.PDBReporter(self.sim_output, self.sim_steps/10))
        simulation.reporters.append(simtk.openmm.app.StateDataReporter(stdout, reportInterval = int(self.sim_steps/10), step = True, potentialEnergy = True, temperature = True))
        simulation.step(self.sim_steps)
        command = "rm -rf " + str(self.sim_output)
        os.system(command)
####################################################################################################################################################################################
