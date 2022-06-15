"""
Functions that involve manipulations of files.
"""

import os

import modules.torsion_drive_outputs as torsion_outputs

def change_names(inpcrd_file, prmtop_file, pdb_file):
    # TODO: replace with three shutil.copy calls
    command = "cp -r " + inpcrd_file + " system_qmmmrebind.inpcrd"
    os.system(command)
    command = "cp -r " + prmtop_file + " system_qmmmrebind.prmtop"
    os.system(command)
    command = "cp -r " + pdb_file + " system_qmmmrebind.pdb"
    os.system(command)


# TODO: unnecessary function
def copy_file(source, destination):

    """
    Copies a file from a source to the destination.
    """
    shutil.copy(source, destination)

def search_in_file(file: str, word: str) -> list:

    """
    Search for the given string in file and return lines
    containing that string along with line numbers.

    Parameters
    ----------
    file : str
        Input file.

    word : str
        Search word.

    Returns
    -------
    list_of_results : list
        List of lists with each element representing the
        line number and the line contents.

    """
    line_number = 0
    list_of_results = []
    with open(file, "r") as f:
        for line in f:
            line_number += 1
            if word in line:
                list_of_results.append((line_number, line.rstrip()))
    return list_of_results

def remove_mm_files(qm_scan_file):
    """
    Delete all generated PDB files.

    Parameters
    ----------
    qm_scan_file : str
        Output scan file containing torsiondrive scans.

    """
    mm_pdb_list = []
    for i in torsion_outputs.get_dihedrals(qm_scan_file):
        if i > 0:
            pdb_file = "plus_" + str(abs(i)) + ".pdb"
        if i < 0:
            pdb_file = "minus_" + str(abs(i)) + ".pdb"
        mm_pdb_list.append(pdb_file)
    for i in mm_pdb_list:
        command = "rm -rf  " + i
        os.system(command)
        command = "rm -rf  " + i[:-4] + ".inpcrd"
        os.system(command)
        command = "rm -rf  " + i[:-4] + ".prmtop"
        os.system(command)

def move_qmmmmrebind_files(
    prmtopfile="system_qmmmrebind.prmtop",
    inpcrdfile="system_qmmmrebind.inpcrd",
    pdbfile="system_qmmmrebind.pdb",
):

    """
    Moves QMMMReBind generated topology and parameter files
    to a new directory .

    Parameters
    ----------
    prmtopfile: str
       QMMMReBind generated prmtop file.

    inpcrdfile: str
        QMMMReBind generated inpcrd file.

    pdbfile: str
        QMMMReBind generated PDB file.

    """
    current_pwd = os.getcwd()
    command = "rm -rf reparameterized_files"
    os.system(command)
    command = "mkdir reparameterized_files"
    os.system(command)
    shutil.copy(
        current_pwd + "/" + prmtopfile,
        current_pwd + "/" + "reparameterized_files" + "/" + prmtopfile,
    )
    shutil.copy(
        current_pwd + "/" + inpcrdfile,
        current_pwd + "/" + "reparameterized_files" + "/" + inpcrdfile,
    )
    shutil.copy(
        current_pwd + "/" + pdbfile,
        current_pwd + "/" + "reparameterized_files" + "/" + pdbfile,
    )

def move_qm_files():

    """
    Moves QM engine generated files to a new directory .

    """
    current_pwd = os.getcwd()
    command = "rm -rf qm_data"
    os.system(command)
    command = "mkdir qm_data"
    os.system(command)
    command = "cp -r " + "*.com* " + current_pwd + "/" + "qm_data"
    os.system(command)
    command = "cp -r " + "*.log* " + current_pwd + "/" + "qm_data"
    os.system(command)
    command = "cp -r " + "*.chk* " + current_pwd + "/" + "qm_data"
    os.system(command)
    command = "cp -r " + "*.fchk* " + current_pwd + "/" + "qm_data"
    os.system(command)


def move_qmmmrebind_files():

    """
    Moves all QMMMREBind files to a new directory.

    """
    current_pwd = os.getcwd()
    command = "rm -rf qmmmrebind_data"
    os.system(command)
    command = "mkdir qmmmrebind_data"
    os.system(command)
    command = "mv " + "*.sdf* " + current_pwd + "/" + "qmmmrebind_data"
    os.system(command)
    command = "mv " + "*.txt* " + current_pwd + "/" + "qmmmrebind_data"
    os.system(command)
    command = "mv " + "*.pdb* " + current_pwd + "/" + "qmmmrebind_data"
    os.system(command)
    command = "mv " + "*.xml* " + current_pwd + "/" + "qmmmrebind_data"
    os.system(command)
    command = "mv " + "*.chk* " + current_pwd + "/" + "qmmmrebind_data"
    os.system(command)
    command = "mv " + "*.fchk* " + current_pwd + "/" + "qmmmrebind_data"
    os.system(command)
    command = "mv " + "*.com* " + current_pwd + "/" + "qmmmrebind_data"
    os.system(command)
    command = "mv " + "*.log* " + current_pwd + "/" + "qmmmrebind_data"
    os.system(command)
    command = "mv " + "*.inpcrd* " + current_pwd + "/" + "qmmmrebind_data"
    os.system(command)
    command = "mv " + "*.prmtop* " + current_pwd + "/" + "qmmmrebind_data"
    os.system(command)
    command = "mv " + "*.parm7* " + current_pwd + "/" + "qmmmrebind_data"
    os.system(command)
    command = "mv " + "*.out* " + current_pwd + "/" + "qmmmrebind_data"
    os.system(command)
    command = "mv " + "*run_command* " + current_pwd + "/" + "qmmmrebind_data"
    os.system(command)
    command = "mv " + "*.dat* " + current_pwd + "/" + "qmmmrebind_data"
    os.system(command)
    command = "mv " + "*.xyz* " + current_pwd + "/" + "qmmmrebind_data"
    os.system(command)