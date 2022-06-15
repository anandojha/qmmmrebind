"""
If a file needs to be modified after it's generated
"""

import re
import os

from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.toolkit.topology import Molecule, Topology
from biopandas.pdb import PandasPdb
import pandas as pd
import parmed
import simtk

import modules.base as base

def delete_guest_angle_params(guest_qm_params_file="guest_qm_params.txt"):
    """
    
    """
    # TODO: change this function to leave out angle parameters in a more
    #  efficient way - without storing a readlines() in memory
    f_params = open(guest_qm_params_file, "r")
    lines_params = f_params.readlines()
    for i in range(len(lines_params)):
        if "Begin writing the Angle Parameters" in lines_params[i]:
            to_begin = int(i)
        if "Finish writing the Angle Parameters" in lines_params[i]:
            to_end = int(i)
    lines_selected = lines_params[:to_begin] + lines_params[to_end + 1 :]
    with open(guest_qm_params_file, "w") as f_:
        f_.write("".join(lines_selected))
    return


def remove_bad_angle_params(
        guest_qm_params_file="guest_qm_params.txt", angle=1.00, k_angle=500):
    with open(guest_qm_params_file, "r") as f_params:
        lines_params = f_params.readlines()
    for i in range(len(lines_params)):
        if "Begin writing the Angle Parameters" in lines_params[i]:
            to_begin = int(i)
        if "Finish writing the Angle Parameters" in lines_params[i]:
            to_end = int(i)
    angle_params = lines_params[to_begin + 1 : to_end]
    lines_to_omit = []
    for i in angle_params:
        if float(re.findall(r"[-+]?\d+[.]?\d*", i)[0]) < float(angle) or float(
            re.findall(r"[-+]?\d+[.]?\d*", i)[1]
        ) > float(k_angle):
            lines_to_omit.append(i)
    for b in lines_to_omit:
        lines_params.remove(b)
    with open(guest_qm_params_file, "w") as file:
        for j in lines_params:
            file.write(j)
            
def xyz_to_pdb(xyz_file, coords_file, template_pdb, system_pdb):
    """
    Converts a XYZ file to a PDB file.

    Parameters
    ----------
    xyz_file : str
        XYZ file containing the coordinates of the system.

    coords_file : str
        A text file containing the coordinates part of XYZ file.

    template_pdb : str
        A pdb file to be used as a template for the required PDB.

    system_pdb : str
        Output PDB file with the coordinates updated in the
        template pdb using XYZ file.

    """
    with open(xyz_file, "r") as f:
        lines = f.readlines()
    needed_lines = lines[2:]
    with open(coords_file, "w") as f:
        for i in needed_lines:
            f.write(i)
    df = pd.read_csv(coords_file, header=None, delimiter=r"\s+")
    df.columns = ["atom", "x", "y", "z"]
    ppdb = PandasPdb()
    ppdb.read_pdb(template_pdb)
    ppdb.df["ATOM"]["x_coord"] = df["x"]
    ppdb.df["ATOM"]["y_coord"] = df["y"]
    ppdb.df["ATOM"]["z_coord"] = df["z"]
    ppdb.to_pdb(system_pdb)

def generate_xml_from_pdb_sdf(system_pdb, system_sdf, system_xml):
    """
    Generates an openforcefield xml file from the pdb file.

    Parameters
    ----------
    system_pdb : str
        Input PDB file.

    system_sdf : str
        SDF file of the system.

    system_xml : str
        XML force field file generated using PDB and SDF files.

    """
    # command = "babel -ipdb " + system_pdb + " -osdf " + system_sdf
    command = "obabel -ipdb " + system_pdb + " -osdf -O " + system_sdf
    os.system(command)
    # off_molecule = openforcefield.topology.Molecule(system_sdf)
    off_molecule = Molecule(system_sdf)
    # force_field = openforcefield.typing.engines.smirnoff.ForceField("openff_unconstrained-1.0.0.offxml")
    force_field = ForceField("openff_unconstrained-1.0.0.offxml")
    system = force_field.create_openmm_system(off_molecule.to_topology())
    pdbfile = simtk.openmm.app.PDBFile(system_pdb)
    structure = parmed.openmm.load_topology(
        pdbfile.topology, system, xyz=pdbfile.positions
    )
    with open(system_xml, "w") as f:
        f.write(simtk.openmm.XmlSerializer.serialize(system))
        
def generate_xml_from_charged_pdb_sdf(
    system_pdb,
    system_init_sdf,
    system_sdf,
    num_charge_atoms,
    index_charge_atom_1,
    charge_atom_1,
    system_xml,
):
    """
    Generates an openforcefield xml file from the pdb
    file via SDF file and openforcefield.

    Parameters
    ----------
    system_pdb : str
        Input PDB file.

    system_init_sdf : str
        SDF file for the system excluding charge information.

    system_sdf : str
        SDF file of the system.

    num_charge_atoms : int
        Total number of charged atoms in the PDB.

    index_charge_atom_1 : int
        Index of the first charged atom.

    charge_atom_1 : float
        Charge on first charged atom.
    system_xml : str
        XML force field file generated using PDB and SDF files.

    """
    # command = "babel -ipdb " + system_pdb + " -osdf " + system_init_sdf
    command = "obabel -ipdb " + system_pdb + " -osdf -O " + system_init_sdf
    os.system(command)
    with open(system_init_sdf, "r") as f1:
        filedata = f1.readlines()
        filedata = filedata[:-2]
    with open(system_sdf, "w+") as out:
        for i in filedata:
            out.write(i)
        line_1 = (
            "M  CHG  "
            + str(num_charge_atoms)
            + "   "
            + str(index_charge_atom_1)
            + "   "
            + str(charge_atom_1)
            + "\n"
        )
        line_2 = "M  END" + "\n"
        line_3 = "$$$$"
        out.write(line_1)
        out.write(line_2)
        out.write(line_3)
    # off_molecule = openforcefield.topology.Molecule(system_sdf)
    off_molecule = Molecule(system_sdf)
    # force_field = openforcefield.typing.engines.smirnoff.ForceField("openff_unconstrained-1.0.0.offxml")
    force_field = ForceField("openff_unconstrained-1.0.0.offxml")
    system = force_field.create_openmm_system(off_molecule.to_topology())
    pdbfile = simtk.openmm.app.PDBFile(system_pdb)
    structure = parmed.openmm.load_topology(
        pdbfile.topology, system, xyz=pdbfile.positions
    )
    with open(system_xml, "w") as f:
        f.write(simtk.openmm.XmlSerializer.serialize(system))

# TODO: rename to singular_chainid?
def singular_resid(pdbfile, qmmmrebind_init_file):

    """
    Returns a PDB file with chain ID = A

    Parameters
    ----------
    pdbfile: str
        Input PDB file

    qmmmrebind_init_file: str
        Output PDB file

    """

    ppdb = PandasPdb().read_pdb(pdbfile)
    ppdb.df["HETATM"]["chain_id"] = "A"
    ppdb.df["ATOM"]["chain_id"] = "A"
    ppdb.to_pdb(
        path=qmmmrebind_init_file, records=None, gz=False, append_newline=True
    )

def add_vectors_inpcrd(pdbfile, inpcrdfile):

    """
    Adds periodic box dimensions to the inpcrd file

    Parameters
    ----------
    pdbfile: str
       PDB file containing the periodic box information.

    inpcrdfile: str
       Input coordinate file.

    """

    pdbfilelines = open(pdbfile, "r").readlines()
    for i in pdbfilelines:
        if "CRYST" in i:
            vector_list = re.findall(r"[-+]?\d*\.\d+|\d+", i)
            vector_list = [float(i) for i in vector_list]
            vector_list = vector_list[1 : 1 + 6]
            line_to_add = (
                "  "
                + base.truncate(vector_list[0])
                + "  "
                + base.truncate(vector_list[1])
                + "  "
                + base.truncate(vector_list[2])
                + "  "
                + base.truncate(vector_list[3])
                + "  "
                + base.truncate(vector_list[4])
                + "  "
                + base.truncate(vector_list[5])
            )
            print(line_to_add)
    with open(inpcrdfile, "a+") as f:
        f.write(line_to_add)


def add_dim_prmtop(pdbfile, prmtopfile):

    """
    Adds periodic box dimensions flag in the prmtop file.

    Parameters
    ----------
    prmtopfile: str
       Input prmtop file.

    pdbfile: str
       PDB file containing the periodic box information.

    """
    pdbfilelines = open(pdbfile, "r").readlines()
    for i in pdbfilelines:
        if "CRYST" in i:
            vector_list = re.findall(r"[-+]?\d*\.\d+|\d+", i)
            vector_list = [float(i) for i in vector_list]
            vector_list = vector_list[1 : 1 + 6]
            vector_list = [i / 10 for i in vector_list]
            vector_list = [base.truncate(i) for i in vector_list]
            vector_list = [i + "E+01" for i in vector_list]
            line3 = (
                "  "
                + vector_list[3]
                + "  "
                + vector_list[0]
                + "  "
                + vector_list[1]
                + "  "
                + vector_list[2]
            )
            print(line3)
    line1 = "%FLAG BOX_DIMENSIONS"
    line2 = "%FORMAT(5E16.8)"
    with open(prmtopfile) as f1, open("intermediate.prmtop", "w") as f2:
        for line in f1:
            if line.startswith("%FLAG RADIUS_SET"):
                line = line1 + "\n" + line2 + "\n" + line3 + "\n" + line
            f2.write(line)
    command = "rm -rf " + prmtopfile
    os.system(command)
    command = "mv  intermediate.prmtop " + prmtopfile
    os.system(command)


def add_period_prmtop(parm_file, ifbox):

    """
    Changes the value of IFBOX if needed for the prmtop / parm file.
    Set to 1 if standard periodic box and 2 when truncated octahedral.
    """
    with open(parm_file) as f:
        parm_lines = f.readlines()
    lines_contain = []
    for i in range(len(parm_lines)):
        if parm_lines[i].startswith("%FLAG POINTERS"):
            lines_contain.append(i + 4)
    line = parm_lines[lines_contain[0]]
    line_new = "%8s  %6s  %6s  %6s  %6s  %6s  %6s  %6s  %6s  %6s" % (
        re.findall(r"\d+", line)[0],
        re.findall(r"\d+", line)[1],
        re.findall(r"\d+", line)[2],
        re.findall(r"\d+", line)[3],
        re.findall(r"\d+", line)[4],
        re.findall(r"\d+", line)[5],
        re.findall(r"\d+", line)[6],
        str(ifbox),
        re.findall(r"\d+", line)[8],
        re.findall(r"\d+", line)[9],
    )
    parm_lines[lines_contain[0]] = line_new + "\n"
    with open(parm_file, "w") as f:
        for i in parm_lines:
            f.write(i)

def add_solvent_pointers_prmtop(non_reparams_file, reparams_file):

    """
    Adds the flag solvent pointers to the topology file.
    """
    f_non_params = open(non_reparams_file, "r")
    lines_non_params = f_non_params.readlines()
    for i in range(len(lines_non_params)):
        if "FLAG SOLVENT_POINTERS" in lines_non_params[i]:
            to_begin = int(i)
    solvent_pointers = lines_non_params[to_begin : to_begin + 3]
    file = open(reparams_file, "a") 
    for i in solvent_pointers:
        file.write(i)

def prmtop_calibration(
    prmtopfile="system_qmmmrebind.prmtop",
    inpcrdfile="system_qmmmrebind.inpcrd",
):

    """
    Standardizes the topology files

    Parameters
    ----------

    prmtopfile: str
       Input prmtop file.

    inpcrdfile: str
       Input coordinate file.

    """
    parm = parmed.load_file(prmtopfile, inpcrdfile)
    parm_1 = parmed.tools.actions.changeRadii(parm, "mbondi3")
    parm_1.execute()
    parm_2 = parmed.tools.actions.setMolecules(parm)
    parm_2.execute()
    parm.save(prmtopfile, overwrite=True)

