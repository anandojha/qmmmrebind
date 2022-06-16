"""

"""

import os
import re

from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.toolkit.topology import Molecule, Topology
from biopandas.pdb import PandasPdb
import pandas as pd
import parmed
import simtk

import modules.constants as const
import modules.base as base
import modules.file_utilities as file_utilities

class GuestAmberXMLAmber:

    """
    A class used to generate a template force field XML file for the ligand
    in order regenerate the reparameterised forcefield XML file.

    This class contain methods to generate a template XML force field through
    openforcefield. XML template generation can be obtained through different
    file formats such as PDB, SDF, and SMI. Methods support charged ligands as
    well. Re-parameterized XML force field files are then generated from the
    template files. Different energy components such as the bond, angle,
    torsional and non-bonded energies are computed for the non-reparametrized
    and the reparameterized force fields. Difference between the
    non-reparameterized and reparameterized force field energies can then be
    analyzed.
    ...

    Attributes
    ----------
    charge : int
        Charge of the ligand.

    num_charge_atoms: int, optional
        Number of charged atoms in the molecule.

    charge_atom_1: int, optional
        Charge on the first charged atom.

    index_charge_atom_1: int, optional
        Index of the first charged atom.

    system_pdb: str, optional
        Ligand PDB file with atom numbers beginning from 1.

    system_mol2: str, optional
        Ligand Mol2 file obtained from PDB file.

    system_in: str, optional
        Prepi file as required by antechamber.

    system_frcmod: str, optional
        FRCMOD file as required by antechamber.

    prmtop_system : str, optional
        Topology file obtained from the ligand PDB.

    inpcrd_system : str, optional
        Coordinate file obtained from the ligand PDB using the
        command saveamberparm.

    system_leap : str, optional
        Amber generated leap file for generating and saving topology
        and coordinate files.

    system_xml: str, optional
        Serialized XML force field file of the ligand.

    system_smi: str, optional
        Ligand SMILES format file.

    system_sdf: str, optional
        Ligand SDF (structure-data) format file.

    system_init_sdf: str, optional
        Ligand SDF (structure-data) format file. This file will be
        generated only if the ligand is charged.

    index_charge_atom_2: int, optional
        Index of the second charged atom of the ligand.

    charge_atom_2: int, optional
        Charge on the second charged atom of the ligand.

    charge_parameter_file: str, optional
        File containing the charges of ligand atoms and their corresponding
        atoms.

    system_qm_pdb: str, optional
        Ligand PDB file with atom numbers beginning from 1.

    bond_parameter_file: str, optional
        Text file containing the bond parameters for the ligand.

    angle_parameter_file: str, optional
        Text file containing the angle parameters of the ligand.

    system_qm_params_file: str, optional
        A text file containing the QM obtained parameters for the
        ligand.

    reparameterised_intermediate_system_xml_file: str, optional
        XML foce field file with bond and angle parameter lines replaced by
        corresponding values obtained from the QM calculations.

    system_xml_non_bonded_file: str, optional
        Text file to write the NonBondedForce Charge Parameters from
        the non-parameterised system XML file.

    system_xml_non_bonded_reparams_file: str, optional
        Text file containing the non-bonded parameters parsed from the
        XML force field file.

    reparameterised_system_xml_file: str, optional
        Reparameterized force field XML file obtained using
        openforcefield.

    non_reparameterised_system_xml_file: str, optional
        Non-reparameterized force field XML file obtained using
        openforcefield.

    prmtop_system_non_params: str, optional
        Amber generated topology file saved from the non-reparameterized
        force field XML file for the ligand.

    inpcrd_system_non_params: str, optional
        Amber generated coordinate file saved from the non-reparameterized
        force field XML file for the ligand.

    prmtop_system_params: str, optional
        Amber generated topology file saved from the reparameterized
        force field XML file for the ligand.

    inpcrd_system_params: str, optional
        Amber generated coordinate file saved from the reparameterized
        force field XML file for the ligand.

    load_topology: str, optional
        Argument to specify how to load the topology. Can either be "openmm"
        or "parmed".

    """

    def __init__(
        self,
        charge=0,
        # TODO: some of these variables are ints, and shouldn't be initialized as strings
        num_charge_atoms="",
        charge_atom_1="",
        index_charge_atom_1="",
        system_pdb="guest_init_II.pdb",
        system_mol2="guest.mol2",
        system_in="guest.in",
        system_frcmod="guest.frcmod",
        prmtop_system="guest.prmtop",
        inpcrd_system="guest.inpcrd",
        system_leap="guest.leap",
        system_xml="guest_init.xml",
        system_smi="guest.smi",
        system_sdf="guest.sdf",
        system_init_sdf="guest_init.sdf",
        index_charge_atom_2=" ",
        charge_atom_2=" ",
        charge_parameter_file="guest_qm_surround_charges.txt",
        system_qm_pdb="guest_init_II.pdb",
        bond_parameter_file="guest_bonds.txt",
        angle_parameter_file="guest_angles.txt",
        system_qm_params_file="guest_qm_params.txt",
        reparameterised_intermediate_system_xml_file="guest_intermediate_reparameterised.xml",
        system_xml_non_bonded_file="guest_xml_non_bonded.txt",
        system_xml_non_bonded_reparams_file="guest_xml_non_bonded_reparams.txt",
        reparameterised_system_xml_file="guest_reparameterised.xml",
        non_reparameterised_system_xml_file="guest_init.xml",
        prmtop_system_non_params="guest_non_params.prmtop",
        inpcrd_system_non_params="guest_non_params.inpcrd",
        prmtop_system_params="guest_params.prmtop",
        inpcrd_system_params="guest_params.inpcrd",
        load_topology="openmm",
    ):

        self.charge = charge
        self.num_charge_atoms = num_charge_atoms
        self.charge_atom_1 = charge_atom_1
        self.index_charge_atom_1 = index_charge_atom_1
        self.system_pdb = system_pdb
        self.system_mol2 = system_mol2
        self.system_in = system_in
        self.system_frcmod = system_frcmod
        self.prmtop_system = prmtop_system
        self.inpcrd_system = inpcrd_system
        self.system_leap = system_leap
        self.system_xml = system_xml
        self.system_smi = system_smi
        self.system_sdf = system_sdf
        self.system_init_sdf = system_init_sdf
        self.index_charge_atom_2 = index_charge_atom_2
        self.charge_atom_2 = charge_atom_2
        self.charge_parameter_file = charge_parameter_file
        self.system_qm_pdb = system_qm_pdb
        self.bond_parameter_file = bond_parameter_file
        self.angle_parameter_file = angle_parameter_file
        self.system_qm_params_file = system_qm_params_file
        self.reparameterised_intermediate_system_xml_file = (
            reparameterised_intermediate_system_xml_file
        )
        self.system_xml_non_bonded_file = system_xml_non_bonded_file
        self.system_xml_non_bonded_reparams_file = (
            system_xml_non_bonded_reparams_file
        )
        self.reparameterised_system_xml_file = reparameterised_system_xml_file
        self.non_reparameterised_system_xml_file = (
            non_reparameterised_system_xml_file
        )
        self.prmtop_system_non_params = prmtop_system_non_params
        self.inpcrd_system_non_params = inpcrd_system_non_params
        self.prmtop_system_params = prmtop_system_params
        self.inpcrd_system_params = inpcrd_system_params
        self.load_topology = load_topology

    def generate_xml_antechamber(self):
        """
        Generates an XML forcefield file from the PDB file through antechamber.
        """
        command = (
            # "babel -ipdb " + self.system_pdb + " -omol2 " + self.system_mol2
            "obabel -ipdb "
            + self.system_pdb
            + " -omol2 -O "
            + self.system_mol2
        )
        os.system(command)
        command = (
            "antechamber -i "
            + self.system_mol2
            + " -fi mol2 -o "
            + self.system_in
            + " -fo prepi -c bcc -nc "
            + str(self.charge)
        )
        os.system(command)
        command = (
            "parmchk2 -i "
            + self.system_in
            + " -o "
            + self.system_frcmod
            + " -f prepi -a Y"
        )
        os.system(command)
        os.system(
            "rm -rf ANTECHAMBER* leap.log sqm* ATOMTYPE.INF PREP.INF NEWPDB.PDB"
        )
        line_1 = "loadamberprep " + self.system_in
        line_2 = "loadamberparams " + self.system_frcmod
        line_3 = "pdb = loadpdb " + self.system_pdb
        line_4 = (
            "saveamberparm pdb "
            + self.prmtop_system
            + " "
            + self.inpcrd_system
        )
        line_5 = "quit"
        with open(self.system_leap, "w") as f:
            f.write("    " + "\n")
            f.write(line_1 + "\n")
            f.write(line_2 + "\n")
            f.write(line_3 + "\n")
            f.write(line_4 + "\n")
            f.write(line_5 + "\n")
        command = "tleap -f " + self.system_leap
        os.system(command)
        parm = parmed.load_file(self.prmtop_system, self.inpcrd_system)
        system = parm.createSystem()
        with open(self.system_xml, "w") as f:
            f.write(simtk.openmm.XmlSerializer.serialize(system))

    def generate_xml_from_pdb_smi(self):
        """
        Generates an XML forcefield file from the SMILES file through
        openforcefield.
        """
        # off_molecule = openforcefield.topology.Molecule(self.system_smi)
        off_molecule = Molecule(self.system_smi)
        # force_field = openforcefield.typing.engines.smirnoff.ForceField("openff_unconstrained-1.0.0.offxml")
        force_field = ForceField("openff_unconstrained-1.0.0.offxml")
        system = force_field.create_openmm_system(off_molecule.to_topology())
        pdbfile = simtk.openmm.app.PDBFile(self.system_pdb)
        structure = parmed.openmm.load_topology(
            pdbfile.topology, system, xyz=pdbfile.positions
        )
        with open(self.system_xml, "w") as f:
            f.write(simtk.openmm.XmlSerializer.serialize(system))

    def generate_xml_from_pdb_sdf(self):
        """
        Generates an XML forcefield file from the SDF file through
        openforcefield.
        """
        command = (
            # "babel -ipdb " + self.system_pdb + " -osdf " + self.system_sdf
            "obabel -ipdb "
            + self.system_pdb
            + " -osdf -O "
            + self.system_sdf
        )
        os.system(command)
        # off_molecule = openforcefield.topology.Molecule(self.system_sdf)
        off_molecule = Molecule(self.system_sdf)
        # force_field = openforcefield.typing.engines.smirnoff.ForceField("openff_unconstrained-1.0.0.offxml")
        force_field = ForceField("openff_unconstrained-1.0.0.offxml")
        system = force_field.create_openmm_system(off_molecule.to_topology())
        pdbfile = simtk.openmm.app.PDBFile(self.system_pdb)
        structure = parmed.openmm.load_topology(
            pdbfile.topology, system, xyz=pdbfile.positions
        )
        with open(self.system_xml, "w") as f:
            f.write(simtk.openmm.XmlSerializer.serialize(system))

    def generate_xml_from_charged_pdb_sdf(self):
        """
        Generates an XML forcefield file for a singly charged ligand molecule
        from the SDF file through openforcefield.
        """
        command = (
            # "babel -ipdb " + self.system_pdb + " -osdf " + self.system_init_sdf
            "obabel -ipdb "
            + self.system_pdb
            + " -osdf -O "
            + self.system_init_sdf
        )
        os.system(command)
        with open(self.system_init_sdf, "r") as f1:
            filedata = f1.readlines()
            filedata = filedata[:-2]
        with open(self.system_sdf, "w+") as out:
            for i in filedata:
                out.write(i)
            line_1 = (
                "M  CHG  "
                + str(self.num_charge_atoms)
                + "   "
                + str(self.index_charge_atom_1)
                + "   "
                + str(self.charge_atom_1)
                + "\n"
            )
            line_2 = "M  END" + "\n"
            line_3 = "$$$$"
            out.write(line_1)
            out.write(line_2)
            out.write(line_3)
        # off_molecule = openforcefield.topology.Molecule(self.system_sdf)
        off_molecule = Molecule(self.system_sdf)
        # force_field = openforcefield.typing.engines.smirnoff.ForceField("openff_unconstrained-1.0.0.offxml")
        force_field = ForceField("openff_unconstrained-1.0.0.offxml")
        system = force_field.create_openmm_system(off_molecule.to_topology())
        pdbfile = simtk.openmm.app.PDBFile(self.system_pdb)
        structure = parmed.openmm.load_topology(
            pdbfile.topology, system, xyz=pdbfile.positions
        )
        with open(self.system_xml, "w") as f:
            f.write(simtk.openmm.XmlSerializer.serialize(system))

    def generate_xml_from_doubly_charged_pdb_sdf(self):
        """
        Generates an XML forcefield file for a singly charged ligand molecule
        from the SDF file through openforcefield.
        """
        command = (
            # "babel -ipdb " + self.system_pdb + " -osdf " + self.system_init_sdf
            "obabel -ipdb "
            + self.system_pdb
            + " -osdf -O "
            + self.system_init_sdf
        )
        os.system(command)
        with open(self.system_init_sdf, "r") as f1:
            filedata = f1.readlines()
            filedata = filedata[:-2]
        with open(self.system_sdf, "w+") as out:
            for i in filedata:
                out.write(i)
            line_1 = (
                "M  CHG  "
                + str(self.num_charge_atoms)
                + "   "
                + str(self.index_charge_atom_1)
                + "   "
                + str(self.charge_atom_1)
                + "   "
                + str(self.index_charge_atom_2)
                + "   "
                + str(self.charge_atom_2)
                + "\n"
            )
            line_2 = "M  END" + "\n"
            line_3 = "$$$$"
            out.write(line_1)
            out.write(line_2)
            out.write(line_3)
        # off_molecule = openforcefield.topology.Molecule(self.system_sdf)
        off_molecule = Molecule(self.system_sdf)
        # force_field = openforcefield.typing.engines.smirnoff.ForceField("openff_unconstrained-1.0.0.offxml")
        force_field = ForceField("openff_unconstrained-1.0.0.offxml")
        system = force_field.create_openmm_system(off_molecule.to_topology())
        pdbfile = simtk.openmm.app.PDBFile(self.system_pdb)
        structure = parmed.openmm.load_topology(
            pdbfile.topology, system, xyz=pdbfile.positions
        )
        with open(self.system_xml, "w") as f:
            f.write(simtk.openmm.XmlSerializer.serialize(system))

    def write_system_params(self):
        """
        Saves the parameters obtained from the QM log files in a text file.
        """
        # Charges from QM files
        df_charges = pd.read_csv(
            self.charge_parameter_file, header=None, delimiter=r"\s+"
        )
        df_charges.columns = ["atom", "charges"]
        qm_charges = df_charges["charges"].values.tolist()
        qm_charges = [round(num, 6) for num in qm_charges]
        # print(qm_charges)
        # Bond Parameters from QM files
        ppdb = PandasPdb()
        ppdb.read_pdb(self.system_qm_pdb)
        atom_name_list = ppdb.df["ATOM"]["atom_number"].values.tolist()
        atom_name_list = [i - 1 for i in atom_name_list]
        # print(atom_name_list)
        df = pd.read_csv(
            self.bond_parameter_file, header=None, delimiter=r"\s+"
        )
        df.columns = [
            "bond",
            "k_bond",
            "bond_length",
            "bond_1",
            "bond_2",
        ]
        # print(df.head())
        bond_1_list = df["bond_1"].values.tolist()
        bond_1_list = [x - 1 + min(atom_name_list) for x in bond_1_list]
        bond_2_list = df["bond_2"].values.tolist()
        bond_2_list = [x - 1 + min(atom_name_list) for x in bond_2_list]
        # print(bond_1_list)
        # print(bond_2_list)
        k_bond_list = df["k_bond"].values.tolist()
        #k_bond_list = [
        #    i * 418.40 for i in k_bond_list
        #]  # kcal/mol * A^2 to kJ/mol * nm^2
        
        k_bond_list = [
            i * const.KCAL_MOL_PER_KJ_MOL * const.ANGSTROMS_PER_NM**2 for i in k_bond_list
        ]  # kcal/mol * A^2 to kJ/mol * nm^2
        k_bond_list = [round(num, 10) for num in k_bond_list]
        # print(k_bond_list)
        bond_length_list = df["bond_length"].values.tolist()
        # TODO: units here? Anstroms per nm?
        bond_length_list = [i / 10.00 for i in bond_length_list]
        bond_length_list = [round(num, 6) for num in bond_length_list]
        # print(bond_length_list)
        # Angle Parameters from QM files
        ppdb = PandasPdb()
        ppdb.read_pdb(self.system_qm_pdb)
        atom_name_list = ppdb.df["ATOM"]["atom_number"].values.tolist()
        atom_name_list = [i - 1 for i in atom_name_list]
        # print(atom_name_list)
        df = pd.read_csv(
            self.angle_parameter_file, header=None, delimiter=r"\s+"
        )
        df.columns = [
            "angle",
            "k_angle",
            "angle_degrees",
            "angle_1",
            "angle_2",
            "angle_3",
        ]
        # print(df.head())
        angle_1_list = df["angle_1"].values.tolist()
        angle_1_list = [x - 1 + min(atom_name_list) for x in angle_1_list]
        # print(angle_1_list)
        angle_2_list = df["angle_2"].values.tolist()
        angle_2_list = [x - 1 + min(atom_name_list) for x in angle_2_list]
        # print(angle_2_list)
        angle_3_list = df["angle_3"].values.tolist()
        angle_3_list = [x - 1 + min(atom_name_list) for x in angle_3_list]
        # print(angle_3_list)
        k_angle_list = df["k_angle"].values.tolist()
        k_angle_list = [
            i * const.KCAL_MOL_PER_KJ_MOL for i in k_angle_list
        ]  # kcal/mol * radian^2 to kJ/mol * radian^2
        k_angle_list = [round(num, 6) for num in k_angle_list]
        # print(k_angle_list)
        angle_list = df["angle_degrees"].values.tolist()
        angle_list = [i * const.RADIANS_PER_DEGREE for i in angle_list]
        angle_list = [round(num, 6) for num in angle_list]
        # print(angle_list)
        xml = open(self.system_qm_params_file, "w")
        xml.write("Begin writing the Bond Parameters" + "\n")
        # TODO: These should use string formatting to become more concise
        for i in range(len(k_bond_list)):
            xml.write(
                "                                "
                + "<Bond"
                + " "
                + "d="
                + '"'
                + str(bond_length_list[i])
                + '"'
                + " "
                + "k="
                + '"'
                + str(k_bond_list[i])
                + '"'
                + " "
                + "p1="
                + '"'
                + str(bond_1_list[i])
                + '"'
                + " "
                + "p2="
                + '"'
                + str(bond_2_list[i])
                + '"'
                + "/>"
                + "\n"
            )
        xml.write("Finish writing the Bond Parameters" + "\n")
        xml.write("Begin writing the Angle Parameters" + "\n")
        for i in range(len(k_angle_list)):
            xml.write(
                "                                "
                + "<Angle"
                + " "
                + "a="
                + '"'
                + str(angle_list[i])
                + '"'
                + " "
                + "k="
                + '"'
                + str(k_angle_list[i])
                + '"'
                + " "
                + "p1="
                + '"'
                + str(angle_1_list[i])
                + '"'
                + " "
                + "p2="
                + '"'
                + str(angle_2_list[i])
                + '"'
                + " "
                + "p3="
                + '"'
                + str(angle_3_list[i])
                + '"'
                + "/>"
                + "\n"
            )
        xml.write("Finish writing the Angle Parameters" + "\n")
        xml.write("Begin writing the Charge Parameters" + "\n")
        for i in range(len(qm_charges)):
            xml.write(
                "<Particle"
                + " "
                + "q="
                + '"'
                + str(qm_charges[i])
                + '"'
                + " "
                + "eps="
                + '"'
                + str(0.00)
                + '"'
                + " "
                + "sig="
                + '"'
                + str(0.00)
                + '"'
                + " "
                + "atom="
                + '"'
                + str(atom_name_list[i])
                + '"'
                + "/>"
                + "\n"
            )
        xml.write("Finish writing the Charge Parameters" + "\n")
        xml.close()

    def write_intermediate_reparameterised_system_xml(self):
        """
        Writes a reparameterised XML force field file for
        ligand but without the QM obtained charges.
        """
        # Bond Parameters
        f_params = open(self.system_qm_params_file, "r")
        lines_params = f_params.readlines()
        # Bond Parameters
        for i in range(len(lines_params)):
            if "Begin writing the Bond Parameters" in lines_params[i]:
                to_begin = int(i)
            if "Finish writing the Bond Parameters" in lines_params[i]:
                to_end = int(i)
        bond_params = lines_params[to_begin + 1 : to_end]
        index_search_replace_bond = []
        for i in bond_params:
            bond_line_to_replace = i
            # print(bond_line_to_replace)
            atom_number_list = [
                re.findall("\d*\.?\d+", i)[3],
                re.findall("\d*\.?\d+", i)[5],
            ]
            # print(atom_number_list)
            comb_1 = (
                "p1="
                + '"'
                + atom_number_list[0]
                + '"'
                + " "
                + "p2="
                + '"'
                + atom_number_list[1]
                + '"'
                + "/>"
            )
            comb_2 = (
                "p1="
                + '"'
                + atom_number_list[1]
                + '"'
                + " "
                + "p2="
                + '"'
                + atom_number_list[0]
                + '"'
                + "/>"
            )
            comb_list_bond = [comb_1, comb_2]
            # print(comb_list_bond)
            list_search_bond = [
                file_utilities.search_in_file(file=self.system_xml, word=comb_1),
                file_utilities.search_in_file(file=self.system_xml, word=comb_2),
            ]
            # print(list_search_bond)
            for j in range(len(list_search_bond)):
                if list_search_bond[j] != []:
                    to_add = (list_search_bond[j], i)
                    # print(to_add)
                    index_search_replace_bond.append(to_add)
        # Angle Parameters
        for i in range(len(lines_params)):
            if "Begin writing the Angle Parameters" in lines_params[i]:
                to_begin = int(i)
            if "Finish writing the Angle Parameters" in lines_params[i]:
                to_end = int(i)
        angle_params = lines_params[to_begin + 1 : to_end]
        index_search_replace_angle = []
        for i in angle_params:
            angle_line_to_replace = i
            # print(angle_line_to_replace)
        index_search_replace_angle = []
        for i in angle_params:
            angle_line_to_replace = i
            # print(angle_line_to_replace)
            atom_number_list = [
                re.findall("\d*\.?\d+", i)[3],
                re.findall("\d*\.?\d+", i)[5],
                re.findall("\d*\.?\d+", i)[7],
            ]
            # print(atom_number_list)
            comb_1 = (
                "p1="
                + '"'
                + atom_number_list[0]
                + '"'
                + " "
                + "p2="
                + '"'
                + atom_number_list[1]
                + '"'
                + " "
                + "p3="
                + '"'
                + atom_number_list[2]
                + '"'
                + "/>"
            )
            comb_2 = (
                "p1="
                + '"'
                + atom_number_list[0]
                + '"'
                + " "
                + "p2="
                + '"'
                + atom_number_list[2]
                + '"'
                + " "
                + "p3="
                + '"'
                + atom_number_list[1]
                + '"'
                + "/>"
            )
            comb_3 = (
                "p1="
                + '"'
                + atom_number_list[1]
                + '"'
                + " "
                + "p2="
                + '"'
                + atom_number_list[0]
                + '"'
                + " "
                + "p3="
                + '"'
                + atom_number_list[2]
                + '"'
                + "/>"
            )
            comb_4 = (
                "p1="
                + '"'
                + atom_number_list[1]
                + '"'
                + " "
                + "p2="
                + '"'
                + atom_number_list[2]
                + '"'
                + " "
                + "p3="
                + '"'
                + atom_number_list[0]
                + '"'
                + "/>"
            )
            comb_5 = (
                "p1="
                + '"'
                + atom_number_list[2]
                + '"'
                + " "
                + "p2="
                + '"'
                + atom_number_list[0]
                + '"'
                + " "
                + "p3="
                + '"'
                + atom_number_list[1]
                + '"'
                + "/>"
            )
            comb_6 = (
                "p1="
                + '"'
                + atom_number_list[2]
                + '"'
                + " "
                + "p2="
                + '"'
                + atom_number_list[1]
                + '"'
                + " "
                + "p3="
                + '"'
                + atom_number_list[0]
                + '"'
                + "/>"
            )
            comb_list_angle = [
                comb_1,
                comb_2,
                comb_3,
                comb_4,
                comb_5,
                comb_6,
            ]
            # print(comb_list_angle)
            list_search_angle = [
                file_utilities.search_in_file(file=self.system_xml, word=comb_1),
                file_utilities.search_in_file(file=self.system_xml, word=comb_2),
                file_utilities.search_in_file(file=self.system_xml, word=comb_3),
                file_utilities.search_in_file(file=self.system_xml, word=comb_4),
                file_utilities.search_in_file(file=self.system_xml, word=comb_5),
                file_utilities.search_in_file(file=self.system_xml, word=comb_6),
            ]
            # print(list_search_angle)
            for j in range(len(list_search_angle)):
                if list_search_angle[j] != []:
                    to_add = (list_search_angle[j], i)
                    # print(to_add)
                    index_search_replace_angle.append(to_add)
        f_org = open(self.system_xml)
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
        f_cop = open(self.reparameterised_intermediate_system_xml_file, "w")
        for i in lines:
            f_cop.write(i)
        f_cop.close()

    def write_reparameterised_system_xml(self):
        """
        Writes a reparameterised XML force field file for the ligand.
        """
        # Bond Parameters
        f_params = open(self.system_qm_params_file, "r")
        lines_params = f_params.readlines()
        # Bond Parameters
        for i in range(len(lines_params)):
            if "Begin writing the Bond Parameters" in lines_params[i]:
                to_begin = int(i)
            if "Finish writing the Bond Parameters" in lines_params[i]:
                to_end = int(i)
        bond_params = lines_params[to_begin + 1 : to_end]
        index_search_replace_bond = []
        # TODO: These should use string formatting to become more concise
        for i in bond_params:
            bond_line_to_replace = i
            # print(bond_line_to_replace)
            atom_number_list = [
                re.findall("\d*\.?\d+", i)[3],
                re.findall("\d*\.?\d+", i)[5],
            ]
            # print(atom_number_list)
            comb_1 = (
                "p1="
                + '"'
                + atom_number_list[0]
                + '"'
                + " "
                + "p2="
                + '"'
                + atom_number_list[1]
                + '"'
                + "/>"
            )
            comb_2 = (
                "p1="
                + '"'
                + atom_number_list[1]
                + '"'
                + " "
                + "p2="
                + '"'
                + atom_number_list[0]
                + '"'
                + "/>"
            )
            comb_list_bond = [comb_1, comb_2]
            # print(comb_list_bond)
            list_search_bond = [
                file_utilities.search_in_file(file=self.system_xml, word=comb_1),
                file_utilities.search_in_file(file=self.system_xml, word=comb_2),
            ]
            # print(list_search_bond)
            for j in range(len(list_search_bond)):
                if list_search_bond[j] != []:
                    to_add = (list_search_bond[j], i)
                    # print(to_add)
                    index_search_replace_bond.append(to_add)
        # Angle Parameters
        for i in range(len(lines_params)):
            if "Begin writing the Angle Parameters" in lines_params[i]:
                to_begin = int(i)
            if "Finish writing the Angle Parameters" in lines_params[i]:
                to_end = int(i)
        angle_params = lines_params[to_begin + 1 : to_end]
        index_search_replace_angle = []
        for i in angle_params:
            angle_line_to_replace = i
            # print(angle_line_to_replace)
        index_search_replace_angle = []
        for i in angle_params:
            angle_line_to_replace = i
            # print(angle_line_to_replace)
            atom_number_list = [
                re.findall("\d*\.?\d+", i)[3],
                re.findall("\d*\.?\d+", i)[5],
                re.findall("\d*\.?\d+", i)[7],
            ]
            # print(atom_number_list)
            comb_1 = (
                "p1="
                + '"'
                + atom_number_list[0]
                + '"'
                + " "
                + "p2="
                + '"'
                + atom_number_list[1]
                + '"'
                + " "
                + "p3="
                + '"'
                + atom_number_list[2]
                + '"'
                + "/>"
            )
            comb_2 = (
                "p1="
                + '"'
                + atom_number_list[0]
                + '"'
                + " "
                + "p2="
                + '"'
                + atom_number_list[2]
                + '"'
                + " "
                + "p3="
                + '"'
                + atom_number_list[1]
                + '"'
                + "/>"
            )
            comb_3 = (
                "p1="
                + '"'
                + atom_number_list[1]
                + '"'
                + " "
                + "p2="
                + '"'
                + atom_number_list[0]
                + '"'
                + " "
                + "p3="
                + '"'
                + atom_number_list[2]
                + '"'
                + "/>"
            )
            comb_4 = (
                "p1="
                + '"'
                + atom_number_list[1]
                + '"'
                + " "
                + "p2="
                + '"'
                + atom_number_list[2]
                + '"'
                + " "
                + "p3="
                + '"'
                + atom_number_list[0]
                + '"'
                + "/>"
            )
            comb_5 = (
                "p1="
                + '"'
                + atom_number_list[2]
                + '"'
                + " "
                + "p2="
                + '"'
                + atom_number_list[0]
                + '"'
                + " "
                + "p3="
                + '"'
                + atom_number_list[1]
                + '"'
                + "/>"
            )
            comb_6 = (
                "p1="
                + '"'
                + atom_number_list[2]
                + '"'
                + " "
                + "p2="
                + '"'
                + atom_number_list[1]
                + '"'
                + " "
                + "p3="
                + '"'
                + atom_number_list[0]
                + '"'
                + "/>"
            )
            comb_list_angle = [
                comb_1,
                comb_2,
                comb_3,
                comb_4,
                comb_5,
                comb_6,
            ]
            # print(comb_list_angle)
            list_search_angle = [
                file_utilities.search_in_file(file=self.system_xml, word=comb_1),
                file_utilities.search_in_file(file=self.system_xml, word=comb_2),
                file_utilities.search_in_file(file=self.system_xml, word=comb_3),
                file_utilities.search_in_file(file=self.system_xml, word=comb_4),
                file_utilities.search_in_file(file=self.system_xml, word=comb_5),
                file_utilities.search_in_file(file=self.system_xml, word=comb_6),
            ]
            # print(list_search_angle)
            for j in range(len(list_search_angle)):
                if list_search_angle[j] != []:
                    to_add = (list_search_angle[j], i)
                    # print(to_add)
                    index_search_replace_angle.append(to_add)
        f_org = open(self.system_xml)
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
        f_cop = open(self.reparameterised_intermediate_system_xml_file, "w")
        for i in lines:
            f_cop.write(i)
        f_cop.close()

        f_params = open(self.system_qm_params_file)
        lines_params = f_params.readlines()
        # Charge Parameters
        for i in range(len(lines_params)):
            if "Begin writing the Charge Parameters" in lines_params[i]:
                to_begin = int(i)
            if "Finish writing the Charge Parameters" in lines_params[i]:
                to_end = int(i)
        charge_params = lines_params[to_begin + 1 : to_end]
        non_bonded_index = []
        for k in charge_params:
            non_bonded_index.append(int(re.findall("[-+]?\d*\.\d+|\d+", k)[3]))
        charge_for_index = []
        for k in charge_params:
            charge_for_index.append(
                float(re.findall("[-+]?\d*\.\d+|\d+", k)[0])
            )

        xml_off = open(self.system_xml)
        xml_off_lines = xml_off.readlines()
        for i in range(len(xml_off_lines)):
            if "<GlobalParameters/>" in xml_off_lines[i]:
                to_begin = int(i)
            if "<Exceptions>" in xml_off_lines[i]:
                to_end = int(i)
        nonbond_params = xml_off_lines[to_begin + 4 : to_end - 1]
        # print(len(nonbond_params))
        f_non_bonded = open(self.system_xml_non_bonded_file, "w")
        for x in nonbond_params:
            f_non_bonded.write(x)
        f_non_bonded = open(self.system_xml_non_bonded_file)
        lines_non_bonded = f_non_bonded.readlines()
        # print(len(lines_non_bonded))
        lines_non_bonded_to_write = []
        for i in range(len(non_bonded_index)):
            line_ = lines_non_bonded[non_bonded_index[i]]
            # print(line_)
            eps = float(re.findall("[-+]?\d*\.\d+|\d+", line_)[0])
            sig = float(re.findall("[-+]?\d*\.\d+|\d+", line_)[2])
            line_to_replace = (
                "                                "
                + "<Particle "
                + "eps="
                + '"'
                + str(eps)
                + '"'
                + " "
                + "q="
                + '"'
                + str(charge_for_index[i])
                + '"'
                + " "
                + "sig="
                + '"'
                + str(sig)
                + '"'
                + "/>"
            )
            lines_non_bonded_to_write.append(line_to_replace)
        data_ = list(zip(non_bonded_index, lines_non_bonded_to_write))

        df_non_bonded_params = pd.DataFrame(
            data_, columns=["line_index", "line"]
        )
        # print(df_non_bonded_params.head())
        f_non_bonded_ = open(self.system_xml_non_bonded_file)
        lines_non_bonded_ = f_non_bonded_.readlines()
        for i in range(len(lines_non_bonded_)):
            if i in non_bonded_index:
                lines_non_bonded_[i] = (
                    df_non_bonded_params.loc[
                        df_non_bonded_params.line_index == i, "line"
                    ].values[0]
                ) + "\n"
        # print(len(lines_non_bonded_))
        f_write_non_bonded_reparams = open(
            self.system_xml_non_bonded_reparams_file, "w"
        )
        for p in range(len(lines_non_bonded_)):
            f_write_non_bonded_reparams.write(lines_non_bonded_[p])
        f_write_non_bonded_reparams.close()
        f_ = open(self.system_xml_non_bonded_reparams_file)
        lines_ = f_.readlines()
        print(len(lines_) == len(lines_non_bonded))

        xml_off = open(self.reparameterised_intermediate_system_xml_file)
        # TODO: implement function(s) to read certain types of files. DRY principle
        xml_off_lines = xml_off.readlines()
        for i in range(len(xml_off_lines)):
            if "<GlobalParameters/>" in xml_off_lines[i]:
                to_begin = int(i)
            if "<Exceptions>" in xml_off_lines[i]:
                to_end = int(i)

        lines_before_params = xml_off_lines[: to_begin + 4]
        f__ = open(self.system_xml_non_bonded_reparams_file)
        lines_params_non_bonded = f__.readlines()
        lines_after_params = xml_off_lines[to_end - 1 :]
        f_reparams_xml = open(self.reparameterised_system_xml_file, "w")
        for x in lines_before_params:
            f_reparams_xml.write(x)
        for x in lines_params_non_bonded:
            f_reparams_xml.write(x)
        for x in lines_after_params:
            f_reparams_xml.write(x)
        f_reparams_xml.close()

    def save_amber_params_non_qm_charges(self):
        """
        Saves amber generated topology files for the ligand
        without the QM charges.
        """
        if self.load_topology == "parmed":
            openmm_system = parmed.openmm.load_topology(
                parmed.load_file(self.system_pdb, structure=True).topology,
                parmed.load_file(self.non_reparameterised_system_xml_file),
            )
        if self.load_topology == "openmm":
            openmm_system = parmed.openmm.load_topology(
                simtk.openmm.app.PDBFile(self.system_pdb).topology,
                parmed.load_file(self.non_reparameterised_system_xml_file),
            )
        openmm_system.save(self.prmtop_system_non_params, overwrite=True)
        openmm_system.coordinates = parmed.load_file(
            self.system_pdb, structure=True
        ).coordinates
        openmm_system.save(self.inpcrd_system_non_params, overwrite=True)
        parm = parmed.load_file(
            self.prmtop_system_non_params, self.inpcrd_system_non_params,
        )

        xml_energy_decomposition = parmed.openmm.energy_decomposition_system(
            openmm_system,
            parmed.load_file(self.non_reparameterised_system_xml_file),
        )
        xml_energy_decomposition_value = [
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in xml_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("HarmonicBondForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in xml_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("HarmonicAngleForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in xml_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("PeriodicTorsionForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in xml_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("NonbondedForce"),
        ]
        xml_energy_decomposition_list = [
            "HarmonicBondForce",
            "HarmonicAngleForce",
            "PeriodicTorsionForce",
            "NonbondedForce",
        ]
        df_energy_xml = pd.DataFrame(
            list(
                zip(
                    xml_energy_decomposition_list,
                    xml_energy_decomposition_value,
                )
            ),
            columns=["Energy_term", "Energy_xml_non_params"],
        )
        df_energy_xml = df_energy_xml.set_index("Energy_term")
        prmtop_energy_decomposition = parmed.openmm.energy_decomposition_system(
            parm, parm.createSystem()
        )
        prmtop_energy_decomposition_value = [
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in prmtop_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("HarmonicBondForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in prmtop_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("HarmonicAngleForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in prmtop_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("PeriodicTorsionForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in prmtop_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("NonbondedForce"),
        ]
        prmtop_energy_decomposition_list = [
            "HarmonicBondForce",
            "HarmonicAngleForce",
            "PeriodicTorsionForce",
            "NonbondedForce",
        ]
        df_energy_prmtop = pd.DataFrame(
            list(
                zip(
                    prmtop_energy_decomposition_list,
                    prmtop_energy_decomposition_value,
                )
            ),
            columns=["Energy_term", "Energy_prmtop_non_params"],
        )
        df_energy_prmtop = df_energy_prmtop.set_index("Energy_term")
        df_compare = pd.concat([df_energy_xml, df_energy_prmtop], axis=1)
        print(df_compare)
        if self.load_topology == "parmed":
            openmm_system = parmed.openmm.load_topology(
                parmed.load_file(self.system_pdb, structure=True).topology,
                parmed.load_file(
                    self.reparameterised_intermediate_system_xml_file
                ),
            )
        if self.load_topology == "openmm":
            openmm_system = parmed.openmm.load_topology(
                simtk.openmm.app.PDBFile(self.system_pdb).topology,
                parmed.load_file(
                    self.reparameterised_intermediate_system_xml_file
                ),
            )
        openmm_system.save(self.prmtop_system_params, overwrite=True)
        openmm_system.coordinates = parmed.load_file(
            self.system_pdb, structure=True
        ).coordinates
        openmm_system.save(self.inpcrd_system_params, overwrite=True)
        parm = parmed.load_file(
            self.prmtop_system_params, self.inpcrd_system_params
        )

        xml_energy_decomposition = parmed.openmm.energy_decomposition_system(
            openmm_system,
            parmed.load_file(
                self.reparameterised_intermediate_system_xml_file
            ),
        )
        xml_energy_decomposition_value = [
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in xml_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("HarmonicBondForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in xml_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("HarmonicAngleForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in xml_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("PeriodicTorsionForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in xml_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("NonbondedForce"),
        ]
        xml_energy_decomposition_list = [
            "HarmonicBondForce",
            "HarmonicAngleForce",
            "PeriodicTorsionForce",
            "NonbondedForce",
        ]
        df_energy_xml = pd.DataFrame(
            list(
                zip(
                    xml_energy_decomposition_list,
                    xml_energy_decomposition_value,
                )
            ),
            columns=["Energy_term", "Energy_xml_params"],
        )
        df_energy_xml = df_energy_xml.set_index("Energy_term")

        prmtop_energy_decomposition = parmed.openmm.energy_decomposition_system(
            parm, parm.createSystem()
        )
        prmtop_energy_decomposition_value = [
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in prmtop_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("HarmonicBondForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in prmtop_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("HarmonicAngleForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in prmtop_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("PeriodicTorsionForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in prmtop_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("NonbondedForce"),
        ]
        prmtop_energy_decomposition_list = [
            "HarmonicBondForce",
            "HarmonicAngleForce",
            "PeriodicTorsionForce",
            "NonbondedForce",
        ]
        df_energy_prmtop = pd.DataFrame(
            list(
                zip(
                    prmtop_energy_decomposition_list,
                    prmtop_energy_decomposition_value,
                )
            ),
            columns=["Energy_term", "Energy_prmtop_params"],
        )
        df_energy_prmtop = df_energy_prmtop.set_index("Energy_term")
        df_compare = pd.concat([df_energy_xml, df_energy_prmtop], axis=1)
        print(df_compare)

    def save_amber_params(self):
        """
        Saves amber generated topology files for the ligand.
        """
        if self.load_topology == "parmed":
            openmm_system = parmed.openmm.load_topology(
                parmed.load_file(self.system_pdb, structure=True).topology,
                parmed.load_file(self.non_reparameterised_system_xml_file),
            )
        if self.load_topology == "openmm":
            openmm_system = parmed.openmm.load_topology(
                simtk.openmm.app.PDBFile(self.system_pdb).topology,
                parmed.load_file(self.non_reparameterised_system_xml_file),
            )
        openmm_system.save(self.prmtop_system_non_params, overwrite=True)
        openmm_system.coordinates = parmed.load_file(
            self.system_pdb, structure=True
        ).coordinates
        openmm_system.save(self.inpcrd_system_non_params, overwrite=True)
        parm = parmed.load_file(
            self.prmtop_system_non_params, self.inpcrd_system_non_params,
        )

        xml_energy_decomposition = parmed.openmm.energy_decomposition_system(
            openmm_system,
            parmed.load_file(self.non_reparameterised_system_xml_file),
        )
        xml_energy_decomposition_value = [
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in xml_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("HarmonicBondForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in xml_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("HarmonicAngleForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in xml_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("PeriodicTorsionForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in xml_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("NonbondedForce"),
        ]
        xml_energy_decomposition_list = [
            "HarmonicBondForce",
            "HarmonicAngleForce",
            "PeriodicTorsionForce",
            "NonbondedForce",
        ]
        df_energy_xml = pd.DataFrame(
            list(
                zip(
                    xml_energy_decomposition_list,
                    xml_energy_decomposition_value,
                )
            ),
            columns=["Energy_term", "Energy_xml_non_params"],
        )
        df_energy_xml = df_energy_xml.set_index("Energy_term")

        prmtop_energy_decomposition = parmed.openmm.energy_decomposition_system(
            parm, parm.createSystem()
        )
        prmtop_energy_decomposition_value = [
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in prmtop_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("HarmonicBondForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in prmtop_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("HarmonicAngleForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in prmtop_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("PeriodicTorsionForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in prmtop_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("NonbondedForce"),
        ]
        prmtop_energy_decomposition_list = [
            "HarmonicBondForce",
            "HarmonicAngleForce",
            "PeriodicTorsionForce",
            "NonbondedForce",
        ]
        df_energy_prmtop = pd.DataFrame(
            list(
                zip(
                    prmtop_energy_decomposition_list,
                    prmtop_energy_decomposition_value,
                )
            ),
            columns=["Energy_term", "Energy_prmtop_non_params"],
        )
        df_energy_prmtop = df_energy_prmtop.set_index("Energy_term")

        df_compare = pd.concat([df_energy_xml, df_energy_prmtop], axis=1)
        print(df_compare)

        if self.load_topology == "parmed":
            openmm_system = parmed.openmm.load_topology(
                parmed.load_file(self.system_pdb, structure=True).topology,
                parmed.load_file(self.reparameterised_system_xml_file),
            )
        if self.load_topology == "openmm":
            openmm_system = parmed.openmm.load_topology(
                simtk.openmm.app.PDBFile(self.system_pdb).topology,
                parmed.load_file(self.reparameterised_system_xml_file),
            )
        openmm_system.save(self.prmtop_system_params, overwrite=True)
        openmm_system.coordinates = parmed.load_file(
            self.system_pdb, structure=True
        ).coordinates
        openmm_system.save(self.inpcrd_system_params, overwrite=True)
        parm = parmed.load_file(
            self.prmtop_system_params, self.inpcrd_system_params
        )

        xml_energy_decomposition = parmed.openmm.energy_decomposition_system(
            openmm_system,
            parmed.load_file(self.reparameterised_system_xml_file),
        )
        xml_energy_decomposition_value = [
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in xml_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("HarmonicBondForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in xml_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("HarmonicAngleForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in xml_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("PeriodicTorsionForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in xml_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("NonbondedForce"),
        ]
        xml_energy_decomposition_list = [
            "HarmonicBondForce",
            "HarmonicAngleForce",
            "PeriodicTorsionForce",
            "NonbondedForce",
        ]
        df_energy_xml = pd.DataFrame(
            list(
                zip(
                    xml_energy_decomposition_list,
                    xml_energy_decomposition_value,
                )
            ),
            columns=["Energy_term", "Energy_xml_params"],
        )
        df_energy_xml = df_energy_xml.set_index("Energy_term")

        prmtop_energy_decomposition = parmed.openmm.energy_decomposition_system(
            parm, parm.createSystem()
        )
        prmtop_energy_decomposition_value = [
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in prmtop_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("HarmonicBondForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in prmtop_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("HarmonicAngleForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in prmtop_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("PeriodicTorsionForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in prmtop_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("NonbondedForce"),
        ]
        prmtop_energy_decomposition_list = [
            "HarmonicBondForce",
            "HarmonicAngleForce",
            "PeriodicTorsionForce",
            "NonbondedForce",
        ]
        df_energy_prmtop = pd.DataFrame(
            list(
                zip(
                    prmtop_energy_decomposition_list,
                    prmtop_energy_decomposition_value,
                )
            ),
            columns=["Energy_term", "Energy_prmtop_params"],
        )
        df_energy_prmtop = df_energy_prmtop.set_index("Energy_term")

        df_compare = pd.concat([df_energy_xml, df_energy_prmtop], axis=1)
        print(df_compare)

    def analyze_diff_energies(self):
        """
        Compares the energies of the ligand obtained from the non-parameterized
        and the parameterized force field files.
        """
        parm_non_params = parmed.load_file(
            self.prmtop_system_non_params, self.inpcrd_system_non_params,
        )
        prmtop_energy_decomposition_non_params = parmed.openmm.energy_decomposition_system(
            parm_non_params, parm_non_params.createSystem()
        )
        prmtop_energy_decomposition_non_params_value = [
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem)
                        for elem in prmtop_energy_decomposition_non_params
                    ]
                    for item in sublist
                ]
            ).get("HarmonicBondForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem)
                        for elem in prmtop_energy_decomposition_non_params
                    ]
                    for item in sublist
                ]
            ).get("HarmonicAngleForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem)
                        for elem in prmtop_energy_decomposition_non_params
                    ]
                    for item in sublist
                ]
            ).get("PeriodicTorsionForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem)
                        for elem in prmtop_energy_decomposition_non_params
                    ]
                    for item in sublist
                ]
            ).get("NonbondedForce"),
        ]
        prmtop_energy_decomposition_non_params_list = [
            "HarmonicBondForce",
            "HarmonicAngleForce",
            "PeriodicTorsionForce",
            "NonbondedForce",
        ]
        df_energy_non_params = pd.DataFrame(
            list(
                zip(
                    prmtop_energy_decomposition_non_params_list,
                    prmtop_energy_decomposition_non_params_value,
                )
            ),
            columns=["Energy_term", "Energy_parm_non_params"],
        )
        df_energy_non_params = df_energy_non_params.set_index("Energy_term")
        # print(df_energy_non_params)
        parm_params = parmed.load_file(
            self.prmtop_system_params, self.inpcrd_system_params
        )
        prmtop_energy_decomposition_params = parmed.openmm.energy_decomposition_system(
            parm_params, parm_params.createSystem()
        )
        prmtop_energy_decomposition_params_value = [
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem)
                        for elem in prmtop_energy_decomposition_params
                    ]
                    for item in sublist
                ]
            ).get("HarmonicBondForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem)
                        for elem in prmtop_energy_decomposition_params
                    ]
                    for item in sublist
                ]
            ).get("HarmonicAngleForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem)
                        for elem in prmtop_energy_decomposition_params
                    ]
                    for item in sublist
                ]
            ).get("PeriodicTorsionForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem)
                        for elem in prmtop_energy_decomposition_params
                    ]
                    for item in sublist
                ]
            ).get("NonbondedForce"),
        ]
        prmtop_energy_decomposition_params_list = [
            "HarmonicBondForce",
            "HarmonicAngleForce",
            "PeriodicTorsionForce",
            "NonbondedForce",
        ]
        df_energy_params = pd.DataFrame(
            list(
                zip(
                    prmtop_energy_decomposition_params_list,
                    prmtop_energy_decomposition_params_value,
                )
            ),
            columns=["Energy_term", "Energy_parm_params"],
        )
        df_energy_params = df_energy_params.set_index("Energy_term")
        # print(df_energy_params)
        df_compare = pd.concat(
            [df_energy_non_params, df_energy_params], axis=1
        )
        df_compare["Energy_difference"] = df_compare[
            "Energy_parm_non_params"
        ].sub(df_compare["Energy_parm_params"], axis=0)
        print(df_compare)


class HostAmberXMLAmber:

    """
    A class used to generate a template force field XML file for the receptor
    in order regenerate the reparameterised forcefield XML file.

    This class contain methods to generate a template XML force field through
    openforcefield. Re-parameterized XML force field files are then
    generated from the template files. Different energy components such as
    bond, angle, torsional and non-bonded energies are computed for the
    non-reparametrized and the reparameterized force fields. Difference
    between the non-reparameterized and reparameterized force field energies
    can then be analyzed.
    ...

    Attributes
    ----------

    system_pdb: str, optional
        Receptor PDB file with atom numbers beginning from 1.

    system_sdf: str, optional
        Receptor SDF (structure-data) format file.

    charge : int
        Charge of the ligand.

    system_mol2: str, optional
        Receptor Mol2 file obtained from PDB file.

    system_in: str, optional
        Prepi file as required by antechamber.

    system_frcmod: str, optional
        FRCMOD file as required by antechamber.

    prmtop_system : str, optional
        Topology file obtained from the receptor PDB.

    inpcrd_system : str, optional
        Coordinate file obtained from the receptor PDB using the
        command saveamberparm.

    system_leap : str, optional
        Amber generated leap file for generating and saving topology
        and coordinate files.

    system_xml: str, optional
        Serilazed XML force field file of the receptor.

    sim_output: str, optional
        PDB file containing the trajectory coordinates for the OpenMM
        simulation.

    sim_steps: str, optional
        Number of steps in the OpenMM MD simulation.

    charge_parameter_file: str, optional
        File containing the charges of receptor atoms and their
        corresponding atoms.

    system_qm_pdb: str, optional
        Receptor QM region's PDB file with atom numbers beginning from 1.

    bond_parameter_file: str, optional
        Text file containing the bond parameters for the receptor.

    angle_parameter_file: str, optional
        Text file containing the angle parameters of the receptor.

    system_qm_params_file: str, optional
        A text file containing the QM obtained parameters for the
        receptor.

    reparameterised_intermediate_system_xml_file: str, optional
        XML force field file with bond and angle parameter lines replaced by
        corresponding values obtained from the QM calculations.

    system_xml_non_bonded_file: str, optional
        Text file to write the NonBondedForce Charge Parameters from
        the non-parameterised system XML file.

    system_xml_non_bonded_reparams_file: str, optional
        Text file containing the non-bonded parameters parsed from the
        XML force field file.

    reparameterised_system_xml_file: str, optional
        Reparameterized force field XML file obtained using
        openforcefield.

    non_reparameterised_system_xml_file: str, optional
        Non-reparameterized force field XML file obtained using
        openforcefield.

    prmtop_system_non_params: str, optional
        Amber generated topology file saved from the non-reparameterized
        force field XML file for the receptor.

    inpcrd_system_non_params: str, optional
        Amber generated coordinate file saved from the non-reparameterized
        force field XML file for the receptor.

    prmtop_system_params: str, optional
        Amber generated topology file saved from the reparameterized
        force field XML file for the receptor.

    inpcrd_system_params: str, optional
        Amber generated coordinate file saved from the reparameterized
        force field XML file for the receptor.

    load_topology: str, optional
        Argument to specify how to load the topology. Can either be "openmm"
        or "parmed".

    """

    def __init__(
        self,
        system_pdb="host.pdb",
        system_sdf="host.sdf",
        charge=0,
        system_mol2="host.mol2",
        system_in="host.in",
        system_frcmod="host.frcmod",
        prmtop_system="host.prmtop",
        inpcrd_system="host.inpcrd",
        system_leap="host.leap",
        system_xml="host.xml",
        sim_output="sim_output.pdb",
        sim_steps=1000,
        charge_parameter_file="host_qm_surround_charges.txt",
        system_qm_pdb="host_qm.pdb",
        bond_parameter_file="host_qm_bonds.txt",
        angle_parameter_file="host_qm_angles.txt",
        system_qm_params_file="host_qm_params.txt",
        reparameterised_intermediate_system_xml_file="host_intermediate_reparameterised.xml",
        system_xml_non_bonded_file="host_xml_non_bonded.txt",
        system_xml_non_bonded_reparams_file="host_xml_non_bonded_reparams.txt",
        reparameterised_system_xml_file="host_reparameterised.xml",
        non_reparameterised_system_xml_file="host.xml",
        prmtop_system_non_params="host_non_params.prmtop",
        inpcrd_system_non_params="host_non_params.inpcrd",
        prmtop_system_params="host_params.prmtop",
        inpcrd_system_params="host_params.inpcrd",
        load_topology="openmm",
    ):
        self.system_pdb = system_pdb
        self.system_sdf = system_sdf
        self.charge = charge
        self.system_mol2 = system_mol2
        self.system_in = system_in
        self.system_frcmod = system_frcmod
        self.prmtop_system = prmtop_system
        self.inpcrd_system = inpcrd_system
        self.system_leap = system_leap
        self.system_xml = system_xml
        self.sim_output = sim_output
        self.sim_steps = sim_steps
        self.charge_parameter_file = charge_parameter_file
        self.system_qm_pdb = system_qm_pdb
        self.bond_parameter_file = bond_parameter_file
        self.angle_parameter_file = angle_parameter_file
        self.system_qm_params_file = system_qm_params_file
        self.reparameterised_intermediate_system_xml_file = (
            reparameterised_intermediate_system_xml_file
        )
        self.system_xml_non_bonded_file = system_xml_non_bonded_file
        self.system_xml_non_bonded_reparams_file = (
            system_xml_non_bonded_reparams_file
        )
        self.reparameterised_system_xml_file = reparameterised_system_xml_file
        self.non_reparameterised_system_xml_file = (
            non_reparameterised_system_xml_file
        )
        self.prmtop_system_non_params = prmtop_system_non_params
        self.inpcrd_system_non_params = inpcrd_system_non_params
        self.prmtop_system_params = prmtop_system_params
        self.inpcrd_system_params = inpcrd_system_params
        self.load_topology = load_topology

    def generate_xml_from_pdb_sdf(self):
        """
        Generates an XML forcefield file from the SDF file through
        openforcefield.
        """
        command = (
            # "babel -ipdb " + self.system_pdb + " -osdf " + self.system_sdf
            "obabel -ipdb "
            + self.system_pdb
            + " -osdf -O "
            + self.system_sdf
        )
        os.system(command)
        # off_molecule = openforcefield.topology.Molecule(self.system_sdf)
        off_molecule = Molecule(self.system_sdf)
        # force_field = openforcefield.typing.engines.smirnoff.ForceField("openff_unconstrained-1.0.0.offxml")
        force_field = ForceField("openff_unconstrained-1.0.0.offxml")
        system = force_field.create_openmm_system(off_molecule.to_topology())
        pdbfile = simtk.openmm.app.PDBFile(self.system_pdb)
        structure = parmed.openmm.load_topology(
            pdbfile.topology, system, xyz=pdbfile.positions
        )
        with open(self.system_xml, "w") as f:
            f.write(simtk.openmm.XmlSerializer.serialize(system))

    def generate_xml_antechamber(self):
        """
        Generates an XML forcefield file from the PDB file through antechamber.
        """
        command = (
            # "babel -ipdb " + self.system_pdb + " -omol2 " + self.system_mol2
            "obabel -ipdb "
            + self.system_pdb
            + " -omol2 -O "
            + self.system_mol2
        )
        os.system(command)
        command = (
            "antechamber -i "
            + self.system_mol2
            + " -fi mol2 -o "
            + self.system_in
            + " -fo prepi -c bcc -nc "
            + str(self.charge)
        )
        os.system(command)
        command = (
            "parmchk2 -i "
            + self.system_in
            + " -o "
            + self.system_frcmod
            + " -f prepi -a Y"
        )
        os.system(command)
        os.system(
            "rm -rf ANTECHAMBER* leap.log sqm* ATOMTYPE.INF PREP.INF NEWPDB.PDB"
        )
        line_1 = "loadamberprep " + self.system_in
        line_2 = "loadamberparams " + self.system_frcmod
        line_3 = "pdb = loadpdb " + self.system_pdb
        line_4 = (
            "saveamberparm pdb "
            + self.prmtop_system
            + " "
            + self.inpcrd_system
        )
        line_5 = "quit"
        with open(self.system_leap, "w") as f:
            f.write("    " + "\n")
            f.write(line_1 + "\n")
            f.write(line_2 + "\n")
            f.write(line_3 + "\n")
            f.write(line_4 + "\n")
            f.write(line_5 + "\n")
        command = "tleap -f " + self.system_leap
        os.system(command)
        parm = parmed.load_file(self.prmtop_system, self.inpcrd_system)
        system = parm.createSystem()
        with open(self.system_xml, "w") as f:
            f.write(simtk.openmm.XmlSerializer.serialize(system))

    def serialize_system(self):
        pdb = simtk.openmm.app.PDBFile(self.system_pdb)
        forcefield = simtk.openmm.app.ForceField("amber14-all.xml")
        system = forcefield.createSystem(pdb.topology)
        integrator = simtk.openmm.LangevinIntegrator(
            300 * simtk.unit.kelvin,
            1 / simtk.unit.picosecond,
            0.002 * simtk.unit.picoseconds,
        )
        simulation = simtk.openmm.app.Simulation(
            pdb.topology, system, integrator
        )
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy(maxIterations=100000)
        state = simulation.context.getState(getEnergy=True)
        energy = state.getPotentialEnergy()
        print(energy)
        simulation.reporters.append(
            simtk.openmm.app.PDBReporter(self.sim_output, self.sim_steps / 10)
        )
        simulation.reporters.append(
            simtk.openmm.app.StateDataReporter(
                stdout,
                reportInterval=int(self.sim_steps / 10),
                step=True,
                potentialEnergy=True,
                temperature=True,
            )
        )
        simulation.step(self.sim_steps)
        command = "rm -rf " + self.sim_output
        os.system(command)
        with open(self.system_xml, "w") as f:
            f.write(simtk.openmm.XmlSerializer.serialize(system))

    def write_system_params(self):
        """
        Saves the parameters obtained from the QM log files in a text file.
        """
        # Charges from QM files
        df_charges = pd.read_csv(
            self.charge_parameter_file, header=None, delimiter=r"\s+"
        )
        df_charges.columns = ["atom", "charges"]
        qm_charges = df_charges["charges"].values.tolist()
        qm_charges = [round(num, 6) for num in qm_charges]
        # print(qm_charges)
        # Bond Parameters from QM files
        ppdb = PandasPdb()
        ppdb.read_pdb(self.system_qm_pdb)
        atom_name_list = ppdb.df["ATOM"]["atom_number"].values.tolist()
        atom_name_list = [i - 1 for i in atom_name_list]
        # print(atom_name_list)
        df = pd.read_csv(
            self.bond_parameter_file, header=None, delimiter=r"\s+"
        )
        df.columns = [
            "bond",
            "k_bond",
            "bond_length",
            "bond_1",
            "bond_2",
        ]
        # print(df.head())
        bond_1_list = df["bond_1"].values.tolist()
        bond_1_list = [x - 1 + min(atom_name_list) for x in bond_1_list]
        bond_2_list = df["bond_2"].values.tolist()
        bond_2_list = [x - 1 + min(atom_name_list) for x in bond_2_list]
        # print(bond_1_list)
        # print(bond_2_list)
        k_bond_list = df["k_bond"].values.tolist()
        k_bond_list = [
            i * const.KCAL_MOL_PER_KJ_MOL * const.ANGSTROMS_PER_NM**2 for i in k_bond_list
        ]  # kcal/mol * A^2 to kJ/mol * nm^2
        k_bond_list = [round(num, 10) for num in k_bond_list]
        # print(k_bond_list)
        bond_length_list = df["bond_length"].values.tolist()
        bond_length_list = [i / 10.00 for i in bond_length_list]
        bond_length_list = [round(num, 6) for num in bond_length_list]
        # print(bond_length_list)
        # Angle Parameters from QM files
        ppdb = PandasPdb()
        ppdb.read_pdb(self.system_qm_pdb)
        atom_name_list = ppdb.df["ATOM"]["atom_number"].values.tolist()
        atom_name_list = [i - 1 for i in atom_name_list]
        # print(atom_name_list)
        df = pd.read_csv(
            self.angle_parameter_file, header=None, delimiter=r"\s+"
        )
        df.columns = [
            "angle",
            "k_angle",
            "angle_degrees",
            "angle_1",
            "angle_2",
            "angle_3",
        ]
        # print(df.head())
        angle_1_list = df["angle_1"].values.tolist()
        angle_1_list = [x - 1 + min(atom_name_list) for x in angle_1_list]
        # print(angle_1_list)
        angle_2_list = df["angle_2"].values.tolist()
        angle_2_list = [x - 1 + min(atom_name_list) for x in angle_2_list]
        # print(angle_2_list)
        angle_3_list = df["angle_3"].values.tolist()
        angle_3_list = [x - 1 + min(atom_name_list) for x in angle_3_list]
        # print(angle_3_list)
        k_angle_list = df["k_angle"].values.tolist()
        k_angle_list = [
            i * 4.184 for i in k_angle_list
        ]  # kcal/mol * radian^2 to kJ/mol * radian^2
        k_angle_list = [round(num, 6) for num in k_angle_list]
        # print(k_angle_list)
        angle_list = df["angle_degrees"].values.tolist()
        angle_list = [(i * math.pi) / 180.00 for i in angle_list]
        angle_list = [round(num, 6) for num in angle_list]
        # print(angle_list)
        xml = open(self.system_qm_params_file, "w")
        xml.write("Begin writing the Bond Parameters" + "\n")
        for i in range(len(k_bond_list)):
            xml.write(
                "                                "
                + "<Bond"
                + " "
                + "d="
                + '"'
                + str(bond_length_list[i])
                + '"'
                + " "
                + "k="
                + '"'
                + str(k_bond_list[i])
                + '"'
                + " "
                + "p1="
                + '"'
                + str(bond_1_list[i])
                + '"'
                + " "
                + "p2="
                + '"'
                + str(bond_2_list[i])
                + '"'
                + "/>"
                + "\n"
            )
        xml.write("Finish writing the Bond Parameters" + "\n")
        xml.write("Begin writing the Angle Parameters" + "\n")
        for i in range(len(k_angle_list)):
            xml.write(
                "                                "
                + "<Angle"
                + " "
                + "a="
                + '"'
                + str(angle_list[i])
                + '"'
                + " "
                + "k="
                + '"'
                + str(k_angle_list[i])
                + '"'
                + " "
                + "p1="
                + '"'
                + str(angle_1_list[i])
                + '"'
                + " "
                + "p2="
                + '"'
                + str(angle_2_list[i])
                + '"'
                + " "
                + "p3="
                + '"'
                + str(angle_3_list[i])
                + '"'
                + "/>"
                + "\n"
            )
        xml.write("Finish writing the Angle Parameters" + "\n")
        xml.write("Begin writing the Charge Parameters" + "\n")
        for i in range(len(qm_charges)):
            xml.write(
                "<Particle"
                + " "
                + "q="
                + '"'
                + str(qm_charges[i])
                + '"'
                + " "
                + "eps="
                + '"'
                + str(0.00)
                + '"'
                + " "
                + "sig="
                + '"'
                + str(0.00)
                + '"'
                + " "
                + "atom="
                + '"'
                + str(atom_name_list[i])
                + '"'
                + "/>"
                + "\n"
            )
        xml.write("Finish writing the Charge Parameters" + "\n")
        xml.close()

    def write_reparameterised_system_xml(self):
        """
        Writes a reparameterised XML force field file for the ligand.
        """
        # Bond Parameters
        f_params = open(self.system_qm_params_file, "r")
        lines_params = f_params.readlines()
        # Bond Parameters
        for i in range(len(lines_params)):
            if "Begin writing the Bond Parameters" in lines_params[i]:
                to_begin = int(i)
            if "Finish writing the Bond Parameters" in lines_params[i]:
                to_end = int(i)
        bond_params = lines_params[to_begin + 1 : to_end]
        index_search_replace_bond = []
        for i in bond_params:
            bond_line_to_replace = i
            # print(bond_line_to_replace)
            atom_number_list = [
                re.findall("\d*\.?\d+", i)[3],
                re.findall("\d*\.?\d+", i)[5],
            ]
            # print(atom_number_list)
            comb_1 = (
                "p1="
                + '"'
                + atom_number_list[0]
                + '"'
                + " "
                + "p2="
                + '"'
                + atom_number_list[1]
                + '"'
                + "/>"
            )
            comb_2 = (
                "p1="
                + '"'
                + atom_number_list[1]
                + '"'
                + " "
                + "p2="
                + '"'
                + atom_number_list[0]
                + '"'
                + "/>"
            )
            comb_list_bond = [comb_1, comb_2]
            # print(comb_list_bond)
            list_search_bond = [
                file_utilities.search_in_file(file=self.system_xml, word=comb_1),
                file_utilities.search_in_file(file=self.system_xml, word=comb_2),
            ]
            # print(list_search_bond)
            for j in range(len(list_search_bond)):
                if list_search_bond[j] != []:
                    to_add = (list_search_bond[j], i)
                    # print(to_add)
                    index_search_replace_bond.append(to_add)
        # Angle Parameters
        for i in range(len(lines_params)):
            if "Begin writing the Angle Parameters" in lines_params[i]:
                to_begin = int(i)
            if "Finish writing the Angle Parameters" in lines_params[i]:
                to_end = int(i)
        angle_params = lines_params[to_begin + 1 : to_end]
        index_search_replace_angle = []
        for i in angle_params:
            angle_line_to_replace = i
            # print(angle_line_to_replace)
        index_search_replace_angle = []
        for i in angle_params:
            angle_line_to_replace = i
            # print(angle_line_to_replace)
            atom_number_list = [
                re.findall("\d*\.?\d+", i)[3],
                re.findall("\d*\.?\d+", i)[5],
                re.findall("\d*\.?\d+", i)[7],
            ]
            # print(atom_number_list)
            comb_1 = (
                "p1="
                + '"'
                + atom_number_list[0]
                + '"'
                + " "
                + "p2="
                + '"'
                + atom_number_list[1]
                + '"'
                + " "
                + "p3="
                + '"'
                + atom_number_list[2]
                + '"'
                + "/>"
            )
            comb_2 = (
                "p1="
                + '"'
                + atom_number_list[0]
                + '"'
                + " "
                + "p2="
                + '"'
                + atom_number_list[2]
                + '"'
                + " "
                + "p3="
                + '"'
                + atom_number_list[1]
                + '"'
                + "/>"
            )
            comb_3 = (
                "p1="
                + '"'
                + atom_number_list[1]
                + '"'
                + " "
                + "p2="
                + '"'
                + atom_number_list[0]
                + '"'
                + " "
                + "p3="
                + '"'
                + atom_number_list[2]
                + '"'
                + "/>"
            )
            comb_4 = (
                "p1="
                + '"'
                + atom_number_list[1]
                + '"'
                + " "
                + "p2="
                + '"'
                + atom_number_list[2]
                + '"'
                + " "
                + "p3="
                + '"'
                + atom_number_list[0]
                + '"'
                + "/>"
            )
            comb_5 = (
                "p1="
                + '"'
                + atom_number_list[2]
                + '"'
                + " "
                + "p2="
                + '"'
                + atom_number_list[0]
                + '"'
                + " "
                + "p3="
                + '"'
                + atom_number_list[1]
                + '"'
                + "/>"
            )
            comb_6 = (
                "p1="
                + '"'
                + atom_number_list[2]
                + '"'
                + " "
                + "p2="
                + '"'
                + atom_number_list[1]
                + '"'
                + " "
                + "p3="
                + '"'
                + atom_number_list[0]
                + '"'
                + "/>"
            )
            comb_list_angle = [
                comb_1,
                comb_2,
                comb_3,
                comb_4,
                comb_5,
                comb_6,
            ]
            # print(comb_list_angle)
            list_search_angle = [
                file_utilities.search_in_file(file=self.system_xml, word=comb_1),
                file_utilities.search_in_file(file=self.system_xml, word=comb_2),
                file_utilities.search_in_file(file=self.system_xml, word=comb_3),
                file_utilities.search_in_file(file=self.system_xml, word=comb_4),
                file_utilities.search_in_file(file=self.system_xml, word=comb_5),
                file_utilities.search_in_file(file=self.system_xml, word=comb_6),
            ]
            # print(list_search_angle)
            for j in range(len(list_search_angle)):
                if list_search_angle[j] != []:
                    to_add = (list_search_angle[j], i)
                    # print(to_add)
                    index_search_replace_angle.append(to_add)
        f_org = open(self.system_xml)
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
        f_cop = open(self.reparameterised_intermediate_system_xml_file, "w")
        for i in lines:
            f_cop.write(i)
        f_cop.close()

        f_params = open(self.system_qm_params_file)
        lines_params = f_params.readlines()
        # Charge Parameters
        for i in range(len(lines_params)):
            if "Begin writing the Charge Parameters" in lines_params[i]:
                to_begin = int(i)
            if "Finish writing the Charge Parameters" in lines_params[i]:
                to_end = int(i)
        charge_params = lines_params[to_begin + 1 : to_end]
        non_bonded_index = []
        for k in charge_params:
            non_bonded_index.append(int(re.findall("[-+]?\d*\.\d+|\d+", k)[3]))
        charge_for_index = []
        for k in charge_params:
            charge_for_index.append(
                float(re.findall("[-+]?\d*\.\d+|\d+", k)[0])
            )

        xml_off = open(self.system_xml)
        xml_off_lines = xml_off.readlines()
        for i in range(len(xml_off_lines)):
            if "<GlobalParameters/>" in xml_off_lines[i]:
                to_begin = int(i)
            if "<Exceptions>" in xml_off_lines[i]:
                to_end = int(i)
        nonbond_params = xml_off_lines[to_begin + 4 : to_end - 1]
        # print(len(nonbond_params))
        f_non_bonded = open(self.system_xml_non_bonded_file, "w")
        for x in nonbond_params:
            f_non_bonded.write(x)
        f_non_bonded = open(self.system_xml_non_bonded_file)
        lines_non_bonded = f_non_bonded.readlines()
        # print(len(lines_non_bonded))
        lines_non_bonded_to_write = []
        for i in range(len(non_bonded_index)):
            line_ = lines_non_bonded[non_bonded_index[i]]
            # print(line_)
            eps = float(re.findall("[-+]?\d*\.\d+|\d+", line_)[0])
            sig = float(re.findall("[-+]?\d*\.\d+|\d+", line_)[2])
            line_to_replace = (
                "                                "
                + "<Particle "
                + "eps="
                + '"'
                + str(eps)
                + '"'
                + " "
                + "q="
                + '"'
                + str(charge_for_index[i])
                + '"'
                + " "
                + "sig="
                + '"'
                + str(sig)
                + '"'
                + "/>"
            )
            lines_non_bonded_to_write.append(line_to_replace)
        data_ = list(zip(non_bonded_index, lines_non_bonded_to_write))

        df_non_bonded_params = pd.DataFrame(
            data_, columns=["line_index", "line"]
        )
        # print(df_non_bonded_params.head())
        f_non_bonded_ = open(self.system_xml_non_bonded_file)
        lines_non_bonded_ = f_non_bonded_.readlines()
        for i in range(len(lines_non_bonded_)):
            if i in non_bonded_index:
                lines_non_bonded_[i] = (
                    df_non_bonded_params.loc[
                        df_non_bonded_params.line_index == i, "line"
                    ].values[0]
                ) + "\n"
        # print(len(lines_non_bonded_))
        f_write_non_bonded_reparams = open(
            self.system_xml_non_bonded_reparams_file, "w"
        )
        for p in range(len(lines_non_bonded_)):
            f_write_non_bonded_reparams.write(lines_non_bonded_[p])
        f_write_non_bonded_reparams.close()
        f_ = open(self.system_xml_non_bonded_reparams_file)
        lines_ = f_.readlines()
        print(len(lines_) == len(lines_non_bonded))

        xml_off = open(self.reparameterised_intermediate_system_xml_file)
        xml_off_lines = xml_off.readlines()
        for i in range(len(xml_off_lines)):
            if "<GlobalParameters/>" in xml_off_lines[i]:
                to_begin = int(i)
            if "<Exceptions>" in xml_off_lines[i]:
                to_end = int(i)

        lines_before_params = xml_off_lines[: to_begin + 4]
        f__ = open(self.system_xml_non_bonded_reparams_file)
        lines_params_non_bonded = f__.readlines()
        lines_after_params = xml_off_lines[to_end - 1 :]
        f_reparams_xml = open(self.reparameterised_system_xml_file, "w")
        for x in lines_before_params:
            f_reparams_xml.write(x)
        for x in lines_params_non_bonded:
            f_reparams_xml.write(x)
        for x in lines_after_params:
            f_reparams_xml.write(x)
        f_reparams_xml.close()

    def save_amber_params(self):
        """
        Saves amber generated topology files for the ligand.
        """
        if self.load_topology == "parmed":
            openmm_system = parmed.openmm.load_topology(
                parmed.load_file(self.system_pdb, structure=True).topology,
                parmed.load_file(self.non_reparameterised_system_xml_file),
            )
        if self.load_topology == "openmm":
            openmm_system = parmed.openmm.load_topology(
                simtk.openmm.app.PDBFile(self.system_pdb).topology,
                parmed.load_file(self.non_reparameterised_system_xml_file),
            )
        openmm_system.save(self.prmtop_system_non_params, overwrite=True)
        openmm_system.coordinates = parmed.load_file(
            self.system_pdb, structure=True
        ).coordinates
        openmm_system.save(self.inpcrd_system_non_params, overwrite=True)
        parm = parmed.load_file(
            self.prmtop_system_non_params, self.inpcrd_system_non_params,
        )

        xml_energy_decomposition = parmed.openmm.energy_decomposition_system(
            openmm_system,
            parmed.load_file(self.non_reparameterised_system_xml_file),
        )
        xml_energy_decomposition_value = [
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in xml_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("HarmonicBondForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in xml_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("HarmonicAngleForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in xml_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("PeriodicTorsionForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in xml_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("NonbondedForce"),
        ]
        xml_energy_decomposition_list = [
            "HarmonicBondForce",
            "HarmonicAngleForce",
            "PeriodicTorsionForce",
            "NonbondedForce",
        ]
        df_energy_xml = pd.DataFrame(
            list(
                zip(
                    xml_energy_decomposition_list,
                    xml_energy_decomposition_value,
                )
            ),
            columns=["Energy_term", "Energy_xml_non_params"],
        )
        df_energy_xml = df_energy_xml.set_index("Energy_term")

        prmtop_energy_decomposition = parmed.openmm.energy_decomposition_system(
            parm, parm.createSystem()
        )
        prmtop_energy_decomposition_value = [
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in prmtop_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("HarmonicBondForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in prmtop_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("HarmonicAngleForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in prmtop_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("PeriodicTorsionForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in prmtop_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("NonbondedForce"),
        ]
        prmtop_energy_decomposition_list = [
            "HarmonicBondForce",
            "HarmonicAngleForce",
            "PeriodicTorsionForce",
            "NonbondedForce",
        ]
        df_energy_prmtop = pd.DataFrame(
            list(
                zip(
                    prmtop_energy_decomposition_list,
                    prmtop_energy_decomposition_value,
                )
            ),
            columns=["Energy_term", "Energy_prmtop_non_params"],
        )
        df_energy_prmtop = df_energy_prmtop.set_index("Energy_term")

        df_compare = pd.concat([df_energy_xml, df_energy_prmtop], axis=1)
        print(df_compare)

        if self.load_topology == "parmed":
            openmm_system = parmed.openmm.load_topology(
                parmed.load_file(self.system_pdb, structure=True).topology,
                parmed.load_file(self.reparameterised_system_xml_file),
            )
        if self.load_topology == "openmm":
            openmm_system = parmed.openmm.load_topology(
                simtk.openmm.app.PDBFile(self.system_pdb).topology,
                parmed.load_file(self.reparameterised_system_xml_file),
            )
        openmm_system.save(self.prmtop_system_params, overwrite=True)
        openmm_system.coordinates = parmed.load_file(
            self.system_pdb, structure=True
        ).coordinates
        openmm_system.save(self.inpcrd_system_params, overwrite=True)
        parm = parmed.load_file(
            self.prmtop_system_params, self.inpcrd_system_params
        )

        xml_energy_decomposition = parmed.openmm.energy_decomposition_system(
            openmm_system,
            parmed.load_file(self.reparameterised_system_xml_file),
        )
        xml_energy_decomposition_value = [
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in xml_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("HarmonicBondForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in xml_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("HarmonicAngleForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in xml_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("PeriodicTorsionForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in xml_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("NonbondedForce"),
        ]
        xml_energy_decomposition_list = [
            "HarmonicBondForce",
            "HarmonicAngleForce",
            "PeriodicTorsionForce",
            "NonbondedForce",
        ]
        df_energy_xml = pd.DataFrame(
            list(
                zip(
                    xml_energy_decomposition_list,
                    xml_energy_decomposition_value,
                )
            ),
            columns=["Energy_term", "Energy_xml_params"],
        )
        df_energy_xml = df_energy_xml.set_index("Energy_term")

        prmtop_energy_decomposition = parmed.openmm.energy_decomposition_system(
            parm, parm.createSystem()
        )
        prmtop_energy_decomposition_value = [
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in prmtop_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("HarmonicBondForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in prmtop_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("HarmonicAngleForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in prmtop_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("PeriodicTorsionForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem) for elem in prmtop_energy_decomposition
                    ]
                    for item in sublist
                ]
            ).get("NonbondedForce"),
        ]
        prmtop_energy_decomposition_list = [
            "HarmonicBondForce",
            "HarmonicAngleForce",
            "PeriodicTorsionForce",
            "NonbondedForce",
        ]
        df_energy_prmtop = pd.DataFrame(
            list(
                zip(
                    prmtop_energy_decomposition_list,
                    prmtop_energy_decomposition_value,
                )
            ),
            columns=["Energy_term", "Energy_prmtop_params"],
        )
        df_energy_prmtop = df_energy_prmtop.set_index("Energy_term")

        df_compare = pd.concat([df_energy_xml, df_energy_prmtop], axis=1)
        print(df_compare)

    def analyze_diff_energies(self):
        """
        Compares the energies of the ligand obtained from the non-parameterized
        and the parameterized force field files.
        """
        parm_non_params = parmed.load_file(
            self.prmtop_system_non_params, self.inpcrd_system_non_params,
        )
        prmtop_energy_decomposition_non_params = parmed.openmm.energy_decomposition_system(
            parm_non_params, parm_non_params.createSystem()
        )
        prmtop_energy_decomposition_non_params_value = [
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem)
                        for elem in prmtop_energy_decomposition_non_params
                    ]
                    for item in sublist
                ]
            ).get("HarmonicBondForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem)
                        for elem in prmtop_energy_decomposition_non_params
                    ]
                    for item in sublist
                ]
            ).get("HarmonicAngleForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem)
                        for elem in prmtop_energy_decomposition_non_params
                    ]
                    for item in sublist
                ]
            ).get("PeriodicTorsionForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem)
                        for elem in prmtop_energy_decomposition_non_params
                    ]
                    for item in sublist
                ]
            ).get("NonbondedForce"),
        ]
        prmtop_energy_decomposition_non_params_list = [
            "HarmonicBondForce",
            "HarmonicAngleForce",
            "PeriodicTorsionForce",
            "NonbondedForce",
        ]
        df_energy_non_params = pd.DataFrame(
            list(
                zip(
                    prmtop_energy_decomposition_non_params_list,
                    prmtop_energy_decomposition_non_params_value,
                )
            ),
            columns=["Energy_term", "Energy_parm_non_params"],
        )
        df_energy_non_params = df_energy_non_params.set_index("Energy_term")
        # print(df_energy_non_params)
        parm_params = parmed.load_file(
            self.prmtop_system_params, self.inpcrd_system_params
        )
        prmtop_energy_decomposition_params = parmed.openmm.energy_decomposition_system(
            parm_params, parm_params.createSystem()
        )
        prmtop_energy_decomposition_params_value = [
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem)
                        for elem in prmtop_energy_decomposition_params
                    ]
                    for item in sublist
                ]
            ).get("HarmonicBondForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem)
                        for elem in prmtop_energy_decomposition_params
                    ]
                    for item in sublist
                ]
            ).get("HarmonicAngleForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem)
                        for elem in prmtop_energy_decomposition_params
                    ]
                    for item in sublist
                ]
            ).get("PeriodicTorsionForce"),
            base.list_to_dict(
                [
                    item
                    for sublist in [
                        list(elem)
                        for elem in prmtop_energy_decomposition_params
                    ]
                    for item in sublist
                ]
            ).get("NonbondedForce"),
        ]
        prmtop_energy_decomposition_params_list = [
            "HarmonicBondForce",
            "HarmonicAngleForce",
            "PeriodicTorsionForce",
            "NonbondedForce",
        ]
        df_energy_params = pd.DataFrame(
            list(
                zip(
                    prmtop_energy_decomposition_params_list,
                    prmtop_energy_decomposition_params_value,
                )
            ),
            columns=["Energy_term", "Energy_parm_params"],
        )
        df_energy_params = df_energy_params.set_index("Energy_term")
        # print(df_energy_params)
        df_compare = pd.concat(
            [df_energy_non_params, df_energy_params], axis=1
        )
        df_compare["Energy_difference"] = df_compare[
            "Energy_parm_non_params"
        ].sub(df_compare["Energy_parm_params"], axis=0)
        print(df_compare)
        
class MergeHostGuestTopology:

    """
    A class used to merge the host and guest topology and coordinate
    files.

    ...

    Attributes
    ----------
    host_prmtop : str
        Topology file of the receptor.

    guest_prmtop : str
        Topology file of the ligand.

    host_inpcrd : str
        Coordinate file of the receptor.

    guest_inpcrd : str
        Coordinate file of the ligand.

    system_prmtop : str
        Topology file of the receptor - ligand complex.

    system_inpcrd : str
        Coordinate file of the receptor - ligand complex.

    """

    def __init__(
        self,
        host_prmtop,
        guest_prmtop,
        host_inpcrd,
        guest_inpcrd,
        system_prmtop,
        system_inpcrd,
    ):

        self.host_prmtop = host_prmtop
        self.guest_prmtop = guest_prmtop
        self.host_inpcrd = host_inpcrd
        self.guest_inpcrd = guest_inpcrd
        self.system_prmtop = system_prmtop
        self.system_inpcrd = system_inpcrd

    def merge_topology_files(self):
        """
        Merge the host and guest topology and coordinate files.
        """
        print(
            "Merging the "
            + self.host_prmtop
            + " "
            + self.guest_prmtop
            + " files"
        )
        print(
            "Merging the "
            + self.host_inpcrd
            + " "
            + self.guest_inpcrd
            + " files"
        )
        host_system = parmed.load_file(self.host_prmtop, xyz=self.host_inpcrd)
        guest_system = parmed.load_file(
            self.guest_prmtop, xyz=self.guest_inpcrd
        )
        system = host_system + guest_system
        system.save(self.system_prmtop, overwrite=True)
        system.save(self.system_inpcrd, overwrite=True)

