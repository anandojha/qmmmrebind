"""

"""

import os
import re

from biopandas.pdb import PandasPdb
import pandas as pd
import parmed
import simtk

import modules.constants as const
import modules.base as base
import modules.file_utilities as file_utilities

class PrepareSolvatedParams:

    """
    A class used to integrate the parameterized topology
    files of the receptor - ligand complex and the solvent.

    This class contain methods to concatanate the solvent (and
    ions ) and the receptor - ligand complex in a single
    parameterized topology file (prmtop and inpcrd).

    ...

    Attributes
    ----------

    init_pdb : str
        Initial PDB file containing the receptor-ligand complex with
        solvent, ions, etc.

    intermediate_pdb : str, optional
        An intermediate PDB file formed during pdb4amber processing.

    solvent_pdb : str, optional
        PDB file containing the water, ions, etc.

    solvent_prmtop : str, optional
        Solvent topology file.

    solvent_inpcrd : str, optional
        Solvent coordinate file.

    solvent_amber_pdb : str, optional
        Solvent PDB file saved from Amber's tleap.

    solvent_leap : str, optional
        Solvent tleap file for parameterizing the solvent.

    system_prmtop : str, optional
        Topology file of the receptor - ligand complex.

    system_inpcrd : str, optional
        Coordinate file of the receptor - ligand complex.

    system_output: str, optional
        PDB file containing the trajectory coordinates for
        the OpenMM simulation.

    sim_steps: str, optional
        Number of steps in the OpenMM MD simulation.

    system_solvent_prmtop : str, optional
        Topology file of the receptor - ligand complex and
        the solvent.

    system_solvent_inpcrd : str, optional
        Coordinate file of the receptor - ligand complex and
        the solvent.

    system_solvent_pdb : str, optional
        PDB file of the receptor - ligand complex and
        the solvent.

    """

    def __init__(
        self,
        init_pdb,
        intermediate_pdb="intermediate.pdb",
        solvent_pdb="solvent.pdb",
        solvent_prmtop="solvent.prmtop",
        solvent_inpcrd="solvent.inpcrd",
        solvent_amber_pdb="solvent_amber.pdb",
        solvent_leap="solvent.leap",
        system_prmtop="system_torsional_params.prmtop",
        system_inpcrd="system_torsional_params.inpcrd",
        system_output="sim_output.pdb",
        sim_steps=1000,
        system_solvent_prmtop="system_qmmmrebind.prmtop",
        system_solvent_inpcrd="system_qmmmrebind.inpcrd",
        system_solvent_pdb="system_qmmmrebind.pdb",
    ):

        self.init_pdb = init_pdb
        self.intermediate_pdb = intermediate_pdb
        self.solvent_pdb = solvent_pdb
        self.solvent_prmtop = solvent_prmtop
        self.solvent_inpcrd = solvent_inpcrd
        self.solvent_amber_pdb = solvent_amber_pdb
        self.solvent_leap = solvent_leap
        self.system_prmtop = system_prmtop
        self.system_inpcrd = system_inpcrd
        self.system_output = system_output
        self.sim_steps = sim_steps
        self.system_solvent_prmtop = system_solvent_prmtop
        self.system_solvent_inpcrd = system_solvent_inpcrd
        self.system_solvent_pdb = system_solvent_pdb

    def create_solvent_pdb(self):
        """
        Generates a PDB file containing the solvent and the ions.
        """
        water_variables = ["HOH", "WAT"]
        ions = [
            "Na+",
            "Cs+",
            "K+",
            "Li+",
            "Rb+",
            "Cl-",
            "Br-",
            "F-",
            "I-",
            "Ca2",
        ]
        pdb_variables = ["END", "CRYST"]
        with open(self.init_pdb) as f1, open(self.intermediate_pdb, "w") as f2:
            for line in f1:
                if (
                    any(
                        water_variable in line
                        for water_variable in water_variables
                    )
                    or any(
                        pdb_variable in line for pdb_variable in pdb_variables
                    )
                    or any(ion in line for ion in ions)
                ):
                    f2.write(line)
        command = (
            "pdb4amber -i " + self.intermediate_pdb + " -o " + self.solvent_pdb
        )
        os.system(command)
        command = (
            "rm -rf "
            + self.solvent_pdb[:-4]
            + "_nonprot.pdb "
            + self.solvent_pdb[:-4]
            + "_renum.txt "
            + self.solvent_pdb[:-4]
            + "_sslink"
        )
        os.system(command)
        command = "rm -rf " + self.intermediate_pdb
        os.system(command)

    def parameterize_solvent_pdb(self):
        """
        Generates a topology file (prmtop) and a coordinate
        file (inpcrd) for the solvent system.
        """
        line_0 = " "
        line_1 = "source leaprc.protein.ff14SB"
        line_2 = "source leaprc.water.tip3p"
        line_3 = "loadAmberParams frcmod.ionsjc_tip3p"
        line_4 = "pdb = loadpdb " + self.solvent_pdb
        line_5 = (
            "saveamberparm pdb "
            + self.solvent_prmtop
            + " "
            + self.solvent_inpcrd
        )
        line_6 = "savepdb pdb " + self.solvent_amber_pdb
        line_7 = "quit"
        with open(self.solvent_leap, "w") as f:
            f.write(line_0 + "\n")
            f.write(line_1 + "\n")
            f.write(line_2 + "\n")
            f.write(line_3 + "\n")
            f.write(line_4 + "\n")
            f.write(line_5 + "\n")
            f.write(line_6 + "\n")
            f.write(line_7 + "\n")
        command = "tleap -f " + self.solvent_leap
        os.system(command)
        command = "rm -rf leap.log " + self.solvent_leap
        os.system(command)

    def run_openmm_solvent_prmtop_inpcrd(self):
        """
        Runs OpenMM MD simulation with prmtop and inpcrd file
        for the solvent.
        """
        print(
            "Running OpenMM simulation for "
            + self.solvent_prmtop
            + " and "
            + self.solvent_inpcrd
        )
        prmtop = simtk.openmm.app.AmberPrmtopFile(self.solvent_prmtop)
        inpcrd = simtk.openmm.app.AmberInpcrdFile(self.solvent_inpcrd)
        system = prmtop.createSystem()
        integrator = simtk.openmm.LangevinIntegrator(
            300 * simtk.unit.kelvin,
            1 / simtk.unit.picosecond,
            0.002 * simtk.unit.picoseconds,
        )
        simulation = simtk.openmm.app.Simulation(
            prmtop.topology, system, integrator
        )
        simulation.context.setPositions(inpcrd.positions)
        if inpcrd.boxVectors is not None:
            simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
        simulation.minimizeEnergy(maxIterations=100000)
        simulation.reporters.append(
            simtk.openmm.app.PDBReporter(
                self.system_output, self.sim_steps / 10
            )
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
        command = "rm -rf " + self.system_output
        os.system(command)

    def run_openmm_solvent_prmtop_pdb(self):
        """
        Runs OpenMM MD simulation with prmtop and PDB file
        for the solvent.
        """
        print(
            "Running OpenMM simulation for "
            + self.solvent_prmtop
            + " and "
            + self.solvent_amber_pdb
        )
        pdb = simtk.openmm.app.PDBFile(self.solvent_amber_pdb)
        prmtop = simtk.openmm.app.AmberPrmtopFile(self.solvent_prmtop)
        system = prmtop.createSystem()
        integrator = simtk.openmm.LangevinIntegrator(
            300 * simtk.unit.kelvin,
            1 / simtk.unit.picosecond,
            0.002 * simtk.unit.picoseconds,
        )
        simulation = simtk.openmm.app.Simulation(
            prmtop.topology, system, integrator
        )
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy(maxIterations=100000)
        simulation.reporters.append(
            simtk.openmm.app.PDBReporter(
                self.system_output, self.sim_steps / 10
            )
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
        command = "rm -rf " + self.system_output
        os.system(command)

    def merge_topology_files_system_solvent(self):
        """
        Merge the system and solvent topology and coordinate
        files.
        """
        print(
            "Merging the "
            + self.system_prmtop
            + " "
            + self.solvent_prmtop
            + " files"
        )
        print(
            "Merging the "
            + self.system_inpcrd
            + " "
            + self.solvent_inpcrd
            + " files"
        )
        system = parmed.load_file(self.system_prmtop, xyz=self.system_inpcrd)
        solvent = parmed.load_file(
            self.solvent_prmtop, xyz=self.solvent_inpcrd
        )
        system_solvent = system + solvent
        system_solvent.save(self.system_solvent_prmtop, overwrite=True)
        system_solvent.save(self.system_solvent_inpcrd, overwrite=True)
        system_solvent.save(self.system_solvent_pdb, overwrite=True)

    def run_openmm_system_solvent_prmtop_inpcrd(self):
        """
        Runs OpenMM MD simulation with prmtop and inpcrd file
        for the solvent - system complex.
        """
        print(
            "Running OpenMM simulation for "
            + self.system_solvent_prmtop
            + " and "
            + self.system_solvent_inpcrd
        )
        prmtop = simtk.openmm.app.AmberPrmtopFile(self.system_solvent_prmtop)
        inpcrd = simtk.openmm.app.AmberInpcrdFile(self.system_solvent_inpcrd)
        system = prmtop.createSystem()
        integrator = simtk.openmm.LangevinIntegrator(
            300 * simtk.unit.kelvin,
            1 / simtk.unit.picosecond,
            0.002 * simtk.unit.picoseconds,
        )
        simulation = simtk.openmm.app.Simulation(
            prmtop.topology, system, integrator
        )
        simulation.context.setPositions(inpcrd.positions)
        if inpcrd.boxVectors is not None:
            simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
        simulation.minimizeEnergy(maxIterations=100000)
        simulation.reporters.append(
            simtk.openmm.app.PDBReporter(
                self.system_output, self.sim_steps / 10
            )
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
        command = "rm -rf " + self.system_output
        os.system(command)

    def run_openmm_system_solvent_prmtop_pdb(self):
        """
        Runs OpenMM MD simulation with prmtop and PDB file
        for the solvent - system complex.
        """
        print(
            "Running OpenMM simulation for "
            + self.system_solvent_prmtop
            + " and "
            + self.system_solvent_pdb
        )
        pdb = simtk.openmm.app.PDBFile(self.system_solvent_pdb)
        prmtop = simtk.openmm.app.AmberPrmtopFile(self.system_solvent_prmtop)
        system = prmtop.createSystem()
        integrator = simtk.openmm.LangevinIntegrator(
            300 * simtk.unit.kelvin,
            1 / simtk.unit.picosecond,
            0.002 * simtk.unit.picoseconds,
        )
        simulation = simtk.openmm.app.Simulation(
            prmtop.topology, system, integrator
        )
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy(maxIterations=100000)
        simulation.reporters.append(
            simtk.openmm.app.PDBReporter(
                self.system_output, self.sim_steps / 10
            )
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
        command = "rm -rf " + self.system_output
        os.system(command)


class SystemAmberSystem:

    """
    A class used to generate a force field XML file for the system
    from the given amber forcefield topology files and
    regenerate the reparameterised forcefield XML file.

    This class contain methods to generate a XML force field through
    parmed if the amber forcefield topology files are given.
    Re-parameterized XML force field files are then generated from
    these XML focefield files. Different energy components such as
    bond, angle, torsional and non-bonded energies are computed for the
    non-reparametrized and the reparameterized force fields. Difference
    between the non-reparameterized and reparameterized force field energies
    can then be analyzed.
    ...

     Attributes
    ----------

    host_pdb: str, optional
        PDB file for the host.

    system_pdb: str, optional
        PDB file for the system (host, guest and solvent).

    prmtop_system: str, optional
        Topology file for the system (host, guest and solvent).

    system_xml: str, optional
        Serialised XML forcefield file generated by parmed.

    charge_parameter_file_guest: str, optional
        Receptor PDB file with atom numbers beginning from 1.

    guest_qm_pdb: str, optional
        Ligand PDB file with atom numbers beginning from 1.

    bond_parameter_file_guest: str, optional
        Text file containing the bond parameters for the ligand.

    angle_parameter_file_guest: str, optional
        Text file containing the angle parameters of the ligand.

    guest_qm_params_file: str, optional
        Text file containing QM obtained parameters for the ligand.

    charge_parameter_file_host: str, optional
        File containing the charges of receptor atoms and their
        corresponding atoms.

    bond_parameter_file_host: str, optional
        Text file containing the bond parameters for the receptor.

    host_qm_pdb: str, optional
        Receptor QM region's PDB file with atom numbers beginning from 1.

    angle_parameter_file_host: str, optional
        Text file containing the angle parameters of the receptor.

    host_qm_params_file: str, optional
        Text file containing QM obtained parameters for the receptor.

    host_guest_qm_params_file: str, optional
        Text file containing QM obtained parameters for the system.

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

    reparameterized_torsional_params_file : str, optional
        Text file containing the forcefield parameters for the
        ligand previously obtained without torsional reparameterization.

    reparameterised_intermediate_torsional_system_xml_file : str, optional
        XML force field file for the system (without the QM charges) obtained
        with torsional reparamaterization.

    reparameterised_torsional_system_xml_file : str, optional
        XML force field file for the system obtained with
        torsional reparamaterization.

    load_topology: str, optional
        Argument to specify how to load the topology. Can either be "openmm"
        or "parmed".

    non_reparameterised_system_xml_file: str, optional
        Non-reparameterized force field XML file.

    prmtop_system_non_params: str, optional
        Non-reparameterized topology file.

    inpcrd_system_non_params: str, optional
        Non-reparameterized INPCRD file.

    prmtop_system_intermediate_params: str, optional
        Reparameterized topology file but without the QM charges.

    inpcrd_system_intermediate_params: str, optional
        Reparameterized INPCRD file but without the QM charges.

    prmtop_system_params: str, optional
        Reparameterized topology file.

    inpcrd_system_params: str, optional
        Reparameterized INPCRD file.

    """

    def __init__(
        self,
        host_pdb="host.pdb",
        system_pdb="",
        prmtop_system="hostguest.parm7",
        system_xml="hostguest.xml",
        charge_parameter_file_guest="guest_qm_surround_charges.txt",
        guest_qm_pdb="guest_init_II.pdb",
        bond_parameter_file_guest="guest_bonds.txt",
        angle_parameter_file_guest="guest_angles.txt",
        guest_qm_params_file="guest_qm_params.txt",
        charge_parameter_file_host="host_qm_surround_charges.txt",
        bond_parameter_file_host="host_qm_bonds.txt",
        host_qm_pdb="host_qm.pdb",
        angle_parameter_file_host="host_qm_angles.txt",
        host_qm_params_file="host_qm_params.txt",
        host_guest_qm_params_file="host_guest_qm_params.txt",
        reparameterised_intermediate_system_xml_file="hostguest_intermediate.xml",
        system_xml_non_bonded_file="hostguest_non_bonded.txt",
        system_xml_non_bonded_reparams_file="hostguest_non_bonded_reparams.txt",
        reparameterised_system_xml_file="hostguest_reparameterised.xml",
        reparameterized_torsional_params_file="reparameterized_torsional_params.txt",
        reparameterised_intermediate_torsional_system_xml_file="reparameterized_torsional_params.txt",
        reparameterised_torsional_system_xml_file="hostguest_torsional_reparameterised.xml",
        load_topology="openmm",
        non_reparameterised_system_xml_file="hostguest.xml",
        prmtop_system_non_params="hostguest.parm7",
        inpcrd_system_non_params="hostguest_non_params.pdb",
        prmtop_system_intermediate_params="hostguest_intermediate.prmtop",
        inpcrd_system_intermediate_params="hostguest_intermediate.inpcrd",
        prmtop_system_params="hostguest_params.prmtop",
        inpcrd_system_params="hostguest_params.inpcrd",
    ):

        self.host_pdb = host_pdb
        self.system_pdb = system_pdb
        self.prmtop_system = prmtop_system
        self.system_xml = system_xml
        self.charge_parameter_file_guest = charge_parameter_file_guest
        self.guest_qm_pdb = guest_qm_pdb
        self.bond_parameter_file_guest = bond_parameter_file_guest
        self.angle_parameter_file_guest = angle_parameter_file_guest
        self.guest_qm_params_file = guest_qm_params_file
        self.charge_parameter_file_host = charge_parameter_file_host
        self.bond_parameter_file_host = bond_parameter_file_host
        self.host_qm_pdb = host_qm_pdb
        self.angle_parameter_file_host = angle_parameter_file_host
        self.host_qm_params_file = host_qm_params_file
        self.host_guest_qm_params_file = host_guest_qm_params_file
        self.reparameterised_intermediate_system_xml_file = (
            reparameterised_intermediate_system_xml_file
        )
        self.system_xml_non_bonded_file = system_xml_non_bonded_file
        self.system_xml_non_bonded_reparams_file = (
            system_xml_non_bonded_reparams_file
        )
        self.reparameterised_system_xml_file = reparameterised_system_xml_file
        self.reparameterized_torsional_params_file = (
            reparameterized_torsional_params_file
        )
        self.reparameterised_intermediate_torsional_system_xml_file = (
            reparameterised_intermediate_torsional_system_xml_file
        )
        self.reparameterised_torsional_system_xml_file = (
            reparameterised_torsional_system_xml_file
        )
        self.load_topology = load_topology
        self.non_reparameterised_system_xml_file = (
            non_reparameterised_system_xml_file
        )
        self.prmtop_system_non_params = prmtop_system_non_params
        self.inpcrd_system_non_params = inpcrd_system_non_params
        self.prmtop_system_intermediate_params = (
            prmtop_system_intermediate_params
        )
        self.inpcrd_system_intermediate_params = (
            inpcrd_system_intermediate_params
        )
        self.prmtop_system_params = prmtop_system_params
        self.inpcrd_system_params = inpcrd_system_params

    def generate_xml_from_prmtop(self):

        """
        Generates a serialsed XML forcefield file through parmed, given
        the PDB file and its corresponding topology file.
        """

        parm = parmed.load_file(self.prmtop_system, self.system_pdb)
        system = parm.createSystem()
        with open(self.system_xml, "w") as f:
            f.write(simtk.openmm.XmlSerializer.serialize(system))

    def write_guest_params_non_zero(self):

        """
        Saves the parameters of the ligand obtained from the QM log files
        in a text file starting from non-zero ( indexing begins from the
        index of the last atom of the receptor ).
        """

        # Charges from QM files
        df_charges = pd.read_csv(
            self.charge_parameter_file_guest, header=None, delimiter=r"\s+"
        )
        df_charges.columns = ["atom", "charges"]
        qm_charges = df_charges["charges"].values.tolist()
        qm_charges = [round(num, 6) for num in qm_charges]
        # print(qm_charges)
        # Bond Parameters from QM files
        ppdb = PandasPdb()
        ppdb.read_pdb(self.guest_qm_pdb)
        atom_name_list = ppdb.df["ATOM"]["atom_number"].values.tolist()
        # atom_name_list = [i - 1 for i in atom_name_list]
        no_host_atoms = base.get_num_host_atoms(self.host_pdb)
        atom_name_list = [i - 1 + no_host_atoms for i in atom_name_list]
        # print(atom_name_list)
        df = pd.read_csv(
            self.bond_parameter_file_guest, header=None, delimiter=r"\s+"
        )
        df.columns = ["bond", "k_bond", "bond_length", "bond_1", "bond_2"]
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
        ppdb.read_pdb(self.guest_qm_pdb)
        atom_name_list = ppdb.df["ATOM"]["atom_number"].values.tolist()
        # atom_name_list = [i - 1 for i in atom_name_list]
        no_host_atoms = base.get_num_host_atoms(self.host_pdb)
        atom_name_list = [i - 1 + no_host_atoms for i in atom_name_list]
        # print(atom_name_list)
        df = pd.read_csv(
            self.angle_parameter_file_guest, header=None, delimiter=r"\s+"
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
        xml = open(self.guest_qm_params_file, "w")
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

    def write_host_params(self):

        """
        Saves the parameters obtained from the QM log files of the
        receptor in a text file.
        """

        # Charges from QM files
        df_charges = pd.read_csv(
            self.charge_parameter_file_host, header=None, delimiter=r"\s+"
        )
        df_charges.columns = ["atom", "charges"]
        qm_charges = df_charges["charges"].values.tolist()
        qm_charges = [round(num, 6) for num in qm_charges]
        # print(qm_charges)
        # Bond Parameters from QM files
        ppdb = PandasPdb()
        ppdb.read_pdb(self.host_qm_pdb)
        atom_name_list = ppdb.df["ATOM"]["atom_number"].values.tolist()
        atom_name_list = [i - 1 for i in atom_name_list]
        # print(atom_name_list)
        df = pd.read_csv(
            self.bond_parameter_file_host, header=None, delimiter=r"\s+"
        )
        df.columns = ["bond", "k_bond", "bond_length", "bond_1", "bond_2"]
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
        ppdb.read_pdb(self.host_qm_pdb)
        atom_name_list = ppdb.df["ATOM"]["atom_number"].values.tolist()
        atom_name_list = [i - 1 for i in atom_name_list]
        # print(atom_name_list)
        df = pd.read_csv(
            self.angle_parameter_file_host, header=None, delimiter=r"\s+"
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
        xml = open(self.host_qm_params_file, "w")
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

    def merge_qm_params(self):

        """
        Saves the parameters of the ligand obtained from the QM log files
        in a text file starting from non-zero ( indexing begins from the
        index of the last atom of the receptor ).
        """

        # Bond Parameters Host
        f_params_host = open(self.host_qm_params_file, "r")
        lines_params_host = f_params_host.readlines()
        # Bond Parameters Host
        for i in range(len(lines_params_host)):
            if "Begin writing the Bond Parameters" in lines_params_host[i]:
                to_begin = int(i)
            if "Finish writing the Bond Parameters" in lines_params_host[i]:
                to_end = int(i)
        bond_params_host = lines_params_host[to_begin + 1 : to_end]
        # Bond Parameters Guest
        f_params_guest = open(self.guest_qm_params_file, "r")
        lines_params_guest = f_params_guest.readlines()
        # Bond Parameters Guest
        for i in range(len(lines_params_guest)):
            if "Begin writing the Bond Parameters" in lines_params_guest[i]:
                to_begin = int(i)
            if "Finish writing the Bond Parameters" in lines_params_guest[i]:
                to_end = int(i)
        bond_params_guest = lines_params_guest[to_begin + 1 : to_end]
        bond_systems_params = bond_params_host + bond_params_guest
        # Angle Parameters Host
        f_params_host = open(self.host_qm_params_file, "r")
        lines_params_host = f_params_host.readlines()
        # Angle Parameters Host
        for i in range(len(lines_params_host)):
            if "Begin writing the Angle Parameters" in lines_params_host[i]:
                to_begin = int(i)
            if "Finish writing the Angle Parameters" in lines_params_host[i]:
                to_end = int(i)
        angle_params_host = lines_params_host[to_begin + 1 : to_end]
        # Angle Parameters Guest
        f_params_guest = open(self.guest_qm_params_file, "r")
        lines_params_guest = f_params_guest.readlines()
        # Angle Parameters Guest
        for i in range(len(lines_params_guest)):
            if "Begin writing the Angle Parameters" in lines_params_guest[i]:
                to_begin = int(i)
            if "Finish writing the Angle Parameters" in lines_params_guest[i]:
                to_end = int(i)
        angle_params_guest = lines_params_guest[to_begin + 1 : to_end]
        angle_systems_params = angle_params_host + angle_params_guest
        # Charge Parameters Host
        f_params_host = open(self.host_qm_params_file, "r")
        lines_params_host = f_params_host.readlines()
        # Charge Parameters Host
        for i in range(len(lines_params_host)):
            if "Begin writing the Charge Parameters" in lines_params_host[i]:
                to_begin = int(i)
            if "Finish writing the Charge Parameters" in lines_params_host[i]:
                to_end = int(i)
        charge_params_host = lines_params_host[to_begin + 1 : to_end]
        # Charge Parameters Guest
        f_params_guest = open(self.guest_qm_params_file, "r")
        lines_params_guest = f_params_guest.readlines()
        # Charge Parameters Guest
        for i in range(len(lines_params_guest)):
            if "Begin writing the Charge Parameters" in lines_params_guest[i]:
                to_begin = int(i)
            if "Finish writing the Charge Parameters" in lines_params_guest[i]:
                to_end = int(i)
        charge_params_guest = lines_params_guest[to_begin + 1 : to_end]
        charge_systems_params = charge_params_host + charge_params_guest
        system_params = open(self.host_guest_qm_params_file, "w")
        system_params.write("Begin writing the Bond Parameters" + "\n")
        for i in range(len(bond_systems_params)):
            system_params.write(bond_systems_params[i])
        system_params.write("Finish writing the Bond Parameters" + "\n")
        system_params.write("Begin writing the Angle Parameters" + "\n")
        for i in range(len(angle_systems_params)):
            system_params.write(angle_systems_params[i])
        system_params.write("Finish writing the Angle Parameters" + "\n")
        system_params.write("Begin writing the Charge Parameters" + "\n")
        for i in range(len(charge_systems_params)):
            system_params.write(charge_systems_params[i])
        system_params.write("Finish writing the Charge Parameters" + "\n")
        system_params.close()

    def write_intermediate_reparameterised_system_xml(self):

        """
        Writes a reparameterised XML force field file for the
        system but without the QM obtained charges.
        """

        # Bond Parameters
        f_params = open(self.host_guest_qm_params_file, "r")
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
            comb_list_angle = [comb_1, comb_2, comb_3, comb_4, comb_5, comb_6]
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
        Writes a reparameterised XML force field file for the system.
        """

        # Bond Parameters
        f_params = open(self.host_guest_qm_params_file, "r")
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
            comb_list_angle = [comb_1, comb_2, comb_3, comb_4, comb_5, comb_6]
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
        f_params = open(self.host_guest_qm_params_file)
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

    def write_torsional_reparams_intermediate(self):
        """
        Generates a XML force field file for the system ( without the
        QM charges ) with reparameterized torsional parameters of the ligand.
        """

        no_host_atoms = base.get_num_host_atoms(self.host_pdb)
        xml_tor = open(self.reparameterized_torsional_params_file, "r")
        xml_tor_lines = xml_tor.readlines()
        xml_tor_lines_renum = []
        for i in xml_tor_lines:
            i = i.replace(
                "p1=" + '"' + str(int(re.findall("\d*\.?\d+", i)[2])) + '"',
                "p1="
                + '"'
                + str(int(int(re.findall("\d*\.?\d+", i)[2]) + no_host_atoms))
                + '"',
            )
            i = i.replace(
                "p2=" + '"' + str(int(re.findall("\d*\.?\d+", i)[4])) + '"',
                "p2="
                + '"'
                + str(int(int(re.findall("\d*\.?\d+", i)[4]) + no_host_atoms))
                + '"',
            )
            i = i.replace(
                "p3=" + '"' + str(int(re.findall("\d*\.?\d+", i)[6])) + '"',
                "p3="
                + '"'
                + str(int(int(re.findall("\d*\.?\d+", i)[6]) + no_host_atoms))
                + '"',
            )
            i = i.replace(
                "p4=" + '"' + str(int(re.findall("\d*\.?\d+", i)[8])) + '"',
                "p4="
                + '"'
                + str(int(int(re.findall("\d*\.?\d+", i)[8]) + no_host_atoms))
                + '"',
            )
            xml_tor_lines_renum.append(i)

        non_zero_k_tor = []
        for i in xml_tor_lines_renum:
            to_find = "k=" + '"' + "0.0" + '"'
            if to_find not in i:
                non_zero_k_tor.append(i)
        # print(non_zero_k_tor)
        p1 = []
        for i in range(len(non_zero_k_tor)):
            p1.append(int(re.findall("\d*\.?\d+", non_zero_k_tor[i])[2]))
        # print(p1)
        p2 = []
        for i in range(len(non_zero_k_tor)):
            p2.append(int(re.findall("\d*\.?\d+", non_zero_k_tor[i])[4]))
        # print(p2)
        p3 = []
        for i in range(len(non_zero_k_tor)):
            p3.append(int(re.findall("\d*\.?\d+", non_zero_k_tor[i])[6]))
        # print(p3)
        p4 = []
        for i in range(len(non_zero_k_tor)):
            p4.append(int(re.findall("\d*\.?\d+", non_zero_k_tor[i])[8]))
        # print(p4)
        periodicity = []
        for i in range(len(non_zero_k_tor)):
            periodicity.append(
                int(re.findall("\d*\.?\d+", non_zero_k_tor[i])[9])
            )
        # print(periodicity)
        xml_tor_reparams = open(
            self.reparameterised_intermediate_system_xml_file, "r"
        )
        xml_tor_reparams_lines = xml_tor_reparams.readlines()
        for j in range(len(xml_tor_reparams_lines)):
            for i in range(len(non_zero_k_tor)):
                to_find_tor = (
                    "p1="
                    + '"'
                    + str(p1[i])
                    + '"'
                    + " "
                    + "p2="
                    + '"'
                    + str(p2[i])
                    + '"'
                    + " "
                    + "p3="
                    + '"'
                    + str(p3[i])
                    + '"'
                    + " "
                    + "p4="
                    + '"'
                    + str(p4[i])
                    + '"'
                    + " "
                    + "periodicity="
                    + '"'
                    + str(periodicity[i])
                    + '"'
                )
                if to_find_tor in xml_tor_reparams_lines[j]:
                    print(xml_tor_reparams_lines[j])
                    xml_tor_reparams_lines[j] = non_zero_k_tor[i]
        with open(
            self.reparameterised_intermediate_torsional_system_xml_file, "w"
        ) as f:
            for i in xml_tor_reparams_lines:
                f.write(i)

    def write_torsional_reparams(self):
        """
        Generates a XML force field file for the system with reparameterized
        torsional parameters of the ligand.
        """

        no_host_atoms = base.get_num_host_atoms(self.host_pdb)
        xml_tor = open(self.reparameterized_torsional_params_file, "r")
        xml_tor_lines = xml_tor.readlines()
        xml_tor_lines_renum = []
        for i in xml_tor_lines:
            i = i.replace(
                "p1=" + '"' + str(int(re.findall("\d*\.?\d+", i)[2])) + '"',
                "p1="
                + '"'
                + str(int(int(re.findall("\d*\.?\d+", i)[2]) + no_host_atoms))
                + '"',
            )
            i = i.replace(
                "p2=" + '"' + str(int(re.findall("\d*\.?\d+", i)[4])) + '"',
                "p2="
                + '"'
                + str(int(int(re.findall("\d*\.?\d+", i)[4]) + no_host_atoms))
                + '"',
            )
            i = i.replace(
                "p3=" + '"' + str(int(re.findall("\d*\.?\d+", i)[6])) + '"',
                "p3="
                + '"'
                + str(int(int(re.findall("\d*\.?\d+", i)[6]) + no_host_atoms))
                + '"',
            )
            i = i.replace(
                "p4=" + '"' + str(int(re.findall("\d*\.?\d+", i)[8])) + '"',
                "p4="
                + '"'
                + str(int(int(re.findall("\d*\.?\d+", i)[8]) + no_host_atoms))
                + '"',
            )
            xml_tor_lines_renum.append(i)

        non_zero_k_tor = []
        for i in xml_tor_lines_renum:
            to_find = "k=" + '"' + "0.0" + '"'
            if to_find not in i:
                non_zero_k_tor.append(i)
        # print(non_zero_k_tor)
        p1 = []
        for i in range(len(non_zero_k_tor)):
            p1.append(int(re.findall("\d*\.?\d+", non_zero_k_tor[i])[2]))
        # print(p1)
        p2 = []
        for i in range(len(non_zero_k_tor)):
            p2.append(int(re.findall("\d*\.?\d+", non_zero_k_tor[i])[4]))
        # print(p2)
        p3 = []
        for i in range(len(non_zero_k_tor)):
            p3.append(int(re.findall("\d*\.?\d+", non_zero_k_tor[i])[6]))
        # print(p3)
        p4 = []
        for i in range(len(non_zero_k_tor)):
            p4.append(int(re.findall("\d*\.?\d+", non_zero_k_tor[i])[8]))
        # print(p4)
        periodicity = []
        for i in range(len(non_zero_k_tor)):
            periodicity.append(
                int(re.findall("\d*\.?\d+", non_zero_k_tor[i])[9])
            )
        # print(periodicity)
        xml_tor_reparams = open(self.reparameterised_system_xml_file, "r")
        xml_tor_reparams_lines = xml_tor_reparams.readlines()
        for j in range(len(xml_tor_reparams_lines)):
            for i in range(len(non_zero_k_tor)):
                to_find_tor = (
                    "p1="
                    + '"'
                    + str(p1[i])
                    + '"'
                    + " "
                    + "p2="
                    + '"'
                    + str(p2[i])
                    + '"'
                    + " "
                    + "p3="
                    + '"'
                    + str(p3[i])
                    + '"'
                    + " "
                    + "p4="
                    + '"'
                    + str(p4[i])
                    + '"'
                    + " "
                    + "periodicity="
                    + '"'
                    + str(periodicity[i])
                    + '"'
                )
                if to_find_tor in xml_tor_reparams_lines[j]:
                    print(xml_tor_reparams_lines[j])
                    xml_tor_reparams_lines[j] = non_zero_k_tor[i]
        with open(self.reparameterised_torsional_system_xml_file, "w") as f:
            for i in xml_tor_reparams_lines:
                f.write(i)

    def save_amber_params_non_qm_charges(self):

        """
        Saves amber generated topology files for the system
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
            self.prmtop_system_non_params, self.inpcrd_system_non_params
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
                    self.reparameterised_intermediate_torsional_system_xml_file
                ),
            )
        if self.load_topology == "openmm":
            openmm_system = parmed.openmm.load_topology(
                simtk.openmm.app.PDBFile(self.system_pdb).topology,
                parmed.load_file(
                    self.reparameterised_intermediate_torsional_system_xml_file
                ),
            )
        openmm_system.save(
            self.prmtop_system_intermediate_params, overwrite=True
        )
        openmm_system.coordinates = parmed.load_file(
            self.system_pdb, structure=True
        ).coordinates
        openmm_system.save(
            self.inpcrd_system_intermediate_params, overwrite=True
        )
        parm = parmed.load_file(
            self.prmtop_system_intermediate_params,
            self.inpcrd_system_intermediate_params,
        )
        xml_energy_decomposition = parmed.openmm.energy_decomposition_system(
            openmm_system,
            parmed.load_file(
                self.reparameterised_intermediate_torsional_system_xml_file
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
        Saves amber generated topology files for the system.
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
            self.prmtop_system_non_params, self.inpcrd_system_non_params
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
                    self.reparameterised_torsional_system_xml_file
                ),
            )
        if self.load_topology == "openmm":
            openmm_system = parmed.openmm.load_topology(
                simtk.openmm.app.PDBFile(self.system_pdb).topology,
                parmed.load_file(
                    self.reparameterised_torsional_system_xml_file
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
            parmed.load_file(self.reparameterised_torsional_system_xml_file),
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


class SystemGuestAmberSystem:

    """
    A class used to generate a force field XML file for the system
    from the given amber forcefield topology files and
    regenerate the reparameterised forcefield XML file but without
    the host QM parameters.

    This class contain methods to generate a XML force field through
    parmed if the amber forcefield topology files are given.
    Re-parameterized XML force field files are then generated from
    these XML focefield files. Different energy components such as
    bond, angle, torsional and non-bonded energies are computed for the
    non-reparametrized and the reparameterized force fields. Difference
    between the non-reparameterized and reparameterized force field energies
    can then be analyzed.
    ...

     Attributes
    ----------

    host_pdb: str, optional
        PDB file for the host.

    system_pdb: str, optional
        PDB file for the system (host, guest and solvent).

    prmtop_system: str, optional
        Topology file for the system (host, guest and solvent).

    system_xml: str, optional
        Serialised XML forcefield file generated by parmed.

    charge_parameter_file_guest: str, optional
        Receptor PDB file with atom numbers beginning from 1.

    guest_qm_pdb: str, optional
        Ligand PDB file with atom numbers beginning from 1.

    bond_parameter_file_guest: str, optional
        Text file containing the bond parameters for the ligand.

    angle_parameter_file_guest: str, optional
        Text file containing the angle parameters of the ligand.

    guest_qm_params_file: str, optional
        Text file containing QM obtained parameters for the ligand.

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

    reparameterized_torsional_params_file : str, optional
        Text file containing the forcefield parameters for the
        ligand previously obtained without torsional reparameterization.

    reparameterised_intermediate_torsional_system_xml_file : str, optional
        XML force field file for the system (without the QM charges) obtained
        with torsional reparamaterization.

    reparameterised_torsional_system_xml_file : str, optional
        XML force field file for the system obtained with
        torsional reparamaterization.

    load_topology: str, optional
        Argument to specify how to load the topology. Can either be "openmm"
        or "parmed".

    non_reparameterised_system_xml_file: str, optional
        Non-reparameterized force field XML file.

    prmtop_system_non_params: str, optional
        Non-reparameterized topology file.

    inpcrd_system_non_params: str, optional
        Non-reparameterized INPCRD file.

    prmtop_system_intermediate_params: str, optional
        Reparameterized topology file but without the QM charges.

    inpcrd_system_intermediate_params: str, optional
        Reparameterized INPCRD file but without the QM charges.

    prmtop_system_params: str, optional
        Reparameterized topology file.

    inpcrd_system_params: str, optional
        Reparameterized INPCRD file.

    """

    def __init__(
        self,
        host_pdb="host.pdb",
        system_pdb="",
        prmtop_system="hostguest.parm7",
        system_xml="hostguest.xml",
        charge_parameter_file_guest="guest_qm_surround_charges.txt",
        guest_qm_pdb="guest_init_II.pdb",
        bond_parameter_file_guest="guest_bonds.txt",
        angle_parameter_file_guest="guest_angles.txt",
        guest_qm_params_file="guest_qm_params.txt",
        reparameterised_intermediate_system_xml_file="hostguest_intermediate.xml",
        system_xml_non_bonded_file="hostguest_non_bonded.txt",
        system_xml_non_bonded_reparams_file="hostguest_non_bonded_reparams.txt",
        reparameterised_system_xml_file="hostguest_reparameterised.xml",
        reparameterized_torsional_params_file="reparameterized_torsional_params.txt",
        reparameterised_intermediate_torsional_system_xml_file="reparameterized_torsional_params.txt",
        reparameterised_torsional_system_xml_file="hostguest_torsional_reparameterised.xml",
        load_topology="openmm",
        non_reparameterised_system_xml_file="hostguest.xml",
        prmtop_system_non_params="hostguest.parm7",
        inpcrd_system_non_params="hostguest_non_params.pdb",
        prmtop_system_intermediate_params="hostguest_intermediate.prmtop",
        inpcrd_system_intermediate_params="hostguest_intermediate.inpcrd",
        prmtop_system_params="hostguest_params.prmtop",
        inpcrd_system_params="hostguest_params.inpcrd",
    ):

        self.host_pdb = host_pdb
        self.system_pdb = system_pdb
        self.prmtop_system = prmtop_system
        self.system_xml = system_xml
        self.charge_parameter_file_guest = charge_parameter_file_guest
        self.guest_qm_pdb = guest_qm_pdb
        self.bond_parameter_file_guest = bond_parameter_file_guest
        self.angle_parameter_file_guest = angle_parameter_file_guest
        self.guest_qm_params_file = guest_qm_params_file
        self.reparameterised_intermediate_system_xml_file = (
            reparameterised_intermediate_system_xml_file
        )
        self.system_xml_non_bonded_file = system_xml_non_bonded_file
        self.system_xml_non_bonded_reparams_file = (
            system_xml_non_bonded_reparams_file
        )
        self.reparameterised_system_xml_file = reparameterised_system_xml_file
        self.reparameterized_torsional_params_file = (
            reparameterized_torsional_params_file
        )
        self.reparameterised_intermediate_torsional_system_xml_file = (
            reparameterised_intermediate_torsional_system_xml_file
        )
        self.reparameterised_torsional_system_xml_file = (
            reparameterised_torsional_system_xml_file
        )
        self.load_topology = load_topology
        self.non_reparameterised_system_xml_file = (
            non_reparameterised_system_xml_file
        )
        self.prmtop_system_non_params = prmtop_system_non_params
        self.inpcrd_system_non_params = inpcrd_system_non_params
        self.prmtop_system_intermediate_params = (
            prmtop_system_intermediate_params
        )
        self.inpcrd_system_intermediate_params = (
            inpcrd_system_intermediate_params
        )
        self.prmtop_system_params = prmtop_system_params
        self.inpcrd_system_params = inpcrd_system_params

    def generate_xml_from_prmtop(self):

        """
        Generates a serialsed XML forcefield file through parmed, given
        the PDB file and its corresponding topology file.
        """

        parm = parmed.load_file(self.prmtop_system, self.system_pdb)
        system = parm.createSystem()
        with open(self.system_xml, "w") as f:
            f.write(simtk.openmm.XmlSerializer.serialize(system))

    def write_guest_params_non_zero(self):

        """
        Saves the parameters of the ligand obtained from the QM log files
        in a text file starting from non-zero ( indexing begins from the
        index of the last atom of the receptor ).
        """

        # Charges from QM files
        df_charges = pd.read_csv(
            self.charge_parameter_file_guest, header=None, delimiter=r"\s+"
        )
        df_charges.columns = ["atom", "charges"]
        qm_charges = df_charges["charges"].values.tolist()
        qm_charges = [round(num, 6) for num in qm_charges]
        # print(qm_charges)
        # Bond Parameters from QM files
        ppdb = PandasPdb()
        ppdb.read_pdb(self.guest_qm_pdb)
        atom_name_list = ppdb.df["ATOM"]["atom_number"].values.tolist()
        # atom_name_list = [i - 1 for i in atom_name_list]
        no_host_atoms = base.get_num_host_atoms(self.host_pdb)
        atom_name_list = [i - 1 + no_host_atoms for i in atom_name_list]
        # print(atom_name_list)
        df = pd.read_csv(
            self.bond_parameter_file_guest, header=None, delimiter=r"\s+"
        )
        df.columns = ["bond", "k_bond", "bond_length", "bond_1", "bond_2"]
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
        ppdb.read_pdb(self.guest_qm_pdb)
        atom_name_list = ppdb.df["ATOM"]["atom_number"].values.tolist()
        # atom_name_list = [i - 1 for i in atom_name_list]
        no_host_atoms = base.get_num_host_atoms(self.host_pdb)
        atom_name_list = [i - 1 + no_host_atoms for i in atom_name_list]
        # print(atom_name_list)
        df = pd.read_csv(
            self.angle_parameter_file_guest, header=None, delimiter=r"\s+"
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
        xml = open(self.guest_qm_params_file, "w")
        xml.write("Begin writing the Bond Parameters" + "\n")
        # TODO: use string formatting and templates to write these lines
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
        Writes a reparameterised XML force field file for the
        system but without the QM obtained charges.
        """

        # Bond Parameters
        f_params = open(self.guest_qm_params_file, "r")
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
            comb_list_angle = [comb_1, comb_2, comb_3, comb_4, comb_5, comb_6]
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
        Writes a reparameterised XML force field file for the system.
        """

        # Bond Parameters
        with open(self.guest_qm_params_file, "r") as f_params:
            lines_params = f_params.readlines()
        # Bond Parameters
        for i in range(len(lines_params)):
            if "Begin writing the Bond Parameters" in lines_params[i]:
                to_begin = int(i)
            if "Finish writing the Bond Parameters" in lines_params[i]:
                to_end = int(i)
        bond_params = lines_params[to_begin + 1 : to_end]
        index_search_replace_bond = []
        # TODO: again, use string formatting.
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
        # TODO: use string formatting (generalize to function?)
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
            comb_list_angle = [comb_1, comb_2, comb_3, comb_4, comb_5, comb_6]
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
        f_params = open(self.guest_qm_params_file)
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

    def write_torsional_reparams_intermediate(self):
        """
        Generates a XML force field file for the system ( without the
        QM charges ) with reparameterized torsional parameters of the ligand.
        """

        no_host_atoms = base.get_num_host_atoms(self.host_pdb)
        xml_tor = open(self.reparameterized_torsional_params_file, "r")
        xml_tor_lines = xml_tor.readlines()
        xml_tor_lines_renum = []
        for i in xml_tor_lines:
            i = i.replace(
                "p1=" + '"' + str(int(re.findall("\d*\.?\d+", i)[2])) + '"',
                "p1="
                + '"'
                + str(int(int(re.findall("\d*\.?\d+", i)[2]) + no_host_atoms))
                + '"',
            )
            i = i.replace(
                "p2=" + '"' + str(int(re.findall("\d*\.?\d+", i)[4])) + '"',
                "p2="
                + '"'
                + str(int(int(re.findall("\d*\.?\d+", i)[4]) + no_host_atoms))
                + '"',
            )
            i = i.replace(
                "p3=" + '"' + str(int(re.findall("\d*\.?\d+", i)[6])) + '"',
                "p3="
                + '"'
                + str(int(int(re.findall("\d*\.?\d+", i)[6]) + no_host_atoms))
                + '"',
            )
            i = i.replace(
                "p4=" + '"' + str(int(re.findall("\d*\.?\d+", i)[8])) + '"',
                "p4="
                + '"'
                + str(int(int(re.findall("\d*\.?\d+", i)[8]) + no_host_atoms))
                + '"',
            )
            xml_tor_lines_renum.append(i)

        non_zero_k_tor = []
        for i in xml_tor_lines_renum:
            to_find = "k=" + '"' + "0.0" + '"'
            if to_find not in i:
                non_zero_k_tor.append(i)
        # print(non_zero_k_tor)
        p1 = []
        for i in range(len(non_zero_k_tor)):
            p1.append(int(re.findall("\d*\.?\d+", non_zero_k_tor[i])[2]))
        # print(p1)
        p2 = []
        for i in range(len(non_zero_k_tor)):
            p2.append(int(re.findall("\d*\.?\d+", non_zero_k_tor[i])[4]))
        # print(p2)
        p3 = []
        for i in range(len(non_zero_k_tor)):
            p3.append(int(re.findall("\d*\.?\d+", non_zero_k_tor[i])[6]))
        # print(p3)
        p4 = []
        for i in range(len(non_zero_k_tor)):
            p4.append(int(re.findall("\d*\.?\d+", non_zero_k_tor[i])[8]))
        # print(p4)
        periodicity = []
        for i in range(len(non_zero_k_tor)):
            periodicity.append(
                int(re.findall("\d*\.?\d+", non_zero_k_tor[i])[9])
            )
        # print(periodicity)
        xml_tor_reparams = open(
            self.reparameterised_intermediate_system_xml_file, "r"
        )
        xml_tor_reparams_lines = xml_tor_reparams.readlines()
        for j in range(len(xml_tor_reparams_lines)):
            for i in range(len(non_zero_k_tor)):
                to_find_tor = (
                    "p1="
                    + '"'
                    + str(p1[i])
                    + '"'
                    + " "
                    + "p2="
                    + '"'
                    + str(p2[i])
                    + '"'
                    + " "
                    + "p3="
                    + '"'
                    + str(p3[i])
                    + '"'
                    + " "
                    + "p4="
                    + '"'
                    + str(p4[i])
                    + '"'
                    + " "
                    + "periodicity="
                    + '"'
                    + str(periodicity[i])
                    + '"'
                )
                if to_find_tor in xml_tor_reparams_lines[j]:
                    print(xml_tor_reparams_lines[j])
                    xml_tor_reparams_lines[j] = non_zero_k_tor[i]
        with open(
            self.reparameterised_intermediate_torsional_system_xml_file, "w"
        ) as f:
            for i in xml_tor_reparams_lines:
                f.write(i)

    def write_torsional_reparams(self):
        """
        Generates a XML force field file for the system with reparameterized
        torsional parameters of the ligand.
        """

        no_host_atoms = base.get_num_host_atoms(self.host_pdb)
        with open(self.reparameterized_torsional_params_file, "r") as xml_tor:
            xml_tor_lines = xml_tor.readlines()
        xml_tor_lines_renum = []
        # TODO: string formatting and clean up this code to be more concise
        for i in xml_tor_lines:
            i = i.replace(
                "p1=" + '"' + str(int(re.findall("\d*\.?\d+", i)[2])) + '"',
                "p1="
                + '"'
                + str(int(int(re.findall("\d*\.?\d+", i)[2]) + no_host_atoms))
                + '"',
            )
            i = i.replace(
                "p2=" + '"' + str(int(re.findall("\d*\.?\d+", i)[4])) + '"',
                "p2="
                + '"'
                + str(int(int(re.findall("\d*\.?\d+", i)[4]) + no_host_atoms))
                + '"',
            )
            i = i.replace(
                "p3=" + '"' + str(int(re.findall("\d*\.?\d+", i)[6])) + '"',
                "p3="
                + '"'
                + str(int(int(re.findall("\d*\.?\d+", i)[6]) + no_host_atoms))
                + '"',
            )
            i = i.replace(
                "p4=" + '"' + str(int(re.findall("\d*\.?\d+", i)[8])) + '"',
                "p4="
                + '"'
                + str(int(int(re.findall("\d*\.?\d+", i)[8]) + no_host_atoms))
                + '"',
            )
            xml_tor_lines_renum.append(i)

        non_zero_k_tor = []
        for i in xml_tor_lines_renum:
            to_find = "k=" + '"' + "0.0" + '"'
            if to_find not in i:
                non_zero_k_tor.append(i)
        # print(non_zero_k_tor)
        p1 = []
        for i in range(len(non_zero_k_tor)):
            p1.append(int(re.findall("\d*\.?\d+", non_zero_k_tor[i])[2]))
        # print(p1)
        p2 = []
        for i in range(len(non_zero_k_tor)):
            p2.append(int(re.findall("\d*\.?\d+", non_zero_k_tor[i])[4]))
        # print(p2)
        p3 = []
        for i in range(len(non_zero_k_tor)):
            p3.append(int(re.findall("\d*\.?\d+", non_zero_k_tor[i])[6]))
        # print(p3)
        p4 = []
        for i in range(len(non_zero_k_tor)):
            p4.append(int(re.findall("\d*\.?\d+", non_zero_k_tor[i])[8]))
        # print(p4)
        periodicity = []
        for i in range(len(non_zero_k_tor)):
            periodicity.append(
                int(re.findall("\d*\.?\d+", non_zero_k_tor[i])[9])
            )
        # print(periodicity)
        xml_tor_reparams = open(self.reparameterised_system_xml_file, "r")
        xml_tor_reparams_lines = xml_tor_reparams.readlines()
        for j in range(len(xml_tor_reparams_lines)):
            for i in range(len(non_zero_k_tor)):
                to_find_tor = (
                    "p1="
                    + '"'
                    + str(p1[i])
                    + '"'
                    + " "
                    + "p2="
                    + '"'
                    + str(p2[i])
                    + '"'
                    + " "
                    + "p3="
                    + '"'
                    + str(p3[i])
                    + '"'
                    + " "
                    + "p4="
                    + '"'
                    + str(p4[i])
                    + '"'
                    + " "
                    + "periodicity="
                    + '"'
                    + str(periodicity[i])
                    + '"'
                )
                if to_find_tor in xml_tor_reparams_lines[j]:
                    print(xml_tor_reparams_lines[j])
                    xml_tor_reparams_lines[j] = non_zero_k_tor[i]
        with open(self.reparameterised_torsional_system_xml_file, "w") as f:
            for i in xml_tor_reparams_lines:
                f.write(i)

    def save_amber_params_non_qm_charges(self):

        """
        Saves amber generated topology files for the system
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
            self.prmtop_system_non_params, self.inpcrd_system_non_params
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
                    self.reparameterised_intermediate_torsional_system_xml_file
                ),
            )
        if self.load_topology == "openmm":
            openmm_system = parmed.openmm.load_topology(
                simtk.openmm.app.PDBFile(self.system_pdb).topology,
                parmed.load_file(
                    self.reparameterised_intermediate_torsional_system_xml_file
                ),
            )
        openmm_system.save(
            self.prmtop_system_intermediate_params, overwrite=True
        )
        openmm_system.coordinates = parmed.load_file(
            self.system_pdb, structure=True
        ).coordinates
        openmm_system.save(
            self.inpcrd_system_intermediate_params, overwrite=True
        )
        parm = parmed.load_file(
            self.prmtop_system_intermediate_params,
            self.inpcrd_system_intermediate_params,
        )
        xml_energy_decomposition = parmed.openmm.energy_decomposition_system(
            openmm_system,
            parmed.load_file(
                self.reparameterised_intermediate_torsional_system_xml_file
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
        Saves amber generated topology files for the system.
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
            self.prmtop_system_non_params, self.inpcrd_system_non_params
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
                    self.reparameterised_torsional_system_xml_file
                ),
            )
        if self.load_topology == "openmm":
            openmm_system = parmed.openmm.load_topology(
                simtk.openmm.app.PDBFile(self.system_pdb).topology,
                parmed.load_file(
                    self.reparameterised_torsional_system_xml_file
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
            parmed.load_file(self.reparameterised_torsional_system_xml_file),
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
