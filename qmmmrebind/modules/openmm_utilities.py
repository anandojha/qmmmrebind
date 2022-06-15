"""
Utilities which use OpenMM.
"""

from sys import stdout

# TODO: resolve and replace these imports
import simtk
from simtk.openmm import app
import parmed

import modules.base as base
import modules.torsion_drive_outputs as torsion_outputs
import modules.file_modify as file_modify

def get_openmm_energies(system_pdb, system_xml):

    """
    Returns decomposed OPENMM energies for the
    system.

    Parameters
    ----------
    system_pdb : str
        Input PDB file

    system_xml : str
        Forcefield file in XML format

    """

    pdb = simtk.openmm.app.PDBFile(system_pdb)
    ff_xml_file = open(system_xml, "r")
    system = simtk.openmm.XmlSerializer.deserialize(ff_xml_file.read())
    integrator = simtk.openmm.LangevinIntegrator(
        300 * simtk.unit.kelvin,
        1 / simtk.unit.picosecond,
        0.002 * simtk.unit.picoseconds,
    )
    simulation = simtk.openmm.app.Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    state = simulation.context.getState(
        getEnergy=True, getParameters=True, getForces=True
    )
    force_group = []
    for i, force in enumerate(system.getForces()):
        force_group.append(force.__class__.__name__)
    forcegroups = {}
    for i in range(system.getNumForces()):
        force = system.getForce(i)
        force.setForceGroup(i)
        forcegroups[force] = i
    energies = {}
    for f, i in forcegroups.items():
        energies[f] = (
            simulation.context.getState(getEnergy=True, groups=2 ** i)
            .getPotentialEnergy()
            ._value
        )
    decomposed_energy = []
    for key, val in energies.items():
        decomposed_energy.append(val)
    df_energy_openmm = pd.DataFrame(
        list(zip(force_group, decomposed_energy)),
        columns=["Energy_term", "Energy_openmm_params"],
    )
    energy_values = [
        list(
            df_energy_openmm.loc[
                df_energy_openmm["Energy_term"] == "HarmonicBondForce"
            ].values[0]
        )[1],
        list(
            df_energy_openmm.loc[
                df_energy_openmm["Energy_term"] == "HarmonicAngleForce"
            ].values[0]
        )[1],
        list(
            df_energy_openmm.loc[
                df_energy_openmm["Energy_term"] == "PeriodicTorsionForce"
            ].values[0]
        )[1],
        list(
            df_energy_openmm.loc[
                df_energy_openmm["Energy_term"] == "NonbondedForce"
            ].values[0]
        )[1],
    ]
    energy_group = [
        "HarmonicBondForce",
        "HarmonicAngleForce",
        "PeriodicTorsionForce",
        "NonbondedForce",
    ]
    df_energy_open_mm = pd.DataFrame(
        list(zip(energy_group, energy_values)),
        columns=["Energy_term", "Energy_openmm_params"],
    )
    df_energy_open_mm = df_energy_open_mm.set_index("Energy_term")
    print(df_energy_open_mm)

def get_non_torsion_mm_energy(system_pdb, load_topology, system_xml):

    """
    Returns sum of all the non-torsional energies (that
    includes HarmonicBondForce, HarmonicAngleForce
    and NonBondedForce) of the system from the PDB
    file given the topology and the forcefield file.

    Parameters
    ----------
    system_pdb : str
        System PDB file to load the openmm system topology
        and coordinates.

    load_topology : {"openmm", "parmed"}
        Argument to specify how to load the topology.

    system_xml : str
        XML force field file for the openmm system.

    Returns
    -------
    Sum of all the non-torsional energies of the system.

    """
    system_prmtop = system_pdb[:-4] + ".prmtop"
    system_inpcrd = system_pdb[:-4] + ".inpcrd"
    if load_topology == "parmed":
        openmm_system = parmed.openmm.load_topology(
            parmed.load_file(system_pdb, structure=True).topology,
            parmed.load_file(system_xml),
        )
    if load_topology == "openmm":
        openmm_system = parmed.openmm.load_topology(
            simtk.openmm.app.PDBFile(system_pdb).topology,
            parmed.load_file(system_xml),
        )
    openmm_system.save(system_prmtop, overwrite=True)
    openmm_system.coordinates = parmed.load_file(
        system_pdb, structure=True
    ).coordinates
    openmm_system.save(system_inpcrd, overwrite=True)
    parm = parmed.load_file(system_prmtop, system_inpcrd)
    prmtop_energy_decomposition = parmed.openmm.energy_decomposition_system(
        parm, parm.createSystem()
    )
    # print(prmtop_energy_decomposition)
    prmtop_energy_decomposition_value_no_torsion = [
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
        ).get("NonbondedForce"),
    ]
    return sum(prmtop_energy_decomposition_value_no_torsion)


def get_mm_potential_energies(qm_scan_file, load_topology, system_xml):

    """
    Returns potential energy of the system from the PDB file
    given the topology and the forcefield file.

    Parameters
    ----------
    qm_scan_file : str
        Output scan file containing torsiondrive scans.

    load_topology : {"openmm", "parmed"}
        Argument to spcify how to load the topology.

    system_xml : str
        XML file to load the openmm system.

    Returns
    -------
    mm_potential_energies : list
        List of all the non torsion mm energies for the
        generated PDB files.

    """
    mm_pdb_list = []
    for i in torsion_outputs.get_dihedrals(qm_scan_file):
        if i > 0:
            pdb_file = "plus_" + str(abs(i)) + ".pdb"
        if i < 0:
            pdb_file = "minus_" + str(abs(i)) + ".pdb"
        mm_pdb_list.append(pdb_file)
    for i in mm_pdb_list:
        mm_pdb_file = i
    mm_potential_energies = []
    for i in mm_pdb_list:
        mm_pdb_file = i
        mm_energy = get_non_torsion_mm_energy(
            system_pdb=i, load_topology=load_topology, system_xml=system_xml,
        )
        mm_potential_energies.append(mm_energy)
    return mm_potential_energies

def relax_init_structure(
    pdbfile,
    prmtopfile,
    qmmmrebindpdb,
    sim_output="output.pdb",
    sim_steps=100000,
):

    """
    Minimizing the initial PDB file with the given topology
    file

    Parameters
    ----------
    pdbfile: str
        Input PDB file.

    prmtopfile : str
        Input prmtop file.

    qmmmrebind_init_file: str
        Output PDB file.

    sim_output: str
        Simulation output trajectory file.

    sim_steps: int
        MD simulation steps.

    """

    prmtop = simtk.openmm.app.AmberPrmtopFile(prmtopfile)
    pdb = simtk.openmm.app.PDBFile(pdbfile)
    system = prmtop.createSystem(
        nonbondedMethod=simtk.openmm.app.PME,
        nonbondedCutoff=1 * simtk.unit.nanometer,
        constraints=simtk.openmm.app.HBonds,
    )
    integrator = simtk.openmm.LangevinIntegrator(
        300 * simtk.unit.kelvin,
        1 / simtk.unit.picosecond,
        0.002 * simtk.unit.picoseconds,
    )
    simulation = simtk.openmm.app.Simulation(
        prmtop.topology, system, integrator
    )
    simulation.context.setPositions(pdb.positions)
    print(simulation.context.getState(getEnergy=True).getPotentialEnergy())
    simulation.minimizeEnergy(maxIterations=10000000)
    print(simulation.context.getState(getEnergy=True).getPotentialEnergy())
    simulation.reporters.append(
        simtk.openmm.app.PDBReporter(sim_output, int(sim_steps / 10))
    )
    simulation.reporters.append(
        simtk.openmm.app.StateDataReporter(
            stdout,
            int(sim_steps / 10),
            step=True,
            potentialEnergy=True,
            temperature=True,
        )
    )
    simulation.reporters.append(
        simtk.openmm.app.PDBReporter(qmmmrebindpdb, sim_steps)
    )
    simulation.step(sim_steps)
    command = "rm -rf " + sim_output
    os.system(command)

def run_openmm_prmtop_inpcrd(
    pdbfile="system_qmmmrebind.pdb",
    prmtopfile="system_qmmmrebind.prmtop",
    inpcrdfile="system_qmmmrebind.inpcrd",
    sim_output="output.pdb",
    sim_steps=10000,
):

    """
    Runs OpenMM simulation with inpcrd and prmtop files.

    Parameters
    ----------
    pdbfile: str
       Input PDB file.

    prmtopfile: str
       Input prmtop file.

    inpcrdfile: str
       Input coordinate file.

    sim_output: str
       Output trajectory file.

    sim_steps: int
       Simulation steps.

    """

    prmtop = simtk.openmm.app.AmberPrmtopFile(prmtopfile)
    inpcrd = simtk.openmm.app.AmberInpcrdFile(inpcrdfile)
    system = prmtop.createSystem(
        nonbondedCutoff=1 * simtk.unit.nanometer,
        constraints=simtk.openmm.app.HBonds,
    )
    integrator = simtk.openmm.LangevinIntegrator(
        300 * simtk.unit.kelvin,
        1 / simtk.unit.picosecond,
        0.002 * simtk.unit.picoseconds,
    )
    simulation = simtk.openmm.app.Simulation(
        prmtop.topology, system, integrator
    )
    if inpcrd.boxVectors is None:
        file_modify.add_vectors_inpcrd(
            pdbfile=pdbfile, inpcrdfile=inpcrdfile,
        )
    if inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
        print(inpcrd.boxVectors)
    simulation.context.setPositions(inpcrd.positions)
    print(simulation.context.getState(getEnergy=True).getPotentialEnergy())
    simulation.minimizeEnergy(maxIterations=1000000)
    print(simulation.context.getState(getEnergy=True).getPotentialEnergy())
    simulation.reporters.append(
        simtk.openmm.app.PDBReporter(sim_output, int(sim_steps / 10))
    )
    simulation.reporters.append(
        simtk.openmm.app.StateDataReporter(
            stdout,
            int(sim_steps / 10),
            step=True,
            potentialEnergy=True,
            temperature=True,
        )
    )
    simulation.step(sim_steps)


def run_openmm_prmtop_pdb(
    pdbfile="system_qmmmrebind.pdb",
    prmtopfile="system_qmmmrebind.prmtop",
    sim_output="output.pdb",
    sim_steps=10000,
):

    """
    Runs OpenMM simulation with pdb and prmtop files.

    Parameters
    ----------
    pdbfile: str
       Input PDB file.

    prmtopfile: str
       Input prmtop file.

    sim_output: str
       Output trajectory file.

    sim_steps: int
       Simulation steps.

    """
    prmtop = simtk.openmm.app.AmberPrmtopFile(prmtopfile)
    pdb = simtk.openmm.app.PDBFile(pdbfile)
    system = prmtop.createSystem(
        nonbondedCutoff=1 * simtk.unit.nanometer,
        constraints=simtk.openmm.app.HBonds,
    )
    integrator = simtk.openmm.LangevinIntegrator(
        300 * simtk.unit.kelvin,
        1 / simtk.unit.picosecond,
        0.002 * simtk.unit.picoseconds,
    )
    simulation = simtk.openmm.app.Simulation(
        prmtop.topology, system, integrator
    )
    simulation.context.setPositions(pdb.positions)
    print(simulation.context.getState(getEnergy=True).getPotentialEnergy())
    simulation.minimizeEnergy(maxIterations=1000000)
    print(simulation.context.getState(getEnergy=True).getPotentialEnergy())
    simulation.reporters.append(
        simtk.openmm.app.PDBReporter(sim_output, int(sim_steps / 10))
    )
    simulation.reporters.append(
        simtk.openmm.app.StateDataReporter(
            stdout,
            int(sim_steps / 10),
            step=True,
            potentialEnergy=True,
            temperature=True,
        )
    )
    simulation.step(sim_steps)

