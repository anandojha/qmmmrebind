"""

"""

from biopandas.pdb import PandasPdb
import pandas as pd

class PrepareGaussianGuest:

    """
    A class used to prepare the QM engine input file (Gaussian)
    for the ligand and run QM calculations with appropriate
    keywords.

    This class contain methods to write an input file (.com extension)
    for the QM engine. It then runs a QM calculation with the given
    basis set and functional. Checkpoint file is then converted to
    a formatted checkpoint file. Output files (.log, .chk, and .fhck)
    will then be used to extract ligand's force field parameters.

    ...

    Attributes
    ----------
    charge : int, optional
        Charge of the ligand.

    multiplicity: int, optional
        Spin Multiplicity (2S+1) of the ligand where S represents
        the total spin of the ligand.

    guest_pdb: str, optional
        Ligand PDB file with atom numbers beginning from 1.

    n_processors : int, optional
        Number of processors to be used for Gaussian program to run and
        set in %NProcShared command of Gaussian.

    memory : int, optional
        Memory (in GB) to be used set in %Mem command of Gaussian.

    functional: str, optional
        Exchange/Correlation or hybrid functional to use in the Gaussian
        QM calculation.

    basis_set: str, optional
        Basis set to use for the Gaussian QM calculation.

    optimisation: str, optional
        set to "OPT" to perform a geometry optimization on the ligand
        specified in the system; else set to an empty string.

    frequency: str, optional
        set to "FREQ" for Gaussian to perform a frequency calculation;
        else set to an empty string.

    add_keywords_I: str, optional
        Specifies the integration grid.

    add_keywords_II: str, optional
        Specifies the QM engine to select one of the methods for
        analyzing the electron density of the system. Methods used
        are based on fitting the molecular electrostatic potential.
        Methods used are : POP=CHELPG (Charges from Electrostatic
        Potentials using a Grid based method) and POP=MK
        (Merz-Singh-Kollman scheme)

    add_keywords_III: str, optional
        Used to include the IOp keyword (to set the internal options to
        specific values) in the Gaussian command.

    gauss_out_file: str, optional
        This file contains the output script obtained after running
        the Gaussian QM calculation.

    fchk_out_file: str, optional
        Formatted checkpoint file obtained from the checkpoint file
        using formchk command.


    """

    def __init__(
        self,
        charge=0,
        multiplicity=1,
        guest_pdb="guest_init_II.pdb",
        n_processors=12,
        memory=50,
        functional="B3LYP",
        basis_set="6-31G",
        optimisation="OPT",
        frequency="FREQ",
        add_keywords_I="INTEGRAL=(GRID=ULTRAFINE)",
        add_keywords_II="POP(MK,READRADII)",
        add_keywords_III="IOP(6/33=2,6/42=6)",
        gauss_out_file="guest.out",
        fchk_out_file="guest_fchk.out",
    ):

        self.charge = charge
        self.multiplicity = multiplicity
        self.guest_pdb = guest_pdb
        self.n_processors = n_processors
        self.memory = memory
        self.functional = functional
        self.basis_set = basis_set
        self.optimisation = optimisation
        self.frequency = frequency
        self.gauss_out_file = gauss_out_file
        self.fchk_out_file = fchk_out_file
        self.add_keywords_I = add_keywords_I
        self.add_keywords_II = add_keywords_II
        self.add_keywords_III = add_keywords_III

    def write_input(self):
        """
        Writes a Gaussian input file for the ligand.
        """

        command_line_1 = "%Chk = " + self.guest_pdb[:-4] + ".chk"
        command_line_2 = "%Mem = " + str(self.memory) + "GB"
        command_line_3 = "%NProcShared = " + str(self.n_processors)
        command_line_4 = (
            "# "
            + self.functional
            + " "
            + self.basis_set
            + " "
            + self.optimisation
            + " "
            + self.frequency
            + " "
            + self.add_keywords_I
            + " "
            + self.add_keywords_II
            + " "
            + self.add_keywords_III
        )
        command_line_5 = " "
        command_line_6 = self.guest_pdb[:-4] + " " + "gaussian input file"
        command_line_7 = " "
        command_line_8 = str(self.charge) + " " + str(self.multiplicity)
        ppdb = PandasPdb()
        ppdb.read_pdb(self.guest_pdb)
        df = ppdb.df["ATOM"]
        df_1 = ppdb.df["ATOM"]["element_symbol"]
        df_1.columns = ["atom"]
        df_2 = df[["x_coord", "y_coord", "z_coord"]]
        df_merged = pd.concat([df_1, df_2], axis=1)
        command_line_9 = df_merged.to_string(header=False, index=False)
        command_line_10 = " "
        command = [
            command_line_1,
            command_line_2,
            command_line_3,
            command_line_4,
            command_line_5,
            command_line_6,
            command_line_7,
            command_line_8,
            command_line_9,
            command_line_10,
        ]
        commands = "\n".join(command)
        with open(self.guest_pdb[:-4] + ".com", "w") as f:
            f.write(commands)

    def run_gaussian(self):
        """
        Runs the Gaussian QM calculation for the ligand locally.
        """
        execute_command = (
            "g16"
            + " < "
            + self.guest_pdb[:-4]
            + ".com"
            + " > "
            + self.guest_pdb[:-4]
            + ".log"
        )
        with open(self.gauss_out_file, "w+") as f:
            sp.run(
                execute_command, shell=True, stdout=f, stderr=sp.STDOUT,
            )

    def get_fchk(self):
        """
        Converts the Gaussian checkpoint file (.chk) to a formatted checkpoint
        file (.fchk).
        """
        execute_command = (
            "formchk"
            + " "
            + self.guest_pdb[:-4]
            + ".chk"
            + " "
            + self.guest_pdb[:-4]
            + ".fchk"
        )
        with open(self.fchk_out_file, "w+") as f:
            sp.run(
                execute_command, shell=True, stdout=f, stderr=sp.STDOUT,
            )

class PrepareGaussianHostGuest:

    """
    A class used to prepare the QM engine input file (Gaussian) for
    the receptor - ligand complex and run the QM calculations with
    the appropriate keywords.

    This class contain methods to write an input file (.com extension)
    for the QM engine for the receptor - ligand complex. It then runs
    a QM calculation with the given basis set and functional. Checkpoint
    file is then converted to a formatted checkpoint file. Output files
    (.log, .chk, and .fhck) will then be used to extract charges for the
    ligand and the receptor.

    ...

    Attributes
    ----------
    charge : int, optional
        Total charge of the receptor - ligand complex.

    multiplicity : int, optional
        Spin Multiplicity (2S+1) of the ligand where S represents
        the total spin of the ligand.

    guest_pdb : str, optional
        Ligand PDB file with atom numbers beginning from 1.

    host_qm_pdb : str, optional
        PDB file for the receptor's QM region.

    n_processors : int, optional
        Number of processors to be used for Gaussian program to run and
        set in %NProcShared command of Gaussian.

    memory : int, optional
        Memory (in GB) to be used set in %Mem command of Gaussian.

    functional: str, optional
        Exchange/Correlation or hybrid functional to use in the Gaussian
        QM calculation.

    basis_set: str, optional
        Basis set to use for the Gaussian QM calculation.

    optimisation: str, optional
        set to "OPT" to perform a geometry optimization on the ligand
        specified in the system; else set to an empty string.

    frequency: str, optional
        set to "FREQ" for Gaussian to perform a frequency calculation;
        else set to an empty string.

    add_keywords_I: str, optional
        Specifies the integration grid.

    add_keywords_II: str, optional
        Specifies the QM engine to select one of the methods for
        analyzing the electron density of the system. Methods used
        are based on fitting the molecular electrostatic potential.
        Methods used are : POP=CHELPG (Charges from Electrostatic
        Potentials using a Grid based method) and POP=MK
        (Merz-Singh-Kollman scheme)

    add_keywords_III: str, optional
        Used to include the IOp keyword (to set the internal options to
        specific values) in the Gaussian command.

    gauss_system_out_file : str, optional
        This file contains the output script obtained after running
        the Gaussian QM calculation.

    fchk_system_out_file : str, optional
        Formatted checkpoint file obtained from the checkpoint file
        using formchk command.

    host_guest_input : str, optional
        Gaussian input file (.com extension) for the receptor - ligand
        QM region.

    qm_guest_charge_parameter_file : str, optional
        File containing the charges of ligand atoms and their corresponding
        atoms. Charge obtained are the polarised charged due to the
        surrounding receptor's region.

    qm_host_charge_parameter_file : str, optional
        File containing the charges of the QM region of the receptor.

    qm_guest_atom_charge_parameter_file : str, optional
        File containing the charges of ligand atoms. Charge obtained
        are the polarised charged due to the surrounding receptor's region.

    """

    def __init__(
        self,
        charge=0,
        multiplicity=1,
        guest_pdb="guest_init_II.pdb",
        host_qm_pdb="host_qm.pdb",
        n_processors=12,
        memory=50,
        functional="B3LYP",
        basis_set="6-31G",
        optimisation="",
        frequency="",
        add_keywords_I="INTEGRAL=(GRID=ULTRAFINE)",
        add_keywords_II="POP(MK,READRADII)",
        add_keywords_III="IOP(6/33=2,6/42=6) SCRF=PCM",
        gauss_system_out_file="system_qm.out",
        fchk_system_out_file="system_qm_fchk.out",
        host_guest_input="host_guest.com",
        qm_guest_charge_parameter_file="guest_qm_surround_charges.txt",
        qm_host_charge_parameter_file="host_qm_surround_charges.txt",
        qm_guest_atom_charge_parameter_file="guest_qm_atom_surround_charges.txt",
    ):

        self.charge = charge
        self.multiplicity = multiplicity
        self.guest_pdb = guest_pdb
        self.host_qm_pdb = host_qm_pdb
        self.n_processors = n_processors
        self.memory = memory
        self.functional = functional
        self.basis_set = basis_set
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
        self.qm_guest_atom_charge_parameter_file = (
            qm_guest_atom_charge_parameter_file
        )

    def write_input(self):
        """
        Writes a Gaussian input file for the receptor - ligand QM region.
        """
        command_line_1 = "%Chk = " + self.host_guest_input[:-4] + ".chk"
        command_line_2 = "%Mem = " + str(self.memory) + "GB"
        command_line_3 = "%NProcShared = " + str(self.n_processors)
        command_line_4 = (
            "# "
            + self.functional
            + " "
            + self.basis_set
            + " "
            + self.optimisation
            + " "
            + self.frequency
            + " "
            + self.add_keywords_I
            + " "
            + self.add_keywords_II
            + " "
            + self.add_keywords_III
        )
        command_line_5 = " "
        command_line_6 = "Gaussian Input File"
        command_line_7 = " "
        command_line_8 = str(self.charge) + " " + str(self.multiplicity)
        ppdb = PandasPdb()
        ppdb.read_pdb(self.guest_pdb)
        df = ppdb.df["ATOM"]
        df_1 = ppdb.df["ATOM"]["element_symbol"]
        df_1.columns = ["atom"]
        df_3 = df[["x_coord", "y_coord", "z_coord"]]
        df_2 = pd.Series(["0"] * len(df), name="decide_freeze")
        df_merged_1 = pd.concat([df_1, df_2, df_3], axis=1)
        ppdb = PandasPdb()
        ppdb.read_pdb(self.host_qm_pdb)
        df = ppdb.df["ATOM"]
        df_1 = ppdb.df["ATOM"]["element_symbol"]
        df_1.columns = ["atom"]
        df_3 = df[["x_coord", "y_coord", "z_coord"]]
        df_2 = pd.Series(["0"] * len(df), name="decide_freeze")
        df_merged_2 = pd.concat([df_1, df_2, df_3], axis=1)
        df_merged = pd.concat([df_merged_1, df_merged_2], axis=0)
        command_line_9 = df_merged.to_string(header=False, index=False)
        command_line_10 = " "
        command = [
            command_line_1,
            command_line_2,
            command_line_3,
            command_line_4,
            command_line_5,
            command_line_6,
            command_line_7,
            command_line_8,
            command_line_9,
            command_line_10,
        ]
        commands = "\n".join(command)

        with open(self.host_guest_input, "w") as f:
            f.write(commands)

    def run_gaussian(self):
        """
        Runs the Gaussian QM calculation for the ligand - receptor region
        locally.
        """
        execute_command = (
            "g16"
            + " < "
            + self.host_guest_input
            + " > "
            + self.host_guest_input[:-4]
            + ".log"
        )
        with open(self.gauss_system_out_file, "w+") as f:
            sp.run(
                execute_command, shell=True, stdout=f, stderr=sp.STDOUT,
            )

    def get_fchk(self):
        """
        Converts the Gaussian checkpoint file (.chk) to a formatted checkpoint
        file (.fchk).
        """
        execute_command = (
            "formchk"
            + " "
            + self.host_guest_input[:-4]
            + ".chk"
            + " "
            + self.host_guest_input[:-4]
            + ".fchk"
        )
        with open(self.fchk_system_out_file, "w+") as f:
            sp.run(
                execute_command, shell=True, stdout=f, stderr=sp.STDOUT,
            )

    def get_qm_host_guest_charges(self):
        """
        Extract charge information for the receptor - ligand QM region.
        """
        log_file = self.host_guest_input[:-4] + ".log"
        with open(log_file, "r") as f:
            lines = f.readlines()
        for i in range(len(lines)):
            if "Fitting point charges to electrostatic potential" in lines[i]:
                to_begin = int(i)
            if " Sum of ESP charges =" in lines[i]:
                to_end = int(i)
        
        # Why + 4?
        charges = lines[to_begin + 4 : to_end]
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
        df_guest = ppdb.df["ATOM"]
        number_guest_atoms = df_guest.shape[0]
        data_tuples = list(zip(atom_list, charge_list_value))
        df_charge = pd.DataFrame(data_tuples, columns=["Atom", "Charge"])
        number_host_atoms = df_charge.shape[0] - number_guest_atoms
        df_charge_guest = df_charge.head(number_guest_atoms)
        df_charge_host = df_charge.tail(number_host_atoms)
        df_charge_only_guest = df_charge_guest["Charge"]
        df_charge_guest.to_csv(
            self.qm_guest_charge_parameter_file,
            index=False,
            header=False,
            sep=" ",
        )
        df_charge_host.to_csv(
            self.qm_host_charge_parameter_file,
            index=False,
            header=False,
            sep=" ",
        )
        df_charge_only_guest.to_csv(
            self.qm_guest_atom_charge_parameter_file,
            index=False,
            header=False,
            sep=" ",
        )

class PrepareGaussianHost:

    """
    A class used to prepare the QM engine input file (Gaussian)
    for the receptor and run QM calculations with appropriate keywords.

    This class contain methods to write an input file (.com extension)
    for the QM engine. It then runs a QM calculation with the given
    basis set and functional. Checkpoint file is then converted to
    a formatted checkpoint file. Output files (.log, .chk, and .fhck)
    will then be used to extract receptors's force field parameters.

    ...

    Attributes
    ----------
    charge : int, optional
        Charge of the receptor.

    multiplicity: int, optional
        Spin Multiplicity (2S+1) of the receptor where S represents
        the total spin of the receptor.

    host_qm_pdb: str, optional
        PDB file of the receptor's QM region with atom numbers
        beginning from 1.

    n_processors : int, optional
        Number of processors to be used for Gaussian program to run and
        set in %NProcShared command of Gaussian.

    memory : int, optional
        Memory (in GB) to be used set in %Mem command of Gaussian.

    functional: str, optional
        Exchange/Correlation or hybrid functional to use in the Gaussian
        QM calculation.

    basis_set: str, optional
        Basis set to use for the Gaussian QM calculation.

    optimisation: str, optional
        set to "OPT" to perform a geometry optimization on the receptor
        specified in the system; else set to an empty string.

    frequency: str, optional
        set to "FREQ" for Gaussian to perform a frequency calculation;
        else set to an empty string.

    add_keywords_I: str, optional
        Specifies the integration grid.

    add_keywords_II: str, optional
        Specifies the QM engine to select one of the methods for
        analyzing the electron density of the system. Methods used
        are based on fitting the molecular electrostatic potential.
        Methods used are : POP=CHELPG (Charges from Electrostatic
        Potentials using a Grid based method) and POP=MK
        (Merz-Singh-Kollman scheme)

    add_keywords_III: str, optional
        Used to include the IOp keyword (to set the internal options to
        specific values) in the Gaussian command.

    gauss_out_file: str, optional
        This file contains the output script obtained after running
        the Gaussian QM calculation.

    fchk_out_file: str, optional
        Formatted checkpoint file obtained from the checkpoint file
        using formchk command.

    """

    def __init__(
        self,
        charge=0,
        multiplicity=1,
        host_qm_pdb="host_qm.pdb",
        n_processors=12,
        memory=50,
        functional="B3LYP",
        basis_set="6-31G",
        optimisation="OPT",
        frequency="FREQ",
        add_keywords_I="INTEGRAL=(GRID=ULTRAFINE) SCF=(maxcycles=4000) SYMMETRY=NONE",
        add_keywords_II="POP(MK,READRADII)",
        add_keywords_III="IOP(6/33=2,6/42=6)",
        gauss_out_file="host_qm.out",
        fchk_out_file="host_qm_fchk.out",
    ):

        self.charge = charge
        self.multiplicity = multiplicity
        self.host_qm_pdb = host_qm_pdb
        self.n_processors = n_processors
        self.memory = memory
        self.functional = functional
        self.basis_set = basis_set
        self.optimisation = optimisation
        self.frequency = frequency
        self.gauss_out_file = gauss_out_file
        self.fchk_out_file = fchk_out_file
        self.add_keywords_I = add_keywords_I
        self.add_keywords_II = add_keywords_II
        self.add_keywords_III = add_keywords_III

    def write_input(self):
        """
        Writes a Gaussian input file for the receptor QM region.
        """
        # TODO: create generic function for Gaussian Input file (DRY principle)
        command_line_1 = "%Chk = " + self.host_qm_pdb[:-4] + ".chk"
        command_line_2 = "%Mem = " + str(self.memory) + "GB"
        command_line_3 = "%NProcShared = " + str(self.n_processors)
        command_line_4 = (
            "# "
            + self.functional
            + " "
            + self.basis_set
            + " "
            + self.optimisation
            + " "
            + self.frequency
            + " "
            + self.add_keywords_I
            + " "
            + self.add_keywords_II
            + " "
            + self.add_keywords_III
        )
        command_line_5 = " "
        command_line_6 = self.host_qm_pdb[:-4] + " " + "gaussian input file"
        command_line_7 = " "
        command_line_8 = str(self.charge) + " " + str(self.multiplicity)
        ppdb = PandasPdb()
        ppdb.read_pdb(self.host_qm_pdb)
        df = ppdb.df["ATOM"]
        df_1 = ppdb.df["ATOM"]["element_symbol"]
        df_1.columns = ["atom"]
        df_2 = df[["x_coord", "y_coord", "z_coord"]]
        df_merged = pd.concat([df_1, df_2], axis=1)
        command_line_9 = df_merged.to_string(header=False, index=False)
        command_line_10 = " "
        command = [
            command_line_1,
            command_line_2,
            command_line_3,
            command_line_4,
            command_line_5,
            command_line_6,
            command_line_7,
            command_line_8,
            command_line_9,
            command_line_10,
        ]
        commands = "\n".join(command)
        with open(self.host_qm_pdb[:-4] + ".com", "w") as f:
            f.write(commands)

    def run_gaussian(self):
        """
        Runs the Gaussian QM calculation for the receptor locally.
        """
        execute_command = (
            "g16"
            + " < "
            + self.host_qm_pdb[:-4]
            + ".com"
            + " > "
            + self.host_qm_pdb[:-4]
            + ".log"
        )
        with open(self.gauss_out_file, "w+") as f:
            sp.run(
                execute_command, shell=True, stdout=f, stderr=sp.STDOUT,
            )

    def get_fchk(self):
        """
        Converts the Gaussian checkpoint file (.chk) to a formatted checkpoint
        file (.fchk).
        """
        execute_command = (
            "formchk"
            + " "
            + self.host_qm_pdb[:-4]
            + ".chk"
            + " "
            + self.host_qm_pdb[:-4]
            + ".fchk"
        )
        with open(self.fchk_out_file, "w+") as f:
            sp.run(
                execute_command, shell=True, stdout=f, stderr=sp.STDOUT,
            )


