"""

"""

class TorsionDriveSims:

    """
    A class used to create a filetree for torsion scan
    using torsionsdrive for the dihedral angles of the ligand.

    This class creates a directory for carrying out torsiondrive
    calculations followed by fitting of torsional parameters. Methods
    in this class are used to run torsiondrive calculations either for
    all of the torsional angles, or for non-hydrogen / heavy atoms
    contributing to the torsional angle.

    ...

    Attributes
    ----------
    charge : int, optional
        Charge of the ligand.

    multiplicity: int, optional
        Spin Multiplicity (2S+1) of the ligand where S represents
        the total spin of the ligand.

    reparameterised_system_xml_file : str, optional
        Reparamaterixed XML force field for the ligand.

    torsion_xml_file : str, optional
        A text file containing torsional parameters from
        reparameterised XML file.

    xyz_file : str, optional
        XYZ file containing the coordinates of the guest molecule.

    psi_input_file : str, optional
        Input file for psi4 QM engine.

    memory : int, optional
        Memory (in GB) to be used.

    basis_set: str, optional
        Basis set to use for the QM engine.

    functional: str, optional
        Exchange/Correlation or hybrid Functional for the QM engine.

    iterations : int, optional
        Maximum number of geometry optimization steps.

    method_torsion_drive : str, optional
        The algorithm/package to use while running the torsiondrive
        scan. Using --native_opt uses QM program native constrained
        optimization algorithm and turns off geomeTRIC package.

    system_bonds_file : str, optional
        Text file containing bond parameters for the ligand.

    tor_dir : str, optional
        Torsiondrive directory containing separate torsiondrive
        folders, each containing files for a separate torsiondrive
        calculation for a particular dihedral angle.

    dihedral_text_file : str, optional
        Dihedral information file for torsiondrive.

    template_pdb : str, optional
        Guest PDB with atoms beginning from 1 to be used as a
        template PDB to retrieve atom indices and symbols.

    torsion_drive_run_file : str, optional
        bash file for torsiondrive calculations.

    dihedral_interval : int, optional
        Grid spacing for dihedral scan, i.e. every n degrees
        (where n is an integer), multiple values will be mapped
        to each dihedral angle.

    engine : str, optional
        Engine for running torsiondrive scan.

    energy_threshold : float, optional
        Only activate grid points if the new optimization is lower than
        the previous lowest energy (in a.u.).

    """

    def __init__(
        self,
        charge=0,
        multiplicity=1,
        reparameterised_system_xml_file="guest_reparameterised.xml",
        torsion_xml_file="guest_torsion_xml.txt",
        xyz_file="guest_coords.xyz",
        psi_input_file="torsion_drive_input.dat",
        memory=50,
        basis_set="6-31G",
        functional="B3LYP",
        iterations=2000,
        method_torsion_drive="native_opt",
        system_bonds_file="guest_bonds.txt",
        tor_dir="torsion_dir",
        dihedral_text_file="dihedrals.txt",
        template_pdb="guest_init_II.pdb",
        torsion_drive_run_file="run_command",
        dihedral_interval=15,
        engine="psi4",
        energy_threshold=0.00001,
    ):

        self.charge = charge
        self.multiplicity = multiplicity
        self.reparameterised_system_xml_file = reparameterised_system_xml_file
        self.torsion_xml_file = torsion_xml_file
        self.xyz_file = xyz_file
        self.psi_input_file = psi_input_file
        self.memory = memory
        self.basis_set = basis_set
        self.functional = functional
        self.iterations = iterations
        self.method_torsion_drive = method_torsion_drive
        self.system_bonds_file = system_bonds_file
        self.tor_dir = tor_dir
        self.dihedral_text_file = dihedral_text_file
        self.template_pdb = template_pdb
        self.torsion_drive_run_file = torsion_drive_run_file
        self.dihedral_interval = dihedral_interval
        self.engine = engine
        self.energy_threshold = energy_threshold

    def write_torsion_drive_run_file(self):
        """
        Saves a bash file for running torsion scans for torsiondrive.
        """
        if self.method_torsion_drive == "geometric":
            torsion_command = (
                "torsiondrive-launch"
                + " "
                + self.psi_input_file
                + " "
                + self.dihedral_text_file
                + " "
                + "-g"
                + " "
                + str(self.dihedral_interval)
                + " "
                + "-e"
                + " "
                + self.engine
                + " "
                + "--energy_thresh"
                + " "
                + str(self.energy_threshold)
                + " "
                + "-v"
            )
        if self.method_torsion_drive == "native_opt":
            torsion_command = (
                "torsiondrive-launch"
                + " "
                + self.psi_input_file
                + " "
                + self.dihedral_text_file
                + " "
                + "-g"
                + " "
                + str(self.dihedral_interval)
                + " "
                + "-e"
                + " "
                + self.engine
                + " "
                + "--energy_thresh"
                + " "
                + str(self.energy_threshold)
                + " "
                + "--"
                + self.method_torsion_drive
                + " "
                + "-v"
            )
        print(torsion_command)
        with open(self.torsion_drive_run_file, "w") as f:
            f.write(torsion_command)

    def write_tor_params_txt(self):
        """
        Saves a text file containing torsional parameters from the reparameterized XML
        force field file.
        """
        xml_off = open(self.reparameterised_system_xml_file, "r")
        xml_off_lines = xml_off.readlines()
        for i in range(len(xml_off_lines)):
            if "<Torsions>" in xml_off_lines[i]:
                to_begin = int(i)
            if "</Torsions>" in xml_off_lines[i]:
                to_end = int(i)
        torsion_params = xml_off_lines[to_begin + 1 : to_end]

        k_list_off = []
        for i in range(len(torsion_params)):
            k_list_off.append(
                float(re.findall("\d*\.?\d+", torsion_params[i])[0])
            )
        k_list_off = [round(num, 10) for num in k_list_off]
        # print(k_list_off)
        p1 = []
        for i in range(len(torsion_params)):
            p1.append(int(re.findall("\d*\.?\d+", torsion_params[i])[2]))
        p1 = [i + 1 for i in p1]
        # print(p1)
        p2 = []
        for i in range(len(torsion_params)):
            p2.append(int(re.findall("\d*\.?\d+", torsion_params[i])[4]))
        p2 = [i + 1 for i in p2]
        # print(p2)
        p3 = []
        for i in range(len(torsion_params)):
            p3.append(int(re.findall("\d*\.?\d+", torsion_params[i])[6]))
        p3 = [i + 1 for i in p3]
        # print(p3)
        p4 = []
        for i in range(len(torsion_params)):
            p4.append(int(re.findall("\d*\.?\d+", torsion_params[i])[8]))
        p4 = [i + 1 for i in p4]
        # print(p4)
        periodicity = []
        for i in range(len(torsion_params)):
            periodicity.append(
                int(re.findall("\d*\.?\d+", torsion_params[i])[9])
            )
        # print(periodicity)
        phase = []
        for i in range(len(torsion_params)):
            phase.append(float(re.findall("\d*\.?\d+", torsion_params[i])[10]))
        phase = [round(num, 8) for num in phase]
        # print(phase)
        data_tuples = list(zip(k_list_off, p1, p2, p3, p4, periodicity, phase))
        df_tor = pd.DataFrame(
            data_tuples,
            columns=["k", "p1", "p2", "p3", "p4", "periodicity", "phase",],
        )
        # print(df_tor.head())
        df_tor.to_csv(
            self.torsion_xml_file, index=False, header=False, sep=" "
        )

    def write_psi4_input(self):
        """
        Writes a psi4 input QM file.
        """
        xyz_lines = open(self.xyz_file, "r").readlines()[2:]
        with open(self.psi_input_file, "w") as f:
            f.write("memory" + " " + str(self.memory) + " " + "GB" + "\n")
            f.write("molecule" + " " + "{" + "\n")
            f.write(str(self.charge) + " " + str(self.multiplicity) + "\n")
            for line in xyz_lines:
                f.write(line)
            f.write("}" + "\n")
            f.write("set" + " " + "{" + "\n")
            f.write("basis" + " " + self.basis_set + "\n")
            if self.method_torsion_drive == "native_opt":
                f.write("GEOM_MAXITER" + " " + str(self.iterations) + "\n")
            f.write("}" + "\n")
            if self.method_torsion_drive == "native_opt":
                f.write(
                    "optimize" + "(" + "'" + self.functional + "'" ")" + "\n"
                )
            if self.method_torsion_drive == "geometric":
                f.write(
                    "gradient" + "(" + "'" + self.functional + "'" ")" + "\n"
                )

    def create_torsion_drive_dir(self):
        """
        Creates a directory for carrying out torsiondrive
        calculations for all the proper dihedral angles.
        """
        df_tor = pd.read_csv(
            self.torsion_xml_file, header=None, delimiter=r"\s+"
        )
        df_tor.columns = [
            "k",
            "p1",
            "p2",
            "p3",
            "p4",
            "periodicity",
            "phase",
        ]
        # print(df_tor.head())
        df_dihedrals = df_tor[["p1", "p2", "p3", "p4"]]
        # print(df_dihedrals.head())
        dihedrals_list_list = []
        for i in range(len(df_dihedrals)):
            dihedrals_list_list.append(df_dihedrals.iloc[i].values.tolist())
        set_list = set()
        unique_dihedrals_list_list = []
        for x in dihedrals_list_list:
            srtd = tuple(sorted(x))
            if srtd not in set_list:
                unique_dihedrals_list_list.append(x)
                set_list.add(srtd)
        # print(unique_dihedrals_list_list)
        os.system("rm -rf " + self.tor_dir)
        os.system("mkdir " + self.tor_dir)
        parent_cwd = os.getcwd()
        shutil.copy(
            parent_cwd + "/" + self.psi_input_file,
            parent_cwd + "/" + self.tor_dir + "/" + self.psi_input_file,
        )
        shutil.copy(
            parent_cwd + "/" + self.template_pdb,
            parent_cwd + "/" + self.tor_dir + "/" + self.template_pdb,
        )
        shutil.copy(
            parent_cwd + "/" + self.torsion_drive_run_file,
            parent_cwd
            + "/"
            + self.tor_dir
            + "/"
            + self.torsion_drive_run_file,
        )
        os.chdir(parent_cwd + "/" + self.tor_dir)
        torsion_drive_dir = os.getcwd()
        for i in range(len(unique_dihedrals_list_list)):
            dir_name = "torsion_drive" + "_" + str(i)
            os.system("rm -rf " + dir_name)
            os.system("mkdir " + dir_name)
            os.chdir(torsion_drive_dir + "/" + dir_name)
            with open(self.dihedral_text_file, "w") as f:
                f.write(
                    "# dihedral definition by atom indices starting from 1"
                    + "\n"
                )
                f.write("# i     j     k     l" + "\n")
                i_ = unique_dihedrals_list_list[i][0]
                j_ = unique_dihedrals_list_list[i][1]
                k_ = unique_dihedrals_list_list[i][2]
                l_ = unique_dihedrals_list_list[i][3]
                f.write(
                    " "
                    + "{:< 6d}".format(i_)
                    + "{:< 6d}".format(j_)
                    + "{:< 6d}".format(k_)
                    + "{:< 6d}".format(l_)
                    + "\n"
                )
                shutil.copy(
                    torsion_drive_dir + "/" + self.psi_input_file,
                    torsion_drive_dir
                    + "/"
                    + dir_name
                    + "/"
                    + self.psi_input_file,
                )
                shutil.copy(
                    torsion_drive_dir + "/" + self.template_pdb,
                    torsion_drive_dir
                    + "/"
                    + dir_name
                    + "/"
                    + self.template_pdb,
                )
                shutil.copy(
                    torsion_drive_dir + "/" + self.torsion_drive_run_file,
                    torsion_drive_dir
                    + "/"
                    + dir_name
                    + "/"
                    + self.torsion_drive_run_file,
                )
                os.chdir(torsion_drive_dir)
        os.system("rm -rf " + self.psi_input_file)
        os.system("rm -rf " + self.template_pdb)
        os.system("rm -rf " + self.torsion_drive_run_file)
        os.chdir(parent_cwd)

    def create_non_H_torsion_drive_dir(self):
        """
        Creates a directory for carrying out torsiondrive
        calculations for all non-hydrogen torsional angles.
        """
        df_tor = pd.read_csv(
            self.torsion_xml_file, header=None, delimiter=r"\s+"
        )
        df_tor.columns = [
            "k",
            "p1",
            "p2",
            "p3",
            "p4",
            "periodicity",
            "phase",
        ]
        # print(df_tor.head())
        ppdb = PandasPdb()
        ppdb.read_pdb(self.template_pdb)
        df_index_symbol = ppdb.df["ATOM"][["atom_number", "element_symbol"]]
        # print(df_index_symbol.head())
        df_dihedrals = df_tor[["p1", "p2", "p3", "p4"]]
        # print(df_dihedrals.head())
        dihedrals_list_list = []
        for i in range(len(df_dihedrals)):
            dihedrals_list_list.append(df_dihedrals.iloc[i].values.tolist())
        set_list = set()
        unique_dihedrals_list_list = []
        for x in dihedrals_list_list:
            srtd = tuple(sorted(x))
            if srtd not in set_list:
                unique_dihedrals_list_list.append(x)
                set_list.add(srtd)
        # print(unique_dihedrals_list_list)
        atom_dihedral_list = []
        for sub_list in unique_dihedrals_list_list:
            atom_dihedral_list.append(
                [
                    df_index_symbol.loc[
                        df_index_symbol["atom_number"] == sub_list[0]
                    ]["element_symbol"].to_list()[0],
                    df_index_symbol.loc[
                        df_index_symbol["atom_number"] == sub_list[1]
                    ]["element_symbol"].to_list()[0],
                    df_index_symbol.loc[
                        df_index_symbol["atom_number"] == sub_list[2]
                    ]["element_symbol"].to_list()[0],
                    df_index_symbol.loc[
                        df_index_symbol["atom_number"] == sub_list[3]
                    ]["element_symbol"].to_list()[0],
                ]
            )
        # print(atom_dihedral_list)
        index_to_include = []
        for i in range(len(atom_dihedral_list)):
            if "H" not in atom_dihedral_list[i]:
                index_to_include.append(i)
        non_H_dihedrals = []
        for i in index_to_include:
            non_H_dihedrals.append(unique_dihedrals_list_list[i])
        # print(non_H_dihedrals)
        non_H_atom_dihedral_list = []
        for sub_list in non_H_dihedrals:
            non_H_atom_dihedral_list.append(
                [
                    df_index_symbol.loc[
                        df_index_symbol["atom_number"] == sub_list[0]
                    ]["element_symbol"].to_list()[0],
                    df_index_symbol.loc[
                        df_index_symbol["atom_number"] == sub_list[1]
                    ]["element_symbol"].to_list()[0],
                    df_index_symbol.loc[
                        df_index_symbol["atom_number"] == sub_list[2]
                    ]["element_symbol"].to_list()[0],
                    df_index_symbol.loc[
                        df_index_symbol["atom_number"] == sub_list[3]
                    ]["element_symbol"].to_list()[0],
                ]
            )
        print(non_H_atom_dihedral_list)
        os.system("rm -rf " + self.tor_dir)
        os.system("mkdir " + self.tor_dir)
        parent_cwd = os.getcwd()
        shutil.copy(
            parent_cwd + "/" + self.psi_input_file,
            parent_cwd + "/" + self.tor_dir + "/" + self.psi_input_file,
        )
        shutil.copy(
            parent_cwd + "/" + self.template_pdb,
            parent_cwd + "/" + self.tor_dir + "/" + self.template_pdb,
        )
        shutil.copy(
            parent_cwd + "/" + self.torsion_drive_run_file,
            parent_cwd
            + "/"
            + self.tor_dir
            + "/"
            + self.torsion_drive_run_file,
        )
        os.chdir(parent_cwd + "/" + self.tor_dir)
        torsion_drive_dir = os.getcwd()
        for i in range(len(non_H_dihedrals)):
            dir_name = "torsion_drive" + "_" + str(i)
            os.system("rm -rf " + dir_name)
            os.system("mkdir " + dir_name)
            os.chdir(torsion_drive_dir + "/" + dir_name)
            with open(self.dihedral_text_file, "w") as f:
                f.write(
                    "# dihedral definition by atom indices starting from 1"
                    + "\n"
                )
                f.write("# i     j     k     l" + "\n")
                i_ = non_H_dihedrals[i][0]
                j_ = non_H_dihedrals[i][1]
                k_ = non_H_dihedrals[i][2]
                l_ = non_H_dihedrals[i][3]
                f.write(
                    " "
                    + "{:< 6d}".format(i_)
                    + "{:< 6d}".format(j_)
                    + "{:< 6d}".format(k_)
                    + "{:< 6d}".format(l_)
                    + "\n"
                )
                shutil.copy(
                    torsion_drive_dir + "/" + self.psi_input_file,
                    torsion_drive_dir
                    + "/"
                    + dir_name
                    + "/"
                    + self.psi_input_file,
                )
                shutil.copy(
                    torsion_drive_dir + "/" + self.template_pdb,
                    torsion_drive_dir
                    + "/"
                    + dir_name
                    + "/"
                    + self.template_pdb,
                )
                shutil.copy(
                    torsion_drive_dir + "/" + self.torsion_drive_run_file,
                    torsion_drive_dir
                    + "/"
                    + dir_name
                    + "/"
                    + self.torsion_drive_run_file,
                )
                os.chdir(torsion_drive_dir)
        os.system("rm -rf " + self.psi_input_file)
        os.system("rm -rf " + self.template_pdb)
        os.system("rm -rf " + self.torsion_drive_run_file)
        os.chdir(parent_cwd)

    def create_non_H_bonded_torsion_drive_dir(self):
        """
        Creates a directory for carrying out torsiondrive
        calculations for all non-hydrogen bonded torsional angles.
        """
        df_tor = pd.read_csv(
            self.torsion_xml_file, header=None, delimiter=r"\s+"
        )
        df_tor.columns = [
            "k",
            "p1",
            "p2",
            "p3",
            "p4",
            "periodicity",
            "phase",
        ]
        # print(df_tor.head())
        ppdb = PandasPdb()
        ppdb.read_pdb(self.template_pdb)
        df_index_symbol = ppdb.df["ATOM"][["atom_number", "element_symbol"]]
        # print(df_index_symbol.head())
        df_dihedrals = df_tor[["p1", "p2", "p3", "p4"]]
        # print(df_dihedrals.head())
        dihedrals_list_list = []
        for i in range(len(df_dihedrals)):
            dihedrals_list_list.append(df_dihedrals.iloc[i].values.tolist())
        set_list = set()
        unique_dihedrals_list_list = []
        for x in dihedrals_list_list:
            srtd = tuple(sorted(x))
            if srtd not in set_list:
                unique_dihedrals_list_list.append(x)
                set_list.add(srtd)
        # print(unique_dihedrals_list_list)
        atom_dihedral_list = []
        for sub_list in unique_dihedrals_list_list:
            atom_dihedral_list.append(
                [
                    df_index_symbol.loc[
                        df_index_symbol["atom_number"] == sub_list[0]
                    ]["element_symbol"].to_list()[0],
                    df_index_symbol.loc[
                        df_index_symbol["atom_number"] == sub_list[1]
                    ]["element_symbol"].to_list()[0],
                    df_index_symbol.loc[
                        df_index_symbol["atom_number"] == sub_list[2]
                    ]["element_symbol"].to_list()[0],
                    df_index_symbol.loc[
                        df_index_symbol["atom_number"] == sub_list[3]
                    ]["element_symbol"].to_list()[0],
                ]
            )
        # print(atom_dihedral_list)
        index_to_include = []
        for i in range(len(atom_dihedral_list)):
            if "H" not in atom_dihedral_list[i]:
                index_to_include.append(i)
        non_H_dihedrals = []
        for i in index_to_include:
            non_H_dihedrals.append(unique_dihedrals_list_list[i])
        # print(non_H_dihedrals)
        non_H_atom_dihedral_list = []
        for sub_list in non_H_dihedrals:
            non_H_atom_dihedral_list.append(
                [
                    df_index_symbol.loc[
                        df_index_symbol["atom_number"] == sub_list[0]
                    ]["element_symbol"].to_list()[0],
                    df_index_symbol.loc[
                        df_index_symbol["atom_number"] == sub_list[1]
                    ]["element_symbol"].to_list()[0],
                    df_index_symbol.loc[
                        df_index_symbol["atom_number"] == sub_list[2]
                    ]["element_symbol"].to_list()[0],
                    df_index_symbol.loc[
                        df_index_symbol["atom_number"] == sub_list[3]
                    ]["element_symbol"].to_list()[0],
                ]
            )
        # print(non_H_atom_dihedral_list)
        df_bonds_all = pd.read_csv(
            self.system_bonds_file, header=None, delimiter=r"\s+"
        )
        df_bonds_all.columns = [
            "bond_names",
            "k",
            "angle",
            "b1",
            "b2",
        ]
        df_bonds = df_bonds_all[["b1", "b2"]]
        bonds_list_list = []
        for i in range(len(df_bonds)):
            bonds_list_list.append(df_bonds.iloc[i].values.tolist())
        # print(bonds_list_list)
        reverse_bond_list_list = []
        for i in bonds_list_list:
            reverse_bond_list_list.append(base.reverse_list(i))
        # print(reverse_bond_list_list)
        bond_list = bonds_list_list + reverse_bond_list_list
        # print(bond_list)
        non_H_dihedral_bonds_list = []
        for i in non_H_dihedrals:
            non_H_dihedral_bonds_list.append(
                [[i[0], i[1]], [i[1], i[2]], [i[2], i[3]]]
            )
        # print(non_H_dihedral_bonds_list)
        bonded_index_to_include = []
        for i in range(len(non_H_dihedral_bonds_list)):
            if [
                non_H_dihedral_bonds_list[i][0] in bond_list,
                non_H_dihedral_bonds_list[i][1] in bond_list,
                non_H_dihedral_bonds_list[i][2] in bond_list,
            ] == [True, True, True]:
                bonded_index_to_include.append(i)
        # print(bonded_index_to_include)
        non_H_bonded_dihedrals = []
        for i in bonded_index_to_include:
            non_H_bonded_dihedrals.append(non_H_dihedrals[i])
        # print(non_H_bonded_dihedrals)
        non_H_bonded_atom_dihedral_list = []
        for sub_list in non_H_bonded_dihedrals:
            non_H_bonded_atom_dihedral_list.append(
                [
                    df_index_symbol.loc[
                        df_index_symbol["atom_number"] == sub_list[0]
                    ]["element_symbol"].to_list()[0],
                    df_index_symbol.loc[
                        df_index_symbol["atom_number"] == sub_list[1]
                    ]["element_symbol"].to_list()[0],
                    df_index_symbol.loc[
                        df_index_symbol["atom_number"] == sub_list[2]
                    ]["element_symbol"].to_list()[0],
                    df_index_symbol.loc[
                        df_index_symbol["atom_number"] == sub_list[3]
                    ]["element_symbol"].to_list()[0],
                ]
            )
        print(non_H_bonded_atom_dihedral_list)
        os.system("rm -rf " + self.tor_dir)
        os.system("mkdir " + self.tor_dir)
        parent_cwd = os.getcwd()
        shutil.copy(
            parent_cwd + "/" + self.psi_input_file,
            parent_cwd + "/" + self.tor_dir + "/" + self.psi_input_file,
        )
        shutil.copy(
            parent_cwd + "/" + self.template_pdb,
            parent_cwd + "/" + self.tor_dir + "/" + self.template_pdb,
        )
        shutil.copy(
            parent_cwd + "/" + self.torsion_drive_run_file,
            parent_cwd
            + "/"
            + self.tor_dir
            + "/"
            + self.torsion_drive_run_file,
        )
        os.chdir(parent_cwd + "/" + self.tor_dir)
        torsion_drive_dir = os.getcwd()
        for i in range(len(non_H_bonded_dihedrals)):
            dir_name = "torsion_drive" + "_" + str(i)
            os.system("rm -rf " + dir_name)
            os.system("mkdir " + dir_name)
            os.chdir(torsion_drive_dir + "/" + dir_name)
            with open(self.dihedral_text_file, "w") as f:
                f.write(
                    "# dihedral definition by atom indices starting from 1"
                    + "\n"
                )
                f.write("# i     j     k     l" + "\n")
                i_ = non_H_bonded_dihedrals[i][0]
                j_ = non_H_bonded_dihedrals[i][1]
                k_ = non_H_bonded_dihedrals[i][2]
                l_ = non_H_bonded_dihedrals[i][3]
                f.write(
                    " "
                    + "{:< 6d}".format(i_)
                    + "{:< 6d}".format(j_)
                    + "{:< 6d}".format(k_)
                    + "{:< 6d}".format(l_)
                    + "\n"
                )
                shutil.copy(
                    torsion_drive_dir + "/" + self.psi_input_file,
                    torsion_drive_dir
                    + "/"
                    + dir_name
                    + "/"
                    + self.psi_input_file,
                )
                shutil.copy(
                    torsion_drive_dir + "/" + self.template_pdb,
                    torsion_drive_dir
                    + "/"
                    + dir_name
                    + "/"
                    + self.template_pdb,
                )
                shutil.copy(
                    torsion_drive_dir + "/" + self.torsion_drive_run_file,
                    torsion_drive_dir
                    + "/"
                    + dir_name
                    + "/"
                    + self.torsion_drive_run_file,
                )
                os.chdir(torsion_drive_dir)
        os.system("rm -rf " + self.psi_input_file)
        os.system("rm -rf " + self.template_pdb)
        os.system("rm -rf " + self.torsion_drive_run_file)
        os.chdir(parent_cwd)

    def run_torsion_sim(self):
        """
        Run torsion scans using torsiondrive locally.
        """
        parent_cwd = os.getcwd()
        target_dir = parent_cwd + "/" + self.tor_dir
        num_folders = 0
        for _, dirnames, filenames in os.walk(target_dir):
            num_folders += len(dirnames)
        for i in range(num_folders):
            dir_ = "torsion_drive" + "_" + str(i)
            os.chdir(parent_cwd + "/" + self.tor_dir + "/" + dir_)
            run_command = "bash" + " " + self.torsion_drive_run_file
            # os.system(run_command)
            print(run_command)
            os.chdir(parent_cwd)
