"""

"""

import os

from biopandas.pdb import PandasPdb
import numpy as np
import statistics
import itertools

class PrepareQMMM:

    """
    A class used to segregate the QM and MM regions.

    This class contain methods to remove the solvent, ions and all
    entities that are exclusive of receptor and the ligand. It also
    defines the Quantum Mechanical (QM) region and the Molecular
    Mechanical (MM) region based upon the distance of the ligand
    from the receptor and the chosen number of receptor residues. It
    is also assumed that the initial PDB file will have the receptor
    followed by the ligand.

    ...

    Attributes
    ----------
    init_pdb : str
        Initial PDB file containing the receptor-ligand complex with
        solvent, ions, etc.

    cleaned_pdb : str
        Formatted PDB file containing only the receptor and the ligand.

    guest_init_pdb : str
        A separate ligand PDB file with atom numbers not beginning from 1.

    host_pdb : str
        A separate receptor PDB file with atom numbers beginning from 1.

    guest_resname : str
        Three letter residue ID for the ligand.

    guest_pdb : str, optional
        Ligand PDB file with atom numbers beginning from 1.

    guest_xyz : str, optional
        A text file of the XYZ coordinates of the ligand.

    distance : float, optional
        The distance required to define the QM region of the receptor.
        This is the distance between the atoms of the ligand and the
        atoms of the receptor.

    residue_list : str, optional
        A text file of the residue numbers of the receptor within the
        proximity (as defined by the distance) from the ligand.

    host_qm_atoms : str, optional
        A text file of the atom numbers of the receptors in the QM
        region.

    host_mm_atoms : str, optional
        A text file of the atom numbers of the receptors in the MM
        region (all atoms except atoms in the QM region)

    host_qm_pdb : str, optional
        PDB file for the receptor's QM region.

    host_mm_pdb : str, optional
        PDB file for the receptor's MM region.

    qm_pdb : str, optional
        PDB file for the QM region (receptor's QM region and the
        ligand).

    mm_pdb : str, optional
        PDB file for the MM region.

    host_mm_region_I_atoms : str, optional
        A text file of the atom numbers of the receptors in the MM
        region preceeding the QM region.

    host_mm_region_II_atoms : str, optional
        A text file of the atom numbers of the receptors in the MM
        region following the QM region.

    host_mm_region_I_pdb : str, optional
        PDB file of the receptor in the MM region preceeding the
        QM region.

    host_mm_region_II_pdb : str, optional
        PDB file of the receptor in the MM region following the
        QM region.

    num_residues : int, optional
        Number of residues required in the QM region of the receptor.
    """

    def __init__(
        self,
        init_pdb,
        distance,
        num_residues,
        guest_resname,
        cleaned_pdb="system.pdb",
        guest_init_pdb="guest_init.pdb",
        host_pdb="host.pdb",
        guest_pdb="guest_init_II.pdb",
        guest_xyz="guest_coord.txt",
        residue_list="residue_list.txt",
        host_qm_atoms="host_qm.txt",
        host_mm_atoms="host_mm.txt",
        host_qm_pdb="host_qm.pdb",
        host_mm_pdb="host_mm.pdb",
        qm_pdb="qm.pdb",
        mm_pdb="mm.pdb",
        host_mm_region_I_atoms="host_mm_region_I.txt",
        host_mm_region_II_atoms="host_mm_region_II.txt",
        host_mm_region_I_pdb="host_mm_region_I.pdb",
        host_mm_region_II_pdb="host_mm_region_II.pdb",
    ):

        self.init_pdb = init_pdb
        self.distance = distance
        self.num_residues = num_residues
        self.guest_resname = guest_resname
        self.cleaned_pdb = cleaned_pdb
        self.guest_init_pdb = guest_init_pdb
        self.host_pdb = host_pdb
        self.guest_pdb = guest_pdb
        self.guest_xyz = guest_xyz
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

    def clean_up(self):
        """
        Reads the given PDB file, removes all entities except the
        receptor and ligand and saves a new pdb file.
        """
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
        intermediate_file_1 = self.cleaned_pdb[:-4] + "_intermediate_1.pdb"
        intermediate_file_2 = self.cleaned_pdb[:-4] + "_intermediate_2.pdb"
        command = (
            "pdb4amber -i "
            + self.init_pdb
            + " -o "
            + intermediate_file_1
            + " --noter --dry"
        )
        os.system(command)
        to_delete = (
            intermediate_file_1[:-4] + "_nonprot.pdb",
            intermediate_file_1[:-4] + "_renum.txt",
            intermediate_file_1[:-4] + "_sslink",
            intermediate_file_1[:-4] + "_water.pdb",
        )
        os.system("rm -rf " + " ".join(to_delete))
        with open(intermediate_file_1) as f1, open(
                intermediate_file_2, "w") as f2:
            for line in f1:
                if not any(ion in line for ion in ions):
                    f2.write(line)
        with open(intermediate_file_2, "r") as f1:
            filedata = f1.read()
        filedata = filedata.replace("HETATM", "ATOM  ")
        with open(self.cleaned_pdb, "w") as f2:
            f2.write(filedata)
        command = "rm -rf " + intermediate_file_1 + " " + intermediate_file_2
        os.system(command)

    def create_host_guest(self):
        """
        Saves separate receptor and ligand PDB files.
        """
        with open(self.cleaned_pdb) as f1, open(self.host_pdb, "w") as f2:
            for line in f1:
                if not self.guest_resname in line and not "CRYST1" in line:
                    f2.write(line)
        with open(self.cleaned_pdb) as f1, open(
            self.guest_init_pdb, "w"
        ) as f2:
            for line in f1:
                if self.guest_resname in line or "END" in line:
                    f2.write(line)

    def realign_guest(self):
        """
        Saves a ligand PDB file with atom numbers beginning from 1.
        """
        ppdb = PandasPdb()
        ppdb.read_pdb(self.guest_init_pdb)
        to_subtract = min(ppdb.df["ATOM"]["atom_number"]) - 1
        ppdb.df["ATOM"]["atom_number"] = (
            ppdb.df["ATOM"]["atom_number"] - to_subtract
        )
        intermediate_file_1 = self.guest_pdb[:-4] + "_intermediate_1.pdb"
        intermediate_file_2 = self.guest_pdb[:-4] + "_intermediate_2.pdb"
        ppdb.to_pdb(path=intermediate_file_1)
        command = (
            "pdb4amber -i "
            + intermediate_file_1
            + " -o "
            + intermediate_file_2
        )
        os.system(command)
        to_delete = (
            intermediate_file_2[:-4] + "_nonprot.pdb",
            intermediate_file_2[:-4] + "_renum.txt",
            intermediate_file_2[:-4] + "_sslink",
        )
        os.system("rm -rf " + " ".join(to_delete))
        with open(intermediate_file_2, "r") as f1:
            filedata = f1.read()
        filedata = filedata.replace("HETATM", "ATOM  ")
        with open(self.guest_pdb, "w") as f2:
            f2.write(filedata)
        command = "rm -rf " + intermediate_file_1 + " " + intermediate_file_2
        os.system(command)

    def get_guest_coord(self):
        """
        Saves a text file of the XYZ coordinates of the ligand.
        """
        ppdb = PandasPdb()
        ppdb.read_pdb(self.guest_pdb)
        xyz = ppdb.df["ATOM"][["x_coord", "y_coord", "z_coord"]]
        xyz_to_list = xyz.values.tolist()
        np.savetxt(self.guest_xyz, xyz_to_list)

    def get_qm_resids(self):
        """
        Saves a text file of the residue numbers of the receptor within the
        proximity (as defined by the distance) from the ligand.
        """
        guest_coord_list = np.loadtxt(self.guest_xyz)
        host_atom_list = []
        for i in range(len(guest_coord_list)):
            reference_point = guest_coord_list[i]
            # TODO: move reads outside of loop
            ppdb = PandasPdb()
            ppdb.read_pdb(self.host_pdb) 
            distances = ppdb.distance(xyz=reference_point, records=("ATOM"))
            all_within_distance = ppdb.df["ATOM"][
                distances < float(self.distance)
            ]
            host_df = all_within_distance["atom_number"]
            host_list = host_df.values.tolist()
            host_atom_list.append(host_list)
        host_atom_list = list(itertools.chain(*host_atom_list))
        host_atom_list = set(host_atom_list)
        host_atom_list = list(host_atom_list)
        host_atom_list.sort()
        ppdb = PandasPdb()
        ppdb.read_pdb(self.host_pdb)
        df = ppdb.df["ATOM"][["atom_number", "residue_number", "residue_name"]]
        index_list = []
        for i in host_atom_list:
            indices = np.where(df["atom_number"] == i)
            indices = list(indices)[0]
            indices = list(indices)
            index_list.append(indices)
        index_list = list(itertools.chain.from_iterable(index_list))
        df1 = df.iloc[
            index_list,
        ]
        # TODO: make it write list of integers
        resid_num = list(df1.residue_number.unique())
        np.savetxt(self.residue_list, resid_num, fmt="%i")

    def get_host_qm_mm_atoms(self):
        """
        Saves a text file of the atom numbers of the receptors in the QM
        region and MM region separately.
        """
        resid_num = np.loadtxt(self.residue_list)
        # approximated_res_list = [int(i) for i in resid_num]
        approximated_res_list = []
        # TODO: what is this doing?
        for i in range(
            int(statistics.median(resid_num))
            - int(int(self.num_residues) / 2),
            int(statistics.median(resid_num))
            + int(int(self.num_residues) / 2),
        ):
            approximated_res_list.append(i)
        ppdb = PandasPdb()
        ppdb.read_pdb(self.host_pdb)
        df = ppdb.df["ATOM"][["atom_number", "residue_number", "residue_name"]]
        host_index_nested_list = []
        for i in approximated_res_list:
            indices = np.where(df["residue_number"] == i)
            #TODO: the program seems to error when this line is removed, which
            # makes no sense.
            indices = list(indices)[0]
            indices = list(indices)
            host_index_nested_list.append(indices)
        host_index_list = list(
            itertools.chain.from_iterable(host_index_nested_list)
        )
        df_atom = df.iloc[host_index_list]
        df_atom_number = df_atom["atom_number"]
        host_atom_list = df_atom_number.values.tolist()
        selected_atoms = []
        selected_atoms.extend(host_atom_list)
        ppdb = PandasPdb()
        ppdb.read_pdb(self.host_pdb)
        len_atoms = []
        for i in range(len(ppdb.df["ATOM"])):
            len_atoms.append(i + 1)
        non_selected_atoms = list(set(len_atoms).difference(selected_atoms))
        assert len(non_selected_atoms) + len(selected_atoms) == len(len_atoms),\
            "Sum of the atoms in the selected and non-selected region "\
            "does not equal the length of list of total atoms."
        np.savetxt(self.host_qm_atoms, selected_atoms, fmt="%i")
        np.savetxt(self.host_mm_atoms, non_selected_atoms, fmt="%i")

    def save_host_pdbs(self):
        """
        Saves a PDB file for the receptor's QM region and MM
        region separately.
        """
        selected_atoms = np.loadtxt(self.host_qm_atoms)
        # TODO: not necessary if savetxt writes in integers
        selected_atoms = [int(i) for i in selected_atoms]
        ppdb = PandasPdb()
        ppdb.read_pdb(self.host_pdb)
        for i in selected_atoms:
            ppdb.df["ATOM"] = ppdb.df["ATOM"][
                ppdb.df["ATOM"]["atom_number"] != i
            ]
        ppdb.to_pdb(
            path=self.host_mm_pdb, records=None, gz=False, append_newline=True,
        )
        non_selected_atoms = np.loadtxt(self.host_mm_atoms)
        non_selected_atoms = [int(i) for i in non_selected_atoms]
        ppdb = PandasPdb()
        ppdb.read_pdb(self.host_pdb)
        for i in non_selected_atoms:
            ppdb.df["ATOM"] = ppdb.df["ATOM"][
                ppdb.df["ATOM"]["atom_number"] != i
            ]
        ppdb.to_pdb(
            path=self.host_qm_pdb, records=None, gz=False, append_newline=True,
        )

    def get_host_mm_region_atoms(self):
        """
        Saves a text file for the atoms of the receptor's MM region
        preceding the QM region and saves another text file for the
        atoms of the receptor's MM region folllowing the QM region.
        """
        resid_num = np.loadtxt(self.residue_list)
        approximated_res_list = []
        for i in range(
            int(statistics.median(resid_num))
            - int(int(self.num_residues) / 2),
            int(statistics.median(resid_num))
            + int(int(self.num_residues) / 2),
        ):
            approximated_res_list.append(i)
        # print(approximated_res_list)
        ppdb = PandasPdb()
        ppdb.read_pdb(self.host_pdb)
        df = ppdb.df["ATOM"][["residue_number"]]
        res_list = list(set(df["residue_number"].to_list()))
        res_mm_list = list(set(res_list).difference(approximated_res_list))
        # print(res_mm_list)
        res_mm_region_I_list = []
        # TODO: This can probably be made into a single loop by comparing i
        # to the maximum value within approximated_res_list
        for i in res_mm_list:
            for j in approximated_res_list:
                if i < j:
                    res_mm_region_I_list.append(i)
        res_mm_region_I_list = list(set(res_mm_region_I_list))
        res_mm_region_II_list = list(
            set(res_mm_list).difference(res_mm_region_I_list)
        )
        # print(res_mm_region_II_list)
        ppdb.read_pdb(self.host_mm_pdb)
        df = ppdb.df["ATOM"][["atom_number", "residue_number", "residue_name"]]
        mm_region_I_index_nested_list = []
        for i in res_mm_region_I_list:
            indices = np.where(df["residue_number"] == i)
            # TODO: again, this is strange code
            indices = list(indices)[0]
            indices = list(indices)
            mm_region_I_index_nested_list.append(indices)
        mm_region_I_index_list = list(
            itertools.chain.from_iterable(mm_region_I_index_nested_list)
        )
        df_atom = df.iloc[mm_region_I_index_list]
        df_atom_number = df_atom["atom_number"]
        mm_region_I_atom_list = df_atom_number.values.tolist()
        mm_region_I_atoms = []
        mm_region_I_atoms.extend(mm_region_I_atom_list)
        mm_region_II_index_nested_list = []
        for i in res_mm_region_II_list:
            indices = np.where(df["residue_number"] == i)
            # TODO: again, this is strange code
            indices = list(indices)[0]
            indices = list(indices)
            mm_region_II_index_nested_list.append(indices)
        mm_region_II_index_list = list(
            itertools.chain.from_iterable(mm_region_II_index_nested_list)
        )
        df_atom = df.iloc[mm_region_II_index_list]
        df_atom_number = df_atom["atom_number"]
        mm_region_II_atom_list = df_atom_number.values.tolist()
        mm_region_II_atoms = []
        mm_region_II_atoms.extend(mm_region_II_atom_list)
        ppdb.read_pdb(self.host_mm_pdb)
        len_atoms = []
        for i in range(len(ppdb.df["ATOM"])):
            len_atoms.append(i + 1)
        assert len(mm_region_I_atoms) + len(mm_region_II_atoms) == len(len_atoms),\
            "Sum of the atoms in the selected and non-selected region "\
            "does not equal the length of list of total atoms."
        np.savetxt(self.host_mm_region_I_atoms, mm_region_I_atoms, fmt="%i")
        np.savetxt(self.host_mm_region_II_atoms, mm_region_II_atoms, fmt="%i")

    def save_host_mm_regions_pdbs(self):
        """
        Saves a PDB file for the receptor's MM region preceding
        the QM region and saves another PDB file for the receptor's
        MM region folllowing the QM region.
        """
        mm_region_I_atoms = np.loadtxt(self.host_mm_region_I_atoms)
        mm_region_I_atoms = [int(i) for i in mm_region_I_atoms]
        mm_region_II_atoms = np.loadtxt(self.host_mm_region_II_atoms)
        mm_region_II_atoms = [int(i) for i in mm_region_II_atoms]
        
        # NOTE: this is a slightly confusing way to define the atoms to 
        # write to a PDB - the members that are *not* in a section, rather
        # than the members that are.
        ppdb = PandasPdb()
        ppdb.read_pdb(self.host_mm_pdb)
        for i in mm_region_II_atoms:
            ppdb.df["ATOM"] = ppdb.df["ATOM"][
                ppdb.df["ATOM"]["atom_number"] != i
            ]
        ppdb.to_pdb(
            path=self.host_mm_region_I_pdb,
            records=None,
            gz=False,
            append_newline=True,
        )
        ppdb = PandasPdb()
        ppdb.read_pdb(self.host_mm_pdb)
        for i in mm_region_I_atoms:
            ppdb.df["ATOM"] = ppdb.df["ATOM"][
                ppdb.df["ATOM"]["atom_number"] != i
            ]
        ppdb.to_pdb(
            path=self.host_mm_region_II_pdb,
            records=None,
            gz=False,
            append_newline=True,
        )

    def get_qm_mm_regions(self):
        """
        Saves separate PDB files for the QM and MM regions.
        QM regions comprise the QM region of the receptor
        and the entire ligand where the MM region comprise
        the non-selected QM regions of the receptor.
        """
        with open(self.host_qm_pdb) as f1, open(self.qm_pdb, "w") as f2:
            for line in f1:
                if "ATOM" in line:
                    f2.write(line)
        with open(self.guest_pdb) as f1, open(self.qm_pdb, "a") as f2:
            for line in f1:
                if "ATOM" in line:
                    f2.write(line)
            f2.write("END")
        with open(self.host_mm_pdb) as f1, open(self.mm_pdb, "w") as f2:
            for line in f1:
                if "ATOM" in line:
                    f2.write(line)
            f2.write("END")



