"""
Define various constants that are used in QMMMRebind
"""

import numpy as np

BOHRS_PER_ANGSTROM = 0.529
HARTREE_PER_KCAL_MOL = 627.509391
#kcal/mol * A^2 to kJ/mol * nm^2
KCAL_MOL_PER_KJ_MOL = 4.184
ANGSTROMS_PER_NM = 10.0
RADIANS_PER_DEGREE = np.pi / 180.0

method_basis_scale_dict = {
    "HF STO-3G": 0.817,
    "HF 3-21G": 0.906,
    "HF 3-21G*": 0.903,
    "HF 6-31G": 0.903,
    "HF 6-31G*": 0.899,
    "HF 6-31G**": 0.903,
    "HF 6-31+G**": 0.904,
    "HF 6-311G*": 0.904,
    "HF 6-311G**": 0.909,
    "HF TZVP": 0.909,
    "HF cc-pVDZ": 0.908,
    "HF cc-pVTZ": 0.91,
    "HF cc-pVQZ": 0.908,
    "HF aug-cc-pVDZ": 0.911,
    "HF aug-cc-pVTZ": 0.91,
    "HF aug-cc-pVQZ": 0.909,
    "HF daug-cc-pVDZ": 0.912,
    "HF daug-cc-pVTZ": 0.905,
    "ROHF 3-21G": 0.907,
    "ROHF 3-21G*": 0.909,
    "ROHF 6-31G": 0.895,
    "ROHF 6-31G*": 0.89,
    "ROHF 6-31G**": 0.855,
    "ROHF 6-31+G**": 0.856,
    "ROHF 6-311G*": 0.856,
    "ROHF 6-311G**": 0.913,
    "ROHF cc-pVDZ": 0.861,
    "ROHF cc-pVTZ": 0.901,
    "LSDA STO-3G": 0.896,
    "LSDA 3-21G": 0.984,
    "LSDA 3-21G*": 0.982,
    "LSDA 6-31G": 0.98,
    "LSDA 6-31G*": 0.981,
    "LSDA 6-31G**": 0.981,
    "LSDA 6-31+G**": 0.985,
    "LSDA 6-311G*": 0.984,
    "LSDA 6-311G**": 0.988,
    "LSDA TZVP": 0.988,
    "LSDA cc-pVDZ": 0.989,
    "LSDA cc-pVTZ": 0.989,
    "LSDA aug-cc-pVDZ": 0.989,
    "LSDA aug-cc-pVTZ": 0.991,
    "BLYP STO-3G": 0.925,
    "BLYP 3-21G": 0.995,
    "BLYP 3-21G*": 0.994,
    "BLYP 6-31G": 0.992,
    "BLYP 6-31G*": 0.992,
    "BLYP 6-31G**": 0.992,
    "BLYP 6-31+G**": 0.995,
    "BLYP 6-311G*": 0.998,
    "BLYP 6-311G**": 0.996,
    "BLYP TZVP": 0.998,
    "BLYP cc-pVDZ": 1.002,
    "BLYP cc-pVTZ": 0.997,
    "BLYP aug-cc-pVDZ": 0.998,
    "BLYP aug-cc-pVTZ": 0.997,
    "B1B95 STO-3G": 0.883,
    "B1B95 3-21G": 0.957,
    "B1B95 3-21G*": 0.955,
    "B1B95 6-31G": 0.954,
    "B1B95 6-31G*": 0.949,
    "B1B95 6-31G**": 0.955,
    "B1B95 6-31+G**": 0.957,
    "B1B95 6-311G*": 0.959,
    "B1B95 6-311G**": 0.96,
    "B1B95 TZVP": 0.957,
    "B1B95 cc-pVDZ": 0.961,
    "B1B95 cc-pVTZ": 0.957,
    "B1B95 aug-cc-pVDZ": 0.958,
    "B1B95 aug-cc-pVTZ": 0.959,
    "B3LYP STO-3G": 0.892,
    "B3LYP 3-21G": 0.965,
    "B3LYP 3-21G*": 0.962,
    "B3LYP 6-31G": 0.962,
    "B3LYP 6-31G*": 0.96,
    "B3LYP 6-31G**": 0.961,
    "B3LYP 6-31+G**": 0.964,
    "B3LYP 6-311G*": 0.966,
    "B3LYP 6-311G**": 0.967,
    "B3LYP TZVP": 0.965,
    "B3LYP cc-pVDZ": 0.97,
    "B3LYP cc-pVTZ": 0.967,
    "B3LYP cc-pVQZ": 0.969,
    "B3LYP aug-cc-pVDZ": 0.97,
    "B3LYP aug-cc-pVTZ": 0.968,
    "B3LYP aug-cc-pVQZ": 0.969,
    "B3PW91 STO-3G": 0.885,
    "B3PW91 3-21G": 0.961,
    "B3PW91 3-21G*": 0.959,
    "B3PW91 6-31G": 0.958,
    "B3PW91 6-31G*": 0.957,
    "B3PW91 6-31G**": 0.958,
    "B3PW91 6-31+G**": 0.96,
    "B3PW91 6-311G*": 0.963,
    "B3PW91 6-311G**": 0.963,
    "B3PW91 TZVP": 0.964,
    "B3PW91 cc-pVDZ": 0.965,
    "B3PW91 cc-pVTZ": 0.962,
    "B3PW91 aug-cc-pVDZ": 0.965,
    "B3PW91 aug-cc-pVTZ": 0.965,
    "mPW1PW91 STO-3G": 0.879,
    "mPW1PW91 3-21G": 0.955,
    "mPW1PW91 3-21G*": 0.95,
    "mPW1PW91 6-31G": 0.947,
    "mPW1PW91 6-31G*": 0.948,
    "mPW1PW91 6-31G**": 0.952,
    "mPW1PW91 6-31+G**": 0.952,
    "mPW1PW91 6-311G*": 0.954,
    "mPW1PW91 6-311G**": 0.957,
    "mPW1PW91 TZVP": 0.954,
    "mPW1PW91 cc-pVDZ": 0.958,
    "mPW1PW91 cc-pVTZ": 0.959,
    "mPW1PW91 aug-cc-pVDZ": 0.958,
    "mPW1PW91 aug-cc-pVTZ": 0.958,
    "PBEPBE STO-3G": 0.914,
    "PBEPBE 3-21G": 0.991,
    "PBEPBE 3-21G*": 0.954,
    "PBEPBE 6-31G": 0.986,
    "PBEPBE 6-31G*": 0.986,
    "PBEPBE 6-31G**": 0.986,
    "PBEPBE 6-31+G**": 0.989,
    "PBEPBE 6-311G*": 0.99,
    "PBEPBE 6-311G**": 0.991,
    "PBEPBE TZVP": 0.989,
    "PBEPBE cc-pVDZ": 0.994,
    "PBEPBE cc-pVTZ": 0.993,
    "PBEPBE aug-cc-pVDZ": 0.994,
    "PBEPBE aug-cc-pVTZ": 0.994,
    "PBE1PBE STO-3G": 0.882,
    "PBE1PBE 3-21G": 0.96,
    "PBE1PBE 3-21G*": 0.96,
    "PBE1PBE 6-31G": 0.956,
    "PBE1PBE 6-31G*": 0.95,
    "PBE1PBE 6-31G**": 0.953,
    "PBE1PBE 6-31+G**": 0.955,
    "PBE1PBE 6-311G*": 0.959,
    "PBE1PBE 6-311G**": 0.959,
    "PBE1PBE TZVP": 0.96,
    "PBE1PBE cc-pVDZ": 0.962,
    "PBE1PBE cc-pVTZ": 0.961,
    "PBE1PBE aug-cc-pVDZ": 0.962,
    "PBE1PBE aug-cc-pVTZ": 0.962,
    "HSEh1PBE STO-3G": 0.883,
    "HSEh1PBE 3-21G": 0.963,
    "HSEh1PBE 3-21G*": 0.96,
    "HSEh1PBE 6-31G": 0.957,
    "HSEh1PBE 6-31G*": 0.951,
    "HSEh1PBE 6-31G**": 0.954,
    "HSEh1PBE 6-31+G**": 0.955,
    "HSEh1PBE 6-311G*": 0.96,
    "HSEh1PBE 6-311G**": 0.96,
    "HSEh1PBE TZVP": 0.96,
    "HSEh1PBE cc-pVDZ": 0.962,
    "HSEh1PBE cc-pVTZ": 0.961,
    "HSEh1PBE aug-cc-pVDZ": 0.962,
    "HSEh1PBE aug-cc-pVTZ": 0.962,
    "TPSSh 3-21G": 0.969,
    "TPSSh 3-21G*": 0.966,
    "TPSSh 6-31G": 0.962,
    "TPSSh 6-31G*": 0.959,
    "TPSSh 6-31G**": 0.959,
    "TPSSh 6-31+G**": 0.963,
    "TPSSh 6-311G*": 0.963,
    "TPSSh TZVP": 0.964,
    "TPSSh cc-pVDZ": 0.972,
    "TPSSh cc-pVTZ": 0.968,
    "TPSSh aug-cc-pVDZ": 0.967,
    "TPSSh aug-cc-pVTZ": 0.965,
    "B97D3 3-21G": 0.983,
    "B97D3 6-31G*": 0.98,
    "B97D3 6-31+G**": 0.983,
    "B97D3 6-311G**": 0.986,
    "B97D3 TZVP": 0.986,
    "B97D3 cc-pVDZ": 0.992,
    "B97D3 cc-pVTZ": 0.986,
    "B97D3 aug-cc-pVTZ": 0.985,
    "MP2 STO-3G": 0.872,
    "MP2 3-21G": 0.955,
    "MP2 3-21G*": 0.951,
    "MP2 6-31G": 0.957,
    "MP2 6-31G*": 0.943,
    "MP2 6-31G**": 0.937,
    "MP2 6-31+G**": 0.941,
    "MP2 6-311G*": 0.95,
    "MP2 6-311G**": 0.95,
    "MP2 TZVP": 0.948,
    "MP2 cc-pVDZ": 0.953,
    "MP2 cc-pVTZ": 0.95,
    "MP2 cc-pVQZ": 0.948,
    "MP2 aug-cc-pVDZ": 0.959,
    "MP2 aug-cc-pVTZ": 0.953,
    "MP2 aug-cc-pVQZ": 0.95,
    "MP2=FULL STO-3G": 0.889,
    "MP2=FULL 3-21G": 0.955,
    "MP2=FULL 3-21G*": 0.948,
    "MP2=FULL 6-31G": 0.95,
    "MP2=FULL 6-31G*": 0.942,
    "MP2=FULL 6-31G**": 0.934,
    "MP2=FULL 6-31+G**": 0.939,
    "MP2=FULL 6-311G*": 0.947,
    "MP2=FULL 6-311G**": 0.949,
    "MP2=FULL TZVP": 0.953,
    "MP2=FULL cc-pVDZ": 0.95,
    "MP2=FULL cc-pVTZ": 0.949,
    "MP2=FULL cc-pVQZ": 0.957,
    "MP2=FULL aug-cc-pVDZ": 0.969,
    "MP2=FULL aug-cc-pVTZ": 0.951,
    "MP2=FULL aug-cc-pVQZ": 0.956,
    "MP3 STO-3G": 0.894,
    "MP3 3-21G": 0.968,
    "MP3 3-21G*": 0.965,
    "MP3 6-31G": 0.966,
    "MP3 6-31G*": 0.939,
    "MP3 6-31G**": 0.935,
    "MP3 6-31+G**": 0.931,
    "MP3 TZVP": 0.935,
    "MP3 cc-pVDZ": 0.948,
    "MP3 cc-pVTZ": 0.945,
    "MP3=FULL 6-31G*": 0.938,
    "MP3=FULL 6-31+G**": 0.932,
    "MP3=FULL TZVP": 0.934,
    "MP3=FULL cc-pVDZ": 0.94,
    "MP3=FULL cc-pVTZ": 0.933,
    "B2PLYP 6-31G*": 0.949,
    "B2PLYP 6-31+G**": 0.952,
    "B2PLYP TZVP": 0.954,
    "B2PLYP cc-pVDZ": 0.958,
    "B2PLYP cc-pVTZ": 0.959,
    "B2PLYP cc-pVQZ": 0.957,
    "B2PLYP aug-cc-pVTZ": 0.961,
    "B2PLYP=FULL 3-21G": 0.952,
    "B2PLYP=FULL 6-31G*": 0.948,
    "B2PLYP=FULL 6-31+G**": 0.951,
    "B2PLYP=FULL TZVP": 0.954,
    "B2PLYP=FULL cc-pVDZ": 0.959,
    "B2PLYP=FULL cc-pVTZ": 0.956,
    "B2PLYP=FULL aug-cc-pVDZ": 0.962,
    "B2PLYP=FULL aug-cc-pVTZ": 0.959,
    "CID 3-21G": 0.932,
    "CID 3-21G*": 0.931,
    "CID 6-31G": 0.935,
    "CID 6-31G*": 0.924,
    "CID 6-31G**": 0.924,
    "CID 6-31+G**": 0.924,
    "CID 6-311G*": 0.929,
    "CID cc-pVDZ": 0.924,
    "CID cc-pVTZ": 0.927,
    "CISD 3-21G": 0.941,
    "CISD 3-21G*": 0.934,
    "CISD 6-31G": 0.938,
    "CISD 6-31G*": 0.926,
    "CISD 6-31G**": 0.918,
    "CISD 6-31+G**": 0.922,
    "CISD 6-311G*": 0.925,
    "CISD cc-pVDZ": 0.922,
    "CISD cc-pVTZ": 0.93,
    "QCISD 3-21G": 0.969,
    "QCISD 3-21G*": 0.961,
    "QCISD 6-31G": 0.964,
    "QCISD 6-31G*": 0.952,
    "QCISD 6-31G**": 0.941,
    "QCISD 6-31+G**": 0.945,
    "QCISD 6-311G*": 0.957,
    "QCISD 6-311G**": 0.954,
    "QCISD TZVP": 0.955,
    "QCISD cc-pVDZ": 0.959,
    "QCISD cc-pVTZ": 0.956,
    "QCISD aug-cc-pVDZ": 0.969,
    "QCISD aug-cc-pVTZ": 0.962,
    "CCD 3-21G": 0.972,
    "CCD 3-21G*": 0.957,
    "CCD 6-31G": 0.96,
    "CCD 6-31G*": 0.947,
    "CCD 6-31G**": 0.938,
    "CCD 6-31+G**": 0.942,
    "CCD 6-311G*": 0.955,
    "CCD 6-311G**": 0.955,
    "CCD TZVP": 0.948,
    "CCD cc-pVDZ": 0.957,
    "CCD cc-pVTZ": 0.934,
    "CCD aug-cc-pVDZ": 0.965,
    "CCD aug-cc-pVTZ": 0.957,
    "CCSD 3-21G": 0.943,
    "CCSD 3-21G*": 0.943,
    "CCSD 6-31G": 0.943,
    "CCSD 6-31G*": 0.944,
    "CCSD 6-31G**": 0.933,
    "CCSD 6-31+G**": 0.934,
    "CCSD 6-311G*": 0.954,
    "CCSD TZVP": 0.954,
    "CCSD cc-pVDZ": 0.947,
    "CCSD cc-pVTZ": 0.941,
    "CCSD cc-pVQZ": 0.951,
    "CCSD aug-cc-pVDZ": 0.963,
    "CCSD aug-cc-pVTZ": 0.956,
    "CCSD aug-cc-pVQZ": 0.953,
    "CCSD=FULL 6-31G*": 0.95,
    "CCSD=FULL TZVP": 0.948,
    "CCSD=FULL cc-pVTZ": 0.948,
    "CCSD=FULL aug-cc-pVTZ": 0.951,
}

element_list = [
    ["1 ", "H ", "Hydrogen"],
    ["2 ", "He", "Helium"],
    ["3 ", "Li", "Lithium"],
    ["4 ", "Be", "Beryllium"],
    ["5 ", "B ", "Boron"],
    ["6 ", "C ", "Carbon"],
    ["7 ", "N ", "Nitrogen"],
    ["8 ", "O ", "Oxygen"],
    ["9 ", "F ", "Fluorine"],
    ["10", "Ne", "Neon"],
    ["11", "Na", "Sodium"],
    ["12", "Mg", "Magnesium"],
    ["13", "Al", "Aluminum"],
    ["14", "Si", "Silicon"],
    ["15", "P ", "Phosphorus"],
    ["16", "S ", "Sulfur"],
    ["17", "Cl", "Chlorine"],
    ["18", "Ar", "Argon"],
    ["19", "K ", "Potassium"],
    ["20", "Ca", "Calcium"],
    ["21", "Sc", "Scandium"],
    ["22", "Ti", "Titanium"],
    ["23", "V ", "Vanadium"],
    ["24", "Cr", "Chromium"],
    ["25", "Mn", "Manganese"],
    ["26", "Fe", "Iron"],
    ["27", "Co", "Cobalt"],
    ["28", "Ni", "Nickel"],
    ["29", "Cu", "Copper"],
    ["30", "Zn", "Zinc"],
    ["31", "Ga", "Gallium"],
    ["32", "Ge", "Germanium"],
    ["33", "As", "Arsenic"],
    ["34", "Se", "Selenium"],
    ["35", "Br", "Bromine"],
    ["36", "Kr", "Krypton"],
    ["37", "Rb", "Rubidium"],
    ["38", "Sr", "Strontium"],
    ["39", "Y ", "Yttrium"],
    ["40", "Zr", "Zirconium"],
    ["41", "Nb", "Niobium"],
    ["42", "Mo", "Molybdenum"],
    ["43", "Tc", "Technetium"],
    ["44", "Ru", "Ruthenium"],
    ["45", "Rh", "Rhodium"],
    ["46", "Pd", "Palladium"],
    ["47", "Ag", "Silver"],
    ["48", "Cd", "Cadmium"],
    ["49", "In", "Indium"],
    ["50", "Sn", "Tin"],
    ["51", "Sb", "Antimony"],
    ["52", "Te", "Tellurium"],
    ["53", "I ", "Iodine"],
    ["54", "Xe", "Xenon"],
    ["55", "Cs", "Cesium"],
    ["56", "Ba", "Barium"],
    ["57", "La", "Lanthanum"],
    ["58", "Ce", "Cerium"],
    ["59", "Pr", "Praseodymium"],
    ["60", "Nd", "Neodymium"],
    ["61", "Pm", "Promethium"],
    ["62", "Sm", "Samarium"],
    ["63", "Eu", "Europium"],
    ["64", "Gd", "Gadolinium"],
    ["65", "Tb", "Terbium"],
    ["66", "Dy", "Dysprosium"],
    ["67", "Ho", "Holmium"],
    ["68", "Er", "Erbium"],
    ["69", "Tm", "Thulium"],
    ["70", "Yb", "Ytterbium"],
    ["71", "Lu", "Lutetium"],
    ["72", "Hf", "Hafnium"],
    ["73", "Ta", "Tantalum"],
    ["74", "W ", "Tungsten"],
    ["75", "Re", "Rhenium"],
    ["76", "Os", "Osmium"],
    ["77", "Ir", "Iridium"],
    ["78", "Pt", "Platinum"],
    ["79", "Au", "Gold"],
    ["80", "Hg", "Mercury"],
    ["81", "Tl", "Thallium"],
    ["82", "Pb", "Lead"],
    ["83", "Bi", "Bismuth"],
    ["84", "Po", "Polonium"],
    ["85", "At", "Astatine"],
    ["86", "Rn", "Radon"],
    ["87", "Fr", "Francium"],
    ["88", "Ra", "Radium"],
    ["89", "Ac", "Actinium"],
    ["90", "Th", "Thorium"],
    ["91", "Pa", "Protactinium"],
    ["92", "U ", "Uranium"],
    ["93", "Np", "Neptunium"],
    ["94", "Pu", "Plutonium"],
    ["95", "Am", "Americium"],
    ["96", "Cm", "Curium"],
    ["97", "Bk", "Berkelium"],
    ["98", "Cf", "Californium"],
    ["99", "Es", "Einsteinium"],
]

def get_vibrational_scaling(functional, basis_set):

    """
    Returns vibrational scaling factor given the functional
    and the basis set for the QM engine.

    Parameters
    ----------
    functional: str
        Functional

    basis_set: str
        Basis set

    Returns
    -------
    vib_scale: float
        Vibrational scaling factor corresponding to the given
        the basis_set and the functional.

    Examples
    --------
    >>> get_vibrational_scaling("QCISD", "6-311G*")
    0.957

    """
    vib_scale = method_basis_scale_dict.get(functional + " " + basis_set)
    return vib_scale

# TODO: replace function with list comprehension using correct constant
def list_kJ_kcal(list_):

    """
    Convert the elements in the list from
    kiloJoules units to kiloCalories units.

    Parameters
    ----------
    list_ : list
        List with elements in units of kJ.

    Returns
    -------
    converted_list : list
        List with elements in units of kcal.

    Examples
    --------
    >>> list_ = [6, 3, 5]
    >>> list_kJ_kcal(list_)
    [1.4340344168260037, 0.7170172084130019, 1.1950286806883366]

    """
    converted_list = [i / 4.184 for i in list_]
    return converted_list

# TODO: replace with list comprehension using constant
def list_hartree_kcal(list_):
    """
    Convert the elements in the list from
    hartree units to kiloCalories units.

    Parameters
    ----------
    list_ : list
        List with elements in units of hartree.

    Returns
    -------
    converted_list : list
        List with elements in units of kcal.

    Examples
    --------
    >>> list_ = [6, 3, 5]
    >>> list_hartree_kcal(list_)
    [3765.0564000000004, 1882.5282000000002, 3137.547]

    """
    converted_list = [i * 627.5094 for i in list_]
    return converted_list

