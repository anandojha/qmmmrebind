"""
Base functions and methods
"""

from biopandas.pdb import PandasPdb
import modules.constants as const

def get_num_host_atoms(host_pdb):

    """
    Reads the host PDB file and returns the
    total number of atoms.
    """

    ppdb = PandasPdb()
    ppdb.read_pdb(host_pdb)
    no_host_atoms = ppdb.df["ATOM"].shape[0]
    return no_host_atoms

# TODO: replace with list[::-1]
def reverse_list(lst):

    """
    Returns the reversed form of a given list.

    Parameters
    ----------
    lst : list
        Input list.

    Returns
    -------
    reversed_list : list
        Reversed input list.

    Examples
    --------
    >>> lst = [5, 4, 7, 2]
    >>> reverse_list(lst)
    [2, 7, 4, 5]

    """
    reversed_list = lst[::-1]
    return reversed_list

# TODO: use a set
def uniq(input_):

    """
    Returns a list with only unique elements from a list
    containing duplicate / repeating elements.

    Parameters
    ----------
    input_ : list
        Input list.

    Returns
    -------
    output : list
        List with only unique elements.

    Examples
    --------
    >>> lst = [2, 4, 2, 9, 10, 35, 10]
    >>> uniq(lst)
    [2, 4, 9, 10, 35]

    """
    output = []
    for x in input_:
        if x not in output:
            output.append(x)
    return output


# TODO: replace instances of function with the list comprehension
def list_to_dict(lst):

    """
    Converts an input list with mapped characters (every
    odd entry is the key of the dictionary and every
    even entry adjacent to the odd entry is its correponding
    value)  to a dictionary.

    Parameters
    ----------
    lst : list
        Input list.

    Returns
    -------
    res_dct : dict
        A dictionary with every element mapped with
        its successive element starting from index 0.

    Examples
    --------
    >>> lst = [5, 9, 3, 6, 2, 7]
    >>> list_to_dict(lst)
    {5: 9, 3: 6, 2: 7}

    """

    res_dct = {lst[i]: lst[i + 1] for i in range(0, len(lst), 2)}
    return res_dct

def scale_list(list_):

    """
    Returns a scaled list with the minimum value
    subtracted from each element of the corresponding list.

    Parameters
    ----------
    list_ : list
        Input list.

    Returns
    -------
    scaled_list : list
        Scaled list.

    Examples
    --------
    >>> list_ = [6, 3, 5, 11, 3, 2, 8, 6]
    >>> scale_list(list_)
    [4, 1, 3, 9, 1, 0, 6, 4]

    """
    scaled_list = [i - min(list_) for i in list_]
    return scaled_list

def list_diff(list_1, list_2):

    """
    Returns the difference between two lists as a list.

    Parameters
    ----------
    list_1 : list
        First list

    list_2 : list
        Second list.

    Returns
    -------
    diff_list : list
        List containing the diferences between the elements of
        the two lists.
    Examples
    --------
    >>> list_1 = [4, 2, 8, 3, 0, 6, 7]
    >>> list_2 = [5, 3, 1, 5, 6, 0, 4]
    >>> list_diff(list_1, list_2)
    [-1, -1, 7, -2, -6, 6, 3]

    """
    diff_list = []
    zipped_list = zip(list_1, list_2)
    for list1_i, list2_i in zipped_list:
        diff_list.append(list1_i - list2_i)
    return diff_list

# TODO: find a better way to do this
def truncate(x):

    """
    Returns a float or an integer with an exact number
    of characters.

    Parameters
    ----------
    x: str
       input value

    """
    if len(str(int(float(x)))) == 1:
        x = format(x, ".8f")
    if len(str(int(float(x)))) == 2:
        x = format(x, ".7f")
    if len(str(int(float(x)))) == 3:
        x = format(x, ".6f")
    if len(str(int(float(x)))) == 4:
        x = format(x, ".5f")
    if len(str(x)) > 10:
        x = round(x, 10)
    return x