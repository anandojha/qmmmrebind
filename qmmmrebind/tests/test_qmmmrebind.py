"""
Unit and regression test for the qmmmrebind package.
"""

# Import package, test suite, and other packages as needed
import qmmmrebind
import pytest
import sys
from .utils import get_data_filename

def test_qmmmrebind_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "qmmmrebind" in sys.modules

def test_len_list_to_dict():
    """Test the list to dict function"""
    test_list = ["a", "b", "c", "d"]
    ret = qmmmrebind.list_to_dict(test_list)
    assert len(ret) == 2

def test_key_list_to_dict():
    """Test the list to dict function"""
    test_list = ["a", "b", "c", "d"]
    ret = qmmmrebind.list_to_dict(test_list)
    assert ret["a"] == "b"

def test_find_word_in_file():
    """Test if a word is in a file"""
    filename = get_data_filename("test_input_file.dat") 
    word = "ucsd"
    ret = qmmmrebind.search_in_file(file=filename, word=word)
    assert ret[0][0] == 3
