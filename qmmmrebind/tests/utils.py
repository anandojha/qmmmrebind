import os

def get_data_filename(filename):

    filename = filename.lower() 

    test_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(test_dir, "data")

    fname = os.path.join(data_dir, filename)

    if not os.path.exists(fname):
        raise OSError("File %s not found!" % fname)

    return fname