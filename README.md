# QMMMReBind - Quantum Mechanics – Molecular Mechanics ( *QMMM* ) forcefield *Re*paramaterisation of the *Bind*ing site for receptor-ligand complexes

| **Status** |[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/anandojha/qmmmrebind.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/anandojha/qmmmrebind/context:python) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black) [![CI](https://github.com/anandojha/qmmmrebind/workflows/CI/badge.svg)](https://github.com/anandojha/qmmmrebind/actions?query=workflow%3ACI)  [![codecov](https://codecov.io/gh/anandojha/QMMMReBind/branch/main/graph/badge.svg)](https://app.codecov.io/gh/anandojha/qmmmrebind)|
| :------ | :------ |

| **Foundation** | [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![python](https://img.shields.io/badge/python-3.8-blue.svg)](https://www.python.org/)|


# <img src="https://github.com/anandojha/qmmmrebind/blob/main/images/qmmmrebind_logo.jpg" width="400">


## Software Requirements :
* Gaussian16
* TorsionDrive
* Psi4

## Installation and Setup Instructions :
* Make sure [anaconda3](https://www.anaconda.com/) is installed on the local machine. Go to this [link](https://www.anaconda.com/products/individual) and download the latest edition of anaconda3 and install. 
* Create a new conda environment with python version 3.8 :
```bash
conda create -n qmmmrebind python=3.8
conda activate qmmmrebind # activate the conda environment
conda install openforcefield # install openforcefield
conda install openbabel -c conda-forge # install openbabel
conda install git # install git
```
* Clone the *QMMMReBind* repository :
```bash
git clone https://github.com/anandojha/qmmmrebind.git
```
* Perform the following steps to get this package installed quickly on a local linux machine (Installation in the home directory is recommended) : 
```bash
cd qmmmrebind
python setup.py install
python setup.py test  # optionally run tests to check for proper installation 
```

Documentation can be found [here](https://qmmmrebind.tiiny.site/)

### Copyright
Copyright (c) 2021, Anupam Anand Ojha
#### Acknowledgements
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.5.
