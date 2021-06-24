# QMMMReBind - Quantum Mechanics – Molecular Mechanics ( *QMMM* ) forcefield *Re*paramaterisation of the *Bind*ing site for receptor-ligand complexes

| **Status** |[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/anandojha/qmmmrebind.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/anandojha/qmmmrebind/context:python) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black) [![CI](https://github.com/anandojha/qmmmrebind/workflows/CI/badge.svg)](https://github.com/anandojha/qmmmrebind/actions?query=workflow%3ACI)  [![codecov](https://codecov.io/gh/anandojha/QMMMReBind/branch/main/graph/badge.svg)](https://app.codecov.io/gh/anandojha/qmmmrebind)|
| :------ | :------ |

| **Foundation** | [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![python](https://img.shields.io/badge/python-3.8-blue.svg)](https://www.python.org/)|
| :------ | :------ |

# <img src="https://github.com/anandojha/qmmmrebind/blob/main/images/qmmmrebind_logo.jpg" width="400">

## Software Requirements :
Make sure to install these packages before running the QMMMReBind:

* Gaussian16
* TorsionDrive & Psi4 : A new conda environment is recommended to run torsiondrive calculations through Psi4 QM engine. Following commands in the terminal can allow us to perform these calculations: 
```bash
conda create -n torsiondrive python=3.6 # create a new conda environment 
conda activate torsiondrive # activate torsiondrive environment
conda install -c conda-forge torsiondrive # install torsiondrive
conda install -c psi4 psi4 dftd3 gcp # install dependent softwares 
```

## Installation and Setup Instructions :
* Make sure [anaconda3](https://www.anaconda.com/) is installed on the local machine. Go to the  [download](https://www.anaconda.com/products/individual) page of anaconda3 and install the latest version of anaconda3.
* Create a new conda environment with python = 3.8 and install the package with the following commands in the terminal: 
```bash
conda create -n qmmmrebind python=3.8
conda activate qmmmrebind # activate the conda environment
conda install -c conda-forge openff-toolkit # install openforcefield toolkit
conda install -c conda-forge openbabel # install openbabel
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

Documentation can be found [here](https://qmmmrebind.readthedocs.io/en/latest/index.html).

## Authors and Contributors
The following people have contributed directly to the coding and validation efforts of SEEKR2 (listed in an alphabetical order of last name). 
The author would like to thank everyone who has helped or will help improve this project by providing feedback, bug reports, or other comments.

* Rommie Amaro, UC San Diego (Principal Investigator)
* Eliseo Marin-Rimoldi, MoLSSI (Project Mentor and Collaborator)
* Anupam Anand Ojha, UC San Diego (Author and Lead Developer)
* Saumya Thakur, IIT Bombay (Documentation)
* Lane Votapka, UC San Diego (Project Mentor and Contributor)

### Copyright
Copyright (c) 2021, Anupam Anand Ojha
#### Acknowledgements
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.5.
