# QMMMReBind - Quantum Mechanics – Molecular Mechanics ( *QMMM* ) forcefield *Re*paramaterisation of the *Bind*ing site for receptor-ligand complexes

[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/anandojha/qmmmrebind/workflows/CI/badge.svg)](https://github.com/anandojha/qmmmrebind/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/anandojha/QMMMReBind/branch/master/graph/badge.svg)](https://codecov.io/gh/anandojha/QMMMReBind/branch/master)

# <img src="https://github.com/anandojha/qmmmrebind/blob/main/images/qmmmrebind_logo.jpg" width="400">


## Software Requirements :
* Gaussian16
* Amber
* TorsionDrive
* Psi4

## Installation and Setup Instructions :
* Make sure [anaconda3](https://www.anaconda.com/) is installed on the local machine. Go to this [link](https://www.anaconda.com/products/individual) and download the latest edition of anaconda3 and install. 
* Create a new conda environment with python version 3.8 :
```bash
conda create -n qmmmrebind python=3.8
```
* Activate the conda environment:
```bash
conda activate qmmmrebind
```
* Install git, openforcefield and openbabel with the following commands :
```bash
conda install openforcefield
```
```bash
conda install -c conda-forge openbabel 
```
```bash
conda install git
```
* Clone the QMMMReBind repository :
```bash
git clone https://github.com/anandojha/qmmmrebind.git
```
* Perform the following steps to get this package installed quickly on a local linux machine. Installation in the home directory is recommended. 
```bash
cd qmmmrebind
```
```bash
python setup.py install
```
* Optionally run tests to check for proper installation. 
```bash
python setup.py test 
```

Documentation can be found [here](https://github.com/anandojha/qmmmrebind/blob/main/qmmmrebind/html/parameterize.html) .

### Copyright
Copyright (c) 2021, Anupam Anand Ojha
#### Acknowledgements
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.5.
