# QMMMReBind - Quantum Mechanics – Molecular Mechanics ( *QMMM* ) forcefield *Re*paramaterisation of the *Bind*ing site for receptor-ligand complexes

[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/anandojha/qmmmrebind/workflows/CI/badge.svg)](https://github.com/anandojha/qmmmrebind/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/anandojha/QMMMReBind/branch/master/graph/badge.svg)](https://codecov.io/gh/anandojha/QMMMReBind/branch/master)

## Software Requirements :
* Gaussian16
* Amber
* TorsionDrive
* Psi4

## Installation and Set up Instructions :
* conda create -n qmmmrebind python=3.8
* conda activate qmmmrebind
* conda install git openforcefield
* conda install -c conda-forge biopandas mendeleev openmm parmed openbabel 
* git clone https://github.com/anandojha/qmmmrebind.git
* cd qmmmrebind
* pytest 

### Copyright

Copyright (c) 2021, Anupam Anand Ojha

#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.5.
