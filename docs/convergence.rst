Tips to get QMMMReBind working:
===============

########################
Gaussian 16  
########################

When selecting a large QM region for the receptor, there mey be convergence failures. To avoid convergence failures, following adjustements are recommended:

* GEOM=ALLCHECKPOINT : Reads the molecular geometry, charge, multiplicity, and title from the checkpoint file. This is often used to start a second calculation at a different level of theory.

* GUESS=READ : Reads the initial guess from the checkpoint file. If the basis set specified is different from the basis set used in the job which generated the checkpoint file, then the wave function will be projected from one basis to the other. This is an efficient way to switch from one basis to another. When a calculation is started using information from a checkpoint file, calculation results will be placed in the exact same checkpoint file, overwriting the original checkpoint file. Thus it is always a good idea to make a backup copy of the checkpoint file.

* To continue a job that has failed or interrupted, change the # line to include "OPT=RESTART"

* For sequential QM calculation, use "GEOM=CHECK GUESS=READ". For a B3LYP/6-31G geometry optimisation and frequency calculation of a large QM region, the following lines of code in the # line is recommended:

1. # PM3 OPT INTEGRAL=(GRID=ULTRAFINE) SCF=(maxcycles=4000) SYMMETRY=NONE
2. # HF STO-3G OPT GUESS=READ INTEGRAL=(GRID=ULTRAFINE) SCF=(maxcycles=4000) SYMMETRY=NONE
3. # HF 6-31G OPT FREQ GUESS=READ INTEGRAL=(GRID=ULTRAFINE) SCF=(maxcycles=4000) SYMMETRY=NONE POP(MK,READRADII) IOP(6/33=2,6/42=6)
4. # BLYP 3-21G OPT FREQ GUESS=READ INTEGRAL=(GRID=ULTRAFINE) SCF=(maxcycles=4000) SYMMETRY=NONE POP(MK,READRADII) IOP(6/33=2,6/42=6)
5. # BLYP 6-31G OPT FREQ GUESS=READ INTEGRAL=(GRID=ULTRAFINE) SCF=(maxcycles=4000) SYMMETRY=NONE POP(MK,READRADII) IOP(6/33=2,6/42=6)
6. # B3LYP 6-31G OPT FREQ GUESS=READ INTEGRAL=(GRID=ULTRAFINE) SCF=(maxcycles=4000) SYMMETRY=NONE POP(MK,READRADII) IOP(6/33=2,6/42=6)6



