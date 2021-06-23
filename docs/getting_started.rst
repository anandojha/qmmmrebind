Getting Started
===============

This page details how to get started with QMMMReBind. 

########################
Introduction 
########################

Understanding the interaction of metabolites or drugs with biomolecules can help improve our understanding of drug discovery and development. Accurate computational prediction of kinetic rates or residence times of an intermolecular encounter can help us identify lead drugs. Kinetic binding parameters are often correlated with efficacy rather than binding affinities [1,2]. The binding kinetic profile of a drug molecule is characterized by the bimolecular association rate constant (k\ :sub:`on`\)  and the dissociation rate constant (k\ :sub:`off`\).  Drug-Target residence time (1/ k\ :sub:`off`\)  has received the most attention from the drug discovery community since drugs with long-lasting target occupancy are often correlated with greater in-vivo efficacy [3]. k\ :sub:`on`\  is also useful for predicting in-vivo efficacy, particularly related to understanding drug rebinding. Both k\ :sub:`on`\  and k\ :sub:`off`\  can be used to compute the binding free energy. A complete kinetic profile (k\ :sub:`on`\  and k\ :sub:`off`\), in addition to binding free energy, is extremely desirable for the prediction of in-vivo efficacy and further optimization of lead compounds. 


Our laboratory has recently developed a multiscale milestoning simulation approach to estimate receptor-ligand binding kinetics computationally [4,5]. This tool, called “Simulation Enabled Estimation of Kinetic Rates” (SEEKR), incorporates the multiscale and parallel implementation of molecular dynamics (MD) and Brownian dynamics (BD) simulations using the milestoning approach to calculate the association and dissociation rates of receptor-ligand complexes. This approach requires orders of magnitude less simulation time than classical MD simulations and comparable or less simulation time than other enhanced sampling techniques. SEEKR has demonstrated successes for calculating receptor-ligand binding association (k\ :sub:`on`\)  and dissociation rates (k\ :sub:`off`\)  for multiple systems (such as the well-studied system of the protease, trypsin, with the noncovalent binder, benzamidine) as well as for rank-ordering a series of small molecules by dissociation rates and binding free energies [6]. SEEKR is among the few simulation approaches that can obtain an entirely computational estimate of the binding kinetic (k\ :sub:`on`\  and k\ :sub:`off`\)  and thermodynamic (G\ :sub:`bind`\) profiles and shows good agreement with experiment, often using less simulation time than other approaches and requiring no biasing or reweighting of its simulations. 


Significant challenges for pharmaceutically relevant receptor-ligand systems include the size and flexibility of the ligands, large-scale conformational rearrangements, and the need for extensive sampling associated with these events. Timescales of these rearrangements are much longer than can be sampled adequately with classical MD simulations. SEEKR calculations using the original implementation also struggle to simulate these more complex cases. As a consequence of limited MD sampling of rare transitions between these states, which are critical for describing the binding and unbinding event, long-timescale MD simulations are required to sample distributions on each milestone. In an effort to increase the efficiency and accuracy of kinetics calculations, Markovian Milestoning with Voronoi Tessellations (MMVT) has been implemented in SEEKR [7]. In this milestoning scheme, milestones can not only subdivide the distance of the ligand from the binding site, but also the other slow degrees of freedom in the system, such as ligand orientation or protein loop and hinge motions. In addition to aiding in the sampling of rare events with the placement of additional milestones, MMVT reduces the simulation time needed for SEEKR calculations as it overcomes the sampling bottleneck associated with the previous implementation, obtaining an equilibrium distribution on each milestone. Kinetics can then be obtained directly from short, parallel simulations within each Voronoi cell.


MMVT-SEEKR holds much potential for additional improvements in sampling, reduction of simulation times, and accuracy. MMVT-SEEKR currently incorporates MD simulations through NAMD and BD simulations through Browndye [8]. OpenMM is an increasingly popular and effective MD engine and is well-suited for running MD calculations using graphical processing units (GPUs) which are significantly faster than single-core CPU implementations, although serial and multithread CPU computations are also possible [9]. OpenMM offers all of the most common MD simulation capabilities, including a wide variety of integration schemes, compatibility with AMBER and CHARMM forcefields, and a high degree of customizability with forces and constraints within the simulation system. The Amaro lab is in the final stages of developing “SEEKR2,” which is a plugin for the molecular dynamics toolkit OpenMM to perform MMVT simulations. Implementing SEEKR within OpenMM not only significantly improved performance benchmarks (due to GPU speedups) but also allowed us to seamlessly connect with other tools within the MolSSI network.


Our laboratory previously demonstrated the effectiveness of SEEKR in predicting and ranking a series of seven small-molecule compounds for the model system , \beta\-cyclodextrin [6]. This ranking was based on estimating  k\ :sub:`on`\  and k\ :sub:`off`\   of seven host-guest systems. Although results were in good agreement with the previously conducted long timescale MD simulations for the same set of ligands with the same forcefield (GAFF and Q4MD), both methods failed in determining the correct orders for the k\ :sub:`on`\ 's [10]. Predicted k\ :sub:`off`\'s  also had deviations from experimental values although the rankings for the k\ :sub:`off`\'s  were accurate. 


A current limitation of SEEKR is that it relies on fixed point charge force fields, even in the bound state (where polarization may be an issue). We hypothesize that the deviations of k\ :sub:`off`\  from experimental values and incorrect prediction of k\ :sub:`on`\   for the above-described set of host-guest systems can mostly be attributed to the less accurate forcefield parameters for these systems. Highly accurate atomistic force fields are essential to achieve precise k\ :sub:`on`\  and k\ :sub:`off`\   as statistics collected within the milestones depends heavily on the forcefield parameters.


Thus we propose to further the multiscale nature of SEEKR by adding a quantum mechanically re-parameterized QM region to the inner-most milestone (bound state). This additional step will enable the development of forcefield parameters for the ligand and specific protein residues within the vicinity of the ligand through quantum mechanical (QM) calculations, thereby eliminating the limitation of polarization effects in the bound state. This is achieved by integrating QM engines such as Gaussian and psi4 for geometry optimization followed by hessian matrix calculation at the optimized geometry [11-13] and extracting the QM obtained charges. Dihedral and torsional angle potential energy terms form an essential part of the forcefields.  Although quantum mechanical optimizations are carried out relaxing the orthogonal degrees of freedom while fixing the target torsion angles on a grid of values, however for complex potential energy surfaces (PES), results may not be satisfactory in turn leading to inaccurate forcefield descriptions. To avoid such cases, we have incorporated another post-processing tool named TorsionDrive [14]. This tool generates energetically minimized structures on a grid of torsion constraints through recursive wavefront propagation algorithm, thus eliminating the drawbacks of traditional scanning approaches. We can expect highly accurate forcefield parameters through this approach.


With the already existing data from our previous estimates, it is straightforward to compare and interpret the kinetics with the newly parameterized force fields.  We expect to achieve higher accuracy in predicting k\ :sub:`on`\  and k\ :sub:`off`\  simultaneously through our revised approach. We thereby propose the idea of development and automation of quantum mechanical forcefield plugin to SEEKR2 package which would incorporate QM engines such as Gaussian as well as packages such as TorsionDrive for the calculation of torsional degrees of freedom to better estimate the k\ :sub:`on`\  and k\ :sub:`off`\  for the host-guest systems. This package, named, QMMReBind, is a standalone package and its incorporation as a plugin to SEEKR2 would add a considerable advantage in a user's flexibility to select forcefield parameters depending upon the systems of interest.


References
**********************

1. Ganotra G, Wade R. Prediction of Drug–Target Binding Kinetics by Comparative Binding Energy Analysis. ACS Medicinal Chemistry Letters. 2018;9(11):1134-1139.

2. Bernetti M, Cavalli A, Mollica L. Protein–ligand (un)binding kinetics as a new paradigm for drug discovery at the crossroad between experiments and modelling. MedChemComm. 2017;8(3):534-550.

3. Guan H, Lamb M, Peng B, Huang S, DeGrace N, Read J et al. Discovery of novel Jak2–Stat pathway inhibitors with extended residence time on target. Bioorganic & Medicinal Chemistry Letters. 2013;23(10):3105-3110.

4. Votapka L, Jagger B, Heyneman A, Amaro R. SEEKR: Simulation Enabled Estimation of Kinetic Rates, A Computational Tool to Estimate Molecular Kinetics and Its Application to Trypsin–Benzamidine Binding. The Journal of Physical Chemistry B. 2017;121(15):3597-3606.

5. Jagger, B., Votapka, L., Amaro, R. (2018). SEEKR: Simulation Enabled Estimation of Kinetic Rates, A Multiscale Approach for the Calculation of Protein-Ligand Association and Dissociation Kinetics. Biophysical Journal, 114(3), 42a. doi: 10.1016/j.bpj.2017.11.281

6. Jagger B, Lee C, Amaro R. Quantitative Ranking of Ligand Binding Kinetics with a Multiscale Milestoning Simulation Approach. The Journal of Physical Chemistry Letters. 2018;9(17):4941-4948.

7. Jagger B, Ojha A, Amaro R. Predicting Ligand Binding Kinetics Using a Markovian Milestoning with Voronoi Tessellations Multiscale Approach. Journal of Chemical Theory and Computation. 2020;16(8):5348-5357.

8. Huber G, McCammon J. Browndye: A software package for Brownian dynamics. Computer Physics Communications. 2010;181(11):1896-1905.

9. Eastman P, Swails J, Chodera J, McGibbon R, Zhao Y, Beauchamp K et al. OpenMM 7: Rapid development of high performance algorithms for molecular dynamics. PLOS Computational Biology. 2017;13(7):e1005659.

10. Tang Z, Chang C. Binding Thermodynamics and Kinetics Calculations Using Chemical Host and Guest: A Comprehensive Picture of Molecular Recognition. Journal of Chemical Theory and Computation. 2017;14(1):303-318.

11. Hagler A. Quantum Derivative Fitting and Biomolecular Force Fields: Functional Form, Coupling Terms, Charge Flux, Nonbond Anharmonicity, and Individual Dihedral Potentials. Journal of Chemical Theory and Computation. 2015;11(12):5555-5572.

12. Turney J, Simmonett A, Parrish R, Hohenstein E, Evangelista F, Fermann J et al. Psi4: an open-source ab initio electronic structure program. Wiley Interdisciplinary Reviews: Computational Molecular Science. 2011;2(4):556-565.

13. Neese, F., Wennmohs, F., Becker, U. and Riplinger, C., 2020. The ORCA quantum chemistry program package. The Journal of Chemical Physics, 152(22), p.224108.

14. Qiu Y, Smith D, Stern C, Feng M, Jang H, Wang L. Driving torsion scans with wavefront propagation. The Journal of Chemical Physics. 2020;152(24):244116.


########################
Software Requirements
########################

Make sure to install these packages before running the QMMMReBind:

* Gaussian16
* TorsionDrive
* Psi4


########################
Installation and Setup Instructions
########################

* Make sure `anaconda3 <https://www.anaconda.com/>`_ is installed on the local machine. 
* Go to the `download <https://www.anaconda.com/products/individual>`_  page of anaconda3 and install the latest version of anaconda3. 
* Create a new conda environment with python = 3.8 and install the package with the following commands in the terminal: 

.. code-block:: python

    conda create -n qmmmrebind python=3.8 # Create a new conda environment

.. code-block:: python

    conda activate qmmmrebind # Activate the conda environment

.. code-block:: python

    conda install openforcefield # Install openforcefield

.. code-block:: python

    conda install openbabel -c conda-forge # Install openbabel

.. code-block:: python

    conda install git # Install git

* Clone the *QMMMReBind* repository :

.. code-block:: python

    git clone https://github.com/anandojha/qmmmrebind.git

* Perform the following steps to get this package installed quickly on a local linux machine (Installation in the home directory is recommended) : 


.. code-block:: python

    cd qmmmrebind

.. code-block:: python

    python setup.py install

.. code-block:: python

    python setup.py test  # Optionally run tests to check for proper installation 

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
6. # B3LYP 6-31G OPT FREQ GUESS=READ INTEGRAL=(GRID=ULTRAFINE) SCF=(maxcycles=4000) SYMMETRY=NONE POP(MK,READRADII) IOP(6/33=2,6/42=6)

