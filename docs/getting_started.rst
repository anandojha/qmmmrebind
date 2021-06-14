Getting Started
===============

This page details how to get started with QMMMReBind. 

########################
Introduction and Background Information
########################


Understanding the interaction of metabolites or drugs with biomolecules can help improve our understanding of drug discovery and development. Accurate computational prediction of kinetic rates or residence times of an intermolecular encounter can help us identify lead drugs. Kinetic binding parameters are often correlated with efficacy rather than binding affinities [1,2]. The binding kinetic profile of a drug molecule is characterized by the bimolecular association rate constant (k\ :sub:`on`\)  and the dissociation rate constant (k\ :sub:`off`\).  Drug-Target residence time (1/ k\ :sub:`off`\)  has received the most attention from the drug discovery community since drugs with long-lasting target occupancy are often correlated with greater in-vivo efficacy [3]. k\ :sub:`on`\  is also useful for predicting in-vivo efficacy, particularly related to understanding drug rebinding. Both k\ :sub:`on`\  and k\ :sub:`off`\  can be used to compute the binding free energy. A complete kinetic profile (k\ :sub:`on`\  and k\ :sub:`off`\), in addition to binding free energy, is extremely desirable for the prediction of in-vivo efficacy and further optimization of lead compounds. 


Our laboratory has recently developed a multiscale milestoning simulation approach to estimate receptor-ligand binding kinetics computationally [4,5]. This tool, called “Simulation Enabled Estimation of Kinetic Rates” (SEEKR), incorporates the multiscale and parallel implementation of molecular dynamics (MD) and Brownian dynamics (BD) simulations using the milestoning approach to calculate the association and dissociation rates of receptor-ligand complexes. This approach requires orders of magnitude less simulation time than classical MD simulations and comparable or less simulation time than other enhanced sampling techniques. SEEKR has demonstrated successes for calculating receptor-ligand binding association (k\ :sub:`on`\)  and dissociation rates (k\ :sub:`off`\)  for multiple systems (such as the well-studied system of the protease, trypsin, with the noncovalent binder, benzamidine) as well as for rank-ordering a series of small molecules by dissociation rates and binding free energies [6]. SEEKR is among the few simulation approaches that can obtain an entirely computational estimate of the binding kinetic (k\ :sub:`on`\  and k\ :sub:`off`\)  and thermodynamic (G\ :sub:`bind`\) profiles and shows good agreement with experiment, often using less simulation time than other approaches and requiring no biasing or reweighting of its simulations. 

Significant challenges for pharmaceutically relevant receptor-ligand systems include the size and flexibility of the ligands, large-scale conformational rearrangements, and the need for extensive sampling associated with these events. Timescales of these rearrangements are much longer than can be sampled adequately with classical MD simulations. SEEKR calculations using the original implementation also struggle to simulate these more complex cases. As a consequence of limited MD sampling of rare transitions between these states, which are critical for describing the binding and unbinding event, long-timescale MD simulations are required to sample distributions on each milestone. In an effort to increase the efficiency and accuracy of kinetics calculations, Markovian Milestoning with Voronoi Tessellations (MMVT) has been implemented in SEEKR [7]. In this milestoning scheme, milestones can not only subdivide the distance of the ligand from the binding site, but also the other slow degrees of freedom in the system, such as ligand orientation or protein loop and hinge motions. In addition to aiding in the sampling of rare events with the placement of additional milestones, MMVT reduces the simulation time needed for SEEKR calculations as it overcomes the sampling bottleneck associated with the previous implementation, obtaining an equilibrium distribution on each milestone. Kinetics can then be obtained directly from short, parallel simulations within each Voronoi cell.


MMVT-SEEKR holds much potential for additional improvements in sampling, reduction of simulation times, and accuracy. MMVT-SEEKR currently incorporates MD simulations through NAMD and BD simulations through Browndye [8]. OpenMM is an increasingly popular and effective MD engine and is well-suited for running MD calculations using graphical processing units (GPUs) which are significantly faster than single-core CPU implementations, although serial and multithread CPU computations are also possible [9]. OpenMM offers all of the most common MD simulation capabilities, including a wide variety of integration schemes, compatibility with AMBER and CHARMM forcefields, and a high degree of customizability with forces and constraints within the simulation system. The Amaro lab is in the final stages of developing “SEEKR2,” which is a plugin for the molecular dynamics toolkit OpenMM to perform MMVT simulations. Implementing SEEKR within OpenMM not only significantly improved performance benchmarks (due to GPU speedups) but also allowed us to seamlessly connect with other tools within the MolSSI network.


Our laboratory previously demonstrated the effectiveness of SEEKR in predicting and ranking a series of seven small-molecule compounds for the model system , \beta\-cyclodextrin [6]. This ranking was based on estimating  k\ :sub:`on`\  and k\ :sub:`off`\   of seven host-guest systems. Although results were in good agreement with the previously conducted long timescale MD simulations for the same set of ligands with the same forcefield (GAFF and Q4MD), both methods failed in determining the correct orders for the k\ :sub:`on`\ 's [10]. Predicted k\ :sub:`off`\'s  also had deviations from experimental values although the rankings for the k\ :sub:`off`\'s  were accurate. 

########################
Software Requirements
########################

Please make sure to install these packages before running the QMMMReBind:

* Gaussian16
* TorsionDrive
* Psi4


########################
Installation and Setup Instructions
########################

* Make sure `anaconda3 <https://www.anaconda.com/>`_ is installed on the local machine. 
* Go to the `download <https://www.anaconda.com/products/individual>`_  page of anaconda3 and install the latest version of anaconda3. 
* Create a new conda environment with python = 3.8 with the following commands in the terminal: 

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

