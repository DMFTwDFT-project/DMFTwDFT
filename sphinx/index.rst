.. DMFTwDFT documentation master file, created by
   sphinx-quickstart on Mon Aug 12 16:32:02 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

DMFTwDFT documentation
======================

DMFTwDFT is an open-source, user-friendly framework to calculate electronic, vibrational and elastic properties in strongly correlated materials (SCM) using beyond-DFT methods such as DFT+U, DFT+Hybrids and DFT+DMFT (Dynamical Mean Field Theory) with a variety of different DFT codes. Currently VASP, Siesta and Quantum Espresso (through Aiida) are supported. 

========
Overview
========

DMFTwDFT offers the following:

1. DMFT has been one of the most successful methods treating many-body ﬂuctuations, by including dynamical but local correlations beyond the static DFT exchange-correlation functional. We provide a DMFT package, DMFTwDFT, interfaced with the Wannier90 code for its efficient extension to various free-licensed DFT codes.

2. We provide the library mode to link the module for computing a DMFT density matrix and updating a charge density within the DFT loops without modifying any DFT source codes signiﬁcantly.

3. We provide a ﬂexible Python-based interface that does not rely much on extensive user experience or speciﬁc parameters to perform the DFT+DMFT calculations of strongly correlated materials.

4. To test and check our implementation, we computed the density of states and the band structure of well-known correlated materials, namely :math:`LaNiO_{3}`, :math:`SrVO_{3}`, and :math:`NiO`. The obtained results are compared to those obtained from other DFT+DMFT implementations. 

5. In the next release, we envision implementing force calculations which will help us to perform phonon calculations in strongly correlated materials, and also implementing the ab-initio Hubbard U calculation using the linear response approach developed by Cococcioni et al.

DMFTwDFT consists of two main scripts to perform the DFT+DMFT calculations.

1. DMFT.py     - Performs the DFT and DMFT calculations. 
2. postDMFT.py - Performs post-processing including analytic contiuation, density of states and band structures.

The ``/scripts`` directory contains several utility scripts. 


=================
Table of Contents
=================

.. toctree::
   :maxdepth: 2
 
   installation 
   tutorials
   library
   developers
   cite



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
