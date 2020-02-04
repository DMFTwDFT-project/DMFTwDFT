.. DMFTwDFT documentation master file, created by
   sphinx-quickstart on Mon Aug 12 16:32:02 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

DMFTwDFT documentation
======================
DMFTwDFT is an open-source, user-friendly framework to calculate electronic, vibrational and elastic properties in strongly correlated materials (SCM) using beyond-DFT methods such as DFT+U, DFT+Hybrids and DFT+DMFT (Dynamical Mean Field Theory) with a variety of different DFT codes. Currently VASP and Siesta are supported. 

DMFTwDFT consists of two main segments.

1. DMFT.py     - Performs the DFT and DMFT calculations. 
2. postDMFT.py - Performs post-processing including analytic contiuation, density of states and band structures.

The ``/scripts`` directory contains several utility scripts. 



Library mode
------------

DMFTwDFT consists of a library mode that can be used to interface other DFT codes to the framework. This calls Fortran subroutines to obtain information to update charge density within the DFT+DMFT loop. Specifically, one can pass the :math:`k`-points information within DFT to the subroutine ``Compute_DMFT()`` and can obtain DMFT weight and Unitary matrix at each :math:`k`-point for computing the charge density :math:`\rho(r)`. The example of using this library can be found in library-mode-test directory. Instructions on how to compile this mode is provided in the installation section. 

For a detailed explanation of DMFTwDFT please refer to:



.. toctree::
   :maxdepth: 2
 
   installation 
   tutorials
   contributors
   cite



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
