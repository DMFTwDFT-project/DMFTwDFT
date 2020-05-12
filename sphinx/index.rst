.. DMFTwDFT documentation master file, created by
   sphinx-quickstart on Mon Aug 12 16:32:02 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

DMFTwDFT documentation
======================

DMFTwDFT is an open-source, user-friendly framework to calculate electronic, vibrational and elastic properties in strongly correlated materials (SCM) using beyond-DFT methods such as DFT+U, DFT+Hybrids and DFT+DMFT (Dynamical Mean Field Theory) with a variety of different DFT codes. Currently VASP, Siesta and Quantum Espresso (through Aiida) are supported. 

DMFTwDFT consists of two main segments.

1. DMFT.py     - Performs the DFT and DMFT calculations. 
2. postDMFT.py - Performs post-processing including analytic contiuation, density of states and band structures.

The ``/scripts`` directory contains several utility scripts. 



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
