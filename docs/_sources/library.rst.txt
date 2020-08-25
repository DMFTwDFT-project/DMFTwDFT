.. _labellibrary:

Library mode
============

DMFTwDFT consists of a library mode that can be used to interface other DFT codes to the framework to enable full charge self-consistent DFT+DMFT calculations. This calls Fortran subroutines to obtain information to update charge density within the DFT+DMFT loop. Specifically, one can pass the :math:`k`-points information within DFT to the subroutine ``dmft_nkij()`` and can obtain the DMFT Occupancy matrix matrix at each :math:`k`-point, :math:`n_{kij}` for computing the charge density :math:`\rho(r)`. The example of using this library can be found in library-mode-test directory. Instructions on how to compile this mode is provided in the installation section. 

For a detailed explanation of DMFTwDFT please refer to the documents in the ``manuals`` directory and the articles `DMFTwDFT <https://arxiv.org/abs/2002.00068>`_ and `PhysRevB.90.235103 <https://journals.aps.org/prb/pdf/10.1103/PhysRevB.90.235103>`_. 