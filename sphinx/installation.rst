Installation
============

DMFTwDFT is currently supports Python 2.x. It will be transitioned to Python 3.x in the future. 
The Python 2.x version does not support ionic relaxation, therefore use an already relaxed structure as input. 

Please install the following dependencies prior to installing DMFTwDFT. 

* matplotlib 
* numpy 
* scipy
* weave
* pychemia
* mpi4py


The structure of the directories is as follows. ::

	DMFTwDFT
		├── bin
		├── docs
		├── examples
		├── reference 
		├── scripts
		├── sphinx
		└── sources
			├── src
			├── dmft_ksum
			├── fort_kpt_tools
			└── CSC-mods

The following section describes the procedure to compile the different componenets required to run DMFTwDFT. Once the executables and libraries have been compiled, they should be inside the ``bin`` directory.

Compiling sources
-----------------

From the ``sources`` directory execute::

	make all

This should compile the follwing executables and libraries and copy them to the ``bin`` directory.

* dmft.x - This module is for achieving DMFT self-consistency. It performs the :math:`k`-point sum and and computes the local Green's function (G_loc.out) and hybridization function (Delta.inp).
* dmft_dos.x - Performs DOS calculation. 
* dmft_ksum_band - Performs band structure calculation. 
* dmft_ksum_partial_band - Performs projected band structure calculation. 
* fort_kpt_tools.so - Fortran based k-points tool.


Compiling library mode
----------------------

The above compilation also generates ``libdmft.a`` which can be used to link DMFTwDFT to DFT codes to enable full charge self-consistent DFT+DMFT calculations. In the case of VASP, add the location of this file to ``LLIBS`` in VASP's makefile.include. The modified VASP files are in the ``sources/CSC-mods`` directory. Copy these to the VASP source directory and recompile VASP. 


External libraries and executables
----------------------------------

DMFTwDFT uses the CTQMC impurity solver and Max-entropy routines developed by Kristjan Haule at Rutgers University and are available in the eDMFT package.
Follow the provided instructions on the `web page <http://hauleweb.rutgers.edu/tutorials/index.html>`_ to download and compile it. 

Once compiled copy the following executables and libraries to the DMFTwDFT ``bin`` directory.

* ctqmc
* gaunt.so
* gutils.so
* skrams
* maxent_routines.so

If the automated compilation is successful, these are found in the eDMFT ``bin`` directory. Otherwise compile them manually inside the ``src/impurity`` directories. gaunt.so and gutils.so are inside the 

Wannier90 library
-----------------

DMFTwDFT requires ``wannier90.x`` and ``w90chk2chk.x`` to be in the ``bin`` directory. You can get them from `<http://www.wannier.org/>`_. VASP should be recompiled with the Wannier90 library.

PATH variables
--------------

Finally, the location of the DMFTwDFT ``bin`` directory should be added to the ``$PATH`` and ``$PYTHONPATH`` environmental variables in your ``.bashrc``.

