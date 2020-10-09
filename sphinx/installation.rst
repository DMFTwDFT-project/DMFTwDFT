Installation
============

DMFTwDFT is currently supported by Python 2.x. It will be transitioned to Python 3.x in the future. 
The Python 2.x version does not support ionic relaxation, therefore use an already relaxed structure as input. 

Please install the following dependencies prior to installing DMFTwDFT.

Python:

* matplotlib 
* numpy 
* scipy
* weave
* mpi4py
* pyprocar
* pychemia

To install execute::

	pip install matplotlib numpy scipy weave mpi4py pyprocar pychemia

NOTE:

PyProcar is used as a support package to find the wannier window of the correlated orbitals. It's not necessary to be installed for DMFTwDFT.
PyChemia is required for the Python 3 version of the code which is still under development. 
For the Python 2 version this is not required. 

Libraries:

* GSL

If not using intel compilers the following libraries should be installed:

* LAPACK
* BLAS
* FFTW


The structure of the directories is as follows. ::

	DMFTwDFT
		├── bin
		├── config 
		├── docs
		├── examples
		├── manuals 
		├── scripts
		├── sphinx
		├── support_packages
		└── sources
			├── src
			├── dmft_ksum
			├── fort_kpt_tools
			└── CSC-mods

The following section describes the procedure to compile the different componenets required to run DMFTwDFT. Once the executables and libraries have been compiled, they should be inside the ``bin`` directory.

Compiling sources
-----------------

First copy a desired ``Makefile.in`` version in the **config** directory to the root directory based on the compiler you wish to use. You may have to specify the locations of the gsl, lapack etc. libraries (default: /usr/local/bin).
Similarly, in the **sources** directory modify intel.make.inc or gfortran.make.inc accordingly. 

To compile, execute::

	python setup.py <compiler> 

where compiler = {intel, gfortran}

This should compile the follwing executables and libraries and copy them to the ``bin`` directory.

* dmft.x - This module is for achieving DMFT self-consistency. It performs the :math:`k`-point sum and and computes the local Green's function (G_loc.out) and hybridization function (Delta.inp).
* dmft_dos.x - Performs DOS calculation. 
* dmft_ksum_band - Performs band structure calculation. 
* dmft_ksum_partial_band - Performs projected band structure calculation. 
* fort_kpt_tools.so - Fortran based k-points calculation module.

External libraries and executables
----------------------------------

DMFTwDFT uses the CTQMC impurity solver and Max-entropy routines developed by Professor Kristjan Haule at Rutgers University available with the `EDMFTF <http://hauleweb.rutgers.edu/tutorials/index.html>`_ package. The following libraries and programs are used:

* ctqmc
* gaunt.so
* gutils.so
* skrams
* maxent_routines.so

If the automated compilation with setup.py is successful, these are found in ``bin`` directory. 

Wannier90 library
-----------------

DMFTwDFT requires ``wannier90.x`` and ``w90chk2chk.x`` to be in the ``bin`` directory. You can get them from `<http://www.wannier.org/>`_. VASP should be recompiled with the Wannier90 library.

PATH variables
--------------

Finally, the location of the DMFTwDFT ``bin`` directory should be added to the ``$PATH`` and ``$PYTHONPATH`` environmental variables in your ``.bashrc``.

Compiling library mode for full charge-self conistent DFT+DMFT calculations (OPTIONAL)
--------------------------------------------------------------------------------------

The above compilation also generates ``libdmft.a`` which can be used to link DMFTwDFT to DFT codes to enable full charge self-consistent DFT+DMFT calculations. Otherwise, the calculations can only be run for self-consistency within DMFT. In the case of VASP please follow the below steps to compile for self-consistency.

1. Generate libdmft.a by the Compiling sources code that will  be needed to link DMFTwDFT to DFT codes to enable full charge self-consistent DFT+DMFT calculations. 

2. Change the VASP makefile.include file. Specify libraries and/or objects to be linked against, in the usual ways::

	LLIBS += -Lparser -lparser -lstdc++ /home/uthpala/wannier90/wannier90-1.2/libwannier.a
         /home/uthpala/Dropbox/git/DMFTwDFT/sources/libdmft.a

3. Before modifying the source code of VASP, we first need to install the VASP as it is. The user should follow the VASP installation instructions from the VASP web site.

4. First copy the modified mlwf.F VASP file from the sources/CSC-mods directory to the VASP source directory and install the VASP again. This step will create some dependenices for next step.

5. Next, copy the other modified/required VASP files such as charge.F,  electron.F,  main.F, and us.F from the sources/CSC-mods directory to the VASP source directory.

6. Finally, recompile VASP. Then rename this vasp executable to vaspDMFT and copy it to teh DMFTwDFT/bin directory.

More information on the library mode can be found in the :ref:`labellibrary` section.

