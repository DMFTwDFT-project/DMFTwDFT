Installation
============
DMFTwDFT is currently supported by Python 2.x. Eventually, we will move on to Python 3.x.
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
		├── scripts
		├── sphinx
		└── sources
			├── src
			├── src_atomd
			├── src_ctqmc
			└── src_post-tools
				├── ksum
				├── maxent_source
				└── skrams

The following section describes the procesdure to compile the different componenets required to run DMFTwDFT. Once the executables and libraries have been compiled, they should be copied into the ``/bin`` directory.

Compiling library mode
----------------------

Compiling the Fortran file ``dmft_lib.F90`` available in the ``/sources/src/`` directory generates ``libdmft.a`` which can be used to link to DFT codes. In the case of VASP, add the location of this file to ``LLIBS`` in makefile.include. To compile: ::

	make lib

Compiling dmft.x
----------------

This module is for achieving DMFT self-consistency.
More details can found in the following article (Appendix B):

`<https://journals.aps.org/prb/pdf/10.1103/PhysRevB.90.235103>`_

It performs the :math:`k`-point sum and and computes the local Green's function (G_loc.out) and hybridization function (Delta.inp). ``atomd.py`` computes local atomic interactions and generates inputs required for the ``ctqmc`` run. These parts are developed by Kristjan Haule at Rutgers University.
The sources to compile dmft.x is located in ``/sources/``. To compile: ::

	make 

Compiling dmft_dos.x
--------------------

This module is used for calculating DMFT partial density of states for the correlated orbitals. The sources are located in ``/source/``. To compile: ::

	make dos

Compiling atomd
---------------

This module computes local atomic interactions and generates inputs required for the ctqmc impurity solver. The sources are located in ``/sources/src_atomd/``. To compile: ::

	make all

Compiling ctqmc
---------------

This is the continuous time quantum monte carlo code to solve the impurity problem. The sources are located in ``/sources/src_ctqmc``. To compile: ::

	make all

This will generate the Fortran libraries ``gaunt.so`` and ``gutils.so``.

Compiling post-tools
--------------------

This is a set of modules to perform DMFT post-processing with the ``postDMFT.py`` script. 
The sources are located under the directory ``/src_post-tools/``, namely ``ksum``, ``maxent_source`` and ``skrams``. To compile, cd into each directory and do: ::

	make all

Once the above executables and libraries are compiled, don't forget to copy them into the ``/bin`` directory.

DMFTwDFT requires ``wannier90.x`` to be in the ``/bin`` directory as well. You can get it from `<http://www.wannier.org/>`_.

Finally, the location of the ``/bin`` directory should be added to the ``$PATH`` and ``$PYTHONPATH`` environmental variables in your ``.bashrc``.

