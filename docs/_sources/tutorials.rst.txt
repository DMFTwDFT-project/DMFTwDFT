Tutorials
=========

The following set of tutorials explain usage of DMFTwDFT. Example files required to run these calculations are available in the /examples directory in the github repo. 
To perform a DFT+DMFT calculation, the following files should be present.

* INPUT.py : contains the input parameters that govern the DMFT calculation. 
* DFT_mu.out : This is a guess for the DFT Fermi energy. 
* para_com.dat : The number of processors used for the DMFT calculation. E.g. mpirun -np 32
* para_com_dft.dat : The number of processors used for the DFT calculation. If not specified it will use the same as para_com.dat. 
* DFT files: VASP = {POSCAR, KPOINTS, POTCAR, INCAR}, Siesta = {.fdf, .psf}	

Before you start remember to add the ``bin`` directory path in ``INPUT.py`` as the value for the key ``path_bin``.
Eg.::

	"path_bin":"~/Dropbox/Research/Projects/DMFTwDFT/bin/"


DMFT calculation
----------------

This performs the DFT + DMFT calculations through the script ``DMFT.py``.
To get a description of the options, run: ::
	
	DMFT.py -h

This script has the following options.

* -dft:
	The choice of DFT code. Currently, VASP and Siesta are supported.

* -relax:
	This flag turns on DFT convergence testing. If the forces are not converged a convergence calculation is attempted and if it fails the user is asked to modify convergence parameters. This is useful for vacancy and defect calculations where a force convergence is required after the vacancy or defect is created in order to obtain a relaxed structure to perform DFT+DMFT with. Currently supported for VASP. This uses PyChemia to check for convergence. The relaxation occurs inside a  ``DFT_relax`` directory. 

* -structurename:
	DFT codes such as Siesta uses input files that contain the name of the system e.g. :math:`SrVO_3.fdf`. Therefore when performing DFT+DMFT calculations with Siesta this flag is required.

* -dmft:
	This flag performs the DMFT calculation using the results from the DFT calculation if a previous DMFT calculation in the same directory is incomplete. 

* -hf:
	This flag performs the Hartree-Fock (HF) calculation to the correlated orbitals specified in INPUT.py if a previous HF calculation in the same directory is incomplete. 

* -force: 
	This flag forces a DMFT or HF calculation even if a previous calculation has been completed. The option to check for completeness is helpful when running many DMFT/HF jobs on a cluster.

* -kmeshtol:
	This controls the tolerance of two k-points belonging to the the same shell in the wannier90 calculation. 	

The calculations are performed in an automatically generated ``DMFT`` or ``HF`` directory where the script was run from. 

E.g.: ::

	$DMFT.py -dft vasp -relax -dmft
	$DMFT.py -dft siesta -structurename SrVO3 -dmft
	
DMFT Post-processing
--------------------

DMFT post-processing is performed with the script ``postDMFT.py``. 
To get a description of the options, run: ::
	
	postDMFT.py -h

This script performs anaytical continuation, density of states and band structure calculations on the DMFT/HF data. Once the DMFT/HF calculations are complete, this script should be initiated within the ``DMFT`` or ``HF`` directories. This script has the following options.

* ac:
	This function performs the Analytic Continuation to obtain the Self Energies on the real axis. It has the option ``-siglistindx`` to specify the last number of Self Energy files (sig.inp) to average for the calculation. 

* dos:
	This function performs the partial density of states of the correlated orbitals. It has the following options:

	- -emin : Minimum energy value for interpolation
	- -emax : Maximum energy value for interpolation
	- -rom : Number of Matsubara Frequencey (:math:`\omega`) points
	- -broaden : Broadening of the dos
	- -show : Display the density of states 
	- -elim : The energy range to plot

* bands:
	This function performs the DMFT band structure calculations. It has the following options:
	
	- -emin : Minimum energy value for interpolation
	- -emax : Maximum energy value for interpolation
	- -rom : Number of Matsubara Frequencey ($\omega$) points
	- -kpband : Number of k-points for band structure calculation
	- -kn : A list of labels for k-points
	- -kp : A list of k-points corresponding to the the k-point labels
	- -plotplain : Flag to plot a plain band structure
	- -plotpartial : Flag to plot a projected band structure
	- -wo : List of Wannier orbitals to project onto the band structure
	- -vlim : Spectral intensity range
	- -show : Display the bands

The projected bands are especially helpful in determining the contribution to bands from different orbitals. 
The calculations are stored in directories ac, dos and bands, respectively. 
The following are some example commands to perform post-processing.

E.g.: ::

	$postDMFT.py ac -siglistindx 4
	$postDMFT.py dos -show
	$postDMFT.py bands -plotplain
	$postDMFT.py bands -plotpartial -wo 4 5 6