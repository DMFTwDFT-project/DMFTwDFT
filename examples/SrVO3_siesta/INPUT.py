import os

################   Input parameters for DFT+DMFT calculations   ##################

###### Main loop parameters ###########

p = {
    "Niter": 1,  # Number of DFT+DMFT iterations
    "Ndft": 1,  # Number of DFT iterations
    "Nit": 2,  # Number of DMFT iterations
    "n_tot": 19,  # Number of total electrons
    "nf": 1.0,  # Number of target Nd
    "nspin": 1,  # Number of total spins
    "atomnames": ["V", "O"],  # The name of atoms
    "orbs": ["d", "p"],  # The name  of orbitals
    "L_rot": [0, 0],  # Whether rotate local axis or not
    "cor_at": [["V1"]],  # Correlated atoms, put symmetry atoms in the same list
    "cor_orb": [
        [["d_z2", "d_x2y2"], ["d_xz", "d_yz", "d_xy"]]
    ],  # DMFT orbitals, other orbitals are treated by HF
    "U": [5.0],  # Intra-U for each cor_at
    "J": [1.0],  # Hund's coupling
    "alpha": [0.2],  # Double counting parameter
    "mix_sig": 0.2,  # Mixing parameter for Sigma
    "q": [24, 24, 24],  # Number of k-points for Wannier k-sum
    "ewin": [-8, 6],  # Energy Window with respect to Fermi energy
    "noms": 400,  # Number of Matsubara frequencies for k-sum
    "dc_type": 1,  # Vdc type
    "mu_iter": 100,  # Steps for the chemical potential convergence
    "Nd_qmc": 1,  # 0: Use Nd_latt, 1: Use Nd_imp
    "path_bin": "~/Dropbox/git/DMFTwDFT/bin/",  # Path to bin files
}


##### CTQMC parameters #########
pC = {
    "exe": ["ctqmc", "# Name of impurity solver"],
    "beta": [30.0, " # Inverse temperature"],
    "M": [5000000, " # Number of Monte Carlo steps"],
    "nom": [100, " # number of Matsubara frequency points to sample"],
    "svd_lmax": [30, "# number of SVD functions to project the solution"],
    "tsample": [200, "# how often to record the measurements"],
    "aom": [1, " # number of frequency points to determin high frequency tail"],
    "warmup": [250000, "  # Warmup"],
    "GlobalFlip": [1000000, "  # Global flip"],
    "Ncout": [500000, " # Ncout"],
    "Naver": [80000000, "  # Naver"],
}


############## CIX parameters ###########
pD = {
    "para=": 1,
    "l=": 2,
    "n=": [6, 7, 8],
    "OCA_G=": False,
    "CoulombF=": "Ising",
}
