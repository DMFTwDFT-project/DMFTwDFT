#!/usr/bin/env python
"""
This checks for Hermitivity of the DMFT Occupancy matrix n_kij
obtained from the dmft_nkij subroutine in the dmft.F90 test module
in DMFTwDFT. n_kij is read from the file dmft-nkij.dat.

shape of n_kij : [iband, jband, kpoints]

The way nkij is saved in dmft-nkij.dat is with jband being the
fastest index, iband being the next fastest and kpoints being
the slowest index.

Usage:

hermitiancheck.py numberofbands numberofkpoints

"""
import numpy as np
import sys

numberofbands = int(sys.argv[1])
numberofkpoints = int(sys.argv[2])

fi = open("dmft-nkij.dat", "r")
header = fi.readline()
data = fi.readlines()
fi.close()

n_kij = np.zeros((numberofbands, numberofbands, numberofkpoints), dtype="complex")

row_counter = 0
for k in range(n_kij.shape[2]):
    for i in range(n_kij.shape[0]):
        for j in range(n_kij.shape[1]):
            n_kij[i, j, k] = float(data[row_counter].split()[4]) + 1j * float(
                data[row_counter].split()[-1].split(")")[0]
            )
            row_counter += 1

# Check for Hermitivity for the (iband, jband) matrix for each
# k-point.
hermitiancheck = []
diagonalcheck = []
for kk in range(n_kij.shape[2]):
    hermitiancheck.append(np.allclose(n_kij[:, :, kk], np.conj(n_kij[:, :, kk]).T))

    # Check if diagonal of each (iband, jband) matrix for each k-points is real.
    # Just a secondary check. Checking for Hermitian already satsifies this.
    diagonalcheck_perkpoint = []
    for jj in range(n_kij.shape[0]):
        diagonalcheck_perkpoint.append(np.isclose(n_kij[jj, jj, kk].imag, 0))
    if all(diagonalcheck_perkpoint):
        diagonalcheck.append(True)
    else:
        diagonalcheck.append(False)

if all(diagonalcheck):
    print("\nThe diagonal of the  (i,j) matrix for each k-point is real.")
else:
    print("\nThe diagonal of the (i,j) matrix for each k-point is NOT real!")
    print(
        "Complex diagonal arrays at k-points: %s \n "
        % [index for index, value in enumerate(diagonalcheck) if value is False]
    )

if all(hermitiancheck):
    print("Hermitian condition of the (i,j) matrix for all k-points is satisfied.")

else:
    print("Hermitian condition of the (i,j) matrix for all k-points is NOT satisfied!")
    print(
        "Non-Hermitian arrays at k-points: %s "
        % [index for index, value in enumerate(hermitiancheck) if value is False]
    )