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

hermitiancheck.py numberofbands numberofkpoints <optional: filename>

"""
import numpy as np
import sys
import scipy.linalg as la

numberofbands = int(sys.argv[1])
numberofkpoints = int(sys.argv[2])

if len(sys.argv) > 3:
    fname = sys.argv[3]
else:
    fname = "dmft-nkij.dat"

fi = open(fname, "r")
header = fi.readline()
data = fi.readlines()
fi.close()

n_kij = np.zeros((numberofbands, numberofbands, numberofkpoints), dtype="complex")

k_weights = []

row_counter = 0
for k in range(n_kij.shape[2]):
    k_weights.append(float(data[row_counter].split()[3]))

    for i in range(n_kij.shape[0]):
        for j in range(n_kij.shape[1]):
            n_kij[i, j, k] = float(data[row_counter].split()[5]) + 1j * float(
                data[row_counter].split()[-1].split(")")[0]
            )
            row_counter += 1

# printing the eigen values for the (iband, jband) array at each k-point

# print("\nEigenvalues for each k-point : ")

fp = open("dmft-nkij_eig_scipy.dat", "w")

fp.write("k, wk, band, n_kij_eig\n")
for ik in range(n_kij.shape[2]):
    eigenvalues = la.eig(n_kij[:, :, ik])
    # print("\nkpoint : %s" % (ik + 1))
    # print(eigenvalues[0].real)

    for ib in range(numberofbands):
        fp.write(
            "{:d} \t {:.4f} \t {:3d} \t {:6.4f}\n".format(
                ik + 1, k_weights[ib], ib + 1, eigenvalues[0][ib].real
            )
        )
fp.close()


# Check for Hermitivity for the (iband, jband) matrix for each
# k-point.
hermitiancheck = []
diagonalcheck = []
diagonal_values = np.zeros((numberofkpoints, numberofbands), dtype="complex")

for kk in range(n_kij.shape[2]):
    hermitiancheck.append(np.allclose(n_kij[:, :, kk], np.conj(n_kij[:, :, kk]).T))

    # Check if diagonal of each (iband, jband) matrix for each k-points is real.
    # Just a secondary check. Checking for Hermitian already satsifies this.
    diagonalcheck_perkpoint = []
    for jj in range(n_kij.shape[0]):
        diagonalcheck_perkpoint.append(np.isclose(n_kij[jj, jj, kk].imag, 0))
        diagonal_values[kk, jj] = n_kij[jj, jj, kk]
    if all(diagonalcheck_perkpoint):
        diagonalcheck.append(True)
    else:
        diagonalcheck.append(False)

fp = open("dmft-nkij_diagonal.dat", "w")
fp.write("k, band, n_kij_diagonal\n")
for ik in range(n_kij.shape[2]):
    for ib in range(numberofbands):
        fp.write(
            "{:d} \t {:3d} \t {:6.4f}\n".format(
                ik + 1, ib + 1, diagonal_values[ik, ib].real
            )
        )
fp.close()


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
