#!/usr/bin/env python

""" Effective mass and Z calculator for DMFTwDFT.

Author: Dr. Hyowon Park

This program calculates the effective mass (m*/m) and quasi-particle residue (Z)
from the self energy calculated from DMFTwDFT.
The order of the columns in sig.inp are:

|---------------------| |-----------------------|----------------------------|
| Matsubara Frequency | | Real part of SE       | Imaginary part SE          |  x Repeats for each group
|---------------------| |-----------------------|----------------------------|

of "cor_orb".

Usage:

    Z.py -siglistindex <# of self energy files to avergae>
         -cor_orb_index <index of cor_orb>

"""


from scipy import *
import scipy.optimize
import sys
import argparse
from argparse import RawTextHelpFormatter
import glob


def fparab(par, x, data):
    chi2 = 0
    for i in range(len(x)):
        chi2 += ((par[0] + x[i] * par[1] + x[i] ** 2 * par[2]) - data[i]) ** 2
    return chi2


def fquad(par, x, data):
    chi2 = 0
    for i in range(len(x)):
        chi2 += (
            (
                par[0]
                + x[i] * par[1]
                + x[i] ** 2 * par[2]
                + x[i] ** 3 * par[3]
                + x[i] ** 4 * par[4]
            )
            - data[i]
        ) ** 2
    return chi2


def calculate_z(args):
    """Calculates Z and m*/m using imaginary self energy
    column of the self energy."""

    # copying the last few self-energies from the DMFT run in the directory above
    siglist = sorted(
        sorted(glob.glob("sig.inp.*"))[1:], key=lambda x: (len(x), float(x[8:]))
    )[-args.siglistindex :]

    Z_avg = 0
    ifit = 5
    ReS = 0
    comments = 5
    for i, fname in enumerate(siglist):
        fi = open(fname, "r")
        line = fi.readline()
        xfit = []
        yfit_r = []
        yfit_i = []
        for i in range(comments):
            fi.readline()
        for i in range(ifit):
            line = fi.readline()

            # Matsubara frequency
            xfit.append(float(line.split()[0]))

            # real column
            yfit_r.append(float(line.split()[2 * args.cor_orb_index - 1]))

            # imaginary column
            yfit_i.append(float(line.split()[2 * args.cor_orb_index]))

        # print xfit, yfit_r, yfit_i
        expan_r = scipy.optimize.fmin_powell(
            fquad, [0, 0, 0, 0, 0], args=(xfit, yfit_r), disp=False
        )
        expan_i = scipy.optimize.fmin_powell(
            fquad, [0, 0, 0, 0, 0], args=(xfit, yfit_i), disp=False
        )
        # print 'estimated derivatives (real): ', expan_r
        # print 'estimated derivatives (imag): ', expan_i
        # print "Z=",1/(1-expan_i[1])
        # print "m*/m=",(1-expan_i[1])
        Z_avg += 1 - expan_i[1]
        ReS += expan_r[0]
    print("\nAverage Z for cor_orb {:1d} = {:0.3f}".format(args.cor_orb_index, Z_avg))
    print(
        "Average m*/m for cor_orb {:1d}  = {:0.3f}".format(
            args.cor_orb_index, Z_avg / (len(siglist) - 1)
        )
    )
    print(
        "Average Re(Sigma) for cor_orb {:1d} = {:0.3f}\n".format(
            args.cor_orb_index, ReS / (len(siglist) - 1)
        )
    )


if __name__ == "__main__":
    args = sys.argv[1:]
    if args:

        parser = argparse.ArgumentParser(
            description=__doc__, formatter_class=RawTextHelpFormatter
        )
        parser.add_argument(
            "-siglistindex",
            type=int,
            default=1,
            help="How many last self energy files to average?",
        )
        parser.add_argument(
            "-cor_orb_index",
            type=int,
            default=1,
            help="cor_orb index (starting from 1).",
        )

        args = parser.parse_args()
        calculate_z(args)
    else:
        print("Usage: Z.py -h")
