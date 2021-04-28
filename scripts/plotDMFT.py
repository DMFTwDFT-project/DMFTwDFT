#!/usr/bin/env python

"""Plotter for DMFTwDFT Green's functions and Self Energies.

Author: Uthpala Herath

This script plots the real and imaginary parts of the Green's function (G_loc.out)
the self-energy (sig.inp) and analytically continued self-energy (Sig.out)
calculated from the DMFTwDFT code.

The order of the columns are:

|---------------------| |-----------------------|----------------------------|
| Matsubara Frequency | | Real part of SE       | Imaginary part SE          | x Repeats for each group
|---------------------| |-----------------------|----------------------------|   of "cor_orb".

Run this script inside the DMFT folder.

Usage:

    plotDMFT.py
            -siglistindex <# of self energy files to average>
            -cor_orb_index <List of cor_orb indexes >
            -cor_orb_labels <Names of cor_orb's>

E.g.- for correlated d-orbitals considering eg and t2g i.e.
"cor_orb":   [[['d_z2','d_x2y2'],['d_xz','d_yz','d_xy']]]

plotDMFT.py -siglistindex 5 -cor_orb_index 1 2 -cor_orb_labels '$e_g$' '$t_{2g}$'


NOTE: siglistindex is only used to average the sig.inp files. Not used for Green's functions
(G_loc) or analytically continued self energy (Sig.out).
Sig.out should be present inside the ac directory  to plot the analytically continued self energy.
"""


import matplotlib

matplotlib.use("pdf")
import matplotlib.pyplot as plt
import glob
import os
import sys
import shutil
import argparse
import subprocess
from argparse import RawTextHelpFormatter
import warnings

warnings.filterwarnings("ignore", module="matplotlib\..*")


# Setting up plotting class
plt.rcParams["mathtext.default"] = "regular"  # Roman ['rm', 'cal', 'it', 'tt', 'sf',
#                                                   'bf', 'default', 'bb', 'frak',
#                                                   'circled', 'scr', 'regular']
plt.rcParams["font.family"] = "Arial"
plt.rc("font", size=22)  # controls default text sizes
plt.rc("axes", titlesize=22)  # fontsize of the axes title
plt.rc("axes", labelsize=22)  # fontsize of the x and y labels
plt.rc("xtick", labelsize=22)  # fontsize of the tick labels
plt.rc("ytick", labelsize=22)  # fontsize of the tick labels
# plt.rc('legend', fontsize=22)    # legend fontsize
# plt.rc('figure', titlesize=22)  # fontsize of the figure title


def plot_dmft(args):
    """Function to plot the Green's functions and self
    energies."""

    if not os.path.exists("plots"):
        os.makedirs("plots")

    ########## Plotting Green's functions #################

    # Unlike, self energy it will only take the last G_loc file without averaging
    # multiple files.
    gf_file = sorted(
        sorted(glob.glob("G_loc.out.*")), key=lambda x: (len(x), float(x[10:]))
    )[-1]

    # setting up figures
    fig1 = plt.figure(figsize=(13, 9))
    ax1 = fig1.add_subplot(111)

    fig2 = plt.figure(figsize=(13, 9))
    ax2 = fig2.add_subplot(111)

    with open(gf_file, "r") as f:
        lines = f.readlines()
        for i, ci in enumerate(args.cor_orb_index):

            x = [float(line.split()[0]) for line in lines]
            y_real = [float(line.split()[2 * ci - 1]) for line in lines]
            y_imag = [float(line.split()[2 * ci]) for line in lines]

            label_real = r"Re ({})".format(str(args.cor_orb_labels[i]))
            label_imag = r"Im ({})".format(str(args.cor_orb_labels[i]))

            ax1.plot(x, y_real, label=label_real)
            ax2.plot(x, y_imag, label=label_imag)

    ax1.set_title("Green's function (Re)")
    ax1.set_xlabel("$i{\omega_n}$")
    ax1.set_ylabel("Re $G(i{\omega_n})$")
    ax1.set_xlim(xmin=0)
    ax1.legend()
    fig1.savefig("./plots/gf-real.pdf")

    ax2.set_title("Green's function (Im)")
    ax2.set_xlabel("$i{\omega_n}$")
    ax2.set_ylabel("Im $G(i{\omega_n})$")
    ax2.set_xlim(xmin=0)
    ax2.legend()
    fig2.savefig("./plots/gf-imag.pdf")

    f.close()

    ########## Plotting Self Energies #################

    se_files = sorted(
        sorted(glob.glob("sig.inp.*"))[1:], key=lambda x: (len(x), float(x[8:]))
    )[-args.siglistindex :]

    # Averaging self-energies with sigaver.py
    if not os.path.exists("./plots/tmp"):
        os.makedirs("./plots/tmp")
    else:
        shutil.rmtree("./plots/tmp/")
        os.makedirs("./plots/tmp")
    for fi in se_files:
        shutil.copy(fi, "./plots/tmp/")

    # Averaging self energies
    print("\nAveraging self-energies from: ")
    print(se_files)
    cmd = "cd plots/tmp/ && sigaver.py sig.inp.*"
    out, err = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    ).communicate()
    if err:
        print(err)
        sys.exit()

    # setting up figures
    fig1 = plt.figure(figsize=(13, 9))
    ax1 = fig1.add_subplot(111)

    fig2 = plt.figure(figsize=(13, 9))
    ax2 = fig2.add_subplot(111)

    with open("./plots/tmp/sig.inpx", "r") as f:
        for i in range(2):
            f.readline()
        lines = f.readlines()
        for i, ci in enumerate(args.cor_orb_index):

            x = [float(line.split()[0]) for line in lines]
            y_real = [float(line.split()[2 * ci - 1]) for line in lines]
            y_imag = [float(line.split()[2 * ci]) for line in lines]

            label_real = r"Re ({})".format(str(args.cor_orb_labels[i]))
            label_imag = r"Im ({})".format(str(args.cor_orb_labels[i]))

            ax1.plot(x, y_real, label=label_real)
            ax2.plot(x, y_imag, label=label_imag)

    ax1.set_title("Self energy (Re)")
    ax1.set_xlabel("$i{\omega_n}$")
    ax1.set_ylabel("Re $\Sigma(i{\omega_n})$")
    ax1.set_xlim(xmin=0)
    ax1.legend()
    fig1.savefig("./plots/sigma-real.pdf")

    ax2.set_title("Self energy (Im)")
    ax2.set_xlabel("$i{\omega_n}$")
    ax2.set_ylabel("Im $\Sigma(i{\omega_n})$")
    ax2.set_xlim(xmin=0)
    ax2.legend()
    fig2.savefig("./plots/sigma-imag.pdf")

    f.close()

    ########## Plotting AC Self Energies #################

    if os.path.exists("./ac/Sig.out"):

        # setting up figures
        fig1 = plt.figure(figsize=(13, 9))
        ax1 = fig1.add_subplot(111)

        fig2 = plt.figure(figsize=(13, 9))
        ax2 = fig2.add_subplot(111)

        with open("./ac/Sig.out", "r") as f:
            for i in range(2):
                f.readline()
            lines = f.readlines()
            for i, ci in enumerate(args.cor_orb_index):

                x = [float(line.split()[0]) for line in lines]
                y_real = [float(line.split()[2 * ci - 1]) for line in lines]
                y_imag = [float(line.split()[2 * ci]) for line in lines]

                label_real = r"Re ({})".format(str(args.cor_orb_labels[i]))
                label_imag = r"Im ({})".format(str(args.cor_orb_labels[i]))

                ax1.plot(x, y_real, label=label_real)
                ax2.plot(x, y_imag, label=label_imag)

        ax1.set_title("Analytically continued self energy (Re)")
        ax1.set_xlabel("$i{\omega_n}$")
        ax1.set_ylabel("Re $\Sigma(i{\omega_n})$")
        ax1.legend()
        fig1.savefig("./plots/sigma-ac-real.pdf")

        ax2.set_title("Analytically continued self energy (Im)")
        ax2.set_xlabel("$i{\omega_n}$")
        ax2.set_ylabel("Im $\Sigma(i{\omega_n})$")
        ax2.legend()
        fig2.savefig("./plots/sigma-ac-imag.pdf")

        f.close()
    else:
        print("No analytic continuation detected. Skipping.")


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
            help="How many self energy files (sig.inp.x.x) to average?",
        )
        parser.add_argument(
            "-cor_orb_index",
            type=int,
            nargs="+",
            default=[1],
            help="List of cor_orb indexes (starting from 1).",
        )
        parser.add_argument(
            "-cor_orb_labels",
            type=str,
            nargs="+",
            default=["orbital1"],
            help="List of cor_orb names.",
        )

        args = parser.parse_args()
        plot_dmft(args)
    else:
        print("Usage: plotDMFT.py -h")
