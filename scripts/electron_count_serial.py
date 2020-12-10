#!/usr/bin/env python

""" Wannier Manifold Electron Occupation Calculater.

This is an extension to the original code developed by Xingyu Liao in Dr. Park's
group at UIC.
It calculates the total number of electrons in the Wannier manifold.
This is an initial requirement to perform DMFT calculations with DMFTwDFT.

This requires that the DMFTwDFT/bin directory is added to $PYTHONPATH.
It takes the atomnames, orbs, cor_at and cor_orb from INPUT.py.
Also requires the DFT input files to be present.

"""


from scipy import *
from os import path
from argparse import RawTextHelpFormatter
import argparse
import sys, subprocess, os
import re
from shutil import copy, copyfile
from itertools import groupby
import numpy as np

import Struct
import VASP
from INPUT import *


class ElectronOccupation:
    """This class has the methods to calculate the electrons in the
    correlated subspace."""

    def __init__(self, args):

        self.dft = args.dft
        self.np = args.np
        self.structurename = args.structurename
        self.nowin = args.nowin

        # Update the total number of bands after DFT calculation.
        # If you don't want this set flag nowin.
        self.updatewanbands = True
        self.wanbands = 0

        # import the VASP class. This can be used for other DFT codes as well.
        self.create_DFTmu()
        self.DFT = VASP.VASP_class(dft=self.dft, structurename=self.structurename)

        # Initial guess for Fermi energy
        self.DFT.EFERMI = 7.0

        # Wannier parameters
        self.kmeshtol = 1e-06
        self.dos_kmesh = 64
        self.num_iter = 1000

        # create bands_plot directory
        if not os.path.exists("bands_plot"):
            os.makedirs("bands_plot")

        # Starting the calculation
        if self.dft == "siesta":
            self.fdf_to_poscar()

        self.dft_run()
        self.postw90_run()
        self.calculator()

    def create_DFTmu(self):
        """
        This function creates a DFT_mu.out with an initial
        guess for the chemical potential. It will be updated once
        the DFT calculation is finished.
        """
        mu = 7.0

        if os.path.exists("DFT_mu.out"):
            os.remove("DFT_mu.out")
        f = open("DFT_mu.out", "w")
        f.write("%f" % mu)
        f.close()

    def fdf_to_poscar(self):
        """
	This function converts the siesta .fdf format to POSCAR for further calculations.
	"""
        fname = self.structurename + ".fdf"
        file = open(fname, "r")
        data = file.read()
        file.close()

        atoms = []
        lattice_constant = float(re.findall(r"LatticeConstant\s*([\d.]*)", data)[0])
        lattice_vectors = re.findall(r"LatticeVectors\s*([\d.\s]*)%endblock", data)

        # creating a numpy array with lattice vectors
        lattice_vec = np.array(lattice_vectors[0].split(), dtype="float64")
        lattice_vec = lattice_vec.reshape(3, 3)
        self.cell = lattice_constant * lattice_vec

        atomic_coordinates = re.findall(
            r"AtomicCoordinatesAndAtomicSpecies\s*([\d.\sa-zA-Z]*)%endblock", data
        )
        atomic_coordinates_lines = atomic_coordinates[0].split("\n")

        atm_coord_len = len(atomic_coordinates_lines)
        for i in range(atm_coord_len - 1):
            atoms.append(atomic_coordinates_lines[i].split()[-1])

        self.symbols = atoms

        species = [i[0] for i in groupby(atoms)]
        species_count = [len(list(group)) for key, group in groupby(atoms)]

        self.positions = np.zeros((atm_coord_len - 1, 3), dtype="float64")
        for counter in range(atm_coord_len - 1):
            self.positions[counter, :] = atomic_coordinates_lines[counter].split()[0:3]

        f = open("POSCAR", "w")
        f.write(" ".join(str(x) for x in species))
        f.write("\n%f\n" % lattice_constant)
        f.write("\n".join(str(x) for x in lattice_vectors[0].split("\n")))
        f.write(" ".join(str(x) for x in species))
        f.write("\n")
        f.write(" ".join(str(x) for x in species_count))
        f.write("\nDirect\n")
        for i in range(atm_coord_len - 1):
            f.write(
                " ".join(
                    [
                        " ".join(map(str, atomic_coordinates_lines[i].split()[0:3])),
                        atomic_coordinates_lines[i].split()[-1],
                        "\n",
                    ]
                )
            )
        f.close()

    def gen_win(self):
        """
        This method generates wannier90.win for the initial DFT run.
        """
        # generating wannier90.win
        TB = Struct.TBstructure("POSCAR", p["atomnames"], p["orbs"])
        TB.Compute_cor_idx(p["cor_at"], p["cor_orb"])

        # Read number of bands from DFT input file
        try:
            if self.dft == "vasp":
                fi = open("INCAR", "r")
                data = fi.read()
                fi.close()
                self.DFT.NBANDS = int(
                    re.findall(r"\n\s*NBANDS\s*=\s*([\d\s]*)", data)[0]
                )
                print("\nNumber of bands read from INCAR = %d " % self.DFT.NBANDS)

            elif self.dft == "siesta":
                fi = open(self.structurename + ".fdf", "r")
                data = fi.read()
                fi.close()
                self.DFT.NBANDS = int(
                    re.findall(r"\n\s*Siesta2Wannier90.NumberOfBands[\s0-9]*", data)[
                        0
                    ].split()[-1]
                )
                print("\nNumber of bands read from .fdf = %d " % self.DFT.NBANDS)

        except:
            self.DFT.NBANDS = 100
            print("WARNING: Number of bands not set in DFT input file!")

        # Setting num_bands in .win file.
        # If set to False num_bands is set to number of DFT bands.
        if list(p.keys()).count("num_bands_win"):
            if p["num_bands_win"]:
                self.wanbands = p["num_bands_win"]
                self.updatewanbands = False
            else:
                self.wanbands = self.DFT.NBANDS
        else:
            self.wanbands = self.DFT.NBANDS

        self.DFT.Create_win(
            TB,
            p["atomnames"],
            p["orbs"],
            p["L_rot"],
            self.wanbands,
            # Initially DFT.EFERMI is taken from DFT_mu.out but will
            # be updated later once the DFT calculation is complete.
            self.DFT.EFERMI + p["ewin"][0],
            self.DFT.EFERMI + p["ewin"][1],
            self.kmeshtol,
            self.num_iter,
        )

        # If exclude_bands are to be included in the .win file.
        # Appending to current .win file.
        if list(p.keys()).count("exclude_bands"):
            if p["exclude_bands"]:
                f = open("wannier90.win", "a")
                f.write("\nexclude_bands :\t")
                f.write(", ".join(str(x) for x in p["exclude_bands"]))
                f.write("\n")
                f.close()

            else:
                pass
        else:
            pass

        # VASP populates the .win file when running but Siesta
        # does not so we need to create a complete .win file for
        # Siesta runs.

        if self.dft == "siesta":

            # Update wannier90.win file.
            f = open("wannier90.win", "a")
            f.write("\nbegin unit_cell_cart\n")
            np.savetxt(f, self.cell)
            f.write("end unit_cell_cart\n\n")

            # writing the atoms cart block
            f.write("begin atoms_frac\n")
            aT = (np.array([self.symbols])).T
            b = self.positions
            atoms_cart = np.concatenate((aT, b), axis=1)
            np.savetxt(f, atoms_cart, fmt="%s")
            f.write("end atoms_frac\n\n")

            # writing the mp_grid line
            fi = open(self.structurename + ".fdf")
            data = fi.read()
            fi.close()
            grid = re.findall(
                r"%block kgrid_Monkhorst_Pack([\s0-9.]*)%endblock kgrid_Monkhorst_Pack",
                data,
            )[0].split()
            f.write("mp_grid= %s %s %s \n" % (grid[0], grid[5], grid[10]))

            # kpoints
            f.write("\nbegin kpoints\n")
            cmd = "kmesh.pl " + grid[0] + " " + grid[5] + " " + grid[10] + " wannier"
            out, err = subprocess.Popen(
                cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            ).communicate()
            f.write(out.decode("utf-8"))
            if err:
                print(err.decode("utf-8"))
            f.write("end kpoints")
            f.close()

        # Call read_num_wann() to store self.num_wann
        self.read_num_wann()
        print("wannier90.win generated.\n")

    def read_num_wann(self):
        """This reads the number of wannier bands from the generatied wannier.win
        file."""

        fi = open("wannier90.win", "r")
        data = fi.readlines()
        fi.close()

        for line in data:
            if re.match("num_wann", line):
                self.num_wann = line.split()[-1]
        print("Number of Wannier functions = %s" % self.num_wann)

    def update_win(self):
        """
        This updates the wannier90.win file with the number of bands and fermi energy
        from the initial DFT calculation.
	"""
        # Updating wannier90.win with the number of DFT bands
        if self.updatewanbands:
            self.DFT.Read_NBANDS()
            self.DFT.Read_EFERMI()
            self.DFT.Update_win(
                self.DFT.NBANDS,
                self.DFT.EFERMI + p["ewin"][0],
                self.DFT.EFERMI + p["ewin"][1],
            )
        else:
            self.DFT.Read_NBANDS()
            self.DFT.Read_EFERMI()
            self.DFT.Update_win(
                self.wanbands,
                self.DFT.EFERMI + p["ewin"][0],
                self.DFT.EFERMI + p["ewin"][1],
            )

        # Update wannier90.win file with the dos block.
        f = open("wannier90.win", "a")
        f.write("\ndos=true")
        f.write("\ndos_kmesh=%s" % self.dos_kmesh)
        f.write("\ndos_project=%s" % self.num_wann)
        f.close()

        print("wannier90.win updated.")

        # find line number of dis_project in wannier90.win.
        fi = open("wannier90.win", "r")
        data = fi.readlines()
        fi.close()

        for ln, line in enumerate(data):
            if re.match("dos_project", line):
                self.dos_project_line = ln + 1

    def dft_run(self):
        """ DFT runner.
        This method performs the initial DFT calculation and the wannier90 calculation..
        """

        # VASP
        if self.dft == "vasp":
            if not self.nowin:
                self.gen_win()  # Generating .win file
            else:
                self.updatewanbands = False

            # initial VASP run
            print("Running VASP...")
            self.vasp_exec = "vasp_std"
            cmd = (
                "mpirun -np " + str(self.np) + " " + self.vasp_exec
            )  # + " > dft.out 2> dft.error"
            out, err = subprocess.Popen(
                cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            ).communicate()
            if err:
                print("DFT calculation failed! Check dft.error for details.\n")
                errdir = "dft.error"
                f = open(errdir, "wb")
                f.write(err)
                f.close()
                sys.exit()
            else:
                print("DFT calculation complete.\n")
                outdir = "dft.out"
                f = open(outdir, "wb")
                f.write(out)
                f.close()

            # Updating .win file and running wannier90.x.
            self.update_win()
            self.run_wan90()

        # Siesta
        elif self.dft == "siesta":
            if not self.nowin:
                self.gen_win()
                copy("wannier90.win", self.structurename + ".win")
            else:
                self.updatewanbands = False

            # Running Wannier90 pre-processor
            self.run_wan90_pp()

            # Running siesta
            print("Running Siesta ...")
            self.siesta_exec = "siesta"
            cmd = (
                "mpirun -np "
                + str(self.np)
                + " "
                + self.siesta_exec
                + "<"
                + self.structurename
                + ".fdf>"
                + self.structurename
                + ".out"
            )
            out, err = subprocess.Popen(
                cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            ).communicate()

            if os.path.exists(self.structurename + ".out"):

                fi = open(self.structurename + ".out", "r")
                done_word = fi.readlines()[-1]
                fi.close()

                if done_word.split()[0] == "Job":
                    print("DFT calculation complete.\n")

                else:
                    print("DFT calculation failed!\n")
                    sys.exit()

            else:
                print("DFT calculation failed!\n")
                sys.exit()

            # need to rename .eigW to .eig to run wannier90
            copy(self.structurename + ".eigW", self.structurename + ".eig")

            if not self.nowin:
                self.update_win()
                copy("wannier90.win", self.structurename + ".win")
            self.run_wan90(self.structurename)

            # renaming files
            copy(self.structurename + ".eig", "wannier90.eig")
            copy(self.structurename + ".chk", "wannier90.chk")
            copy(self.structurename + ".win", "wannier90.win")
            copy(self.structurename + ".amn", "wannier90.amn")

    def run_wan90_pp(self):
        """
	This function performs the wannier90 pre-processing required by some DFT codes like siesta.
        Outputs a .nnkp file which is required for the DFT calculaiton.
        """
        cmd = "wannier90.x -pp" + " " + self.structurename
        out, err = subprocess.Popen(
            cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        ).communicate()
        if err:
            print(err.decode("utf-8"))
            sys.exit()
        else:
            print(out.decode("utf-8"))

    def run_wan90(self, filename="wannier90"):
        """
        Running wannier90.x to generate .chk file.
        """

        print("\nRunning wannier90.x ...")
        cmd = "wannier90.x" + " " + filename
        out, err = subprocess.Popen(
            cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        ).communicate()
        if err:
            print("wannier90 calculation failed!")
            print(err.decode("utf-8"))
            sys.exit()
        else:
            print("wannier90 calculation complete.")
            print(out.decode("utf-8"))

    def postw90_run(self):
        NUM_OF_ORBT = int(self.num_wann)  ## NUMBER OF WANNIER BAND
        LINE_dos = self.dos_project_line  ### LINENUMBER WHICH DEFINE THE dos_project
        for i in range(NUM_OF_ORBT):
            print("%s th iteration is running" % (i + 1))
            fi = open("wannier90.win", "r")
            WIN = array(fi.readlines())
            fi.close()
            WIN[LINE_dos - 1] = "dos_project=" + str(i + 1) + "\n"
            fi = open("wannier90.win", "w")
            for j, winline in enumerate(WIN):
                fi.write(winline)
            fi.close()
            print("postw90 is running")
            cmd = "postw90.x wannier90"
            out, err = subprocess.Popen(
                cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            ).communicate()
            os.system("cp wannier90-dos.dat bands_plot/dos" + str(i + 1))
            print("-------------------------------------------------")

    def Load_DOS(self, filename, Fermi):
        fi = open(filename, "r")
        dos = fi.readlines()
        xy = zeros((len(dos), 2))
        for i, den in enumerate(dos):
            xy[i] = array(array(den.split()).astype(float))
            xy[i][0] -= Fermi
        fi.close()
        return xy

    def Integration(self, xy):
        total = 0.0
        for i in range(len(xy) - 1):
            total = (
                total + abs(xy[i + 1][0] - xy[i][0]) * (xy[i][1] + xy[i + 1][1]) / 2.0
            )
        return total

    def calculator(self):
        occ = 0.0
        for i in range(int(self.num_wann)):
            if not path.exists("./bands_plot/dos" + str(i + 1)):
                print("File does not exist!")
                exit()
            sepE = 0
            speI = 0
            dos = self.Load_DOS("./bands_plot/dos" + str(i + 1), self.DFT.EFERMI)
            for j, ele in enumerate(dos):
                if ele[0] > sepE:
                    speI = j
                    break
            print("%s : %s " % ((i + 1), self.Integration(dos[:speI, :])))

            occ = self.Integration(dos[:speI, :]) + occ
        print("Total electron occupancy in Wannier manifold : %s" % occ)
        np.savetxt("occupancy.out", [occ], fmt="%10.5f")


if __name__ == "__main__":
    args = sys.argv[1:]
    if args:

        des = (
            "This script counts the total number of electrons in the Wannier sub space."
        )
        parser = argparse.ArgumentParser(
            description=des, formatter_class=RawTextHelpFormatter
        )

        parser.add_argument(
            "-dft",
            default="vasp",
            type=str,
            help="Choice of DFT code.",
            choices=["vasp", "siesta", "qe"],
        )

        parser.add_argument(
            "-structurename",
            type=str,
            help="Name of the structure. Not required for VASP. ",
            default=None,
        )

        parser.add_argument("-np", default=1, type=int, help="Number of processors.")

        parser.add_argument(
            "-nowin",
            action="store_true",
            help="Flag to disable automatic generation of .win file.",
        )

        args = parser.parse_args()
        ElectronOccupation(args)

    else:
        print("Usage: electron_count.py -h")
