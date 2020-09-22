#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import argparse
import os
import re
import shutil
import subprocess
import sys
from argparse import RawTextHelpFormatter
from itertools import groupby
from shutil import copyfile

import numpy as np

import Struct
import VASP
from INPUT import *
from splash import welcome


class DMFTLauncher:
    """ DMFTwDFT Initialization and Calculation (Python 2.x version).

	This class contains methods to run the initial DFT and wannier90 calculations
	to generate inputs for the DMFT calculation and performs it.

        This version does not support ionic covergence. Please use an optimized structure for the calculation.

	Run with:
	DMFT.py <options>

	<options>:
	-dft <vasp,siesta,qe>
	-dmft : This flag performs dmft calculation
	-hf : This flag performs Hartree-Fock calcualtion
	-force : This flag forces a dmft or hf calculation even if it has been completed
	-kmeshtol : k-mesh tolerance for wannier90
        -aiida : Flag for aiida calculations
        -nowin : Flag to disable automatic generation of .win file
        -v : Enable verbosity

	"""

    def __init__(self, args):
        """
	Contains common functions for all methods.
	This launches the dmft calculation as well.
	"""

        print("Initializing calculation...\n")

        if os.path.exists("para_com.dat"):
            fipa = open("para_com.dat", "r")
            self.para_com = str(fipa.readline())[:-1]
            fipa.close()
        else:
            self.para_com = "mpirun -np 40"

        if os.path.exists("para_com_dft.dat"):
            fid = open("para_com_dft.dat")
            self.para_com_dft = str(fid.readline())[:-1]
            fid.close()
        else:
            self.para_com_dft = self.para_com

        # initial chemical potential
        self.create_DFTmu()

        # import the VASP class. This can be used for other DFT codes as well.
        self.DFT = VASP.VASP_class(
            dft=args.dft, aiida=args.aiida, structurename=args.structurename
        )

        self.dir = os.getcwd()  # dft running directory (current directory)
        self.structurename = (
            args.structurename
        )  # name of structure. Required for siesta (structurename.fdf).
        self.kmeshtol = args.kmeshtol  # kmesh tolerence for wannier mesh
        self.force = args.force  # force dmft calculation True of False
        self.v = args.v  # Verbosity
        self.dft = args.dft  # DFT type
        self.aiida = args.aiida  # Flag for aiida calculation
        self.lowdin = args.lowdin  # Flag for Siesta Lowdin
        self.nowin = args.nowin  # Flag for .win generation

        ####### DFT and wannier90 initialization #######

        # wannier90 executable
        self.wannier90_exec = "wannier90.x"
        self.wanbands = 0
        self.updatewanbands = True

        # VASP calculation
        if self.dft == "vasp":
            self.vasp_exec = "vasp_std"  # vasp executable

        # Siesta calculation
        elif self.dft == "siesta":
            self.siesta_exec = "siesta"  # siesta executable
            self.fdf_to_poscar()

        # QE calculation
        elif self.dft == "qe" and self.aiida is False:
            # Nothing to do here for now since we only have
            # QE aiida calculations.
            self.qe_exec = "pw.x"
            self.qe_to_poscar()

        # aiida calculation
        if self.aiida:
            self.win_to_poscar(input="aiida.win")

        ####### DMFT initialization #######

        # Generate initial self energy.
        self.gen_sig()

        # Setting  DMFT type
        if args.dmft:
            self.type = "DMFT"

        elif args.hf:
            self.type = "HF"

        # Launch DMFT calculation
        self.run_dmft()

    # ------------------------------- UTILITIES -------------------------------------------------

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

    def gen_sig(self):
        """
        This method generates the initial self energy file sig.inp.
        """
        cmd = "sigzero.py"
        out, err = subprocess.Popen(
            cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        ).communicate()
        if err:
            print(err.decode("utf-8"))
            sys.exit()
        else:
            print("Initial self-energy file generated.")

    def fdf_to_poscar(self):
        """
	This function converts the siesta .fdf format to POSCAR for further calculations.
	"""
        # file = pychemia.code.siesta.SiestaInput(self.structurename + ".fdf")
        # self.st = file.get_structure()
        # pychemia.code.vasp.write_poscar(
        #     self.st, filepath="POSCAR", newformat=True, direct=True, comment=None
        # )
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

    def qe_to_poscar(self):
        """Creates a POSCAR from a Quantum Espressso scf input file.
           CELL_PARAMETERS must be defined.
        """
        fname = self.structurename + ".scf.in"
        file = open(fname, "r")
        data = file.read()
        file.close()
        lattice_vec = re.findall(r"CELL_PARAMETERS\s*[a-zA-Z]*([e\d\s.+-]*)", data)
        lattice_vec = [float(x) for x in lattice_vec[0].split()]
        self.cell = np.array((lattice_vec), dtype="float64").reshape(3, 3)
        nat = int(re.findall(r"nat\s*=\s*([\d]*)", data)[0])
        raw_positions = re.findall(
            r"ATOMIC_POSITIONS\s*crystal([e\d.+\sa-zA-Z]*)\n", data
        )
        full_structure = np.array(raw_positions[0].split()).reshape(nat, 4)
        self.positions = np.zeros((nat, 3), dtype="float64")
        self.symbols = []
        for icount, i in enumerate(full_structure):
            self.positions[icount] = i[1:]
            self.symbols.append(i[0])
        species = [i[0] for i in groupby(self.symbols)]
        species_count = [len(list(group)) for key, group in groupby(self.symbols)]

        # Writing to POSCAR
        f = open("POSCAR", "w")
        f.write(" ".join(str(x) for x in species))
        f.write("\n%f\n" % 1.0)
        for i in range(len(self.cell)):
            f.write("%f %f %f\n" % (self.cell[i, 0], self.cell[i, 1], self.cell[i, 2]))
        f.write(" ".join(str(x) for x in species))
        f.write("\n")
        f.write(" ".join(str(x) for x in species_count))
        f.write("\n")
        f.write("Direct\n")

        for i in range(len(full_structure)):
            f.write(
                "%s %s %s %s\n"
                % (
                    full_structure[i, 1],
                    full_structure[i, 2],
                    full_structure[i, 3],
                    full_structure[i, 0],
                )
            )
        f.close()

    def read_poscar(self, fname="POSCAR"):
        "Reads a POSCAR and sets self.cell, self.symbols and self.positions."

        file = open(fname, "r")
        data = file.readlines()
        file.close()

        lattice_constant = float(data[1])
        lattice_vec = np.array(
            (
                [float(x) for x in data[2].split()],
                [float(x) for x in data[3].split()],
                [float(x) for x in data[4].split()],
            )
        )
        self.cell = lattice_constant * lattice_vec
        num_atoms = np.sum([int(x) for x in data[6].split()])
        full_structure = np.array((data[8 : 8 + num_atoms]), dtype="str")
        self.positions = np.zeros((num_atoms, 3), dtype="float64")
        self.symbols = []
        for icount, i in enumerate(full_structure):
            self.positions[icount] = i.split()[0:3]
            self.symbols.append(i.split()[-1])

    def win_to_poscar(self, input="aiida.win"):
        """This function generates a POSCAR from a .win file.
        Used in aiida calculations where the .win file is provided."""

        file = open(input, "r")
        data = file.read()
        file.close()

        # Unit cell
        unit_cell_data = re.findall(
            r"begin\s*unit\s*_cell_cart\s*[a-zA-Z]*([+-e\s\d.]*)end", data
        )
        unit_cell = np.array(unit_cell_data[0].split(), dtype="float64")
        unit_cell = unit_cell.reshape(3, 3)

        # coordinate type (cartesian or direct)
        coord_type = re.findall(r"begin\s*atoms_([a-z]*)", data)[0]

        # atomic coordinates
        atom_coords = np.array(
            re.findall(r"begin\s*atoms_[\sa-z]*([\sa-zA-Z\d.]*)end", data)[0].split()
        )
        atom_coords = atom_coords.reshape(int(len(atom_coords) / 4), 4)

        atoms = atom_coords[:, 0]
        species = [i[0] for i in groupby(atoms)]
        species_count = [len(list(group)) for key, group in groupby(atoms)]

        f = open("POSCAR", "w")
        f.write(" ".join(str(x) for x in species))
        f.write("\n%f\n" % 1.0)
        for i in range(len(unit_cell)):
            f.write("%f %f %f\n" % (unit_cell[i, 0], unit_cell[i, 1], unit_cell[i, 2]))
        f.write(" ".join(str(x) for x in species))
        f.write("\n")
        f.write(" ".join(str(x) for x in species_count))
        f.write("\n")
        if coord_type == "cart":
            f.write("Cartesian\n")
        elif coord_type == "frac":
            f.write("Direct\n")

        for i in range(len(atom_coords)):
            f.write(
                "%s %s %s %s\n"
                % (
                    atom_coords[i, 1],
                    atom_coords[i, 2],
                    atom_coords[i, 3],
                    atom_coords[i, 0],
                )
            )

        f.close()

    def copy_files(self):
        """
        This creates a directory DMFT or HF in the root directory
        and copies all the necessary files.
        """

        # creating directory for DMFT
        if os.path.exists(self.type):
            if os.path.exists(self.type + "/imp.0/"):
                shutil.rmtree(self.type + "/imp.0/")
                # os.makedirs("DMFT")
        else:
            os.makedirs(self.type)

        # copying files into DMFT or HF directory
        if self.structurename != None and self.dft != None:
            cmd = (
                "cd "
                + self.type
                + " && Copy_input.py ../ "
                + "-structurename "
                + self.structurename
                + " -dft "
                + self.dft
            )
        else:
            cmd = "cd " + self.type + " && Copy_input.py ../ "
        out, err = subprocess.Popen(cmd, shell=True).communicate()
        if err:
            print("File copy failed!\n")
            print(err)
            sys.exit()
        else:
            print(out)
            print(
                "\n"
                + self.type
                + " initialization complete. Ready to run calculation.\n"
            )

    def gen_pw2wannier90(self):
        # Generates pw2wannier90 file for QE+wannier90 calculation.
        fi = open(self.structurename + ".pw2wannier90.in", "wt")
        fi.write("&inputpp\n")
        fi.write("\toutdir = './'\n")
        fi.write("\tprefix = '%s'\n" % self.structurename)
        fi.write("\tseedname = '%s'\n" % self.structurename)
        fi.write("\tspin_component = 'none'\n")
        fi.write("\twrite_mmn = .true.\n")
        fi.write("\twrite_amn = .true.\n")
        fi.write("\twrite_unk = .false.\n")
        fi.write("\twrite_dmn = .false.\n")
        fi.write("\twan_mode = 'standalone'\n")
        fi.write("/\n")
        fi.close()

    # ----------------------------- WANNIER90 -------------------------------------------

    def gen_win(self):
        """
        This method generates wannier90.win for the initial DFT run.
        """

        # generating wannier90.win
        TB = Struct.TBstructure("POSCAR", p["atomnames"], p["orbs"])
        TB.Compute_cor_idx(p["cor_at"], p["cor_orb"])
        # print((TB.TB_orbs))

        # Read number of bands from DFT input file
        try:
            if self.dft == "vasp":
                fi = open("INCAR", "r")
                data = fi.read()
                fi.close()
                self.DFT.NBANDS = int(
                    re.findall(r"\n\s*NBANDS\s*=\s*([\d\s]*)", data)[0]
                )
                print("Number of bands read from INCAR = %d " % self.DFT.NBANDS)

            elif self.dft == "siesta":
                fi = open(self.structurename + ".fdf", "r")
                data = fi.read()
                fi.close()
                self.DFT.NBANDS = int(
                    re.findall(r"\n\s*Siesta2Wannier90.NumberOfBands[\s0-9]*", data)[
                        0
                    ].split()[-1]
                )
                print("Number of bands read from .fdf = %d " % self.DFT.NBANDS)

            elif self.dft == "qe":
                fi = open(self.structurename + ".scf.out", "r")
                data = fi.read()
                fi.close()
                self.DFT.NBANDS = int(
                    re.findall(
                        r"\n\s*number\s*of\s*Kohn-Sham\s*states=([\s\d]*)", data
                    )[0]
                )
                print("Number of bands read from .scf.out = %d " % self.DFT.NBANDS)

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
            print("wannier90.win generated.")

        # Similar to Siesta, Quantum Espresso requires the generation
        # of the .win file manually.
        # self.cell, self.symbols and self.positions is set from
        # read_poscar().

        elif self.dft == "qe" and self.aiida == False:

            # reading poscar that should be generated.
            self.read_poscar()

            # Update wannier90.win file then rename it
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
            fi = open(self.structurename + ".scf.in")
            data = fi.read()
            fi.close()
            self.grid = re.findall(r"K_POINTS\s*automatic\s*([\d\s]*)", data,)[
                0
            ].split()[0:3]
            self.grid = [int(x) for x in self.grid]
            f.write("mp_grid= %s %s %s \n" % (self.grid[0], self.grid[1], self.grid[2]))

            # kpoints
            f.write("\nbegin kpoints\n")
            cmd = (
                "kmesh.pl "
                + str(self.grid[0])
                + " "
                + str(self.grid[1])
                + " "
                + str(self.grid[2])
                + " wannier"
            )
            out, err = subprocess.Popen(
                cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            ).communicate()
            f.write(out.decode("utf-8"))
            if err:
                print(err.decode("utf-8"))
            f.write("end kpoints")
            f.close()
            print("wannier90.win generated.")

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

        print("wannier90.win updated.")

        # Updating DFT_mu.out
        np.savetxt("DFT_mu.out", [self.DFT.EFERMI])

    def run_wan90_pp(self):
        """
	This function performs the wannier90 pre-processing required by some DFT codes like siesta.
        Outputs a .nnkp file which is required for the DFT calculaiton.
        """
        cmd = self.wannier90_exec + " -pp" + " " + self.structurename
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
        Running wannier90.x to generate .chk and .eig files.
        """

        print("Running wannier90...")
        cmd = self.para_com + " " + self.wannier90_exec + " " + filename
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

    # -------------------------------- DFT ----------------------------------------------

    def run_dft(self):
        """
        This function  calls the dft calculations and the wannier calculations
        """

        # VASP
        if self.dft == "vasp":
            if not self.nowin:
                self.gen_win()  # Generating .win file
            else:
                self.updatewanbands = False
            self.vasp_run(self.dir)
            self.update_win()
            self.run_wan90()
            self.copy_files()

        # Siesta
        elif self.dft == "siesta":
            if not self.lowdin:  # Generate .win file. Lowdin does this internally.
                if not self.nowin:
                    self.gen_win()
                    shutil.copy("wannier90.win", self.structurename + ".win")
                else:
                    self.updatewanbands = False

            # Running Siesta
            self.siesta_run(self.dir)

            # need to rename .eigW to .eig to run wannier90
            shutil.copy(self.structurename + ".eigW", self.structurename + ".eig")

            if not self.lowdin:
                if not self.nowin:
                    self.update_win()
                    shutil.copy("wannier90.win", self.structurename + ".win")
                self.run_wan90(self.structurename)

            # renaming files
            shutil.copy(self.structurename + ".eig", "wannier90.eig")
            shutil.copy(self.structurename + ".chk", "wannier90.chk")
            shutil.copy(self.structurename + ".win", "wannier90.win")
            shutil.copy(self.structurename + ".amn", "wannier90.amn")
            self.copy_files()

        # Quantum Espresso (Without aiida)
        elif self.dft == "qe" and self.aiida == False:
            if not self.nowin:
                self.gen_win()  # Assume POSCAR is present.
                shutil.copy("wannier90.win", self.structurename + ".win")
            else:
                self.updatewanbands = False

            # Running Quantum Espresso for SCF
            self.qe_run(cal_type="scf", dir=self.dir)

            # Setting up NSCF run.
            shutil.copy(self.structurename + ".scf.in", self.structurename + ".nscf.in")

            # replacing calculation = 'scf' with calculation = 'nscf' in the .nscf.in file.
            fi = open(self.structurename + ".nscf.in", "r")
            data = fi.read()
            fi.close()
            data = data.replace("'scf'", "'nscf'")
            fi = open(self.structurename + ".nscf.in", "w")
            fi.write(data)
            fi.close()

            # Removing the last two lines:
            # K_POINTS automatic
            # <grid>
            fi = open(self.structurename + ".nscf.in", "r")
            data = fi.readlines()
            fi.close()

            fi = open(self.structurename + ".nscf.in", "w")
            for i in data[0:-2]:
                fi.write(i)
            fi.close()

            # Appending the K_POINTS crystal grid using
            # kmesh.pl.
            f = open(self.structurename + ".nscf.in", "a")
            cmd = (
                "kmesh.pl "
                + str(self.grid[0])
                + " "
                + str(self.grid[1])
                + " "
                + str(self.grid[2])
            )
            out, err = subprocess.Popen(
                cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            ).communicate()
            f.write(out.decode("utf-8"))
            if err:
                print(err.decode("utf-8"))
            f.close()

            # Save the Total energy from the SCF calculation.
            # Number of bands and fermi energy are updated when
            # update_win() is called.
            self.DFT.Read_OSZICAR()

            # Updating .win file.
            if not self.nowin:
                self.update_win()
                shutil.copy("wannier90.win", self.structurename + ".win")

            # Running Quantum Espresso for NSCF calculation.
            self.qe_run(cal_type="nscf", dir=self.dir)

            # Run wannier90 pre-processing.
            self.run_wan90_pp()

            # Run Quantum Espresso - wannier90 interface.
            self.gen_pw2wannier90()
            self.qe_run(cal_type="pw2wannier90", dir=self.dir)

            # Running wannier90.x
            self.run_wan90(self.structurename)

            # renaming files
            shutil.copy(self.structurename + ".eig", "wannier90.eig")
            shutil.copy(self.structurename + ".chk", "wannier90.chk")
            shutil.copy(self.structurename + ".win", "wannier90.win")
            shutil.copy(self.structurename + ".amn", "wannier90.amn")
            self.copy_files()

        # aiida
        elif self.aiida:
            # renaming files
            shutil.copy("aiida.eig", "wannier90.eig")
            shutil.copy("aiida.chk", "wannier90.chk")
            shutil.copy("aiida.win", "wannier90.win")
            shutil.copy("aiida.amn", "wannier90.amn")
            self.copy_files()

    def vasp_run(self, dir):
        """
	This method runs the inital VASP calculation.
	"""

        # initial VASP run
        print("Running VASP in %s" % dir)
        cmd = (
            "cd " + dir + " && " + self.para_com_dft + " " + self.vasp_exec
        )  # + " > dft.out 2> dft.error"
        out, err = subprocess.Popen(
            cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        ).communicate()
        if err:
            print("DFT calculation failed! Check dft.error for details.\n")
            errdir = dir + os.sep + "dft.error"
            f = open(errdir, "wb")
            f.write(err)
            f.close()
            sys.exit()
        else:
            print("DFT calculation complete.\n")
            outdir = dir + os.sep + "dft.out"
            f = open(outdir, "wb")
            f.write(out)
            f.close()

    def siesta_run(self, dir):
        """
	This method runs the initial siesta calculation.
        """
        # wannier90 pre-processing
        if not self.lowdin:
            self.run_wan90_pp()

        # Running siesta
        print("Running Siesta in %s" % dir)
        cmd = (
            "cd "
            + dir
            + " && "
            + self.para_com_dft
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

    def qe_run(self, cal_type="scf", dir=None):
        """
        This runs the Quantum Espresso SCF/NSCF calculation with pw.x.

        """

        def run_pwx(finname=None, foutname=None, cal_type="scf"):

            print("\nRunning Quantum Espresso %s in %s" % (cal_type, dir))
            cmd = (
                "cd "
                + dir
                + " && "
                + self.para_com_dft
                + " "
                + self.qe_exec
                + " -in "
                + finname
                + " > "
                + foutname
            )
            out, err = subprocess.Popen(
                cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            ).communicate()

            if os.path.exists(foutname):

                fi = open(foutname, "r")
                done_word = re.search(r"JOB\s*DONE.", fi.read())
                fi.close()

                if done_word:
                    print("Quantum Espresso %s  calculation complete.\n" % cal_type)

                else:
                    print("Quantum Espresso %s calculation failed!\n" % cal_type)
                    sys.exit()

            else:
                print("Quantum Espresso %s calculation failed!\n" % cal_type)
                sys.exit()

        def run_pw2wannier90(finname=None, foutname=None):

            print("Running Quantum Espresso pw2wannier90.x ...")
            cmd = (
                "cd "
                + dir
                + " && "
                + self.para_com_dft
                + " "
                + "pw2wannier90.x"
                + " -in "
                + finname
                + " > "
                + foutname
            )
            out, err = subprocess.Popen(
                cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            ).communicate()

            if os.path.exists(foutname):

                fi = open(foutname, "r")
                done_word = re.search(r"JOB\s*DONE.", fi.read())
                fi.close()

                if done_word:
                    print("Quantum Espresso pw2wannier90.x calculation complete.\n")

                else:
                    print("Quantum Espresso pw2wannier90.x calculation failed!\n")
                    sys.exit()

            else:
                print("Quantum Espresso pw2wannier90.x calculation failed!\n")
                sys.exit()

        # setting input and output file names
        if cal_type == "scf":
            finname = self.structurename + ".scf.in"
            foutname = self.structurename + ".scf.out"

        elif cal_type == "nscf":
            finname = self.structurename + ".nscf.in"
            foutname = self.structurename + ".nscf.out"

        elif cal_type == "pw2wannier90":
            finname = self.structurename + ".pw2wannier90.in"
            foutname = self.structurename + ".pw2wannier90.out"

        # Setting program to run
        if cal_type == "scf" or cal_type == "nscf":
            # Running pw.x
            run_pwx(finname, foutname, cal_type)

        elif cal_type == "pw2wannier90":
            # Running pw2wannier90.x
            run_pw2wannier90(finname, foutname)

    # --------------------- DMFT -----------------------------------------------

    def run_dmft(self):
        """
        This first checks if there is a previous DMFT or HF calculation and runs
        only if that run is incomplete unless forced.
        """

        separator_art = """
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        """

        # Checking for previous DMFT run in the directory
        pathstr = self.type + os.sep + "INFO_TIME"

        if os.path.exists(pathstr):
            fi = open(pathstr, "r")
            done_word = fi.readlines()[-1]
            fi.close()

            if done_word.split()[0] == "Calculation":
                print("Existing " + self.type + " calculation is complete.")
                if self.force:
                    # forcing DMFT calculation
                    print("-force flag enabled. Restarting " + self.type + "...")
                    self.run_dft()
                    print(separator_art)
                    print(
                        "*-*-*-*-*- Starting "
                        + self.type
                        + " calculation -*-*-*-*-* \n"
                    )

                    if self.dft != None and self.structurename != None:
                        if self.type == "HF":
                            cmd = (
                                "cd "
                                + self.type
                                + " && "
                                + "RUNDMFT.py -hf -dft "
                                + self.dft
                                + " -structurename "
                                + self.structurename
                                + " -aiida "
                                + str(self.aiida)
                            )
                        else:
                            cmd = (
                                "cd "
                                + self.type
                                + " && "
                                + "RUNDMFT.py -dft "
                                + self.dft
                                + " -structurename "
                                + self.structurename
                                + " -aiida "
                                + str(self.aiida)
                            )

                    elif self.dft != None:
                        if self.type == "HF":
                            cmd = (
                                "cd "
                                + self.type
                                + " && "
                                + "RUNDMFT.py -hf -dft "
                                + self.dft
                                + " -aiida "
                                + str(self.aiida)
                            )
                        else:
                            cmd = (
                                "cd "
                                + self.type
                                + " && "
                                + "RUNDMFT.py -dft "
                                + self.dft
                                + " -aiida "
                                + str(self.aiida)
                            )

                    else:
                        if self.type == "HF":
                            cmd = "cd " + self.type + " && " + "RUNDMFT.py -hf"
                        else:
                            cmd = "cd " + self.type + " && " + "RUNDMFT.py "

                    if self.v:
                        subprocess.Popen(cmd, shell=True,).communicate()

                    else:
                        out, err = subprocess.Popen(
                            cmd,
                            shell=True,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                        ).communicate()

                        if err:
                            print(
                                self.type
                                + " calculation failed! Check "
                                + self.type
                                + ".error for details.\n"
                            )
                            errdir = self.type + os.sep + self.type + ".error"
                            f = open(errdir, "wb")
                            f.write(err)
                            f.close()
                            sys.exit()
                        else:
                            print(self.type + " calculation complete.\n")
                            outdir = self.type + os.sep + self.type + ".out"
                            f = open(outdir, "wb")
                            f.write(out)
                            f.close()

                else:
                    # exit when exiting DMFT calculation is complete.
                    print("-force flag disabled. Exiting. ")
                    sys.exit()

            else:
                # Incomplete DMFT calculation.
                print("Incomplete " + self.type + " calculation found.")
                self.run_dft()
                print(separator_art)
                print("*-*-*-*-*- Starting " + self.type + " calculation -*-*-*-*-* \n")
                if self.dft != None and self.structurename != None:
                    if self.type == "HF":
                        cmd = (
                            "cd "
                            + self.type
                            + " && "
                            + "RUNDMFT.py -hf -dft "
                            + self.dft
                            + " -structurename "
                            + self.structurename
                        )
                    else:
                        cmd = (
                            "cd "
                            + self.type
                            + " && "
                            + "RUNDMFT.py -dft "
                            + self.dft
                            + " -structurename "
                            + self.structurename
                        )
                elif self.dft != None:
                    if self.type == "HF":
                        cmd = (
                            "cd "
                            + self.type
                            + " && "
                            + "RUNDMFT.py -hf -dft "
                            + self.dft
                            + " -aiida "
                            + str(self.aiida)
                        )
                    else:
                        cmd = (
                            "cd "
                            + self.type
                            + " && "
                            + "RUNDMFT.py -dft "
                            + self.dft
                            + " -aiida "
                            + str(self.aiida)
                        )

                else:
                    if self.type == "HF":
                        cmd = "cd " + self.type + " && " + "RUNDMFT.py -hf"
                    else:
                        cmd = "cd " + self.type + " && " + "RUNDMFT.py "

                if self.v:
                    subprocess.Popen(cmd, shell=True,).communicate()

                else:
                    out, err = subprocess.Popen(
                        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                    ).communicate()

                    if err:
                        print(
                            self.type
                            + " calculation failed! Check "
                            + self.type
                            + ".error for details.\n"
                        )
                        errdir = self.type + os.sep + self.type + ".error"
                        f = open(errdir, "wb")
                        f.write(err)
                        f.close()
                        sys.exit()
                    else:
                        print(self.type + " calculation complete.\n")
                        outdir = self.type + os.sep + self.type + ".out"
                        f = open(outdir, "wb")
                        f.write(out)
                        f.close()
        else:
            # no DMFT/INFO_TIME found
            self.run_dft()
            print(separator_art)
            print("*-*-*-*-*- Starting " + self.type + " calculation -*-*-*-*-* \n")
            if self.dft != None and self.structurename != None:
                if self.type == "HF":
                    cmd = (
                        "cd "
                        + self.type
                        + " && "
                        + "RUNDMFT.py -hf -dft "
                        + self.dft
                        + " -structurename "
                        + self.structurename
                    )
                else:
                    cmd = (
                        "cd "
                        + self.type
                        + " && "
                        + "RUNDMFT.py -dft "
                        + self.dft
                        + " -structurename "
                        + self.structurename
                    )

            elif self.dft != None:
                if self.type == "HF":
                    cmd = (
                        "cd "
                        + self.type
                        + " && "
                        + "RUNDMFT.py -hf -dft "
                        + self.dft
                        + " -aiida "
                        + str(self.aiida)
                    )
                else:
                    cmd = (
                        "cd "
                        + self.type
                        + " && "
                        + "RUNDMFT.py -dft "
                        + self.dft
                        + " -aiida "
                        + str(self.aiida)
                    )

            else:
                if self.type == "HF":
                    cmd = "cd " + self.type + " && " + "RUNDMFT.py -hf"
                else:
                    cmd = "cd " + self.type + " && " + "RUNDMFT.py "

            if self.v:
                subprocess.Popen(cmd, shell=True,).communicate()

            else:
                out, err = subprocess.Popen(
                    cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                ).communicate()

                if err:
                    print(
                        self.type
                        + " calculation failed! Check "
                        + self.type
                        + ".error for details.\n"
                    )
                    errdir = self.type + os.sep + self.type + ".error"
                    f = open(errdir, "wb")
                    f.write(err)
                    f.close()
                    sys.exit()
                else:
                    print(self.type + " calculation complete.\n")
                    outdir = self.type + os.sep + self.type + ".out"
                    f = open(outdir, "wb")
                    f.write(out)
                    f.close()


if __name__ == "__main__":
    args = sys.argv[1:]
    if args:
        welcome()
        des = "This script performs DFT+DMFT calculations through maximally localized Wannier functions.\nFor post-processing, run postDMFT.py."
        parser = argparse.ArgumentParser(
            description=des, formatter_class=RawTextHelpFormatter
        )

        # parser for dft
        parser.add_argument(
            "-dft",
            default="vasp",
            type=str,
            help="Choice of DFT code for the DMFT calculation.",
            choices=["vasp", "siesta", "qe"],
        )
        type_parser = parser.add_mutually_exclusive_group()
        type_parser.add_argument(
            "-dmft",
            action="store_true",
            help="Flag to run DMFT. Checks for a previous DMFT calculation and runs only if it is incomplete.",
        )
        type_parser.add_argument(
            "-hf",
            action="store_true",
            help="Flag to perform Hartree-Fock calculation to the correlated orbitals.",
        )
        parser.add_argument(
            "-force",
            action="store_true",
            help="Flag to force DMFT or HF calculation even if a previous calculation has been completed.",
        )
        parser.add_argument(
            "-structurename",
            type=str,
            help="Name of the structure. Not required for VASP. ",
            default=None,
        )
        parser.add_argument(
            "-aiida", help="Flag for aiida calculation. ", action="store_true"
        )
        parser.add_argument(
            "-kmeshtol",
            default=0.00001,
            type=float,
            help="The tolerance to control if two k-points belong to the same shell in wannier90.",
        )
        parser.add_argument(
            "-lowdin", action="store_true", help="Flag to use Siesta Lowdin version.",
        )

        parser.add_argument(
            "-nowin",
            action="store_true",
            help="Flag to disable automatic generation of .win file.",
        )

        parser.add_argument(
            "-v", action="store_true", help="Enable verbosity.",
        )

        args = parser.parse_args()
        DMFTLauncher(args)

    else:
        print("Usage: DMFT.py -h")
