#!/usr/bin/env python2
import argparse
import glob
import os
import re
import shutil
import sys
from argparse import RawTextHelpFormatter
from os.path import getsize

from scipy import *

if __name__ == "__main__":

    # top level parser
    des = "This tool copies all the required files necessary for the DMFT calculation."
    parser = argparse.ArgumentParser(
        description=des, formatter_class=RawTextHelpFormatter
    )
    parser.add_argument("cpdr", default="../", type=str, help="DFT directory")
    parser.add_argument(
        "-structurename",
        type=str,
        help="Name of the structure. Not required for VASP. ",
    )
    parser.add_argument(
        "-post",
        type=str,
        help="Copy files for dos or band structure calculation?",
        choices=["dos", "bands"],
    )
    args = parser.parse_args()

    cpdr = args.cpdr
    post = args.post
    structurename = args.structurename

    if not os.path.exists(cpdr):
        print("The directory " + cpdr + " does not exists!")
        sys.exit(1)

    if os.path.exists(cpdr + "/seedname.dat"):
        shutil.copy2(cpdr + "/seedname.dat", ".")
        fi = open(cpdr + "/seedname.dat", "r")
        structurename = str(fi.readline())
        fi.close()
    else:
        structurename = None

    if post == "dos":

        if structurename is None:
            dosfiles = [
                "wannier90.chk",
                "wannier90.eig",
                "DFT_mu.out",
                "DMFT_mu.out",
                "dmft_params.dat",
                "INPUT.py",
                "para_com.dat",
            ]
        else:
            dosfiles = [
                str(structurename) + ".chk",
                str(structurename) + ".eig",
                "DFT_mu.out",
                "DMFT_mu.out",
                "dmft_params.dat",
                "INPUT.py",
                "para_com.dat",
            ]

        for i, files in enumerate(dosfiles):
            if os.path.exists(cpdr + "/" + files):
                print("Copying dos file " + files + " to the current directory")
                shutil.copy2(cpdr + "/" + files, ".")
            else:
                print(files + " must exist in a " + cpdr + " directory! Exiting!")
                sys.exit(1)

        exec(compile(open("INPUT.py", "rb").read(), "INPUT.py", "exec"))

    elif post == "bands":

        if structurename is None:

            bandsfiles = [
                "wannier90.chk",
                "wannier90.eig",
                "DMFT_mu.out",
                "dmft_params.dat",
                "INPUT.py",
                "para_com.dat",
            ]
        else:
            bandsfiles = [
                str(structurename) + ".chk",
                str(structurename) + ".eig",
                "DMFT_mu.out",
                "dmft_params.dat",
                "INPUT.py",
                "para_com.dat",
            ]

        for i, files in enumerate(bandsfiles):
            if os.path.exists(cpdr + "/" + files):
                print("Copying bands file " + files + " to the current directory")
                shutil.copy2(cpdr + "/" + files, ".")
            else:
                print(files + " must exist in a " + cpdr + " directory! Exiting!")
                sys.exit(1)

        if structurename:
            # dmft_ksum_band looks for wannier90 files.
            # Renaming is just easier.
            shutil.copy(str(structurename) + ".chk", "wannier90.chk")
            shutil.copy(str(structurename) + ".eig", "wannier90.eig")

        exec(compile(open("INPUT.py", "rb").read(), "INPUT.py", "exec"))

    else:
        DFTfiles = {
            "VASP": [
                "OUTCAR",
                "OSZICAR",
                "POSCAR",
                "POTCAR",
                "KPOINTS",
                "INCAR",
                "WAVECAR",
                "para_com.dat",
            ],
            "Siesta": [
                str(structurename) + ".out",
                "para_com.dat",
                "POSCAR",
                "seedname.dat",
            ],
            "aiida": ["aiida.out", "POSCAR"],
            "Quantum Espresso": [
                str(structurename) + ".scf.out",
                "para_com.dat",
                "POSCAR",
                "seedname.dat",
            ],
        }
        Ldft = False
        for dft in DFTfiles:
            outfile = DFTfiles[dft][0]
            if os.path.exists(cpdr + "/" + outfile):
                print(dft + " results have been found in a " + cpdr + " directory!")
                Ldft = True
                break
            if Ldft == False:
                print(
                    "No "
                    + dft
                    + " results have been found in a "
                    + cpdr
                    + " directory!"
                )  # ; break #sys.exit(1)

        if Ldft == True:
            for i, files in enumerate(DFTfiles[dft][:]):
                # copy files
                if os.path.exists(cpdr + "/" + files):
                    print("Copying DFT file " + files + " to the current directory")
                    shutil.copy2(cpdr + "/" + files, ".")
                else:
                    if i < 2:
                        print(
                            files + " must exist in a " + cpdr + " directory! Exiting!"
                        )
                        sys.exit(1)
                    else:
                        print(files + " does not exist in a " + cpdr + " directory!")
                        print(files + " will be needed for charge updates!")
            if os.path.exists(cpdr + "/DFT_mu.out"):
                print("Copying DFT file DFT_mu.out to the current directory")
                shutil.copy2(cpdr + "/DFT_mu.out", ".")
            else:
                print("Making DFT_mu.out file")
                fi = open(cpdr + "/" + outfile, "r")

                # vasp case
                if dft == "VASP":
                    EFermi = 0
                    for line in fi:
                        if re.search("Fermi energy", line) or re.search(
                            "E-fermi", line
                        ):
                            line_fermi = line
                    # print line_fermi
                    val = re.search(r"(\-?\d+\.?\d*)", line_fermi)
                    # print val
                    EFermi = float(val.group(1))
                    savetxt("./DFT_mu.out", array([EFermi]))
                    fi.close()

                # Siesta case
                elif dft == "Siesta":
                    data = fi.read()
                    Efermi = float(
                        re.findall(r"Fermi\s=[\s0-9+-.]*", data)[0].split()[-1]
                    )
                    savetxt("./DFT_mu.out", array([EFermi]))
                    fi.close()

                # Quantum Espresso case
                elif dft == "Quantum Espresso" or dft == "aiida":
                    data = fi.read()
                    Efermi = float(
                        re.findall(r"the\s*Fermi\s*energy\s*is\s*([\s\d.]*)ev", data)[0]
                    )
                    savetxt("./DFT_mu.out", array([EFermi]))
                    fi.close()

        if structurename is None:
            Wannierfiles = [
                "wannier90.chk",
                "wannier90.eig",
                "wannier90.win",
                "wannier90.amn",
            ]
        else:
            Wannierfiles = [
                str(structurename) + ".chk",
                str(structurename) + ".eig",
                str(structurename) + ".win",
                str(structurename) + ".amn",
            ]

        for i, files in enumerate(Wannierfiles):
            # copy files
            if os.path.exists(cpdr + "/" + files):
                print("Copying Wannier file " + files + " to the current directory")
                shutil.copy2(cpdr + "/" + files, ".")

            else:
                print(files + " must exist in a " + cpdr + " directory! Exiting!")
                sys.exit(1)
                # if i<2:
                #   print files+" must exist in a "+cpdr+" directory! Exiting!"; sys.exit(1)
                # else:
                #   print files+" does not exist in a "+cpdr+" directory!"
                #   print files+" will be needed for charge update!"

        DMFTfiles = ["sig.inp", "DMFT_mu.out", "INPUT.py"]
        if os.path.exists(cpdr + "/sig.inp"):
            print("Copying DMFT file sig.inp to the current directory")
            shutil.copy2(cpdr + "/sig.inp", ".")
        else:
            print("sig.inp file does not exist! Must be generated using sigzero.py")
        if os.path.exists(cpdr + "/INPUT.py"):
            print("Copying DMFT file INPUT.py to the current directory")
            shutil.copy2(cpdr + "/INPUT.py", ".")
        else:
            print("sig.inp file does not exist! Must be generated using sigzero.py")
        if os.path.exists(cpdr + "/DMFT_mu.out"):
            print("Copying DMFT file DMFT_mu.out to the current directory")
            shutil.copy2(cpdr + "/DMFT_mu.out", ".")
        else:
            print("DMFT_mu.out file does not exist! Copying from DFT_mu.out")
            shutil.copy2("./DFT_mu.out", "./DMFT_mu.out")

        if os.path.exists(cpdr + "/para_com_dft.dat"):
            print("Copying DFT file para_com_dft.dat to the current directory")
            shutil.copy2(cpdr + "/para_com_dft.dat", ".")
        else:
            print("para_com_dft.dat file does not exist! Copying from para_com.dat")
            shutil.copy2("./para_com.dat", "./para_com_dft.dat")
