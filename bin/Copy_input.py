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
    parser.add_argument(
        "-dft",
        default="vasp",
        type=str,
        help="Choice of DFT code for the DMFT calculation.",
        choices=["vasp", "siesta", "qe"],
    )
    args = parser.parse_args()

    cpdr = args.cpdr
    post = args.post
    structure = args.structurename
    dft = args.dft

    if not os.path.exists(cpdr):
        print("The directory " + cpdr + " does not exists!")
        sys.exit(1)

    if post == "dos":
        dosfiles = [
            "wannier90.chk",
            "wannier90.eig",
            "DFT_mu.out",
            "DMFT_mu.out",
            "dmft_params.dat",
            "INPUT.py",
            "para_com.dat",
        ]
        for files in dosfiles:
            if os.path.exists(cpdr + "/" + files):
                print("Copying dos file " + files + " to the current directory")
                shutil.copy2(cpdr + "/" + files, ".")
            else:
                print(files + " must exist in a " + cpdr + " directory! Exiting!")
                sys.exit(1)
        exec(compile(open("INPUT.py", "rb").read(), "INPUT.py", "exec"))

    elif post == "bands":
        bandsfiles = [
            "wannier90.chk",
            "wannier90.eig",
            "DMFT_mu.out",
            "dmft_params.dat",
            "INPUT.py",
            "para_com.dat",
            "SigMdc.out",
            "SigMdc_dn.out",
        ]
        for files in bandsfiles:
            if os.path.exists(cpdr + "/" + files):
                print("Copying bands file " + files + " to the current directory")
                shutil.copy2(cpdr + "/" + files, ".")
            else:
                print(files + " must exist in a " + cpdr + " directory! Exiting!")
                sys.exit(1)
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
            "Siesta": [str(structure) + ".out", "para_com.dat", "POSCAR"],
            "aiida": ["aiida.out", "POSCAR"],
            "Quantum Espresso": [str(structure) + ".scf.out", "para_com.dat", "POSCAR"],
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
                if dft == "vasp":
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
                if dft == "siesta":
                    data = fi.read()
                    Efermi = float(
                        re.findall(r"Fermi\s=[\s0-9+-.]*", data)[0].split()[-1]
                    )
                    savetxt("./DFT_mu.out", array([EFermi]))
                    fi.close()

        Wannierfiles = [
            "wannier90.chk",
            "wannier90.eig",
            "wannier90.win",
            "wannier90.amn",
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
