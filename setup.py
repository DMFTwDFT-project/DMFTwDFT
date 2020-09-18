#!/usr/bin/env python
""" DMFTwDFT setup.

Please copy Makefile.in from config directory.
Then run setup.py <compiler>.

E.g. - setup.py intel

All the executables will be copied to the bin directory.
Don't forget to install wannier90 and recompile VASP with wannier90.
Also copy wannier90.x and w90chk2chk.x to the bin directory.

"""
import sys
import os
import shutil
import requests
import subprocess
import tarfile
import glob
import argparse
from argparse import RawTextHelpFormatter

sys.path.insert(1, "./bin")
import splash


def main(args):

    """Installation main function."""

    # call cleanup
    cleanup()

    # print welcome message
    splash.welcome()

    # --------------- COMPILING INTERNAL SOURCES -----------------------------

    compiler = str(args.compiler)

    # Running the Makefile to compile internal sources.
    if compiler == "intel":
        print("Compiler : intel\n")
        shutil.copy("./sources/make.inc.intel", "./sources/make.inc")
    elif compiler == "gfortran":
        print("Compiler : gfortran\n")
        shutil.copy("./sources/make.inc.gfortran", "./sources/make.inc")

    print("Compiling internal sources...\n")
    cmd = "cd sources; make clean; make all > internal.log 2>&1 "
    out, err = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    ).communicate()

    # Checking if all internal sources have been compiled
    file_list = [
        "dmft.x",
        "dmft_dos.x",
        "libdmft.a",
        "./dmft_ksum/dmft_ksum_band",
        "./dmft_ksum/dmft_ksum_partial_band",
        "./fort_kpt_tools/fort_kpt_tools.so",
    ]

    result_array = []

    for fi in file_list:
        result = os.path.exists("./sources/" + fi)
        result_array.append(result)
        print("Compiled file %s exists : %s " % (fi, result))

    if all(result_array):
        print("Internal compilation complete.")
    else:
        print(
            "Internal compilation failed! Check internal.log for details. Make sure gfortran.make.inc or intel.make.inc in sources points to the correct libraries."
        )
        sys.exit()

    # --------------- COMPILING EXTERNAL SOURCES -----------------------------

    # Download from Professor Haule's EDMFT website.
    url = "http://hauleweb.rutgers.edu/downloads/EDMFTF.tgz"
    print("\nCompiling external sources...")
    print("Downloading sources from EDMFTF...")
    r = requests.get(url, stream=True)
    with open("./sources/EDMFTF.tgz", "wb") as edmft_file:
        for chunk in r.iter_content(chunk_size=1024):
            # writing one chunk at a time to pdf file
            if chunk:
                edmft_file.write(chunk)

    # Extracting EDMFTF.tgz file
    print("Extracting...\n")
    tar = tarfile.open("./sources/EDMFTF.tgz")
    tar.extractall("./sources/")
    tar.close()

    # Getting current EDMFTF folder
    EDMFTF_list = []
    for foldername in glob.glob("./sources/EDMFTF*/"):
        EDMFTF_list.append(foldername)
    EDMFTF_folder = EDMFTF_list[-1]

    # Copy Makefile.in to /EDMFTF/src/ directory.
    src_dir = EDMFTF_folder + "/src/"
    shutil.copy("Makefile.in", src_dir)

    # Compiling ctqmc
    ctqmc_dir = EDMFTF_folder + "/src/impurity/ctqmc/"
    print("Compiling ctqmc...")
    cmd = "cd " + ctqmc_dir + "; make clean; make ctqmc > ctqmc.log 2>&1 "
    out, err = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    ).communicate()
    if os.path.exists(ctqmc_dir + "ctqmc"):
        print("Complete.\n")
    else:
        print("ctqmc compilation failed! Check ctqmc.log for details.")
        sys.exit()
    # Copy to bin directory
    shutil.copy(ctqmc_dir + "ctqmc", "./bin/")

    # Compiling atomd (gaunt.so, gutils.so)
    atomd_dir = EDMFTF_folder + "/src/impurity/atomd/"
    print("Compiling atomd : gaunt.so, gutils.so...")
    cmd = "cd " + atomd_dir + "; make clean; make all > atomd.log 2>&1 "
    out, err = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    ).communicate()
    if os.path.exists(atomd_dir + "gaunt.so") and os.path.exists(
        atomd_dir + "gutils.so"
    ):
        print("Complete.\n")
    else:
        print("atomd compilation failed! Check atomd.log for details.")
        sys.exit()
    # Copy to bin directory
    shutil.copy(atomd_dir + "gaunt.so", "./bin/")
    shutil.copy(atomd_dir + "gutils.so", "./bin/")

    # Compiling maxent_routines
    maxent_dir = EDMFTF_folder + "/src/impurity/maxent_source/"
    print("Compiling maxent_routines...")
    cmd = "cd " + maxent_dir + "; make clean; make all > maxent_routines.log 2>&1 "
    out, err = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    ).communicate()
    if os.path.exists(maxent_dir + "maxent_routines.so"):
        print("Complete.\n")
    else:
        print(
            "maxent_routines compilation failed! Check maxent_routines.log for details."
        )
        sys.exit()
    # Copy to bin directory
    shutil.copy(maxent_dir + "maxent_routines.so", "./bin/")

    # Compiling skrams
    skrams_dir = EDMFTF_folder + "/src/impurity/skrams/"
    print("Compiling skrams...")
    cmd = "cd " + skrams_dir + "; make clean; make all > skrams.log 2>&1 "
    out, err = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    ).communicate()
    if os.path.exists(skrams_dir + "skrams"):
        print("Complete.\n")
    else:
        print("skrams compilation failed! Check skrams.log for details.")
        sys.exit()
    # Copy to bin directory
    shutil.copy(skrams_dir + "skrams", "./bin/")

    # Compilation complete
    print("DMFTwDFT compilation complete!")
    print(
        "Please add the bin directory to $PATH and $PYTHONPATH variables in your .bashrc."
    )
    print("Thank you!")


def cleanup():
    """ Cleanup. """
    if os.path.exists("./sources/internal.log"):
        os.remove("./sources/internal.log")
    if os.path.exists("./sources/EDMFTF.tgz"):
        os.remove("./sources/EDMFTF.tgz")
    try:
        for foldername in glob.glob("./sources/EDMFTF*"):
            shutil.rmtree(foldername)
    except (FileNotFoundError, IOError):
        pass

    # Cleaning bin folder
    bin_files = [
        "dmft.x",
        "dmft_dos.x",
        "dmft_ksum_band",
        "dmft_ksum_partial_band",
        "fort_kpt_tools.so",
        "ctqmc",
        "gaunt.so",
        "gutils.so",
        "skrams",
        "maxent_routines.so",
    ]

    for bin_i in bin_files:
        if os.path.exists("./bin/" + bin_i):
            os.remove("./bin/" + bin_i)


if "__main__" == __name__:
    args = sys.argv[1:]
    if args:
        # top level parser
        parser = argparse.ArgumentParser(
            description="DMFTwDFT setup. \n Please copy Makefile.in from config directory.",
            formatter_class=RawTextHelpFormatter,
        )
        parser.add_argument(
            "compiler",
            type=str,
            help="Compiler.",
            choices=["intel", "gfortran"],
            default="intel",
        )
        args = parser.parse_args()
        main(args)
    else:
        print(
            "Usage: setup.py <compiler> \n Please copy Makefile.in from config directory."
        )
        print("OPTIONS : {intel, gfortran}")
