#!/usr/bin/env python

import argparse
import glob
import os
import re
import shutil
import socket
import subprocess
import sys
import time
from argparse import RawTextHelpFormatter
from collections import OrderedDict
from copy import deepcopy
from os.path import getsize

import scipy.interpolate
from scipy import *

import DMFT_MOD
import Fileio
import IMP_SOLVER
import Struct
import VASP

###########################################################################
#### This program executes VASP+DMFT using CTQMC impurity solver ##########
#### This part controls the overall self-consistent steps        ##########
############################################  Dr. Hyowon Park    ##########


def now():
    return time.strftime(" at %c")  #%Y:%H HH:%M MM:%S SS')


def CreateINCAR(params_vasp):
    " Creates input file (INCAR) for vasp"
    f = open("INCAR", "w")
    for p in params_vasp:
        print >> f, p, "  ", params_vasp[p][0], "\t", params_vasp[p][1]
    f.close()


def check_convergence():
    """
    Check if DMFT convergence is reached by checking self-energy of
    last 3 iterations
    """
    # sigdiff tolerance
    if list(p.keys()).count("sig_tol"):
        if p["sig_tol"]:
            sigdifftol = p["sig_tol"]
        else:
            sigdifftol = 1e-03
    else:
        sigdifftol = 1e-03

    sigdiff_array = []

    # file to save self-energy difference
    dSig = open("dSigma.dat", "a+")

    for i in range(1, 4):
        # opening INFO_ITER
        fi = open("INFO_ITER", "r")
        lastline = fi.readlines()[-i]
        fi.close()

        # self energies
        sig1 = float(lastline.split()[4])
        sig2 = float(lastline.split()[5])

        sigdiff = abs(sig1 - sig2)

        if i == 1:
            print ("dSigma = {:.6f}".format(sigdiff))
            dSig.write("{:.6f}\n".format(sigdiff))

        if sigdiff < sigdifftol:
            sigdiff_array.append(True)
        else:
            sigdiff_array.append(False)

    dSig.close()

    if all(sigdiff_array):
        return True
    else:
        return False

def write_iterations_log(itt, it):
    """Writes to iterations.log to preserve iterations history
    through resumed calculations."""

    if itt == 1 and it == 1:
        # add new line to each new DMFT calculation
        with open("iterations.log","a") as fi:
            fi.write(" ")

    # read the last line
    fi = open("iterations.log","r")
    lines = fi.readlines()
    fi.close()

    # now edit the last line of the list of lines
    new_last_line = "{:d} {:d}\n".format(itt, it)
    lines[-1] = new_last_line

    # now write the modified list back out to the file
    with open("iterations.log",'w') as fi:
        fi.writelines(lines)


if __name__ == "__main__":
    # top level parser
    des = "This script performs the  DMFT calculation."
    parser = argparse.ArgumentParser(
        description=des, formatter_class=RawTextHelpFormatter
    )
    parser.add_argument(
        "-structurename",
        type=str,
        help="Name of the structure. Not required for VASP. ",
        default=None,
    )
    parser.add_argument(
        "-dft",
        default="vasp",
        type=str,
        help="Choice of DFT code for the DMFT calculation.",
        choices=["vasp", "siesta", "qe"],
    )
    parser.add_argument(
        "-aiida", help="Flag for aiida calculation. ", type=str, default="False"
    )
    parser.add_argument(
        "-hf",
        action="store_true",
        help="Flag to perform Hartree-Fock calculation to the correlated orbitals.",
    )
    args = parser.parse_args()

    # setting aiida flag from string to boolean.
    if args.aiida == "True":
        args.aiida = True
    elif args.aiida == "False":
        args.aiida = False

    # END OF ARGPARSE SECTION #

    execfile("INPUT.py")  # Read input file

    main_out = open("INFO_TIME", "w")
    main_iter = open("INFO_ITER", "w")
    DMFT_iter = open("INFO_KSUM", "w")
    DM_iter = open("INFO_DM", "w")
    E_iter = open("INFO_ENERGY", "w")
    DFT_iter = open("INFO_DFT_loop", "w")

    if os.path.exists("para_com.dat"):
        fipa = open("para_com.dat", "r")
        para_com = str(fipa.readline())[:-1]
        fipa.close()
    else:
        para_com = "mpirun -np 1"

    if os.path.exists("para_com_dft.dat"):
        fid = open("para_com_dft.dat")
        para_com_dft = str(fid.readline())[:-1]
        fid.close()
    else:
        para_com_dft = para_com

    if os.path.exists("para_com2.dat"):
        fipa = open("para_com2.dat", "r")
        para_com2 = str(fipa.readline())[:-1]
        fipa.close()
    else:
        para_com2 = ""

    # delete dSigma.dat file on start
    if os.path.exists("dSigma.dat"):
        os.remove("dSigma.dat")

    ########### Global variables ###########################

    # Number of wannier bands to be set in the .win file
    # when it is not the same as the DFT bands
    wanbands = 0

    updatewanbands = False

    ############ Initial Preparation ########################

    if p["Nit"] > 0 and p["Niter"] > 1:
        main_out.write(
            "-----------Charge Self-consistent DFT+DMFT calculations-----------"
        )
        print ("\nCalculation type : Charge self-consistent DFT+DMFT calculation\n")

    if p["Nit"] > 0 and p["Niter"] == 1:
        main_out.write(
            "----------Non-charge Self-consistent DFT+DMFT calculations-----------"
        )
        print ("\nCalculation type : Non-charge self-consistent DFT+DMFT calculation\n")
    main_out.write("\n")
    main_out.flush()

    main_out.write("Calculation Starts" + now())
    main_out.write("\n")
    main_out.flush()

    DMFT_iter.write(
        "%12s %12s %12s %12s %12s %12s %12s' "
        % ("# mu", "TOT_ELEC", "Nd[0]", "Nd[-1]", "EKIN", "Sigoo-Vdc", "<d_epsk>")
    )
    DMFT_iter.write("\n")
    DMFT_iter.flush()

    DM_iter.write("%12s" % ("# Occupancy"))
    DM_iter.write("\n")
    DM_iter.flush()

    main_iter.write(
        "%3s %6s %8s %8s %16s %16s %16s %16s"
        % (
            "#",
            "# DMFT",
            "Nd_latt",
            "Nd_imp",
            "(Sigoo-Vdc)_latt",
            "(Sigoo-Vdc)_imp",
            "TOT_E(Tr(SigG))",
            "TOT_E(EPOT_imp)",
        )
    )
    main_iter.write("\n")
    main_iter.flush()

    E_iter.write(
        "%12s %18s %18s %12s"
        % ("# DFT_E", "<epsk>-<epsk>_o", "1/2*Tr(Sigma*G)", "EPOT_imp")
    )
    E_iter.write("\n")
    E_iter.flush()

    DFT_iter.write("%10s %10s %10s %10s" % ("# mu", "TOT_ELEC", "Nd[0]", "EKIN"))
    DFT_iter.write("\n")
    DFT_iter.flush()

    CHGDIFF = 0.0
    #### Read POSCAR ########
    TB = Struct.TBstructure("POSCAR", p["atomnames"], p["orbs"])
    cor_at = p["cor_at"]
    cor_orb = p["cor_orb"]
    TB.Compute_cor_idx(cor_at, cor_orb)

    # ordering dictionary items to print more clearly
    print ("Wannier orbitals in correlated subspace:")
    orb_dic = TB.TB_orbs
    l = list(orb_dic.items())
    l.sort()
    orb_dic = OrderedDict(l)
    for i in orb_dic:
        print ("%s : %s" % (i, orb_dic[i]))

    if TB.LHF == ".TRUE.":
        p["Nd_qmc"] = 0
    U = p["U"]
    J = p["J"]
    T = 1.0 / pC["beta"][0]
    noms = p["noms"]

    ########### Create symmetry index ################
    N_atoms = 0
    for i, ats in enumerate(cor_at):
        N_atoms += len(ats)
    atm_idx = zeros(N_atoms, dtype=int)
    loc_idx = -1
    for i, ats in enumerate(cor_at):
        loc_idx += 1
        for at in ats:
            at_idx = int(at[-1]) - 1
            atm_idx[at_idx] = loc_idx
    # print atm_idx

    sym_idx = []
    loc_idx = 0
    for i, ats in enumerate(cor_at):
        d_orb = TB.TB_orbs[ats[0]]
        # idx=zeros(len(d_orb),dtype=int)
        idx = zeros(len(d_orb) * p["nspin"], dtype=int)
        for ispin in range(p["nspin"]):
            for j, orbs in enumerate(cor_orb[i]):
                loc_idx += 1
                for orb in orbs:
                    # idx[d_orb.index(orb)]=loc_idx
                    idx[d_orb.index(orb) + ispin * len(d_orb)] = loc_idx
        sym_idx.append(idx.tolist())

    DMFT = DMFT_MOD.DMFT_class(p, pC, TB)
    DFT = VASP.VASP_class(
        dft=args.dft, structurename=args.structurename, aiida=args.aiida
    )

    ETOT_old = 0.0
    ETOT2_old = 0.0
    ETOT = 0.0
    ETOT2 = 0.0
    CHGDIFF = 0.0
    CHGDIFF2 = 0.0
    shutil.copy2("sig.inp", "sig.inp.0")

    # ----------------------- Starting DFT+DMFT loop -----------------------------

    for itt in range(p["Niter"]):
        main_out.write("--- Starting charge loop " + str(itt + 1) + now() + "---")
        main_out.write("\n")
        main_out.flush()

        Fileio.Create_dmft_params(p, pC, N_atoms, atm_idx, sym_idx)

        for it in range(p["Nit"]):
            main_out.write("--- Starting DMFT loop " + str(it + 1) + now() + "---")
            print ("\n----- Starting DMFT loop : %s -----\n" % (str(it + 1)))
            main_out.write("\n")
            main_out.flush()

            if itt == 0 and it == 0:
                for i, ats in enumerate(cor_at):
                    if p["nspin"] == 2:
                        pD["para="] = 0
                    else:
                        pD["para="] = 1
                    if p["orbs"][0] == "f":
                        pD["l="] = 3
                    elif p["orbs"][0] == "d":
                        pD["l="] = 2
                    pD["J="] = float(J[i])
                    pD["Eimp="] = zeros(
                        p["nspin"] * len(cor_orb[i])
                    )  # array(ed[i])-ed[i][0]+array(DMFT.sig_st[i])-DMFT.sig_st[i][0]
                    IMP_SOLVER.Create_atomd(pD)
                    IMP_SOLVER.Create_Trans(
                        TB.ncor_orb_list[i], p["nspin"], ats[0], cor_orb[i], TB
                    )
                    cmd = "cp Trans.dat Trans" + str(i + 1) + ".dat"
                    print os.popen(cmd).read()
                    cmd = (
                        p["path_bin"]
                        + "atom_d.py 'atom_d.inp' 2> atom_d.error || { echo 'ATOMD run failed!'; exit 1; }"
                    )
                    # cmd = "./atom_d.py 'atom_d.inp' 2> atom_d.error || { echo 'ATOMD run failed!'; exit 1; }"
                    out, err = subprocess.Popen(
                        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
                    ).communicate()
                    print out  # , err
                    cmd = "cp UC.dat UC" + str(i + 1) + ".dat"
                    print os.popen(cmd).read()

            #

            if it == 0:

                DMFT.EKIN0 = 0
                #print ("Running dmft0.x...")
                print ("Running XHF0.py...")

                # XHF0.py and WAN90.py require wannier90. files instead of
                # structurename. files. Therefore, we need to rename them to
                # wannier90. format when running with Siesta or Quantum Espresso.

                if args.dft == "siesta" or args.dft == "qe":
                    shutil.copy(args.structurename + ".eig", "wannier90.eig")
                    shutil.copy(args.structurename + ".chk", "wannier90.chk")
                    shutil.copy(args.structurename + ".win", "wannier90.win")
                    shutil.copy(args.structurename + ".amn", "wannier90.amn")

                # cmd = (
                #     para_com
                #     + " "
                #     + p["path_bin"]
                #     + "dmft0.x > ksum_output_dmft0.x 2> ksum_error_dmft0.x"
                # )
                cmd = (
                    para_com
                    + " "
                    + p["path_bin"]
                    + "XHF0.py > ksum_output_XHF0 2> ksum_error_XHF0"
                )

                out, err = subprocess.Popen(
                    cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
                ).communicate()

                fiEKIN = open("INFO_KSUM", "r")
                lines = fiEKIN.readlines()
                DMFT.EKIN0 = float(lines[-1].split()[4])
                fiEKIN.close()

            if TB.LHF == ".FALSE.":
                # cmd = para_com+" "+p['path_bin']+"dmft_ksum_sp > ksum_output 2> ksum_error"
                print ("Running dmft.x...")
                cmd = (
                    para_com
                    + " "
                    + p["path_bin"]
                    + "dmft.x > ksum_output_dmft.x 2> ksum_error_dmft.x"
                )
                out, err = subprocess.Popen(
                    cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
                ).communicate()
            else:
                # print ("Running XHF...")
                print ("Running dmft.x for HF...")
                cmd = (
                    para_com
                    + " "
                    + p["path_bin"]
                    + "dmft.x > ksum_output_dmft.x_HF 2> ksum_error_dmft.x_HF"  # should it be XHF0.py?
                )
                out, err = subprocess.Popen(
                    cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
                ).communicate()
            print out, err

            fi = open("DMFT_mu.out", "r")
            DMFT.mu = float(fi.readline())
            fi.close()

            DMFT.Update_Sigoo_and_Vdc(TB, "sig.inp", p["nspin"])
            DMFT.Update_Nlatt(TB, p["nspin"])

            main_out.write("--- DMFT Ksum is finished " + now() + "---")
            main_out.write("\n")
            main_out.flush()

            # Ed_loc = array(loadtxt("Ed.out"))
            Ed_loc = array(loadtxt("Ed.out"), ndmin=1)
            print "\nEd.out:", Ed_loc
            Ed = []
            om_loc, Delta = Fileio.Read_complex_multilines("Delta.out")
            loc_idx = -1
            for i, ats in enumerate(cor_at):
                Delta_loc = zeros(
                    (len(cor_orb[i]) * p["nspin"], len(om_loc)), dtype=complex
                )
                Ed.append([])
                for ispin in range(p["nspin"]):
                    for j, orbs in enumerate(cor_orb[i]):
                        loc_idx += 1
                        print "---i:", i, " ", cor_at[i], " ---loc_idx:", loc_idx
                        Ed[i].append(Ed_loc[loc_idx])
                        Delta_loc[j + ispin * len(cor_orb[i]), :] = deepcopy(
                            Delta[loc_idx, :]
                        )
                Fileio.Print_complex_multilines(
                    Delta_loc, om_loc, "Delta" + str(i + 1) + ".inp"
                )
            Ed = array(Ed)
            print "Ed:", Ed
            if TB.LHF == ".FALSE.":
                IMP_SOLVER.RUN_CTQMC(
                    p, pC, pD, it, itt, para_com, DMFT.mu, Ed, DMFT.Vdc, args.hf
                )
            main_out.write("--- Impurity solver is finished " + now() + "---")
            main_out.write("\n")
            main_out.flush()

            ################# Mix sig.inp ###############################3

            DMFT.Read_Sig(TB, p["nspin"])
            DMFT.Compute_Energy(DFT, TB, Ed)
            DMFT.Compute_Sigoo_and_Vdc(p, TB)
            DMFT.Mix_Sig_and_Print_sig_inp(
                TB, p["Nd_qmc"], p["mix_sig"], "sig.inp", p["nspin"]
            )
            shutil.copy2("sig.inp", "sig.inp." + str(itt + 1) + "." + str(it + 1))
            shutil.copy2("G_loc.out", "G_loc.out." + str(itt + 1) + "." + str(it + 1))

            DMFT.ETOT = DFT.E - DMFT.EKIN0 + DMFT.EKIN + DMFT.EPOT - DMFT.Edc
            DMFT.ETOT2 = DFT.E - DMFT.EKIN0 + DMFT.EKIN + DMFT.EPOT2 - DMFT.Edc_imp
            # Print out results
            main_iter.write(
                "%3d %3d %12.6f %12.6f %14.6f %14.6f %16.6f %16.6f %16.6f"
                % (
                    itt + 1,
                    it + 1,
                    DMFT.Nd_latt[0],
                    DMFT.Nd_imp[0],
                    DMFT.Sigoo[0][0] - DMFT.Vdc[0],
                    DMFT.Sigoo_imp[0][0] - DMFT.Vdc_imp[0],
                    DMFT.ETOT,
                    DMFT.ETOT2,
                    CHGDIFF,
                )
            )
            main_iter.write("\n")
            main_iter.flush()

            # write to iteratinos.log
            write_iterations_log(itt+1, it+1)

            E_iter.write(
                "%14.6f %14.6f %14.6f %14.6f"
                % (
                    DFT.E,
                    DMFT.EKIN - DMFT.EKIN0,
                    DMFT.EPOT - DMFT.Edc,
                    DMFT.EPOT2 - DMFT.Edc_imp,
                )
            )
            E_iter.write("\n")
            E_iter.flush()

            if itt > 2 or it > 2:
                # Check if convergence is reached
                sigdiff = check_convergence()

                if sigdiff:
                    print ("\nConvergence reached.")
                    print ("\nCalculation complete.")
                    main_out.write("\nConvergence reached." + now())
                    main_out.write("\nCalculation Ends" + now())
                    main_out.write("\n")
                    main_out.flush()
                    sys.exit()
                else:
                    pass

        # ------------------------- Full charge self-consistent DFT+DMFT calculation -------------------------------------

        if itt < p["Niter"] - 1:

            if args.dft == "vasp":

                main_out.write("--- Running vaspDMFT " + now() + "---")
                main_out.write("\n")
                main_out.flush()
                print ("\n----- Starting DFT loop : %s -----\n" % (str(itt + 2)))
                print ("\n--- Running vaspDMFT ---\n")

                if itt == 0:
                    f = open("INCAR", "a")
                    print >> f, "NELM= " + str(p["Ndft"])
                    f.close()

                cmd = (
                    para_com_dft
                    + " "
                    + p["path_bin"]
                    + "vaspDMFT > vasp.out 2> vasp.error || { echo 'Parallel run failed!'; exit 1; }"
                )
                out, err = subprocess.Popen(
                    cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
                ).communicate()
                shutil.copy2("CHGCAR", "CHGCAR." + str(itt))
                if itt > 0:
                    DFT.Read_NELECT()
                    CHGDIFF = DFT.Diff_CHGCAR(
                        "CHGCAR." + str(itt - 1), "CHGCAR." + str(itt)
                    )
                    print ("Charge difference = %f" % CHGDIFF)

                DFT.Read_NBANDS()
                DFT.Read_EFERMI()

                # Setting num_bands in .win file.
                # If set to False num_bands is set to number of DFT bands.
                if list(p.keys()).count("num_bands_win"):
                    if p["num_bands_win"]:
                        wanbands = p["num_bands_win"]
                        updatewanbands = False
                    else:
                        wanbands = DFT.NBANDS
                else:
                    wanbands = DFT.NBANDS

                # Updating .win
                DFT.Update_win(
                    wanbands, DFT.EFERMI + p["ewin"][0], DFT.EFERMI + p["ewin"][1]
                )

                print os.popen("rm wannier90.chk").read()
                print os.popen("rm wannier90.chk.fmt").read()
                main_out.write(
                    "-------------- Running wannier 90 "
                    + str(itt + 1)
                    + "----------------"
                )
                print ("Running wannier90...")
                main_out.write("\n")
                main_out.flush()
                # parallel support is available for wannier90 versions above v1.2
                cmd = para_com + " " + p["path_bin"] + "wannier90.x wannier90"
                # cmd = p["path_bin"] + "wannier90.x wannier90"
                out, err = subprocess.Popen(
                    cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
                ).communicate()
                if os.path.isfile("wannier90.chk"):
                    print ("wannier90 calculation complete.")
                    print out  # , err
                else:
                    print ("wannier90 calculation failed! Exiting.")
                    sys.exit()

            elif args.dft == "siesta":

                main_out.write("--- Running siesta " + now() + "---")
                main_out.write("\n")
                main_out.flush()
                print ("\n----- Starting DFT loop : %s -----\n" % (str(itt + 2)))
                print ("\n--- Running siesta ---\n")

                # Renaming wannier90.win files to siesta files.
                # shutil.copy("wannier90.eig", args.structurename + ".eig")
                # shutil.copy("wannier90.amn", args.structurename + ".amn")
                # shutil.copy("wannier90.chk", args.structurename + ".chk")
                # shutil.copy("wannier90.win", args.structurename + ".win")

                # Copying the .psf and .fdf files from top directory
                fdfsource = "../" + args.structurename + ".fdf"
                shutil.copy(fdfsource, ".")
                for file in glob.glob(r"../*.psf"):
                    shutil.copy(file, ".")

                # Running wannier90 pre-processor to obtain .nnkp file
                cmd = p["path_bin"] + "wannier90.x" + " -pp" + " " + args.structurename
                out, err = subprocess.Popen(
                    cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
                ).communicate()
                print (out)

                # Running siesta
                cmd = (
                    para_com_dft
                    + " "
                    + p["path_bin"]
                    + "siesta"
                    + "<"
                    + args.structurename
                    + ".fdf>"
                    + args.structurename
                    + ".out"
                    + " 2>"
                    + "siesta-dmft.error"
                )
                out, err = subprocess.Popen(
                    cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
                ).communicate()

                if os.path.exists(args.structurename + ".out"):
                    fi = open(args.structurename + ".out", "r")
                    done_word = fi.readlines()[-1]
                    fi.close()
                    if done_word.split()[0] == "Job":
                        print ("siesta calculation complete.\n")
                    else:
                        print ("siesta calculation failed!\n")
                        sys.exit()
                else:
                    print ("siesta calculation failed!\n")
                    sys.exit()

                # need to rename .eigW to .eig to run wannier90
                shutil.copy(args.structurename + ".eigW", args.structurename + ".eig")

                # Update Bands and Fermi energy
                DFT.Read_NBANDS()
                DFT.Read_EFERMI()

                # Setting num_bands in .win file.
                # If set to False num_bands is set to number of DFT bands.
                if list(p.keys()).count("num_bands_win"):
                    if p["num_bands_win"]:
                        wanbands = p["num_bands_win"]
                        updatewanbands = False
                    else:
                        wanbands = DFT.NBANDS
                else:
                    wanbands = DFT.NBANDS

                # Copying seedname.win to wannier90.win to update bands and fermi.
                shutil.copy(args.structurename + ".win", "wannier90.win")
                # Updating .win
                DFT.Update_win(
                    wanbands, DFT.EFERMI + p["ewin"][0], DFT.EFERMI + p["ewin"][1]
                )
                # Copying back to seedname.win
                shutil.copy("wannier90.win", args.structurename + ".win")

                # Running wannier90
                chk_seedname_rm = "rm " + args.structurename + ".chk"
                chk_seedname_fmt_rm = "rm " + args.structurename + ".chk.fmt"

                print os.popen(chk_seedname_rm).read()
                print os.popen(chk_seedname_fmt_rm).read()
                main_out.write(
                    "-------------- Running wannier 90 "
                    + str(itt + 1)
                    + "----------------"
                )
                print ("Running wannier90...")
                main_out.write("\n")
                main_out.flush()
                # parallel support is available for wannier90 versions above v1.2
                cmd = (
                    para_com
                    + " "
                    + p["path_bin"]
                    + "wannier90.x"
                    + " "
                    + args.structurename
                )
                # cmd = p["path_bin"] + "wannier90.x wannier90"
                out, err = subprocess.Popen(
                    cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
                ).communicate()
                if os.path.isfile(args.structurename + ".chk"):
                    print ("wannier90 calculation complete.")
                    print out  # , err
                else:
                    print ("wannier90 calculation failed! Exiting.")
                    sys.exit()

                # Renaming siesta files to wannier90 files
                # shutil.copy(args.structurename + ".eig", "wannier90.eig")
                # shutil.copy(args.structurename + ".chk", "wannier90.chk")
                # shutil.copy(args.structurename + ".win", "wannier90.win")
                # shutil.copy(args.structurename + ".amn", "wannier90.amn")

    main_out.write("Calculation Ends" + now())
    print ("\nCalculation Ends.")
    main_out.write("\n")
    main_out.flush()
