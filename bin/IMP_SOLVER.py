#!/usr/bin/env python2

import copy
import os
import re
import shutil
import subprocess
import sys
import time
# import generate_cix
# import generate_5orb
from time import gmtime, strftime

from scipy import *

import Fileio


def CreateInputFile(params_ctqmc):
    " Creates input file (PARAMS) for CT-QMC solver"
    f = open("PARAMS", "w")
    print >> f, "# Input file for continuous time quantum Monte Carlo"
    for p in params_ctqmc.keys():
        print >> f, p, "\t", params_ctqmc[p][0], "\t", params_ctqmc[p][1]
    f.close()


def Create_atomd(params_atomd):
    " Creates input file for cix file"
    f = open("atom_d.inp", "w")
    for p in params_atomd.keys():
        if (
            p == "n="
            or p == "Ewindow="
            or p == "Ex="
            or p == "Ncentral="
            or p == "Ep="
            or p == "Eimp="
        ):
            print >> f, p, "[",
            for elem in params_atomd[p]:
                print >> f, str(elem) + ",",
            print >> f, "]"
        elif p == "CoulombF=":
            print >> f, p, '"' + params_atomd[p] + '"'
        else:
            print >> f, p, params_atomd[p]
    f.close()


def Create_Trans(norb, nspin, at, cor_orbs, TB):

    d_orb = TB.TB_orbs[at]
    orb_Trans = {}
    orb_Trans["f1"] = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0]
    orb_Trans["f2"] = [0.0, 0.0, 1j / sqrt(2), 0.0, 1j / sqrt(2), 0.0, 0.0]
    orb_Trans["f3"] = [0.0, 0.0, 1 / sqrt(2), 0.0, -1 / sqrt(2), 0.0, 0.0]
    orb_Trans["f4"] = [0.0, 1j / sqrt(2), 0.0, 0.0, 0.0, -1j / sqrt(2), 0.0]
    orb_Trans["f5"] = [0.0, 1 / sqrt(2), 0.0, 0.0, 0.0, 1 / sqrt(2), 0.0]
    orb_Trans["f6"] = [1j / sqrt(2), 0.0, 0.0, 0.0, 0.0, 0.0, 1j / sqrt(2)]
    orb_Trans["f7"] = [1 / sqrt(2), 0.0, 0.0, 0.0, 0.0, 0.0, -1 / sqrt(2)]
    orb_Trans["d_z2"] = [0.0, 0.0, 1.0, 0.0, 0.0]
    orb_Trans["d_xz"] = [0.0, 1 / sqrt(2), 0.0, -1 / sqrt(2), 0.0]
    orb_Trans["d_yz"] = [0.0, -1j / sqrt(2), 0.0, -1j / sqrt(2), 0.0]
    orb_Trans["d_x2y2"] = [1 / sqrt(2), 0.0, 0.0, 0.0, 1 / sqrt(2)]
    orb_Trans["d_xy"] = [-1j / sqrt(2), 0.0, 0.0, 0.0, 1j / sqrt(2)]
    f = open("Trans.dat", "w")
    print >> f, norb * nspin, norb * nspin, " #  size of Sigind and CF"
    print >> f, " #---- Sigind follows "
    mat_row = zeros((norb * nspin, norb * nspin), dtype=int)
    loc_idx = 1
    for ispin in range(nspin):
        for j, orbs in enumerate(cor_orbs):
            for orb in orbs:
                idx = TB.idx[at][orb]
                idx2 = d_orb.index(orb)
                if TB.cor_idx[idx] > 1:
                    mat_row[idx2 + ispin * norb, idx2 + ispin * norb] = loc_idx
            loc_idx += 1
    for i in range(norb * nspin):
        print >> f, " ",
        for j in range(norb * nspin):
            print >> f, mat_row[i, j], "  ",
        print >> f, ""
    print >> f, "#---- CF follows"
    for i, orb in enumerate(d_orb):
        print >> f, " ",
        Trans = orb_Trans[orb]
        for elem in Trans:
            print >> f, "%10.8f %10.8f  " % (elem.real, elem.imag),
        print >> f, ""


def RUN_CTQMC(p, pC, pD, it, itt, para_com, mu, ed, vdc, hf):
    U = p["U"]
    J = p["J"]
    cor_at = p["cor_at"]
    cor_orb = p["cor_orb"]
    # mu=loadtxt('DMFT_mu.out')
    # ed=[]
    # fi=open('ED.out','r')
    # for line in fi.readlines():
    #   ed.append(map(float,line.split()))
    # fi.close()
    # sig_st=[]
    # fi=open('sig_st.out','r')
    # for line in fi.readlines():
    #   sig_st.append(map(float,line.split()))
    # fi.close()
    for i in range(len(cor_at)):
        if len(cor_orb[i]) > 0:
            dir_name = "imp." + str(i) + "/"
            if it == 0 and itt == 0:
                cmd = "mkdir " + dir_name
                print os.popen(cmd).read()
                print os.popen(
                    "cp Trans" + str(i + 1) + ".dat " + dir_name + "Trans.dat"
                ).read()
                # cmd = 'rm status.*'
                # print os.popen(cmd).read()
                if os.path.exists("status" + str(i) + ".tar.gz"):
                    print os.popen("mv status" + str(i) + ".tar.gz " + dir_name).read()

            shutil.copy2("Delta" + str(i + 1) + ".inp", dir_name + "Delta.inp")

            os.chdir(dir_name)

            if it == 0 and itt == 0:
                if os.path.exists("status" + str(i) + ".tar.gz"):
                    print os.popen("tar xzvf status" + str(i) + ".tar.gz").read()
                    print os.popen("rm status" + str(i) + ".tar.gz").read()

            if p["orbs"][i] == "f":
                pD["l="] = 3
            elif p["orbs"][i] == "d":
                pD["l="] = 2
            pD["J="] = float(J[i])

            if hf:
                print ("\nSetting up CTQMC for Hartree-Fock...")
                pD["Eimp="] = array(ed[i]) - ed[i][0] - array([0, 100.0])
                print ("Eimp = %s \n" % pD["Eimp="])
            else:
                print ("\nSetting up CTQMC for DMFT...")
                pD["Eimp="] = array(ed[i]) - ed[i][0]
                print ("Eimp = %s \n" % pD["Eimp="])

            # changed ed[i] to ed[0] in sp version.

            Create_atomd(pD)
            # IMP_SOLVER.Create_Trans(TB.ncor_orb_list[i],p['nspin'],ats[0],cor_orb[i],TB)
            # cmd = 'cp Trans'+str(i+1)+'.dat Trans.dat'
            # print os.popen(cmd).read()
            cmd = (
                p["path_bin"]
                + "atom_d.py 'atom_d.inp' 2> atom_d.error || { echo 'ATOMD run failed!'; exit 1; }"
            )
            out, err = subprocess.Popen(
                cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            ).communicate()
            print out  # , err
            # cmd = 'cp actqmc.cix actqmc'+str(i+1)+'.cix'
            # print os.popen(cmd).read()
            # pD['J=']=float(J[i])
            # pD['Eimp=']=array(ed[i])-ed[i][0]
            # if it==0:
            #   #generate_cix.Print_cix(J[i],ed[i])
            #   #generate_5orb.Print_cix(J[i],ed[i][0],ed[i][1])
            #   Create_atomd(pD)
            #   cmd = p['path_bin']+"atom_d.py 'atom_d.inp' 2> atom_d.error || { echo 'ATOMD run failed!'; exit 1; }"
            #   out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
            #   print out #, err
            #   cmd = 'cp actqmc.cix actqmc'+str(i+1)+'.cix'
            #   print os.popen(cmd).read()
            # if it==0:
            #   generate_cix.Print_cix(J[i],ed[i])
            #   cmd = 'cp impurity.cix impurity'+str(i+1)+'.cix'
            #   print os.popen(cmd).read()
            # params_ctqmc['cix'][0]='impurity'+str(i+1)+'.cix'

            pC["mu"] = [
                float(mu - ed[i][0] + vdc[i]),
                "# Chemical potential",
            ]  ### Impurity mu=mu_lattice-Ed
            #         pC['mu']=[float(mu-ed[i]+vdc[i]),]
            pC["U"] = [float(U[i]), "# Coulomb repulsion (F0)"]
            pC["J"] = [float(J[i]), "# Hund's coupling"]
            pC["cix"] = ["actqmc.cix", "# cix file"]
            # pC['Delta'][0]='Delta'+str(i+1)+'.inp'
            pC["Delta"] = ["Delta.inp", "# Delta file"]
            CreateInputFile(pC)
            # cmd = 'cp '+dir_name+'status.* .'
            # print os.popen(cmd).read()
            # Running ctqmc
            print "--- Running qmc for atom", i, "iteration: ", str(it + 1), " ---"
            strftime("%a, %d %b %Y %H:%M:%S", gmtime())

            cmd = (
                para_com
                + " "
                + p["path_bin"]
                + pC["exe"][0]
                + " > nohup.out  2> error"
                + str(i + 1)
                + ".out || { echo 'Parallel run failed!'; exit 1; }"
            )
            out, err = subprocess.Popen(
                cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            ).communicate()
            print out  # , err

            # Some copying to store data obtained so far (at each iteration)
            # cmd = 'cp Gf.out ../Gf'+str(i+1)+'.out'
            # print os.popen(cmd).read() # copying Gf
            shutil.copy2("ctqmc.log", "ctqmc.log." + str(itt) + "." + str(it))
            shutil.copy2("Sig.out", "Sig.out." + str(itt) + "." + str(it))
            shutil.copy2("Gf.out", "Gf.out." + str(itt) + "." + str(it))
            # cmd = 'cp Gf.out Gf.out.'+str(itt)+'.'+str(it)
            # print os.popen(cmd).read() # copying Gf
            ##cmd = 'cp Delta'+str(i+1)+'.inp Delta'+str(i+1)+'.inp.'+str(itt)+'.'+str(it)
            ##print os.popen(cmd).read() # copying Gf
            # cmd = 'cp Sig.out ../Sig'+str(i+1)+'.out'
            # print os.popen(cmd).read() # copying Gf
            # cmd = 'cp Sig.out Sig.out.'+str(itt)+'.'+str(it)
            # print os.popen(cmd).read() # copying Gf
            ##cmd = 'cp Gtau.dat Gtau'+str(i+1)+'.dat'
            ##print os.popen(cmd).read() # copying Gf
            # cmd = 'cp nohup.out nohup.out.'+str(itt)+'.'+str(it)
            # print os.popen(cmd).read() # copying Gf
            # cmd = 'cp Probability.dat Probability.dat.'+str(itt)+'.'+str(it)
            # print os.popen(cmd).read() # copying Gf
            # cmd = 'cp histogram.dat histogram.dat.'+str(itt)+'.'+str(it)
            # print os.popen(cmd).read()

            os.chdir("..")


if __name__ == "__main__":
    execfile("INPUT.py")
    pD["Eimp="] = [1, 2, 3]
    Create_atomd(pD)
