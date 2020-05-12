#!/usr/bin/env python2
import copy
import re
import sys

from scipy import *


def Read_complex_multilines(D_name, skipline=0):
    """ This function reads Sig.out file """
    print "Reading a file ", D_name
    fi = open(D_name, "r")
    for i in range(skipline):
        fi.readline()
    lines = fi.readlines()
    fi.close()
    m = len(lines[0].split())
    Nom = len(lines)
    om_data = zeros(Nom, dtype=float)
    Data = zeros(((m - 1) / 2, Nom), dtype=complex)
    for iom, line in enumerate(lines):
        l = line.split()
        om_data[iom] = float(l[0])
        for ib in range((m - 1) / 2):
            Data[ib, iom] = complex(float(l[1 + ib * 2]), float(l[2 + ib * 2]))
    return om_data, Data


def Read_float_multilines(D_name):
    """ This function reads Sig.out file """
    print "Reading a file ", D_name
    fi = open(D_name, "r")
    lines = fi.readlines()
    fi.close()
    m = len(lines[0].split())
    Nom = len(lines)
    # om_data=zeros(Nom,dtype=float)
    Data = zeros((m, Nom))
    for iom, line in enumerate(lines):
        l = line.split()
        for ib in range(m):
            Data[ib, iom] = float(l[ib])
    return Data


def Read_float(D_name):
    """ This function reads Sig.out file """
    print "Reading a file ", D_name
    fi = open(D_name, "r")
    lines = fi.readlines()
    fi.close()
    Data = array(map(float, lines[0].split()))
    return Data


def Read_complex_Data(D_name):
    """ This function reads Sig.out file """
    print "Reading a file ", D_name
    fi = open(D_name, "r")
    line = fi.readline()
    val = re.search(r"TrSigmaG=(\-?\d+\.?\d*)", line)
    TrSigmaG = float(val.group(1))
    val = re.search(r"mu=(\-?\d+\.?\d*)", line)
    mu = float(val.group(1))
    val = re.search(r"Ekin=(\-?\d+\.?\d*)", line)
    Ekin = float(val.group(1))
    val = re.search(r"Epot=(\-?\d+\.?\d*)", line)
    Epot = float(val.group(1))
    val = re.search(r"nf=(\d+\.?\d*)", line)
    nf_q = float(val.group(1))
    mom = eval(line.split()[-1][4:])
    lines = fi.readlines()
    Nom = len(lines)
    m = len(lines[0].split())
    om_data = zeros(Nom, dtype=float)
    Data = zeros(((m - 1) / 2, Nom), dtype=complex)
    iom = 0
    for line in lines:
        l = line.split()
        om_data[iom] = float(l[0])
        for ib in range((m - 1) / 2):
            Data[ib, iom] = complex(float(l[1 + ib * 2]), float(l[2 + ib * 2]))
        iom = iom + 1
    if iom != Nom:
        print "Something is wrong!"
        sys.exit(1)
    if len(mom) != (m - 1) / 2:
        print "Something is wrong!"
        sys.exit(1)

    return (om_data, Data, TrSigmaG, Epot, nf_q, mom, Ekin, mu)


def Print_complex(data, mesh, filename):
    n1 = len(mesh)
    fi = open(filename, "w")
    for i in range(n1):
        print >> fi, "%.14f " % (mesh[i]),
        print >> fi, "%.14f %.14f " % (data[i].real, data[i].imag)


def Print_complex_multilines(data, mesh, filename, headers=[]):
    n0 = len(data)
    n1 = len(mesh)
    fi = open(filename, "w")
    for header in headers:
        print >> fi, header
    for i in range(n1):
        for j in range(n0):
            if j == 0:
                print >> fi, "%20.15f " % (mesh[i]),
            print >> fi, "%20.15f %20.15f " % (data[j, i].real, data[j, i].imag),
        print >> fi, ""


def Print_float(data, filename):
    n0 = len(data)
    fi = open(filename, "w")
    for j in range(n0):
        print >> fi, "%.14f " % (data[j]),
    print >> fi, ""


def Read_float(filename):
    fi = open(filename, "r")
    Data = array(map(float, fi.readline().split()))
    return Data


def Create_dmft_params(p, pC, N_atoms, atm_idx, sym_idx):
    f = open("dmft_params.dat", "w")
    print >> f, "# Number of k-points in Wannier basis="
    print >> f, p["q"][0], p["q"][1], p["q"][2]
    print >> f, "# Total number of electrons="
    print >> f, p["n_tot"]
    #   print >> f, "# Temperature [eV]="
    #   print >> f, 1.0/pC['beta'][0]
    print >> f, "# Number of om points for k-sum"
    print >> f, p["noms"]
    print >> f, "# Number of iterations for mu"
    print >> f, p["mu_iter"]
    print >> f, "# Number of total spin"
    print >> f, p["nspin"]
    print >> f, "# Number of total correlated atoms"
    print >> f, N_atoms
    print >> f, "# Number of correlated orbitals per atom"
    #   print >> f, len(sym_idx[atm_idx[0]])
    print >> f, len(sym_idx[atm_idx[0]]) / p["nspin"]
    print >> f, "# Orbital index for the self-energy at each atom"
    for i in range(N_atoms):
        for j in range(len(sym_idx[atm_idx[i]])):
            print >> f, sym_idx[atm_idx[i]][j],
        print >> f, ""
    f.close()


def Create_INPUT(p, pC, TB, T_high, noms_high, LFORCE=".FALSE."):
    atm_idx = []
    idx = 1
    for ats in p["cor_at"]:
        for at in ats:
            atm_idx.append(idx)
        idx += 1
    f = open("VASP.input", "w")
    print >> f, LFORCE
    print >> f, TB.LHF
    print >> f, p["n_tot"]
    print >> f, p["nspin"]
    print >> f, p["nfine"]
    print >> f, TB.ncor_orb
    print >> f, TB.max_cor_orb
    if TB.LHF == ".TRUE.":
        print >> f, "1"
    else:
        print >> f, p["noms"]
    if TB.LHF == ".TRUE.":
        print >> f, "1"
    else:
        print >> f, noms_high
    if TB.LHF == ".TRUE.":
        print >> f, "1"
    else:
        print >> f, p["noms"] + p["nomlog"]
    print >> f, 1.0 / pC["beta"][0]
    print >> f, T_high
    for i in range(len(atm_idx)):
        print >> f, atm_idx[i],
    print >> f, ""
    for i in range(len(atm_idx)):
        print >> f, p["U"][atm_idx[i] - 1],
    print >> f, ""
    for i in range(len(atm_idx)):
        print >> f, p["J"][atm_idx[i] - 1],
    print >> f, ""
    for i in range(len(atm_idx)):
        print >> f, p["alpha"][atm_idx[i] - 1],
    print >> f, ""
    f.close()
    if TB.LHF == ".FALSE.":
        f = open("ksum.input", "w")
        print >> f, p["q"][0], p["q"][1], p["q"][2]
        print >> f, p["noms"], p["noms"] + p["nomlog"]
        print >> f, p["nspin"]
        print >> f, TB.ncor_orb
        print >> f, TB.max_cor_orb
        for i in range(len(atm_idx)):
            print >> f, atm_idx[i],
        print >> f, ""
        print >> f, 1.0 / pC["beta"][0]
        #      print >> f, p['Nd_f']
        print >> f, p["n_tot"]
        print >> f, p["mu_iter"]
        print >> f, p["mix_sig"]
        for i in range(len(atm_idx)):
            print >> f, p["U"][atm_idx[i] - 1],
        print >> f, ""
        for i in range(len(atm_idx)):
            print >> f, p["alpha"][atm_idx[i] - 1],
        print >> f, ""
        for i in range(len(atm_idx)):
            print >> f, p["J"][atm_idx[i] - 1],
        print >> f, ""
        f.close()
    else:
        f = open("ksum.input", "w")
        print >> f, p["nspin"]
        print >> f, TB.ncor_orb
        print >> f, TB.max_cor_orb
        print >> f, p["n_tot"]
        print >> f, p["mu_iter"]
        for i in range(len(atm_idx)):
            print >> f, atm_idx[i],
        print >> f, ""
        for i in range(len(atm_idx)):
            print >> f, p["U"][atm_idx[i] - 1],
        print >> f, ""
        for i in range(len(atm_idx)):
            print >> f, p["alpha"][atm_idx[i] - 1],
        print >> f, ""
        for i in range(len(atm_idx)):
            print >> f, p["J"][atm_idx[i] - 1],
        print >> f, ""
        f.close()
