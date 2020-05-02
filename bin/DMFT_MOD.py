#!/usr/bin/env python2

import copy
import os
import sys

import scipy.interpolate
from scipy import *

import Fileio


def Interpolate(ommesh_input, data_input, ommesh_output, sorder):

    data_output = zeros((shape(data_input)[0], len(ommesh_output)), dtype=complex)
    data_input = array(data_input)
    for i in range(len(data_input)):
        SSpline = scipy.interpolate.splrep(
            ommesh_input, data_input[i].real, k=sorder, s=0
        )
        data_output[i, :] += scipy.interpolate.splev(ommesh_output, SSpline)
        SSpline = scipy.interpolate.splrep(
            ommesh_input, data_input[i].imag, k=sorder, s=0
        )
        data_output[i, :] += 1j * scipy.interpolate.splev(ommesh_output, SSpline)
    return data_output


def Fermi(eps):
    if eps > 100.0:
        return 0.0
    elif eps < -100.0:
        return 1.0
    else:
        return 1.0 / (1.0 + exp(eps))


def Compute_OCC(Np, Nm, d_orb):
    OCC_loc = array(
        list((Np[: len(d_orb)] + Nm[: len(d_orb)]) / 2)
        + list((Np[: len(d_orb)] - Nm[: len(d_orb)]) / 2)
    )
    return OCC_loc


def Compute_TrG1G2(om, G1, G2, T, small=1e-3):
    TrG1G2 = 0.0
    if abs(G1[-1]) < small or abs(G2[-1]) < small:
        TrG1G2 += 2 * T * sum(G1 * G2).real
    else:
        A = 1.0 / ((1.0 / G1[-1]).imag / om[-1])
        B = (1.0 / G1[-1]).real * A
        G1_tail = zeros(len(om), dtype=complex)
        G1_s = zeros(len(om), dtype=complex)
        G1_tail = A / (1j * om + B)
        G1_s = G1 - A / (1j * om + B)
        Ag = 1.0 / ((1.0 / G2[-1]).imag / om[-1])
        Bg = (1.0 / G2[-1]).real * Ag
        G2_tail = zeros(len(om), dtype=complex)
        G2_s = zeros(len(om), dtype=complex)
        G2_tail = Ag / (1j * om + Bg)
        G2_s = G2 - Ag / (1j * om + Bg)
        TrG1G2 += 2 * T * sum(G1_tail * G2_s).real
        TrG1G2 += 2 * T * sum(G1_s * G2_tail).real
        TrG1G2 += 2 * T * sum(G1_s * G2_s).real
        if abs(B - Bg) < small:
            TrG1G2 += exp(Bg / T) / (1 + exp(Bg / T)) ** 2 / (2 * T)
            TrG1G2 += exp(B / T) / (1 + exp(B / T)) ** 2 / (2 * T)
        else:
            TrG1G2 += A * Ag * (Fermi(B / T) - Fermi(Bg / T)) / (B - Bg)

    return TrG1G2


class DMFT_class:
    """
   A general class to manipulate DMFT modules
   """

    def __init__(self, p, pC, TB):

        cor_at = p["cor_at"]
        cor_orb = p["cor_orb"]
        self.cor_at = cor_at
        self.cor_orb = cor_orb
        if os.path.exists("DMFT_mu.out"):
            self.mu = float(loadtxt("DMFT_mu.out"))
        else:
            print "DMFT_mu.out file is missing"
            exit()

        # self.Ed=array([loadtxt('Ed.out')])
        # self.Natom=zeros(len(self.cor_at),dtype=int)
        self.Nd_latt = zeros(len(self.cor_at), dtype=float)
        self.VDC = zeros(len(self.cor_at), dtype=float)
        self.T = 1.0 / pC["beta"][0]

        self.N_latt = zeros((len(self.cor_at), TB.max_cor_orb), dtype=float)
        self.MOM = zeros((len(self.cor_at), TB.max_cor_orb), dtype=float)
        # for i in range(len(cor_at)): self.VDC[i]=(p['J'][i]-p['U'][i])/2.

        # self.nom=p['nom']
        # noms=p['noms']
        # nomlog=p['nomlog']

        # self.EqLogMesh(p['noms'],p['nomlog'],p['nom'],self.T)

        # for i,ats in enumerate(cor_at):
        #   self.Natom[i]=len(cor_at[i])

    def Update_Sigoo_and_Vdc(self, TB, sig_file, nspin):

        self.Sigoo = zeros((len(self.cor_at), 2 * TB.max_cor_orb), dtype=float)
        self.Vdc = zeros(len(self.cor_at), dtype=float)
        fi = open(sig_file, "r")
        self.nom, self.ncor_orb = map(int, fi.readline()[16:].split())
        fi.readline()
        fi.readline()
        Sigoo = eval(fi.readline()[8:])
        Vdc = eval(fi.readline()[7:])
        fi.close()
        idx = 0
        for i, ats in enumerate(self.cor_at):
            for j, orbs in enumerate(self.cor_orb[i]):
                if len(orbs) > 0:
                    idx += 1
        if self.ncor_orb != nspin * idx:
            print "The number of correlated orbitals in sig.inp and INPUT.py is not consistent"
            exit()

        idx2 = 0
        for i, ats in enumerate(self.cor_at):
            d_orb = TB.TB_orbs[ats[0]]
            self.Vdc[i] = Vdc[idx2]
            # self.Sigoo.append([]);#self.Vdc.append([])
            for j, orbs in enumerate(self.cor_orb[i]):
                if len(orbs) > 0:
                    for orb in orbs:
                        self.Sigoo[i, d_orb.index(orb)] = Sigoo[idx2]
                        if nspin == 1:
                            self.Sigoo[i, len(d_orb) + d_orb.index(orb)] = Sigoo[idx2]
                        else:
                            self.Sigoo[i, len(d_orb) + d_orb.index(orb)] = Sigoo[
                                idx + idx2
                            ]
                    idx2 += 1
        # print self.nom,self.ncor_orb
        # print self.Sigoo, self.Vdc

    def Compute_Energy(self, DFT, TB, ed):
        # """This module compute totoal energy using DMFT"""
        self.ETOT = 0.0
        self.ETOT2 = 0.0
        self.EPOT = 0.0
        self.EPOT2 = 0.0
        self.EDC = 0.0
        self.EKIN = 0.0
        # self.EKIN0=0.0
        DFT.Read_OSZICAR("OSZICAR")
        ETOT_imp = 0.0
        TrDeltaG = 0.0
        for i, ats in enumerate(self.cor_at):
            self.EPOT2 += len(ats) * self.TrSigG[i]
            # E_KIN2+=len(ats)*DMFT.Ekin[i]
            om, Delta = Fileio.Read_complex_multilines("imp." + str(i) + "/Delta.inp")
            om, G_loc = Fileio.Read_complex_multilines("G_loc.out")
            # om2,Gf=Fileio.Read_complex_multilines('imp.'+str(i)+'/Gf.out',1)
            # nom_s=min(len(om),len(om2))
            for j, orbs in enumerate(self.cor_orb[i]):
                TrDeltaG += (
                    2
                    * len(ats)
                    * len(orbs)
                    * Compute_TrG1G2(om, Delta[j], G_loc[j], self.T)
                )
                # TrDeltaG+=Compute_TrG1G2(om[:nom_s],Delta[j][:nom_s],Gf[j][:nom_s],self.T)
            ETOT_imp += len(ats) * (
                self.Ekin_imp[i] + self.Epot_imp[i] - self.mu_imp[i] * self.Nd_imp[i]
            )
        fiDMFT = open("INFO_KSUM", "r")
        Eline = fiDMFT.readlines()[-1:]
        self.EKIN = float(Eline[0].split()[4])
        fiDMFT.close()
        om, Sig_loc = Fileio.Read_complex_multilines("sig.inp", 5)
        om, G_loc = Fileio.Read_complex_multilines("G_loc.out")
        # print shape(Sig_loc)
        # print shape(G_loc)
        idx = -1
        for i, ats in enumerate(self.cor_at):
            d_orb = TB.TB_orbs[ats[0]]
            for j, orbs in enumerate(self.cor_orb[i]):
                idx += 1
                self.EPOT += (
                    len(ats)
                    * len(orbs)
                    * Compute_TrG1G2(om, Sig_loc[idx], G_loc[idx], self.T)
                )
                for orb in orbs:
                    # spin X2
                    self.EPOT += (
                        0.5
                        * len(ats)
                        * self.Sigoo[i][j]
                        * self.N_latt[i, d_orb.index(orb)]
                    )
        VdcNd = 0.0
        VdcNd2 = 0.0
        VdcNd3 = 0.0
        for i, ats in enumerate(self.cor_at):
            VdcNd2 += len(ats) * self.Vdc[i] * self.Nd_imp[i]  # *DMFT.Nd_imp[i]**2
            VdcNd3 += len(ats) * self.Vdc[i] * self.Nd_latt[i]
            VdcNd += len(ats) * self.mu * self.Nd_latt[i]
            d_orb = TB.TB_orbs[ats[0]]
            for j, orbs in enumerate(self.cor_orb[i]):
                for orb in orbs:
                    VdcNd -= len(ats) * ed[i][j] * self.N_latt[i, d_orb.index(orb)]
        # self.EPOT2=ETOT_imp-TrDeltaG+VdcNd+VdcNd3
        # print self.EPOT,self.EPOT2

    def Read_Sig(self, TB, nspin):
        self.Sig = []
        self.Nd_imp = zeros(len(self.cor_at), dtype=float)
        self.N_imp = zeros((len(self.cor_at), TB.max_cor_orb), dtype=float)
        self.MOM_imp = zeros((len(self.cor_at), TB.max_cor_orb), dtype=float)
        self.Epot_imp = []
        self.Ekin_imp = []
        self.mu_imp = []
        self.TrSigG = []
        for i, ats in enumerate(self.cor_at):
            d_orb = TB.TB_orbs[ats[0]]
            fileSig = "imp." + str(i) + "/Sig.out"
            if os.path.exists(fileSig):
                (
                    self.ommesh,
                    Sig_file,
                    TrS,
                    Epot,
                    nf_q,
                    mom,
                    Ekin,
                    imp_mu,
                ) = Fileio.Read_complex_Data(fileSig)
                if len(Sig_file) != nspin * len(self.cor_orb[i]):
                    print "The number of correated orbital is not same as Sig file column"
                    exit()
                if len(mom) != nspin * len(self.cor_orb[i]):
                    print "The number of correated orbital is not same as mom list in Sig file"
                    exit()
                # if self.nom>len(ommesh_long): print "nom should be decreased!"; exit()
                self.Nd_imp[i] = nf_q
                for j, orbs in enumerate(self.cor_orb[i]):
                    for orb in orbs:
                        self.N_imp[i, d_orb.index(orb)] = mom[j] / len(orbs)
                        # self.N_imp[TB.idx[at][orb]]=mom[j]/len(orbs)
                        if nspin == 2:
                            self.N_imp[i, d_orb.index(orb)] += mom[
                                j + len(self.cor_orb[i])
                            ] / len(orbs)
                            self.MOM_imp[i, d_orb.index(orb)] = (
                                mom[j] - mom[j + len(self.cor_orb[i])]
                            ) / len(orbs)
                for j in range(nspin * len(self.cor_orb[i])):
                    for iom in range(len(self.ommesh)):
                        if Sig_file[j, iom].imag > 0:
                            Sig_file[j, iom] = Sig_file[j, iom].real + 0.0j
                    self.Sig.append(Sig_file[j])

                self.TrSigG.append(TrS)
                self.Epot_imp.append(Epot)
                self.Ekin_imp.append(Ekin)
                self.mu_imp.append(imp_mu)
        self.Sig = array(self.Sig)

    def Cmp_Sig_highT(self, T_high, noms_high, nspin):
        self.ommesh_highT = zeros(noms_high)
        for i in range(noms_high):
            self.ommesh_highT[i] = pi * T_high * (2 * i + 1)
        self.Sig_highT = Interpolate(self.ommesh, self.Sig, self.ommesh_highT, 3)
        if nspin == 2:
            self.Sig_dn_highT = Interpolate(
                self.ommesh, self.Sig_dn, self.ommesh_highT, 3
            )

    def Compute_Sigoo_and_Vdc(self, p, TB):
        """Compute Energy and Sigma"""
        nspin = p["nspin"]
        U = p["U"]
        J = p["J"]
        dc_type = p["dc_type"]
        alpha = p["alpha"]
        self.Sigoo = zeros((len(self.cor_at), 2 * TB.max_cor_orb), dtype=float)
        self.Vdc = zeros(len(self.cor_at), dtype=float)
        self.Sigoo_imp = zeros((len(self.cor_at), 2 * TB.max_cor_orb), dtype=float)
        self.Vdc_imp = zeros(len(self.cor_at), dtype=float)
        # self.Eoo=0.0#zeros(len(self.cor_at),dtype=float)
        # self.Eoo_imp=0.0#zeros(len(self.cor_at),dtype=float)
        self.Edc = 0.0
        self.Edc_imp = 0.0
        # self.EHF_cor=zeros(len(self.cor_at),dtype=float)
        #      self.SigMdc=zeros(TB.ncor_orb,dtype=float)
        #      if nspin==2: self.SigMdc_dn=zeros(TB.ncor_orb,dtype=float)
        for i, ats in enumerate(self.cor_at):
            d_orb = TB.TB_orbs[ats[0]]
            fi = open("UC" + str(i + 1) + ".dat", "r")
            UC = []
            for line in fi.readlines():
                UC.append(map(float, line.split()))
            if len(UC) != 2 * len(d_orb):
                print "The size of UC is not consistent with orb"
                exit()
            UC = array(UC) + U[i] - diag(ones(2 * len(d_orb)) * U[i])
            OCC = Compute_OCC(self.N_latt[i], self.MOM[i], d_orb)
            self.Sigoo[i, : 2 * len(d_orb)] = dot(UC, OCC)
            OCC = Compute_OCC(self.N_imp[i], self.MOM_imp[i], d_orb)
            self.Sigoo_imp[i, : 2 * len(d_orb)] = dot(UC, OCC)

            ###### Compute VDC & EDC #######
            if dc_type == 1:
                self.Vdc[i] = (U[i] - alpha[i]) * (self.Nd_latt[i] - 0.5) - J[i] / 2 * (
                    self.Nd_latt[i] - 1
                )
                self.Vdc_imp[i] = (U[i] - alpha[i]) * (self.Nd_imp[i] - 0.5) - J[
                    i
                ] / 2 * (self.Nd_imp[i] - 1)
                self.Edc += len(ats) * (
                    (U[i] - alpha[i]) * self.Nd_latt[i] * (self.Nd_latt[i] - 1) / 2.0
                    - J[i] * self.Nd_latt[i] * (self.Nd_latt[i] - 2) / 4.0
                )
                self.Edc_imp += len(ats) * (
                    (U[i] - alpha[i]) * self.Nd_imp[i] * (self.Nd_imp[i] - 1) / 2.0
                    - J[i] * self.Nd_imp[i] * (self.Nd_imp[i] - 2) / 4.0
                )
            elif dc_type == 2:
                self.Vdc[i] = U[i] * (self.Nd_latt[i] - alpha[i] - 0.5) - J[i] / 2 * (
                    self.Nd_latt[i] - alpha[i] - 1
                )
                self.Vdc_imp[i] = U[i] * (self.Nd_imp[i] - alpha[i] - 0.5) - J[
                    i
                ] / 2 * (self.Nd_imp[i] - alpha[i] - 1)
                self.Edc += len(ats) * (
                    U[i]
                    * (self.Nd_latt[i] - alpha[i])
                    * (self.Nd_latt[i] - alpha[i] - 1)
                    / 2.0
                    - J[i]
                    * (self.Nd_latt[i] - alpha[i])
                    * (self.Nd_latt[i] - alpha[i] - 2)
                    / 4.0
                )
                self.Edc_imp += len(ats) * (
                    U[i]
                    * (self.Nd_imp[i] - alpha[i])
                    * (self.Nd_imp[i] - alpha[i] - 1)
                    / 2.0
                    - J[i]
                    * (self.Nd_imp[i] - alpha[i])
                    * (self.Nd_imp[i] - alpha[i] - 2)
                    / 4.0
                )
            elif dc_type == 3:
                self.Vdc[i] = U[i] * (self.Nf[i] - 0.5) - J[i] / 2 * (self.Nf[i] - 1)
                self.Vdc_imp[i] = U[i] * (self.Nf[i] - 0.5) - J[i] / 2 * (
                    self.Nf[i] - 1
                )
            else:
                print "This dc type is not supported!"
                exit()

        # self.EHF+=len(ats)*0.5*dot(OCC,self.VHF[i][:2*len(d_orb)])
        # self.EHF_imp+=len(ats)*0.5*dot(OCC_imp,self.VHF_imp[i][:2*len(d_orb)])
        # self.EHF_latt+=len(ats)*0.5*dot(OCC_latt,self.VHF_latt[i][:2*len(d_orb)])

    #         for at in ats:
    #            for orb in TB.TB_orbs[at]:
    #               idx=TB.idx[at][orb]
    #               self.SigMdc[idx] = self.VHF[i,d_orb.index(orb)]-self.VDC[i]
    #         if nspin==2:
    #            for at in ats:
    #               for orb in TB.TB_orbs[at]:
    #                  idx=TB.idx[at][orb]
    #                  self.SigMdc_dn[idx] = self.VHF[i,d_orb.index(orb)+len(d_orb)]-self.VDC[i]
    #      if (os.path.exists('SigMdc.out')):
    #         self.SigMdc_old=Fileio.Read_float('SigMdc.out')
    #      else:
    #         self.SigMdc_old=copy.deepcopy(self.SigMdc)
    #      if nspin==2:
    #         if (os.path.exists('SigMdc_dn.out')):
    #            self.SigMdc_dn_old=Fileio.Read_float('SigMdc_dn.out')
    #         else:
    #            self.SigMdc_dn_old=copy.deepcopy(self.SigMdc_dn)

    def Symmetrize_orb(self, TB, data, nspin):
        out = []
        for i, ats in enumerate(self.cor_at):
            d_orb = TB.TB_orbs[ats[0]]
            for j, orbs in enumerate(self.cor_orb[i]):
                data_loc = 0.0
                for orb in orbs:
                    data_loc += data[i, d_orb.index(orb)]
                data_loc /= len(orbs)
                out.append(data_loc)
            if nspin == 2:
                for j, orbs in enumerate(self.cor_orb[i]):
                    data_loc = 0.0
                    for orb in orbs:
                        data_loc += data[i, d_orb.index(orb) + len(d_orb)]
                    data_loc /= len(orbs)
                    out.append(data_loc)

        return out

    def Mix_Sig_and_Print_sig_inp(self, TB, Nd_qmc, mix_sig, sig_file, nspin):
        fi = open(sig_file, "r")
        self.nom, self.ncor_orb = map(int, fi.readline()[16:].split())
        fi.readline()
        fi.readline()
        Sigoo_old = eval(fi.readline()[8:])
        Vdc_old = eval(fi.readline()[7:])
        fi.close()
        om, Sig_old = Fileio.Read_complex_multilines(sig_file, 5)
        sym_Sigoo = self.Symmetrize_orb(TB, self.Sigoo, nspin)
        sym_Sigoo_imp = self.Symmetrize_orb(TB, self.Sigoo_imp, nspin)
        sym_Vdc = []
        sym_Vdc_imp = []
        for i, ats in enumerate(self.cor_at):
            for j, orbs in enumerate(self.cor_orb[i]):
                sym_Vdc.append(self.Vdc[i])
                sym_Vdc_imp.append(self.Vdc_imp[i])
            if nspin == 2:
                for j, orbs in enumerate(self.cor_orb[i]):
                    sym_Vdc.append(self.Vdc[i])
                    sym_Vdc_imp.append(self.Vdc_imp[i])

        # print self.Sigoo_imp,sym_Sigoo_imp
        # print array(Sigoo_old),sym_Sigoo_imp
        if Nd_qmc == 1:
            new_Sigoo = array(Sigoo_old) + mix_sig * (
                array(sym_Sigoo_imp) - array(Sigoo_old)
            )
            new_Vdc = array(Vdc_old) + mix_sig * (array(sym_Vdc_imp) - array(Vdc_old))
        else:
            new_Sigoo = array(Sigoo_old) + mix_sig * (
                array(sym_Sigoo) - array(Sigoo_old)
            )
            new_Vdc = array(Vdc_old) + mix_sig * (array(sym_Vdc) - array(Vdc_old))
        # print new_Sigoo,new_Vdc
        self.SigMdc = new_Sigoo - new_Vdc
        # print self.SigMdc
        # print om, self.ommesh
        for i in range(len(self.Sig)):
            self.Sig[i, :] -= sym_Sigoo_imp[i]
        # print self.Sig
        if len(om) != len(self.ommesh):
            Sig_old = Interpolate(om, Sig_old, self.ommesh, 3)
            self.nom = len(self.ommesh)
        self.Sig = Sig_old + mix_sig * (self.Sig - Sig_old)
        # print self.Sig
        # print shape(self.Sig)
        header1 = "# nom,ncor_orb= " + str(self.nom) + " " + str(self.ncor_orb)
        header2 = "# T= %18.15f" % (self.T)  # +str(self.T)
        header3 = "# s_oo-Vdc= "
        for i in range(self.ncor_orb):
            header3 += "%18.15f " % (self.SigMdc[i])
        header4 = "# s_oo= " + str(new_Sigoo.tolist())
        header5 = "# Vdc= " + str(new_Vdc.tolist())
        # print header1
        # print header2
        # print header3
        # print header4
        # print header5
        Fileio.Print_complex_multilines(
            self.Sig,
            self.ommesh,
            "sig.inp",
            [header1, header2, header3, header4, header5],
        )

    #      self.Sig=array(self.Sig_old)+mix_sig*(array(self.Sig)-array(self.Sig_old))
    #      self.SigMdc=array(self.SigMdc_old)+mix_sig*(array(self.SigMdc)-array(self.SigMdc_old))
    #      if nspin==2:
    #         self.Sig_dn=array(self.Sig_dn_old)+mix_sig*(array(self.Sig_dn)-array(self.Sig_dn_old))
    #         self.SigMdc_dn=array(self.SigMdc_dn_old)+mix_sig*(array(self.SigMdc_dn)-array(self.SigMdc_dn_old))

    def Compute_Delta(self, T, nspin, cor_at, cor_orb, TB, nom, delta=0.0):
        #####  Store local Green function and local Self energy with equidistant mesh as a list type ##########
        ommesh, tot_GLOC = Fileio.Read_complex_multilines("G_loc.out")
        ommesh, tot_SIGLOC = Fileio.Read_complex_multilines("SigMoo.out")
        SIGMDC = loadtxt("SigMdc.out")
        if nspin == 2:
            ommesh, tot_GLOC_dn = Fileio.Read_complex_multilines("G_loc_dn.out")
            ommesh, tot_SIGLOC_dn = Fileio.Read_complex_multilines("SigMoo_dn.out")
            SIGMDC_dn = loadtxt("SigMdc_dn.out")
        DMFT_mu = loadtxt("DMFT_mu.out")
        tot_ed = loadtxt("Ed.out")
        #      tot_sig_st=loadtxt('Sig_st.out')
        ed = []
        GLOC = []
        SIGLOC = []
        for i, ats in enumerate(cor_at):
            GLOC.append([])
            SIGLOC.append([])
            ed.append([])
            # sig_st.append([])
            for j, orbs in enumerate(cor_orb[i]):
                Gf_avg = zeros(len(ommesh), dtype=complex)
                Sig_avg = zeros(len(ommesh), dtype=complex)
                ed_avg = 0.0
                # sigst_avg=0.0
                for at in ats:
                    for orb in orbs:
                        idx = TB.idx[at][orb]
                        Gf_avg += tot_GLOC[idx]
                        Sig_avg += tot_SIGLOC[idx] + SIGMDC[idx]
                        ed_avg += tot_ed[idx]
                        # sigst_avg+=DMFT.Sig_st[i][j]#tot_sig_st[idx]
                Gf_avg /= len(ats) * len(orbs)
                Sig_avg /= len(ats) * len(orbs)  # ;Sig_avg-=sig_st[i]
                ed_avg /= len(ats) * len(orbs)
                # sigst_avg/=len(ats)*len(orbs)
                GLOC[i].append(list(Gf_avg))
                SIGLOC[i].append(list(Sig_avg))
                ed[i].append(ed_avg)
                # sig_st[i].append(sigst_avg)
        if nspin == 2:
            GLOC_dn = []
            SIGLOC_dn = []
            for i, ats in enumerate(cor_at):
                GLOC_dn.append([])
                SIGLOC_dn.append([])
                for j, orbs in enumerate(cor_orb[i]):
                    Gf_avg = zeros(len(ommesh), dtype=complex)
                    Sig_avg = zeros(len(ommesh), dtype=complex)
                    for at in ats:
                        for orb in orbs:
                            idx = TB.idx[at][orb]
                            Gf_avg += tot_GLOC_dn[idx]
                            Sig_avg += tot_SIGLOC_dn[idx] + SIGMDC_dn[idx]
                    Gf_avg /= len(ats) * len(orbs)
                    Sig_avg /= len(ats) * len(orbs)  # ;Sig_avg-=sig_st[i]
                    GLOC_dn[i].append(list(Gf_avg))
                    SIGLOC_dn[i].append(list(Sig_avg))
        for i in range(len(GLOC)):
            if len(GLOC[i]) > 0:
                Delta_s = zeros((nspin * len(GLOC[i]), len(ommesh)), dtype=complex)
                for j in range(len(GLOC[i])):
                    Delta_s[j, :] = (
                        1j * ommesh
                        + DMFT_mu
                        - ed[i][j]
                        - array(SIGLOC[i][j])
                        + 1j * delta
                        - 1.0 / array(GLOC[i][j])
                    )
                if nspin == 2:
                    for j in range(len(GLOC[i])):
                        Delta_s[j + len(GLOC[i]), :] = (
                            1j * ommesh
                            + DMFT_mu
                            - ed[i][j]
                            - array(SIGLOC_dn[i][j])
                            + 1j * delta
                            - 1.0 / array(GLOC_dn[i][j])
                        )
            ######  Interpolate Delta ####
            ommesh_new = pi * T * (2 * arange(nom) + 1)

            Delta = Interpolate(ommesh, Delta_s, ommesh_new, 3)
            Fileio.Print_complex_multilines(
                Delta, ommesh_new, "Delta" + str(i + 1) + ".inp"
            )
        return DMFT_mu, ed  # ,sig_st

    def mu_bcast(self, comm):
        self.mu = comm.bcast(self.mu, root=0)

    def Gather_Ksum_HF(self, comm, TOT_NKPTS, NWANN):
        rank = comm.Get_rank()
        size = comm.Get_size()
        self.EKIN = comm.gather(self.EKIN, root=0)
        if rank == 0:
            self.EKIN = sum(self.EKIN) / TOT_NKPTS * 2
        if rank == 0:
            dm = zeros((size, NWANN), dtype=float)
        else:
            dm = None
        comm.Gather(self.dm, dm, root=0)
        if rank == 0:
            self.dm = zeros(NWANN, dtype=float)
            for i in range(NWANN):
                self.dm[i] = sum(dm[:, i]) / TOT_NKPTS * 2

    def Gather_Ksum_DMFT(self, comm, TOT_NKPTS, NWANN):
        rank = comm.Get_rank()
        size = comm.Get_size()
        self.EKIN = comm.gather(self.EKIN, root=0)
        if rank == 0:
            self.EKIN = sum(self.EKIN) / TOT_NKPTS * 2
        if rank == 0:
            dm = zeros((size, NWANN), dtype=float)
        else:
            dm = None
        comm.Gather(self.dm, dm, root=0)
        if rank == 0:
            self.dm = zeros(NWANN, dtype=float)
            for i in range(NWANN):
                self.dm[i] = sum(dm[:, i]) / TOT_NKPTS * 2
        if rank == 0:
            Gloc = zeros((size, self.ncor_orb, len(self.ommesh)), dtype=complex)
        else:
            Gloc = None
        comm.Gather(self.Gloc, Gloc, root=0)
        if rank == 0:
            self.Gloc = zeros((self.ncor_orb, len(self.ommesh)), dtype=complex)
            for i in range(self.ncor_orb):
                for iom in range(len(self.ommesh)):
                    self.Gloc[i, iom] = sum(Gloc[:, i, iom]) / TOT_NKPTS

    def Gather_Ksum_DMFT2(self, comm, NFDIR, TOT_NKPTS, NWANN):
        rank = comm.Get_rank()
        size = comm.Get_size()
        if rank == 0:
            ddm = zeros((size, NFDIR, NWANN), dtype=float)
        else:
            ddm = None
        comm.Gather(self.ddm, ddm, root=0)
        if rank == 0:
            self.ddm = zeros((NFDIR, NWANN), dtype=float)
            for i in range(size):
                self.ddm[:, :] += ddm[i, :, :]
            self.ddm = self.ddm / TOT_NKPTS * 2
        if rank == 0:
            dGloc = zeros((size, NFDIR, self.ncor_orb, len(self.ommesh)), dtype=complex)
        else:
            dGloc = None
        comm.Gather(self.dGloc, dGloc, root=0)
        if rank == 0:
            self.dGloc = zeros((NFDIR, self.ncor_orb, len(self.ommesh)), dtype=complex)
            for i in range(size):
                self.dGloc[:, :, :] += dGloc[i, :, :, :]
            self.dGloc = self.dGloc / TOT_NKPTS

    def EqLogMesh(self, noms, nomlog, nom, T):
        """ This function computes the linear mesh for samll omega and the log mesh for large omega
      """
        self.ommesh = zeros(noms + nomlog, dtype=float)
        for i in range(noms):
            self.ommesh[i] = pi * T * (2 * i + 1)
        logi = log(pi * T * (2 * noms + 1))
        logf = log(pi * T * (2 * (nom - 1) + 1))
        for i in range(noms, noms + nomlog):
            self.ommesh[i] = exp(logi + (logf - logi) * (i - noms) / (nomlog - 1))

    def Update_Nlatt(self, TB, nspin=2):
        self.N_latt = zeros((len(self.cor_at), TB.max_cor_orb), dtype=float)
        self.MOM = zeros((len(self.cor_at), TB.max_cor_orb), dtype=float)
        self.Nd_latt = zeros(len(self.cor_at), dtype=float)
        line = open("INFO_DM", "r").readlines()[-1].split()
        for i, ats in enumerate(self.cor_at):
            d_orb = TB.TB_orbs[ats[0]]
            for at in ats:
                for j, orb in enumerate(d_orb):
                    idx = TB.idx[at][orb]
                    print "index: ", idx, "for atom: ", at, ", orbital: ", orb
                    self.N_latt[i, j] += float(line[idx])
                    if nspin == 2:
                        self.MOM[i, j] += float(line[idx + TB.ncor_orb])
            self.Nd_latt[i] = sum(self.N_latt[i])
            self.N_latt[i, :] /= len(ats)
            if nspin == 2:
                self.MOM[i, :] /= len(ats)
            self.Nd_latt[i] /= len(ats)

    def Compute_VDC_latt(self, Uprime, J, dc_type):
        self.VDC_old = copy.deepcopy(self.VDC)
        self.VDC = zeros(len(self.Nd_latt), dtype=float)
        for i in range(len(self.Nd_latt)):
            if dc_type == 1:
                # self.VDC[i]=(Uprime[i]-2*J[i])*(self.Nd_latt[i]-0.5)-J[i]/2*(self.Nd_latt[i]-3)
                self.VDC[i] = Uprime[i] * (self.Nd_latt[i] - 0.5) - J[i] / 2 * (
                    self.Nd_latt[i] - 1
                )
                # self.VDC[i]=Uprime[i]*(self.N_imp[i]-0.5)-J[i]/2*(self.N_imp[i]-1)
            elif dc_type == 2:
                self.VDC[i] = (Uprime[i] - 2 * J[i]) * (0.9 * self.Nd_latt[i]) - J[
                    i
                ] / 2 * (2 * self.Nd_latt[i] / 5)
            else:
                print "dc type is wrong!"
                exit()

    def Mix_DC(self, Uprime, J, self_dc, mix_dc, Nd_f):
        for i in range(len(self.Nd_latt)):
            if self_dc == True:
                self.VDC[i] = self.VDC_old[i] + mix_dc * (self.VDC[i] - self.VDC_old[i])
            else:
                VDC_f = (Uprime[i] - 2 * J[i]) * (Nd_f[i] - 0.5) - J[i] / 2 * (
                    Nd_f[i] - 3
                )
                self.VDC[i] = self.VDC_old[i] + mix_dc * (VDC_f - self.VDC[i])

    def dm_update(self, cor_at):
        fi = open("dm.out", "r")
        self.mu = float(fi.readline().split()[1])
        dc_array = fi.readline().split()
        if len(dc_array) == len(cor_at) + 1:
            for i in range(len(cor_at)):
                self.VDC[i] = float(dc_array[i + 1])
        else:
            for i in range(len(cor_at)):
                self.VDC[i] = float(dc_array[1])
        for i in range(3):
            fi.readline()
        self.dm = array([float(dmi) for dmi in fi.readline().split()])
        if len(self.VDC) != len(cor_at):
            print "Specify VDC in dm.out as many as correlated atom"
            exit()
        fi.close()

    def Mix_SigMdc(self, mix_sig):
        # self.SigMdc=self.SigMdc_old+mix_sig*(self.SigMdc-self.SigMdc_old)
        self.SigMdc = array(self.SigMdc_old) + mix_sig * (
            array(self.SigMdc) - array(self.SigMdc_old)
        )

    def Compute_FORCE(self, cor_at, cor_orb, TB, nom, NFDIR):
        """ This function computes the FORCE of DMFT using E_pot=1/2*Tr[G_loc*Sig_loc]
          The static part Estat is HF part (omega->infty) and dynamtical part is DMFT part
      """
        self.FORCE = zeros(NFDIR)
        for ifd, ifdir in enumerate(range(NFDIR)):
            self.dGLOC = []
            for i, ats in enumerate(cor_at):
                self.dGLOC.append([])
                for j, orbs in enumerate(cor_orb[i]):
                    dGf_avg = zeros(len(self.ommesh), dtype=complex)
                    for at in ats:
                        for orb in orbs:
                            idx = TB.idx[at][orb]
                            dGf_avg += self.dGloc[ifd][idx]
                    dGf_avg /= len(ats) * len(orbs)
                    self.dGLOC[i].append(list(dGf_avg))
            ######  Interpolate ####
            ommesh_new = pi * self.T * (2 * arange(nom) + 1)
            dGf_loc = []
            Sig_loc = []
            for i in range(len(self.dGLOC)):
                dGf_loc.append([])
                Sig_loc.append([])
                Gf_intpol = Interpolate(self.ommesh, self.dGLOC[i], ommesh_new, 3)
                Sig_intpol = Interpolate(self.ommesh, self.SIGLOC[i], ommesh_new, 3)
                for j in range(len(self.dGLOC[i])):
                    dGf_loc[i].append(list(Gf_intpol[j]))
                    Sig_loc[i].append(list(Sig_intpol[j]))

            d_orb = ["d_z2", "d_x2y2", "d_xz", "d_yz", "d_xy"]
            loc_idx = 0
            Fstat = 0.0
            Fdyna = 0.0
            for i, ats in enumerate(cor_at):
                for at in ats:
                    for orb in d_orb:
                        idx = TB.idx[at][orb]
                        Fstat += self.SigMdc[idx, -1].real * self.ddm[ifd, idx]
                    ####  DMFT part ################
                    # factor 4 = negative omega, spin
                for j, orbs in enumerate(cor_orb[i]):
                    Sig0 = Sig_loc[i][j][-1].real
                    Sig1 = Sig_loc[i][j][-1].imag * ommesh_new[-1]
                    Gf1 = -dGf_loc[i][j][-1].imag * ommesh_new[-1]
                    Fdyna += (
                        len(ats)
                        * len(orbs)
                        * (
                            (
                                sum(
                                    array(dGf_loc[i][j]) * (array(Sig_loc[i][j]) - Sig0)
                                ).real
                                - Gf1 * Sig1 * sum(1 / ommesh_new ** 2)
                            )
                            * self.T
                            * 4
                            + Gf1 * Sig1 / 2 / self.T
                        )
                    )
            self.FORCE[ifd] = Fstat + Fdyna
            # print Fstat,Fdyna

    def Compute_EPOT(self, cor_at, cor_orb, TB, nom):
        """ This function computes the POT energy of DMFT using E_pot=1/2*Tr[G_loc*Sig_loc]
          The static part Estat is HF part (omega->infty) and dynamtical part is DMFT part
      """
        ######  Interpolate ####
        ommesh_new = pi * self.T * (2 * arange(nom) + 1)
        Gf_loc = []
        Sig_loc = []
        for i in range(len(self.GLOC)):
            Gf_loc.append([])
            Sig_loc.append([])
            Gf_intpol = Interpolate(self.ommesh, self.GLOC[i], ommesh_new, 3)
            Sig_intpol = Interpolate(self.ommesh, self.SIGLOC[i], ommesh_new, 3)
            for j in range(len(self.GLOC[i])):
                Gf_loc[i].append(list(Gf_intpol[j]))
                Sig_loc[i].append(list(Sig_intpol[j]))
            # for j in range(len(self.GLOC[i])):
            #   Gf_intpol=zeros(nom,dtype=complex)
            #   Sig_intpol=zeros(nom,dtype=complex)
            #   GfSpline = scipy.interpolate.splrep(ommesh, array(Gf_loc_small[i][j]).real, k=3, s=0)
            #   Gf_intpol += scipy.interpolate.splev(ommesh_new, GfSpline)
            #   GfSpline = scipy.interpolate.splrep(ommesh, array(Gf_loc_small[i][j]).imag, k=3, s=0)
            #   Gf_intpol += 1j*scipy.interpolate.splev(ommesh_new, GfSpline)
            #   Gf_loc[i].append(list(Gf_intpol))
            #   SigSpline = scipy.interpolate.splrep(ommesh, array(Sig_loc_small[i][j]).real, k=3, s=0)
            #   Sig_intpol += scipy.interpolate.splev(ommesh_new, SigSpline)
            #   SigSpline = scipy.interpolate.splrep(ommesh, array(Sig_loc_small[i][j]).imag, k=3, s=0)
            #   Sig_intpol += 1j*scipy.interpolate.splev(ommesh_new, SigSpline)
            #   Sig_loc[i].append(list(Sig_intpol))

        self.EPOT = 0.0
        Estat = 0.0
        Estat2 = 0.0
        Edyna = 0.0
        self.EPOT2 = 0.0
        Edyna2 = 0.0

        d_orb = ["d_z2", "d_x2y2", "d_xz", "d_yz", "d_xy"]
        loc_idx = 0
        for i, ats in enumerate(cor_at):
            Estat += len(ats) * self.EHF[i]
            Estat2 += len(ats) * self.EHF_cor[i]
            # for at in ats:
            # for orb in d_orb:
            # idx=TB.idx[at][orb]
            # Estat+=0.5*(self.SigMdc[idx,-1].real+self.VDC[i])*self.dm[idx]
            ####  DMFT part ################
            # factor 4 = negative omega, spin
            for j, orbs in enumerate(cor_orb[i]):
                Sig0 = Sig_loc[i][j][-1].real
                Sig1 = Sig_loc[i][j][-1].imag * ommesh_new[-1]
                Edyna += (
                    len(ats)
                    * len(orbs)
                    * 0.5
                    * (
                        (
                            sum(
                                array(Gf_loc[i][j]) * (array(Sig_loc[i][j]) - Sig0)
                            ).real
                            - Sig1 * sum(1 / ommesh_new ** 2)
                        )
                        * self.T
                        * 4
                        + Sig1 / 2 / self.T
                    )
                )
            if len(cor_orb[i]) > 0:
                Edyna2 += len(ats) * self.Eimp[i]
        self.EPOT = Estat + Edyna
        self.EPOT2 = Edyna2
        # self.EPOT2=Estat-Estat2+Edyna2

    def Print_Gloc(self, print_at, TB):
        for at in print_at:
            Fileio.Print_complex_multilines(
                array(
                    [
                        self.Gloc[i]
                        for i in range(
                            TB.idx[at][TB.TB_orbs[at][0]],
                            TB.idx[at][TB.TB_orbs[at][-1]] + 1,
                        )
                    ]
                ),
                self.ommesh,
                "G_loc_" + at + ".out",
            )


if __name__ == "__main__":

    import Struct, VASP

    execfile("INPUT.py")  # Read input file
    TB = Struct.TBstructure("POSCAR", p["atomnames"], p["orbs"])
    cor_at = p["cor_at"]
    cor_orb = p["cor_orb"]
    TB.Compute_cor_idx(cor_at, cor_orb)
    print TB.TB_orbs
    DFT = VASP.VASP_class()
    DMFT = DMFT_class(p, pC, TB)
    DMFT.Update_Sigoo_and_Vdc(TB, "sig.inp")
    DMFT.Update_Nlatt(TB, p["nspin"])
    DMFT.Read_Sig(TB, p["nspin"])
    # print DMFT.Sig
    # DMFT.Compute_Sigoo_and_Vdc(Nd_qmc,p,TB)
    ed = array([loadtxt("Ed.out")])
    print ed
    DMFT.Compute_Energy(DFT, TB, ed)
    DMFT.Compute_Sigoo_and_Vdc(p, TB)
    DMFT.Mix_Sig_and_Print_sig_inp(TB, p["Nd_qmc"], p["mix_sig"], "sig.inp")
