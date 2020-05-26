#!/usr/bin/env python2
import os
import re
import subprocess

from scipy import *


class VASP_class:
    def __init__(self, dft="vasp", structurename=None, aiida=False):
        # pass
        if not os.path.exists("DFT_mu.out"):
            print "DFT_mu.out should exist!"
            exit()
        self.EFERMI = float(loadtxt("DFT_mu.out"))
        self.NBANDS = 0
        self.dft = dft
        self.structurename = structurename
        self.aiida = aiida

    def Print_CHGCAR(self, headlines, CHG, augCHG, augIDX):
        fi = open("CHGCAR", "w")
        for line in headlines:
            print >> fi, line,
        Ngrid = len(CHG)
        for i in range(Ngrid / 5):
            print >> fi, "",
            for j in range(i * 5, (i + 1) * 5):
                if CHG[j] >= 0:
                    expon = int(log10(CHG[j])) + 1
                    print >> fi, "%.11fE%+.2d" % (CHG[j] / 10 ** expon, expon),
                else:
                    expon = int(log10(abs(CHG[j]))) + 1
                    num = "%.11f" % (abs(CHG[j]) / 10 ** expon)
                    print >> fi, "-" + num[1:] + "E%+.2d" % expon,
            print >> fi, ""
        if Ngrid % 5 > 0:
            print >> fi, "",
            for i in range((Ngrid / 5) * 5, Ngrid):
                if CHG[i] >= 0:
                    expon = int(log10(CHG[i])) + 1
                    print >> fi, "%.11fE%+.2d" % (CHG[i] / 10 ** expon, expon),
                else:
                    expon = int(log10(abs(CHG[i]))) + 1
                    num = "%.11f" % (abs(CHG[i]) / 10 ** expon)
                    print >> fi, "-" + num[1:] + "E%+.2d" % expon,
            print >> fi, ""
        for i, idx in enumerate(augIDX):
            print >> fi, "augmentation occupancies %3d %3d" % (i + 1, idx)
            for j in range(idx / 5):
                print >> fi, "",
                for k in range(j * 5, (j + 1) * 5):
                    if augCHG[i][k] >= 0:
                        print >> fi, " %.7E" % augCHG[i][k],
                    else:
                        print >> fi, "-%.7E" % abs(augCHG[i][k]),
                print >> fi, ""
            if idx % 5 > 0:
                print >> fi, "",
                for j in range((idx / 5) * 5, idx):
                    if augCHG[i][j] >= 0:
                        print >> fi, " %.7E" % augCHG[i][j],
                    else:
                        print >> fi, "-%.7E" % abs(augCHG[i][j]),
                print >> fi, ""
        fi.close()

    def Read_CHGCAR(self, finame="CHGCAR"):
        fi = open(finame, "r")
        headlines = []
        for i in range(6):
            headlines.append(fi.readline())
        line = fi.readline()
        headlines.append(line)
        num_atom = [int(na) for na in line.split()]
        headlines.append(fi.readline())
        for i in range(sum(num_atom)):
            headlines.append(fi.readline())
        headlines.append(fi.readline())
        line = fi.readline()
        headlines.append(line)
        grid = [int(na) for na in line.split()]
        Ngrid = grid[0] * grid[1] * grid[2]
        CHG = zeros(Ngrid)
        for i in range(Ngrid / 5):
            CHG[i * 5 : (i + 1) * 5] = [float(na) for na in fi.readline().split()]
        if Ngrid % 5 > 0:
            CHG[(Ngrid / 5) * 5 :] = [float(na) for na in fi.readline().split()]
        CHG = array(CHG)
        CHGSUM = sum(CHG) / Ngrid
        print "total charge=", sum(CHG) / Ngrid
        CHG = CHG * self.NELECT / CHGSUM
        print "Normalized charge=", sum(CHG) / Ngrid
        augCHG = []
        augIDX = []
        while True:
            line = fi.readline()
            if len(line.split()) > 0 and line.split()[0] == "augmentation":
                num_proj = int(line.split()[3])
                augIDX.append(num_proj)
                augment = zeros(num_proj)
                for i in range(num_proj / 5):
                    augment[i * 5 : (i + 1) * 5] = [
                        float(na) for na in fi.readline().split()
                    ]
                if num_proj % 5 > 0:
                    augment[(num_proj / 5) * 5 :] = [
                        float(na) for na in fi.readline().split()
                    ]
                augCHG.append(augment)

            # AC=[sum(na) for na in augCHG]
            # print 'AUG charge=',sum(AC)
            else:
                fi.close()
                break
        return headlines, CHG, augCHG, augIDX

    def Diff_CHGCAR(self, fi1name="CHGCAR.0", fi2name="CHGCAR.1"):
        fi = open(fi1name, "r")
        for i in range(6):
            fi.readline()
        num_atom = [int(na) for na in fi.readline().split()]
        fi.readline()
        for i in range(sum(num_atom)):
            fi.readline()
        fi.readline()
        grid = [int(na) for na in fi.readline().split()]
        Ngrid1 = grid[0] * grid[1] * grid[2]
        CHG1 = zeros(Ngrid1)
        for i in range(Ngrid1 / 5):
            CHG1[i * 5 : (i + 1) * 5] = [float(na) for na in fi.readline().split()]
        if Ngrid1 % 5 > 0:
            CHG1[(Ngrid1 / 5) * 5 :] = [float(na) for na in fi.readline().split()]
        CHG1 = array(CHG1)
        # print sum(CHG)/Ngrid
        CHGSUM = sum(CHG1) / Ngrid1
        CHG1 = CHG1 * self.NELECT / CHGSUM
        augCHG1 = []
        for line in fi.readlines():
            if line.split()[0] != "augmentation":
                augCHG1.append([float(na) for na in line.split()])
        AC = [sum(na) for na in augCHG1]
        # print sum(AC)
        fi.close()

        fi = open(fi2name, "r")
        for i in range(6):
            fi.readline()
        num_atom = [int(na) for na in fi.readline().split()]
        fi.readline()
        for i in range(sum(num_atom)):
            fi.readline()
        fi.readline()
        grid = [int(na) for na in fi.readline().split()]
        Ngrid2 = grid[0] * grid[1] * grid[2]
        CHG2 = zeros(Ngrid2)
        for i in range(Ngrid2 / 5):
            CHG2[i * 5 : (i + 1) * 5] = [float(na) for na in fi.readline().split()]
        if Ngrid1 % 5 > 0:
            CHG2[(Ngrid1 / 5) * 5 :] = [float(na) for na in fi.readline().split()]
        CHG2 = array(CHG2)
        # print sum(CHG)/Ngrid
        CHGSUM = sum(CHG2) / Ngrid2
        CHG2 = CHG2 * self.NELECT / CHGSUM
        augCHG2 = []
        for line in fi.readlines():
            if line.split()[0] != "augmentation":
                augCHG2.append([float(na) for na in line.split()])
        AC = [sum(na) for na in augCHG2]
        # print sum(AC)
        fi.close()

        if Ngrid1 != Ngrid2:
            print "Two CHGCAR are not game grid!"
            exit()
        # print sqrt(sum((CHG1-CHG2)**2)/Ngrid1)
        diff_chg = sqrt(sum((CHG1 - CHG2) ** 2) / Ngrid1)
        diff_aug = 0.0
        for i, a1 in enumerate(augCHG1):
            diff_aug += sum((array(a1) - array(augCHG2[i])) ** 2)
        sqrtdiff_aug = sqrt(diff_aug)
        print ("Augmentation difference = %f " % sqrtdiff_aug)
        return diff_chg

    def Read_OSZICAR(self, finame="OSZICAR"):
        # Reading the Total energy from DFT calculations.

        if self.dft == "vasp":
            fi = open(finame, "r")
            self.E = float(fi.readlines()[-1].split()[2])
            fi.close()
            # print ("Total energy read from VASP = %f eV" % self.E)

        elif self.dft == "siesta":
            fi = open(self.structurename + ".out", "r")
            data = fi.read()
            fi.close()
            self.E = float(re.findall(r"Total\s=[\s0-9+-.]*", data)[0].split()[-1])
            # print ("Total energy read from Siesta = %f eV" % self.E)

        elif self.dft == "qe":
            if self.aiida:
                fi = open("aiida.out", "r")
            else:
                fi = open(self.structurename + ".scf.out", "r")
            data = fi.read()
            fi.close()
            total_energy = float(
                re.findall(r"total\senergy\s*=\s*([-\d.\s]*)Ry", data)[-1]
            )
            # in QE this is in Ry. Converting to eV.
            self.E = 13.60569 * total_energy

        print ("Total energy = %f eV" % self.E)

    def Read_NBANDS(self):
        if self.dft == "vasp":
            cmd = "grep NBANDS= OUTCAR | sed 's/.* //g'"
            out, err = subprocess.Popen(
                cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            ).communicate()
            self.NBANDS = int(out[:-1])  # Ignoring \n

        elif self.dft == "siesta":
            # Reading the number of bands from siesta output
            fi = open(self.structurename + ".out", "r")
            data = fi.read()
            fi.close()
            self.NBANDS = int(
                re.findall(r"Siesta2Wannier90.NumberOfBands[\s0-9]*", data)[0].split()[
                    -1
                ]
            )
        elif self.dft == "qe" and self.aiida == False:
            fi = open(self.structurename + ".scf.out", "r")
            data = fi.read()
            fi.close()
            self.NBANDS = int(
                re.findall(r"\n\s*number\s*of\s*Kohn-Sham\s*states=([\s\d]*)", data)[0]
            )

        print ("Number of bands = %d " % self.NBANDS)

    def Read_NELECT(self):
        cmd = "grep NELECT OUTCAR"
        out, err = subprocess.Popen(
            cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        ).communicate()
        self.NELECT = float(out.split()[2])

    def Read_EFERMI(self):
        if self.dft == "vasp":
            fi = open("OUTCAR", "r")
            for line in fi:
                if re.search("Fermi energy", line) or re.search("E-fermi", line):
                    line_fermi = line
            # print line_fermi
            val = re.search(r"(\-?\d+\.?\d*)", line_fermi)
            # print val
            self.EFERMI = float(val.group(1))
            savetxt("DFT_mu.out", array([self.EFERMI]))
            print ("Fermi energy = %f eV" % self.EFERMI)

        elif self.dft == "siesta":
            # Reading the Fermi energy from siesta output
            fi = open(self.structurename + ".out", "r")
            data = fi.read()
            fi.close()
            self.EFERMI = float(re.findall(r"Fermi\s=[\s0-9+-.]*", data)[0].split()[-1])
            savetxt("DFT_mu.out", array([self.EFERMI]))
            print ("Fermi energy = %f eV" % self.EFERMI)

        elif self.dft == "qe":
            # Reading the Fermi energy from QE outputs.
            if self.aiida:
                fi = open("aiida.out", "r")
            else:
                fi = open(self.structurename + ".scf.out", "r")
            data = fi.read()
            fi.close()
            self.EFERMI = float(
                re.findall(r"the\s*Fermi\s*energy\s*is\s*([\s\d.]*)ev", data)[0]
            )
            print ("Fermi energy = %f eV" % self.EFERMI)

            savetxt("DFT_mu.out", array([self.EFERMI]))

    def Read_EBAND(self):
        cmd = "grep EBANDS OUTCAR | tail -n 1"
        out, err = subprocess.Popen(
            cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        ).communicate()
        self.EBAND = float(out.split()[3])

    def Create_win(
        self,
        TB,
        atomnames,
        orbs,
        L_rot,
        nb=60,
        emin=0.0,
        emax=10.0,
        kmesh_tol=0.0001,
        num_iter=500,
    ):
        if len(L_rot) != len(atomnames):
            print "Check INPUT.py for atomnames and L_rot!"
            exit()
        if len(orbs) != len(atomnames):
            print "Check INPUT.py for atomnames and L_rot!"
            exit()
        nwann = 0
        for i, orb in enumerate(orbs):
            nwann += (2 * TB.L[orb] + 1) * TB.num_atoms[atomnames[i]]
        fi = open("wannier90.win", "w")
        if nb == 0:
            pass
        else:
            print >> fi, "num_bands   =  " + str(nb)
        print >> fi, "dis_win_min  = " + str(emin)
        print >> fi, "dis_win_max  = " + str(emax)
        print >> fi, "num_wann   = " + str(nwann)
        print >> fi, "kmesh_tol   = " + str(kmesh_tol)
        print >> fi, "num_iter   = " + str(num_iter)
        print >> fi, ""
        print >> fi, "begin projections"
        for i, at in enumerate(atomnames):
            if L_rot[i] == 1:
                for j in range(TB.num_atoms[at]):
                    cordi, rot_vec = TB.Rot_axis(at + str(j + 1))
                    print >> fi, "f=" + cordi + ":l=" + str(
                        TB.L[orbs[i]]
                    ) + ":" + rot_vec
            elif L_rot[i] == 0:
                print >> fi, at + ":" + str(orbs[i])
            else:
                print "L_rot is wrong"
                exit()
        print >> fi, "end projections"
        fi.close()

    def Update_win(self, nb=60, emin=0.0, emax=10.0):
        fi = open("wannier90.win", "r")
        lines = fi.readlines()
        fi.close()
        fi = open("wannier90.win", "w")
        print >> fi, "num_bands   =  " + str(nb)
        print >> fi, "dis_win_min  = " + str(emin)
        print >> fi, "dis_win_max  = " + str(emax)
        for line in lines:
            lsplit = line.split()
            if len(lsplit) > 0:
                if line.split()[0] == "num_bands":
                    pass
                elif line.split()[0] == "dis_win_min":
                    pass
                elif line.split()[0] == "dis_win_max":
                    pass
                else:
                    print >> fi, line,
            else:
                print >> fi, line,
        fi.close()


if __name__ == "__main__":

    mix_CHG = 0.3
    VASP = VASP_class()
    # VASP.Read_OSZICAR()
    # print VASP.E
    # VASP.Read_NBANDS()
    # print VASP.NBANDS
    # VASP.Read_EFERMI()
    # print VASP.EFERMI
    # VASP.Update_win()
    # VASP.Read_EBAND()
    # print VASP.EBAND
    VASP.Read_NELECT()
    print VASP.NELECT
    headlines, CHG_old, augCHG_old, augIDX_old = VASP.Read_CHGCAR("CHGCAR.0")
    headlines, CHG, augCHG, augIDX = VASP.Read_CHGCAR("CHGCAR.1")
    CHG_new = CHG_old + mix_CHG * (CHG - CHG_old)
    augCHG_new = []
    if augIDX == augIDX_old:
        for idx in range(len(augIDX)):
            augCHG_new.append(
                augCHG_old[idx] + mix_CHG * (augCHG[idx] - augCHG_old[idx])
            )
    else:
        print "Cannot mix CHGCAR!"
        exit()
    ##print augCHG_new
    VASP.Print_CHGCAR(headlines, CHG_new, augCHG_new, augIDX)
    ##for line in headlines:
    ##   print line,
    print VASP.Diff_CHGCAR("CHGCAR.0", "CHGCAR.1")
    # print VASP.Diff_CHGCAR('CHGCAR.0','CHGCAR')
