#!/usr/bin/env python2
import argparse
import copy
import glob
import itertools
import math
import os
import re
import shutil
import math
import subprocess
import sys
from argparse import RawTextHelpFormatter

import matplotlib

matplotlib.use("pdf")
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate
from matplotlib.font_manager import FontProperties, fontManager
from pylab import *
from scipy import *
from scipy import interpolate

import Fileio

# import oreo
# import Re_wt
import Struct
from INPUT import *
from splash import welcome


# Setting up plotting class
plt.rcParams["mathtext.default"] = "regular"  # Roman ['rm', 'cal', 'it', 'tt', 'sf',
#                                                   'bf', 'default', 'bb', 'frak',
#                                                   'circled', 'scr', 'regular']
plt.rcParams["font.family"] = "Georgia"
plt.rc("font", size=22)  # controls default text sizes
plt.rc("axes", titlesize=22)  # fontsize of the axes title
plt.rc("axes", labelsize=22)  # fontsize of the x and y labels
plt.rc("xtick", labelsize=22)  # fontsize of the tick labels
plt.rc("ytick", labelsize=22)  # fontsize of the tick labels
# plt.rc('legend', fontsize=22)    # legend fontsize
# plt.rc('figure', titlesize=22)  # fontsize of the figure title


class PostProcess:
    """DMFTwDFT post processing tool.

	This class contains methods to perform post processing of the DMFT calculations.
	Run inside the DMFT or HF directories.

        Usage:
	postDMFT.py <options>

	-h for help.

	<options>:
	ac    : Performs analytic continuation.
	dos   : Performs density of states calculation.
	bands : Performs band structure calculation.


	"""

    def __init__(self):
        """
        Initializes the following:
        """
        # mpirun
        if os.path.exists("para_com.dat"):
            fipa = open("para_com.dat", "r")
            self.para_com = str(fipa.readline())[:-1]
            fipa.close()
        else:
            self.para_com = ""

        # dmft_bin
        self.path_bin = p["path_bin"]

    def checksig(self):
        """
		Checks if Sig.out is created in the ac directory.
		"""

        if os.path.exists("./ac/Sig.out"):
            return True
        else:
            print("Analytic Continuation incomplete. Sig.out not found.")
            return False

    def interpol(self, emin, emax, rom, broaden, dest_dir, sp=False):
        """
	This performs the interpolation of points on the real axis.
	"""
        print("\nInterpolating points on real axis...")
        headerline = 2
        om, Sig = Fileio.Read_complex_multilines("./ac/Sig.out", headerline)
        s_oo = None
        Vdc = None
        # The exec() function doesn't work properly on Python3 so I had to use a workaround:
        fi = open("./ac/Sig.out", "r")
        line1 = fi.readline()
        s_oo = re.findall(r"\s*([0-9.+-]*)", line1)
        while "" in s_oo:
            s_oo.remove("")
        line2 = fi.readline()
        Vdc = re.findall(r"\s*([0-9.+-]*)", line2)
        while "" in Vdc:
            Vdc.remove("")

            # exec(ar[0])
            # m=re.search('#(.*)',line)
            # exec(m.group(1).strip())
        # s_oo_Vdc=np.array(s_oo)-array(Vdc)
        fi.close()
        s_oo_Vdc = np.array((np.array(s_oo)).astype(np.float)) - np.array(
            (np.array(Vdc)).astype(np.float)
        )

        ommesh = np.linspace(emin, emax, rom)

        # non spin polarized case
        if sp == False:

            Sig_tot = np.zeros((len(Sig), rom), dtype=complex)

            for i in range(len(Sig)):
                SigSpline = interpolate.splrep(om, Sig[i].real, k=1, s=0)
                Sig_tot[i, :] += interpolate.splev(ommesh, SigSpline)
                SigSpline = interpolate.splrep(om, Sig[i].imag, k=1, s=0)
                Sig_tot[i, :] += 1j * interpolate.splev(ommesh, SigSpline)

            header1 = "# nom,ncor_orb= " + str(len(ommesh)) + " " + str(len(Sig_tot))
            # header2='# T= %18.15f'%(1.0/pC['beta'][0])#+str(self.T)
            header2 = "# T= %18.15f" % (broaden)  # +str(self.T)
            header3 = "# s_oo-Vdc= "
            for i in range(len(s_oo_Vdc)):
                header3 += "%18.15f " % (s_oo_Vdc[i])
            header4 = "# s_oo= " + str(s_oo)
            header5 = "# Vdc= " + str(Vdc)
            if dest_dir == "dos":
                Fileio.Print_complex_multilines(
                    Sig_tot,
                    ommesh,
                    "./dos/sig.inp_real",
                    [header1, header2, header3, header4, header5],
                )

            # create sig.inp_real
            if dest_dir == "bands":
                Fileio.Print_complex_multilines(
                    Sig_tot,
                    ommesh,
                    "./bands/sig.inp_real",
                    [header1, header2, header3, header4, header5],
                )

        # Spin polarized calculation
        if sp:

            Sig_tot = np.zeros((int(len(Sig) / 2), rom), dtype=complex)
            Sig_tot_dn = np.zeros((int(len(Sig) / 2), rom), dtype=complex)

            for i in range(int(len(Sig) / 2)):
                # spin
                SigSpline = interpolate.splrep(om, Sig[i].real, k=1, s=0)
                Sig_tot[i, :] += interpolate.splev(ommesh, SigSpline)
                SigSpline = interpolate.splrep(om, Sig[i].imag, k=1, s=0)
                Sig_tot[i, :] += 1j * interpolate.splev(ommesh, SigSpline)

                # spin down
                SigSpline = interpolate.splrep(
                    om, Sig[i + int(len(Sig) / 2)].real, k=1, s=0
                )
                Sig_tot_dn[i, :] += interpolate.splev(ommesh, SigSpline)
                SigSpline = interpolate.splrep(
                    om, Sig[i + int(len(Sig) / 2)].imag, k=1, s=0
                )
                Sig_tot_dn[i, :] += 1j * interpolate.splev(ommesh, SigSpline)

            header1 = "# nom,ncor_orb= " + str(len(ommesh)) + " " + str(len(Sig_tot))
            header2 = "# T= %18.15f" % (broaden)
            header3 = "# s_oo-Vdc= "
            header3_dn = "# s_oo-Vdc= "
            for i in range(int(len(s_oo_Vdc) / 2)):
                header3 += "%18.15f " % (s_oo_Vdc[i])
                header3_dn += "%18.15f " % (s_oo_Vdc[i + int(len(s_oo_Vdc) / 2)])
            header4 = "# s_oo= " + str(s_oo[0 : int(len(s_oo) / 2)])
            header4_dn = "# s_oo= " + str(s_oo[int(len(s_oo) / 2) :])
            header5 = "# Vdc= " + str(Vdc[0 : int(len(Vdc) / 2)])
            header5_dn = "# Vdc= " + str(Vdc[int(len(Vdc) / 2) :])

            if dest_dir == "dos":
                Fileio.Print_complex_multilines(
                    Sig_tot,
                    ommesh,
                    "./dos/sig.inp_real",
                    [header1, header2, header3, header4, header5],
                )
                Fileio.Print_complex_multilines(
                    Sig_tot_dn,
                    ommesh,
                    "./dos/sig.inp_real_dn",
                    [header1, header2, header3_dn, header4_dn, header5_dn],
                )

            # create sig.inp_real
            if dest_dir == "bands":
                Fileio.Print_complex_multilines(
                    Sig_tot,
                    ommesh,
                    "./bands/sig.inp_real",
                    [header1, header2, header3, header4, header5],
                )
                Fileio.Print_complex_multilines(
                    Sig_tot_dn,
                    ommesh,
                    "./bands/sig.inp_real_dn",
                    [header1, header2, header3_dn, header4_dn, header5_dn],
                )

        print("Interpolation complete.\n")

    def anal_cont(self, args):
        """
		This method performs the analytic continuation.
		"""

        siglistindx = args.siglistindx

        # creating directory for ac
        if os.path.exists("ac"):
            shutil.rmtree("ac")
            os.makedirs("ac")
        else:
            os.makedirs("ac")

        # copying the last few self-energies from the DMFT run in the directory above
        siglist = sorted(
            sorted(glob.glob("sig.inp.*"))[1:], key=lambda x: (len(x), float(x[8:]))
        )[-siglistindx:]
        for file in siglist:
            shutil.copy(file, "ac")

        # averaging self energies
        print("\nAveraging self-energies from: ")
        print(siglist)
        cmd = "cd ac && sigaver.py sig.inp.*"
        out, err = subprocess.Popen(
            cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        ).communicate()
        if err:
            print(err)
            print("Averaging self-energies Failed!\n")
            sys.exit()
        else:
            print("Self-energies averaged.\n")

        # copy maxent_params.dat from source if not in DMFT directory
        if os.path.exists("maxent_params.dat"):
            shutil.copyfile("maxent_params.dat", "./ac/maxent_params.dat")
        else:
            src = self.path_bin + os.sep + "maxent_params.dat"
            src_path = os.path.abspath(os.path.expanduser(os.path.expandvars(src)))
            shutil.copyfile(src_path, "./ac/maxent_params.dat")

        # Analytic continuation
        print("Running analytic continuation...")
        cmd = "cd ac && maxent_run.py sig.inpx >ac.out 2>ac.error"
        out, err = subprocess.Popen(
            cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        ).communicate()

        if os.path.exists("./ac/Sig.out"):
            print("Analytic continuation complete.\n")
        else:
            print("Analytic continuation failed! Check ac.error for details.\n")
            sys.exit()

    def Create_kpath(self, KPoints, nk_band):
        print("\nGenerating k-path...")

        def dist(a, b):
            return np.lib.scimath.sqrt(np.sum((np.array(a) - np.array(b)) ** 2))

        # returns the distance of given a and b points
        KPoints = np.array(KPoints)
        path_len = []
        for i in range(len(KPoints) - 1):
            path_len.append(dist(KPoints[i + 1], KPoints[i]))
        path_nk = list(map(int, nk_band * np.array(path_len) / np.sum(path_len)))
        klist = []
        dist_K = [0.0]
        dist_SK = [0.0]
        for i, nkk in enumerate(path_nk):
            for n in range(nkk):
                klist.append(KPoints[i] + (KPoints[i + 1] - KPoints[i]) * n / nkk)
                if len(klist) > 1:
                    dist_K.append(dist_K[-1] + dist(klist[-1], klist[-2]))
            dist_SK.append(dist_SK[-1] + path_len[i])

        # Add the ending point
        klist.append(KPoints[-1])
        dist_K.append(dist_K[-1] + dist(klist[-1], klist[-2]))
        return np.array(klist), np.array(dist_K), np.array(dist_SK)

    def genksum(self, rom, kpband):
        fp = open("./bands/ksum.input", "w")
        fp.write("%d" % kpband)
        fp.write("\n%d" % rom)
        fp.write("\n%d" % p["nspin"])
        fp.write("\n%d" % 5)
        fp.write("\n%d" % 5)
        fp.write("\n%f" % 0.01)
        fp.write("\n%d" % p["n_tot"])
        fp.write("\n%d" % p["mu_iter"])
        fp.close()

    def Make_coor_list(self, ord_ST, nord_ST):
        idx = []
        for ST in ord_ST:
            for i, nord in enumerate(nord_ST):
                for cmps in nord:
                    if ST == cmps:
                        idx.append(i)
                        break
        return idx

    def dos(self, args):
        """
		This method performs the Density of States calculation.
		"""

        dest_dir = "dos"

        # creating directory for dos
        if os.path.exists("dos"):
            print("dos directory already exists.")
        else:
            os.makedirs("dos")

        # copy Sig.out into /dos
        if self.checksig():
            shutil.copyfile("./ac/Sig.out", "./dos/Sig.out")
        else:
            sys.exit()

        # interpolating
        self.interpol(args.emin, args.emax, args.rom, args.broaden, dest_dir, args.sp)

        # copying files from DMFT directory to dos directory
        cmd = "cd dos && Copy_input.py ../ -post dos"
        out, err = subprocess.Popen(cmd, shell=True).communicate()
        if err:
            print("File copy failed!\n")
            print(err)
            sys.exit()
        else:
            print(out)

        if args.sp == False:

            # running dmft_dos.x
            print("Calculating DMFT DOS...")
            cmd = "cd dos && " + self.para_com + " " + "dmft_dos.x"
            out, err = subprocess.Popen(
                cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            ).communicate()
            if err:
                print(err)
                print("DMFT DOS calculation failed!\n")
                sys.exit()
            else:
                print("DMFT DOS calculation complete.\n")

        else:
            # Generating files to plot DOS
            # read dmft_params.dat and set spin = 1
            file = open("./dos/dmft_params.dat", "r")
            data = file.readlines()
            file.close()
            data[9] = "1\n"
            file = open("./dos/dmft_params.dat", "w")
            for ele in data:
                file.write(ele)
            file.close()

            # running dmft_dos.x for spin up
            print("Calculating DMFT DOS for spin up...")
            cmd = "cd dos && " + self.para_com + " " + "dmft_dos.x"
            out, err = subprocess.Popen(
                cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            ).communicate()
            if err:
                print(err)
                print("DMFT DOS calculation failed!\n")
                sys.exit()
            else:
                print("DMFT DOS calculation complete.\n")

            # copying G_loc.out to G_loc_up.out
            shutil.copyfile("./dos/G_loc.out", "./dos/G_loc_up.out")

            # copying sig.inp_real_dn to sig.inp_real
            shutil.copyfile("./dos/sig.inp_real_dn", "./dos/sig.inp_real")

            # running dmft_dos.x for spin down
            print("Calculating DMFT DOS for spin down...")
            cmd = "cd dos && " + self.para_com + " " + "dmft_dos.x"
            out, err = subprocess.Popen(
                cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            ).communicate()
            if err:
                print(err)
                print("DMFT DOS calculation failed!\n")
                sys.exit()
            else:
                print("DMFT DOS calculation complete.\n")

            # copying G_loc.out to G_loc_dn.out
            shutil.copyfile("./dos/G_loc.out", "./dos/G_loc_dn.out")

        # calling dos plotting
        self.plot_dos(args)

    def plot_dos(self, args):
        """
	This method plots the density of states plot.
	"""

        # mapping d and p orbitals with atoms
        d_indices = [i for i, x in enumerate(p["orbs"]) if x == "d"]
        d_atoms = [p["atomnames"][i] for i in d_indices]

        p_indices = [i for i, x in enumerate(p["orbs"]) if x == "p"]
        p_atoms = [p["atomnames"][i] for i in p_indices]

        # Finding the number of columns
        fi = open("./dos/G_loc.out", "r")
        firstline = fi.readline()
        fi.close()
        num_columns = len(firstline.split())

        # index arrays for eg and t2g
        eg_index = []
        t2g_index = []
        icounter_eg = 2
        icounter_t2g = 4

        # if cor_at has more than one atom,
        # it will sum up the total contributions
        for i0 in range(len(p["cor_at"])):
            for i in range(len(p["cor_at"][0])):
                # egj
                eg_index.append(icounter_eg)
                eg_index.append(icounter_eg + 6)
                icounter_eg = icounter_eg + 10

                # t2g
                t2g_index.append(icounter_t2g)
                t2g_index.append(icounter_t2g + 2)
                t2g_index.append(icounter_t2g + 6)
                icounter_t2g = icounter_t2g + 10

        # index for p
        # initialized by last index of t2g index array
        if p_atoms:
            p_index = []
            icounter_p = t2g_index[-1] + 2
            for i in range(icounter_p, num_columns, 2):
                p_index.append(i)

        if args.sp == False:
            print("Plotting DOS...")
            with open("./dos/G_loc.out", "r") as f:
                lines = f.readlines()
                x = [float(line.split()[0]) for line in lines]

                # d-eg
                y_eg_sum = [
                    sum([float(line.split()[k]) for k in eg_index]) for line in lines
                ]

                # d-t2g
                y_t2g_sum = [
                    sum([float(line.split()[k]) for k in t2g_index]) for line in lines
                ]

                # px, px, pz
                if p_atoms:
                    y_p_sum = [
                        sum([float(line.split()[k]) for k in p_index]) for line in lines
                    ]

            if p_atoms:
                y_p = [-1 * count / math.pi for count in y_p_sum]
            y_eg = [-1 * count / math.pi for count in y_eg_sum]
            y_t2g = [-1 * count / math.pi for count in y_t2g_sum]

            # Plotting
            fig = plt.figure(figsize=(13, 9))
            ax = fig.add_subplot(111)

            if len(p["cor_at"]) < 2:
                eg_label = p["atomnames"][0] + " $d-e_g$"
                t2g_label = p["atomnames"][0] + " $d-t_{2g}$"
            else:
                eg_label = " $d-e_g$"
                t2g_label = " $d-t_{2g}$"

            ax.plot(x, y_eg, "r", label=eg_label)
            ax.plot(x, y_t2g, "b", label=t2g_label)
            if p_atoms:
                p_label = p["atomnames"][1] + " $p$"
                ax.plot(x, y_p, "g", label=p_label)
            ax.set_title("DMFT PDOS")
            ax.set_xlabel("Energy (eV)")
            ax.set_ylabel("DOS (states eV/cell)")
            ax.set_xlim(args.elim)
            ax.axvline(x=0, color="gray", linestyle="--")
            plt.legend()
            fig.savefig("./dos/DMFT-PDOS.png")
            if args.show:
                plt.show()
            f.close()

        # Spin polarized case
        if args.sp:
            print("Plotting spin polarized DOS...")
            with open("./dos/G_loc_up.out", "r") as f:
                lines = f.readlines()
                x = [float(line.split()[0]) for line in lines]

                # d-eg spin up
                y_eg_sum = [
                    sum([float(line.split()[k]) for k in eg_index]) for line in lines
                ]

                # d-t2g spin up
                y_t2g_sum = [
                    sum([float(line.split()[k]) for k in t2g_index]) for line in lines
                ]

                # px, px, pz
                if p_atoms:
                    y_p_sum = [
                        sum([float(line.split()[k]) for k in p_index]) for line in lines
                    ]

            if p_atoms:
                y_p = [-1 * count / 3.14 for count in y_p_sum]
            y_eg = [-1 * count / 3.14 for count in y_eg_sum]
            y_t2g = [-1 * count / 3.14 for count in y_t2g_sum]

            # spin down component
            with open("./dos/G_loc_dn.out", "r") as f:
                lines = f.readlines()
                x_dn = [float(line.split()[0]) for line in lines]

            # d-eg spin down
            y_eg_dn_sum = [
                sum([float(line.split()[k]) for k in eg_index]) for line in lines
            ]

            # d-t2g spin down
            y_t2g_dn_sum = [
                sum([float(line.split()[k]) for k in t2g_index]) for line in lines
            ]

            # px, px, pz
            if p_atoms:
                y_p_dn_sum = [
                    sum([float(line.split()[k]) for k in p_index]) for line in lines
                ]

            # x -1 for spin down components
            if p_atoms:
                y_p_dn = [1 * count / 3.14 for count in y_p_dn_sum]
            y_eg_dn = [1 * count / 3.14 for count in y_eg_dn_sum]
            y_t2g_dn = [1 * count / 3.14 for count in y_t2g_dn_sum]

            # Plotting
            fig = plt.figure(figsize=(13, 9))
            ax = fig.add_subplot(111)

            if len(p["cor_at"]) < 2:
                eg_label = p["atomnames"][0] + " $d-e_g$"
                t2g_label = p["atomnames"][0] + " $d-t_{2g}$"
            else:
                eg_label = " $d-e_g$"
                t2g_label = " $d-t_{2g}$"

            ax.plot(x, y_eg, "r", label=eg_label)
            ax.plot(x, y_t2g, "b", label=t2g_label)
            ax.plot(x_dn, y_eg_dn, "r")
            ax.plot(x_dn, y_t2g_dn, "b")
            if p_atoms:
                p_label = p["atomnames"][1] + " $p$"
                ax.plot(x, y_p, "g", label=p_label)
                ax.plot(x_dn, y_p_dn, "g")

            ax.set_title("DMFT PDOS")
            ax.set_xlabel("Energy (eV)")
            ax.set_ylabel("DOS (states eV/cell)")
            ax.axvline(x=0, color="gray", linestyle="--")
            ax.axhline(y=0, color="gray", linestyle="--")
            ax.set_xlim(min(min(x), min(x_dn)), max(max(x), max(x_dn)))
            plt.legend()
            fig.savefig("./dos/DMFT-PDOS_sp.png")
            if args.show:
                plt.show()
            f.close()

    def bands(self, args):
        """
	This method performs the band structure calculations.
	"""
        if args.elim:
            self.emin = args.elim[0]
            self.emax = args.elim[1]
        else:
            self.emin = -6.0
            self.emax = 6.0

        dest_dir = "bands"
        dummy_broaden = 1.0

        # spin polarized calculation?
        sp = args.sp

        # creating directory for bands
        if os.path.exists("bands"):
            print("bands directory already exists.")
        else:
            os.makedirs("bands")

        # copy Sig.out into /bands
        if self.checksig():
            shutil.copyfile("./ac/Sig.out", "./bands/Sig1.out")
        else:
            sys.exit()

        # interpolating
        self.interpol(self.emin, self.emax, args.rom, dummy_broaden, dest_dir, sp)

        #############################Xingu's contribution###################################################################################

        #########################Read POSCAR ###################################################
        TB = Struct.TBstructure("POSCAR", p["atomnames"], p["orbs"])
        cor_at = p["cor_at"]
        cor_orb = p["cor_orb"]
        TB.Compute_cor_idx(cor_at, cor_orb)
        ############################# Generate sequence list ##################################
        sort_atm = sorted([atm for atms in cor_at for atm in atms])
        atm_sqn = self.Make_coor_list(sort_atm, cor_at)
        # print("atm_sqn : %s" % atm_sqn)
        orb_sqn = []
        for i in atm_sqn:
            len_sf = 0
            for j in range(i):
                len_sf += max(Make_coor_list(TB.TB_orbs[cor_at[j][0]], cor_orb[j])) + 1
            orb_idx = self.Make_coor_list(TB.TB_orbs[cor_at[i][0]], cor_orb[i])
            for idx in orb_idx:
                orb_sqn.append(idx + len_sf)
        # print("orb_sqn : %s" % orb_sqn)
        ################################# Read sig.inp_real ####################################

        # non spin-polarized calculation
        if sp == False:

            fi = open("./bands/sig.inp_real", "r")
            [nom, ncor_orb] = [int(ele) for ele in fi.readline().split()[-2:]]
            temperature = float(fi.readline().split()[-1])
            sigmdc_tmp = [float(ele) for ele in fi.readline().split()[-1 * ncor_orb :]]
            sigoo = fi.readline().split()
            vdc = fi.readline().split()
            sig_real = np.array([ele.split() for ele in fi.readlines()])
            fi.close()
            assert np.shape(sig_real)[0] == nom
            ###################### Write SigMoo ##################
            SigMoo = np.zeros((nom, len(orb_sqn) * 2 + 1), dtype=float)
            SigMoo[:, 0] = np.array(sig_real[:, 0])
            for i in range(len(orb_sqn)):
                SigMoo[:, i * 2 + 1] = np.array(sig_real[:, orb_sqn[i] * 2 + 1])
                SigMoo[:, i * 2 + 2] = np.array(sig_real[:, orb_sqn[i] * 2 + 2])
            np.savetxt("./bands/SigMoo_real.out", SigMoo, fmt="%.10f")
            ############################ Srite SigMdc.out ########################
            SigMdc = np.array([sigmdc_tmp[i] for i in orb_sqn])
            np.savetxt("./bands/SigMdc.out", SigMdc[None], fmt="%.12f")

        # spin polarized calculation
        if sp:

            siginpreal_files = ["sig.inp_real", "sig.inp_real_dn"]
            SigMooreal_files = ["SigMoo_real.out", "SigMoo_dn_real.out"]
            SigMdc_files = ["SigMdc.out", "SigMdc_dn.out"]

            for filecounter in range(len(siginpreal_files)):
                filestr = "./bands/" + siginpreal_files[filecounter]
                fi = open(filestr, "r")
                [nom, ncor_orb] = [int(ele) for ele in fi.readline().split()[-2:]]
                temperature = float(fi.readline().split()[-1])
                sigmdc_tmp = [
                    float(ele) for ele in fi.readline().split()[-1 * ncor_orb :]
                ]
                sigoo = fi.readline().split()
                vdc = fi.readline().split()
                sig_real = np.array([ele.split() for ele in fi.readlines()])
                fi.close()
                assert np.shape(sig_real)[0] == nom
                ###################### Write SigMoo ##################
                SigMoo = np.zeros((nom, len(orb_sqn) * 2 + 1), dtype=float)
                SigMoo[:, 0] = np.array(sig_real[:, 0])
                for i in range(len(orb_sqn)):
                    SigMoo[:, i * 2 + 1] = np.array(sig_real[:, orb_sqn[i] * 2 + 1])
                    SigMoo[:, i * 2 + 2] = np.array(sig_real[:, orb_sqn[i] * 2 + 2])
                filestr = "./bands/" + SigMooreal_files[filecounter]
                np.savetxt(filestr, SigMoo, fmt="%.10f")
                ############################ Write SigMdc.out ########################
                SigMdc = np.array([sigmdc_tmp[i] for i in orb_sqn])
                # print(sigmdc_tmp)
                # print(SigMdc)
                filestr = "./bands/" + SigMdc_files[filecounter]
                np.savetxt(filestr, SigMdc[None], fmt="%.12f")

        # ################################################################################################################

        if args.autokp or args.compare:
            args.knames, ticks, discontinuities, args.kplist = self.readKPOINTS(args)

        print("kplist : %s" % args.kplist)
        print("knames : %s" % args.knames)

        # Iterating args.kpband to get correct k-list.
        indexerror = True

        while indexerror:

            try:
                # generating k-path
                klist, dist_K, dist_SK = self.Create_kpath(args.kplist, args.kpband)

                with open("./bands/klist.dat", "w") as fi:

                    # Put within try, exception to select correct kpband value.
                    print("Trying number of k-points (kpband) = %d" % args.kpband)
                    for i in range(args.kpband):
                        kcheck = 0
                        for j, d in enumerate(dist_SK):
                            # print(
                            #     "i:%d dist_K: %s dist_SK: %s abs: %f"
                            #     % (i, str(dist_K[i]), str(d), float(dist_K[i] - d))
                            # )
                            if abs(dist_K[i] - d) < 1e-10:
                                fi.write(
                                    "%.14f  %.14f  %.14f  %.14f  %s \n"
                                    % (
                                        dist_K[i],
                                        klist[i][0],
                                        klist[i][1],
                                        klist[i][2],
                                        args.knames[j],
                                    )
                                )
                                kcheck = 1
                                break
                        if kcheck == 0:
                            fi.write(
                                "%.14f  %.14f  %.14f  %.14f \n"
                                % (dist_K[i], klist[i][0], klist[i][1], klist[i][2])
                            )

                indexerror = False

            except IndexError:
                # Iterating args.kpband by 1.
                args.kpband += 1
                indexerror = True

        print("k-path generated.\n")

        # copying files from DMFT directory to dos directory
        cmd = "cd bands && Copy_input.py ../ -post bands"
        out, err = subprocess.Popen(cmd, shell=True).communicate()
        if err:
            print("File copy failed!\n")
            print(err)
            sys.exit()
        else:
            print(out)

        # generating ksum.input
        self.genksum(args.rom, args.kpband)

        # running dmft_ksum_band

        if args.plotplain:
            print("\nCalculating plain band structure...")
            cmd = "cd bands && " + self.para_com + " " + "dmft_ksum_band"
            out, err = subprocess.Popen(
                cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            ).communicate()
            if err:
                print(err)
                print("Band structure calculation failed!\n")
                sys.exit()
            else:
                print("Band structure calculation complete.\n")
                self.plot_plain_bands(args)

        if args.plotpartial:
            print("\nCalculating projected band structure...")
            cmd = "cd bands && " + self.para_com + " " + "dmft_ksum_partial_band"
            out, err = subprocess.Popen(
                cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            ).communicate()
            if err:
                print(err)
                print("Band structure calculation failed!\n")
                sys.exit()
            else:
                print("Band structure calculation complete.\n")
                self.plot_partial_bands(args)

        if sp:
            print("\nCalculating spin-polarized band structure...")
            cmd = "cd bands && " + self.para_com + " " + "dmft_ksum_band"
            out, err = subprocess.Popen(
                cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            ).communicate()
            if err:
                print(err)
                print("Band structure calculation failed!\n")
                sys.exit()
            else:
                print("Band structure calculation complete.\n")
                self.plot_sp_bands(args)

    def plot_plain_bands(self, args):
        """
        This method plots the regular DMFT band structure.
        """

        print("Plotting plain band structure...")

        nk = 0
        SKP = []
        SKPoints = []
        distk = []
        kpts = []
        fi = open("./bands/klist.dat", "r")
        for line in fi.readlines():
            line = line.split()
            distk.append(float(line[0]))
            kpts.append([float(line[1]), float(line[2]), float(line[3])])
            if len(line) == 5:
                SKP.append(float(line[0]))
                SKPoints.append(line[4])
        fi.close()

        fi = open("./bands/ksum.input", "r")
        numk = int(fi.readline())
        print("numk = %s" % numk)
        nom = int(fi.readline())
        print("nom = %s" % nom)
        fi.close()
        A_k = []
        dist_k = []
        om = []
        kpts = []
        fi = open("./bands/Gk.out", "r")
        for i in range(numk):
            kpts.append(list(map(float, fi.readline().split()[1:])))
            A_k.append([])
            om.append([])
            for j in range(nom):
                line = list(map(float, fi.readline().split()))
                A_k[i].append(-1 * line[2] / 3.14159265)
                om[i].append(line[0])
            A_k[i] = np.array(A_k[i])
        fi.close()

        A_k = np.transpose(A_k)[::-1]

        (ymin, ymax) = (om[0][0], om[0][-1])
        (xmin, xmax) = (distk[0], distk[-1])

        # Setting spectral function ranges.

        if args.vlim is None and args.normalize:
            vmax = max([max(p) for p in A_k])
            vmin = min([min(p) for p in A_k])
            norm = cm.colors.Normalize(vmax=vmax, vmin=vmin)
            print(
                "Normalizing spectral function range to : [%.2f, %.2f]" % (vmin, vmax)
            )

        elif args.vlim and args.normalize:
            vmax = args.vlim[1]
            vmin = args.vlim[0]
            norm = cm.colors.Normalize(vmax=vmax, vmin=vmin)
            print(
                "Normalizing spectral function range to : [%.2f, %.2f]" % (vmin, vmax)
            )

        elif args.vlim and args.normalize is False:
            vmax = args.vlim[1]
            vmin = args.vlim[0]
            norm = None
            print("Setting spectral function range to : [%.2f, %.2f]" % (vmin, vmax))

        else:
            vmax = 10.0
            vmin = 0.0
            norm = None
            print("Setting spectral function range to : [%.2f, %.2f]" % (vmin, vmax))

        # Plotting
        fig = plt.figure(figsize=(13, 9))
        ax = fig.add_subplot(111)

        im = ax.imshow(
            A_k,
            cmap=plt.get_cmap(args.cmap),
            vmin=vmin,
            vmax=vmax,
            norm=norm,
            extent=[xmin, xmax, ymin, ymax],
            aspect="auto",
        )
        cb = fig.colorbar(im, orientation="vertical", pad=0.05, shrink=1.0, ax=ax,)
        cb.ax.tick_params()

        ax.set_xticks(SKP)
        ax.set_xticklabels(SKPoints)
        ax.set_xlabel(r"$k$-path")
        ax.set_ylabel(r"$E-E_F$ [eV]")
        ax.axhline(y=0, color="black", ls="--")

        fig.tight_layout()

        if args.compare:
            fig, ax = self.plotDFTBands(args, fig=fig, ax=ax)

        if args.show:
            plt.show()
        fig.savefig("./bands/A_k.eps", format="eps", dpi=1200)

    def plot_partial_bands(self, args):
        """
        This method plots partial bands for orbitals. The order of the orbitals is the Wannier orbital order.
        """
        print("Wannier orbitals list : %s" % str(args.wanorbs))

        SKP = []
        SKPoints = []
        distk = []
        kpts = []
        fi = open("./bands/klist.dat", "r")
        for line in fi.readlines():
            line = line.split()
            distk.append(float(line[0]))
            kpts.append([float(line[1]), float(line[2]), float(line[3])])
            if len(line) == 5:
                SKP.append(float(line[0]))
                SKPoints.append(line[4])
        fi.close()

        fi = open("./bands/ksum.input", "r")
        numk = int(fi.readline())
        print("numk = %s" % numk)
        nom = int(fi.readline())
        print("nom = %s" % nom)
        fi.close()

        A_k = []
        dist_k = []
        om = []

        fi_gk = open("./bands/Gk.out", "r")
        data = re.findall("k=\s*[0-9E+-.\sorb=\s]*", fi_gk.read())
        fi_gk.close()
        filtered_orbs = []

        for i in range(numk):
            # kpts.append(list(map(float,data[i].split('\n')[0].split()[1:])))

            for orbs in args.wanorbs:
                filtered_orbs.append(
                    data[i].split("\n")[1:][(orbs - 1) * nom + orbs : orbs * nom + orbs]
                )

        for orb_counter in range(numk * len(args.wanorbs)):
            for j in range(nom):
                A_k.append(
                    -1 * float(filtered_orbs[orb_counter][j].split()[2]) / 3.14159265
                )
                om.append(float(filtered_orbs[orb_counter][j].split()[0]))

        A_k = np.array(A_k)
        A_k = A_k.reshape(numk, nom * len(args.wanorbs))

        A_kblend = np.zeros((len(args.wanorbs), numk, nom))
        A_ktotal = np.zeros((numk, nom))

        nom_counter = 0
        for orb in range(len(args.wanorbs)):
            A_kblend[orb, :, :] = A_k[:, nom_counter : nom_counter + nom]
            nom_counter = nom_counter + nom
            A_ktotal = A_ktotal + A_kblend[orb, :, :]

        A_ktotal = np.transpose(A_ktotal)[::-1]

        om = np.array(om)
        om = om.reshape(numk, nom * len(args.wanorbs))

        (ymin, ymax) = (om[0][0], om[0][-1])  # 500x100 energy matrix
        (xmin, xmax) = (distk[0], distk[-1])

        # Setting spectral function ranges.

        if args.vlim is None and args.normalize:
            vmax = max([max(p) for p in A_ktotal])
            vmin = min([min(p) for p in A_ktotal])
            norm = cm.colors.Normalize(vmax=vmax, vmin=vmin)
            print(
                "Normalizing spectral function range to : [%.2f, %.2f]" % (vmin, vmax)
            )

        elif args.vlim and args.normalize:
            vmax = args.vlim[1]
            vmin = args.vlim[0]
            norm = cm.colors.Normalize(vmax=vmax, vmin=vmin)
            print(
                "Normalizing spectral function range to : [%.2f, %.2f]" % (vmin, vmax)
            )

        elif args.vlim and args.normalize is False:
            vmax = args.vlim[1]
            vmin = args.vlim[0]
            norm = None
            print("Setting spectral function range to : [%.2f, %.2f]" % (vmin, vmax))

        else:
            vmax = 10.0
            vmin = 0.0
            norm = None
            print("Setting spectral function range to : [%.2f, %.2f]" % (vmin, vmax))

        # Plotting
        fig = plt.figure(figsize=(13, 9))
        ax = fig.add_subplot(111)

        im = ax.imshow(
            A_ktotal,
            cmap=plt.get_cmap(args.cmap),
            vmin=vmin,
            vmax=vmax,
            norm=norm,
            extent=[xmin, xmax, ymin, ymax],
            aspect="auto",
        )
        cb = fig.colorbar(im, orientation="vertical", pad=0.05, shrink=1.0, ax=ax,)
        cb.ax.tick_params()

        ax.set_xticks(SKP)
        ax.set_xticklabels(SKPoints)
        ax.set_xlabel(r"$k$-path")
        ax.set_ylabel(r"$E-E_F$ [eV]")
        ax.axhline(y=0, color="black", ls="--")

        fig.tight_layout()

        if args.compare:
            fig, ax = self.plotDFTBands(args, fig=fig, ax=ax)

        if args.show:
            plt.show()
        fig.savefig("./bands/A_k_partial.eps", format="eps", dpi=1200)

    def plot_sp_bands(self, args):
        """
	This method plots spin-polarized bands.
	"""

        print("Plotting spin-polarized band structure...")

        nk = 0
        SKP = []
        SKPoints = []
        distk = []
        kpts = []
        fi = open("./bands/klist.dat", "r")
        for line in fi.readlines():
            line = line.split()
            distk.append(float(line[0]))
            kpts.append([float(line[1]), float(line[2]), float(line[3])])
            if len(line) == 5:
                SKP.append(float(line[0]))
                SKPoints.append(line[4])
        fi.close()

        # Spin up dataset
        fi = open("./bands/ksum.input", "r")
        numk = int(fi.readline())
        print("numk = %s" % numk)
        nom = int(fi.readline())
        print("nom = %s" % nom)
        fi.close()
        A_k = []
        dist_k = []
        om = []
        kpts = []
        fi = open("./bands/Gk.out", "r")
        for i in range(numk):
            kpts.append(list(map(float, fi.readline().split()[1:])))
            A_k.append([])
            om.append([])
            for j in range(nom):
                line = list(map(float, fi.readline().split()))
                A_k[i].append(-1 * line[2] / 3.14159265)
                om[i].append(line[0])
            A_k[i] = np.array(A_k[i])
        fi.close()
        A_k = np.transpose(A_k)[::-1]

        # Spin down dataset
        A_k2 = []
        dist_k2 = []
        om2 = []
        kpts2 = []
        fi = open("./bands/Gk_dn.out", "r")
        for i2 in range(numk):
            kpts2.append(list(map(float, fi.readline().split()[1:])))
            A_k2.append([])
            om2.append([])
            for j2 in range(nom):
                line2 = list(map(float, fi.readline().split()))
                A_k2[i2].append(-1 * line2[2] / 3.14159265)
                om2[i2].append(line[0])
            A_k2[i2] = np.array(A_k2[i2])
        fi.close()
        A_k2 = np.transpose(A_k2)[::-1]

        (ymin, ymax) = (om[0][0], om[0][-1])
        (xmin, xmax) = (distk[0], distk[-1])

        # Setting spectral function ranges.

        if args.vlim is None and args.normalize:
            vmax1 = max([max(p) for p in A_k])
            vmin1 = min([min(p) for p in A_k])
            vmax2 = max([max(p) for p in A_k2])
            vmin2 = min([min(p) for p in A_k2])
            vmax = max([vmax1, vmax2])
            vmin = min([vmin1, vmin2])
            norm = cm.colors.Normalize(vmax=vmax, vmin=vmin)
            print(
                "Normalizing spectral function range to : [%.2f, %.2f]" % (vmin, vmax)
            )

        elif args.vlim and args.normalize:
            vmax = args.vlim[1]
            vmin = args.vlim[0]
            norm = cm.colors.Normalize(vmax=vmax, vmin=vmin)
            print(
                "Normalizing spectral function range to : [%.2f, %.2f]" % (vmin, vmax)
            )

        elif args.vlim and args.normalize is False:
            vmax = args.vlim[1]
            vmin = args.vlim[0]
            norm = None
            print("Setting spectral function range to : [%.2f, %.2f]" % (vmin, vmax))

        else:
            vmax = 10.0
            vmin = 0.0
            norm = None
            print("Setting spectral function range to : [%.2f, %.2f]" % (vmin, vmax))

        # Plotting
        fig = plt.figure(figsize=(13, 9))

        if args.sp and args.spinup is False and args.spindown is False:
            # subplots

            # Big subplot
            ax = fig.add_subplot(111)
            # Turn off axis lines and ticks of the big subplot
            ax.spines["top"].set_color("none")
            ax.spines["bottom"].set_color("none")
            ax.spines["left"].set_color("none")
            ax.spines["right"].set_color("none")
            ax.tick_params(
                labelcolor="w", top=False, bottom=False, left=False, right=False
            )

            # spin up
            ax1 = fig.add_subplot(1, 2, 1)
            im1 = ax1.imshow(
                A_k,
                cmap=plt.get_cmap(args.cmap),
                vmin=vmin,
                vmax=vmax,
                norm=norm,
                extent=[xmin, xmax, ymin, ymax],
                aspect="auto",
            )
            cb1 = fig.colorbar(
                im1, orientation="vertical", pad=0.05, shrink=1.0, ax=ax1,
            )
            cb1.ax.tick_params()
            ax1.set_title("Spin Up")
            ax1.set_xticks(SKP)
            ax1.set_xticklabels(SKPoints)
            ax1.axhline(y=0, color="black", ls="--")

            # spin down
            ax2 = fig.add_subplot(1, 2, 2)
            im2 = ax2.imshow(
                A_k2,
                cmap=plt.get_cmap(args.cmap),
                vmin=vmin,
                vmax=vmax,
                norm=norm,
                extent=[xmin, xmax, ymin, ymax],
                aspect="auto",
            )
            cb2 = fig.colorbar(
                im2, orientation="vertical", pad=0.05, shrink=1.0, ax=ax2,
            )
            cb2.ax.tick_params()
            ax2.set_title("Spin Down")
            ax2.set_xticks(SKP)
            ax2.set_xticklabels(SKPoints)
            ax2.axhline(y=0, color="black", ls="--")

            # Set common labels
            ax.set_xlabel(r"$k$-path")
            ax.set_ylabel(r"$E-E_F$")

            # Set common labels
            # fig.text(0.2, 0.04, "$k$-path", ha="center", va="center")
            # fig.text(0.5, 0.2, "$E-E_F$", ha="center", va="center", rotation="vertical")

            # common colorbar
            # fig.subplots_adjust(right=0.8)
            # cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
            # cb = fig.colorbar(
            #     im, orientation="vertical", pad=0.05, shrink=1.0, cax=cbar_ax,
            # )
            # cb.ax.tick_params()

            fig.tight_layout()

            if args.compare:
                fig, ax1 = self.plotDFTBands(args, fig=fig, ax=ax1)
                fig, ax2 = self.plotDFTBands(args, fig=fig, ax=ax2)

            if args.show:
                show()
            fig.savefig("./bands/A_k_sp.eps", format="eps", dpi=1200)

        elif args.sp and args.spinup and args.spindown is False:

            # Plotting only spin up bands
            ax = fig.add_subplot(111)

            im = ax.imshow(
                A_k,
                cmap=plt.get_cmap(args.cmap),
                vmin=vmin,
                vmax=vmax,
                norm=norm,
                extent=[xmin, xmax, ymin, ymax],
                aspect="auto",
            )
            cb = fig.colorbar(im, orientation="vertical", pad=0.05, shrink=1.0, ax=ax,)
            cb.ax.tick_params()

            ax.set_title("Spin Up")
            ax.set_xticks(SKP)
            ax.set_xticklabels(SKPoints)
            ax.set_xlabel(r"$k$-path")
            ax.set_ylabel(r"$E-E_F$ [eV]")
            ax.axhline(y=0, color="black", ls="--")

            fig.tight_layout()

            if args.compare:
                fig, ax = self.plotDFTBands(args, fig=fig, ax=ax)

            if args.show:
                plt.show()
            fig.savefig("./bands/A_k_spinup.eps", format="eps", dpi=1200)

        elif args.sp and args.spindown and args.spinup is False:

            # Plotting only spin down bands
            ax = fig.add_subplot(111)

            im = ax.imshow(
                A_k2,
                cmap=plt.get_cmap(args.cmap),
                vmin=vmin,
                vmax=vmax,
                norm=norm,
                extent=[xmin, xmax, ymin, ymax],
                aspect="auto",
            )
            cb = fig.colorbar(im, orientation="vertical", pad=0.05, shrink=1.0, ax=ax,)
            cb.ax.tick_params()

            ax.set_title("Spin Down")
            ax.set_xticks(SKP)
            ax.set_xticklabels(SKPoints)
            ax.set_xlabel(r"$k$-path")
            ax.set_ylabel(r"$E-E_F$ [eV]")
            ax.axhline(y=0, color="black", ls="--")

            fig.tight_layout()

            if args.compare:
                fig, ax = self.plotDFTBands(args, fig=fig, ax=ax)

            if args.show:
                plt.show()
            fig.savefig("./bands/A_k_spindown.eps", format="eps", dpi=1200)

    def plotDFTBands(self, args, fig=None, ax=None):
        """This function plots DFT bands using the
        EIGENVAL file."""

        # plotting style
        marker = "."
        markersize = 0.01
        color = "lawngreen"
        linewidth = 0.25
        linestyle = "dashed"

        ##### Reading the EIGENVAL file #####

        fi = open("EIGENVAL", "r")
        for i in range(5):
            skip = fi.readline()
        header = fi.readline()
        skip = fi.readline()
        data = fi.readlines()
        fi.close()

        numkpoints = int(header.split()[1])
        numbands = int(header.split()[2])

        kpoints = np.zeros((numkpoints, 3), dtype="float64")
        bands = np.zeros((numkpoints, numbands), dtype="float64")

        kpointscounter = 0
        for i in data:
            if len(i.split()) == 4:
                kpoints[kpointscounter, :] = [float(x) for x in i.split()[0:3]]
                kpointscounter += 1
            if len(i.split()) == 3:
                bands[kpointscounter - 1, int(i.split()[0]) - 1] = float(i.split()[1])

        ##### Reading DFT_mu.out to get Fermi energy #####
        fi = open("DFT_mu.out", "r")
        fermi = float(fi.readline())
        fi.close

        # get knames and kticks from readKPOINTS()
        ticksNames, ticks, discontinuities, kplist = self.readKPOINTS(args)

        ##### Plotting bands #####

        bands = (bands.transpose() - np.array(fermi)).transpose()
        bands = bands.transpose()

        if kpoints is not None:
            xaxis = [0]

            #### MODIFIED FOR DISCONTINOUS BANDS ####
            if ticks:

                # counters for number of discontinuities
                icounter = 1
                ii = 0

                for i in range(1, len(kpoints) - len(discontinuities)):
                    d = kpoints[icounter] - kpoints[icounter - 1]
                    d = np.sqrt(np.dot(d, d))
                    xaxis.append(d + xaxis[-1])
                    icounter += 1
                    ii += 1
                    if ii in discontinuities:
                        icounter += 1
                        ii += 1
                        xaxis.append(xaxis[-1])
                xaxis = np.array(xaxis)

                # plotting
                for i_tick in range(len(ticks) - 1):
                    x = xaxis[ticks[i_tick] : ticks[i_tick + 1] + 1]
                    y = bands.transpose()[ticks[i_tick] : ticks[i_tick + 1] + 1, :]
                    ax.plot(
                        x,
                        y,
                        marker=marker,
                        markersize=markersize,
                        color=color,
                        linewidth=linewidth,
                        linestyle=linestyle,
                    )

            #### END  OF MODIFIED DISCONTINUOUS BANDS ####

        ax.set_xlim(xaxis.min(), xaxis.max())

        # Handling ticks
        if ticks:
            # added for meta-GGA calculations
            if ticks[0] > 0:
                ax.set_xlim(left=xaxis[ticks[0]])
            ticks = [xaxis[x] for x in ticks]
            ax.set_xticks(ticks)
            ax.set_xticklabels(ticksNames)
            ax.set_ylim([self.emin, self.emax])
            # for xc in ticks:
            #    ax.axvline(x=xc, color="k")
        # ax.axhline(color="black", linestyle="--")

        return fig, ax

    def readKPOINTS(self, args):
        """Reads KPOINTS file to get knames and kticks."""

        ##### Reading the KPOINTS file #####

        # Getting the high symmetry point names from KPOINTS file
        f = open("KPOINTS", "r")
        KPread = f.read()
        f.close()

        KPmatrix = re.findall("reciprocal[\s\S]*", KPread)
        tick_labels = np.array(re.findall("!\s*(.*)", KPmatrix[0]))
        knames = []
        knames = [tick_labels[0]]

        # Creating kplist
        splitarray = KPmatrix[0].split("!")
        kplist = []
        for i in splitarray[0:-1]:
            kplist.append([float(x) for x in i.split()[-3:]])
        kplist = list(kplist for kplist, _ in itertools.groupby(kplist))

        # Checking for discontinuities
        discont_indx = []
        icounter = 1
        while icounter < len(tick_labels) - 1:
            if tick_labels[icounter] == tick_labels[icounter + 1]:
                knames.append(tick_labels[icounter])
                icounter = icounter + 2
            else:
                discont_indx.append(icounter)
                knames.append(tick_labels[icounter] + "|" + tick_labels[icounter + 1])
                icounter = icounter + 2
        knames.append(tick_labels[-1])
        discont_indx = list(dict.fromkeys(discont_indx))

        # End of discontinuity check

        # Improve latex rendering
        for i in range(len(knames)):
            if knames[i] == "GAMMA":
                knames[i] = "\Gamma"
            else:
                pass
        knames = [str("$" + latx + "$") for latx in knames]

        # getting the number of grid points from the KPOINTS file
        f2 = open("KPOINTS", "r")
        KPreadlines = f2.readlines()
        f2.close()
        numgridpoints = int(KPreadlines[1].split()[0])

        kticks = [0]
        gridpoint = 0
        for kt in range(len(knames) - 1):
            gridpoint = gridpoint + numgridpoints
            kticks.append(gridpoint - 1)

        # creating an array for discontunuity k-points. These are the indexes
        # of the discontinuity k-points.
        discontinuities = []
        for k in discont_indx:
            discontinuities.append(kticks[int(k / 2) + 1])
        if discontinuities:
            print("discont. list  : %s " % str(discontinuities))

        return knames, kticks, discontinuities, kplist


if __name__ == "__main__":
    args = sys.argv[1:]
    if args:

        welcome()
        des = "This script performs Analytic Contiunation, Density of States and Band structure calculations from DMFTwDFT outputs.\nRun inside DMFT or HF directory."
        parser = argparse.ArgumentParser(
            description=des, formatter_class=RawTextHelpFormatter
        )
        subparsers = parser.add_subparsers(help="sub-command help")

        # parser for ac
        parser_ac = subparsers.add_parser("ac", help="Analytic Continuation")
        parser_ac.add_argument(
            "-siglistindx",
            default=2,
            type=int,
            help="How many last self energy files to average?",
        )
        parser_ac.set_defaults(func=PostProcess().anal_cont)

        # parser for dos
        parser_dos = subparsers.add_parser("dos", help="DMFT Density of States")
        parser_dos.add_argument(
            "-emin", default=-5.0, type=float, help="Minimum value for interpolation"
        )
        parser_dos.add_argument(
            "-emax", default=5.0, type=float, help="Maximum value for interpolation"
        )
        parser_dos.add_argument(
            "-sp", action="store_true", help="Flag to plot spin-polarized DOS"
        )
        parser_dos.add_argument(
            "-rom", default=1000, type=int, help="Matsubara Frequency (omega) points"
        )
        parser_dos.add_argument("-broaden", default=0.03, type=float, help="Broadening")
        parser_dos.add_argument(
            "-show", action="store_true", help="Display the density of states"
        )
        parser_dos.add_argument(
            "-elim", type=float, nargs=2, help="Energy range to plot"
        )
        parser_dos.set_defaults(func=PostProcess().dos)

        # parser for bands
        parser_bands = subparsers.add_parser("bands", help="DMFT Bandstructure")
        parser_bands.add_argument(
            "-elim", type=float, nargs=2, help="Energy range to plot"
        )
        parser_bands.add_argument(
            "-rom", default=1000, type=int, help="Matsubara Frequency (omega) points"
        )
        parser_bands.add_argument(
            "-kpband",
            default=500,
            type=int,
            help="Number of k-points for band structure calculation",
        )
        parser_bands.add_argument(
            "-kn",
            "--knames",
            default=["$\Gamma$", "$X$", "$M$", "$\Gamma$", "$R$"],
            type=str,
            nargs="+",
            help="Names of the k-points",
        )
        parser_bands.add_argument(
            "-kp",
            "--kplist",
            default=[[0, 0, 0], [0.5, 0, 0], [0.5, 0.5, 0], [0, 0, 0], [0.5, 0.5, 0.5]],
            type=int,
            nargs="+",
            action="append",
            help="List of k-points as an array",
        )
        parser_bands.add_argument(
            "-autokp",
            action="store_true",
            help="Flag to use KPOINTS file to obtain k-path.",
        )
        parser_bands.add_argument(
            "-cmap",
            default="hot",
            type=str,
            help="Colormap to use for plotting spectral function.",
        )
        parser_bands.add_argument(
            "-plotplain", action="store_true", help="Flag to plot plain band structure"
        )
        parser_bands.add_argument(
            "-sp",
            action="store_true",
            help="Flag to plot spin-polarized band structure. Default plots both spin up and spin down plots on the same figure.",
        )
        parser_bands.add_argument(
            "-spinup",
            action="store_true",
            help="Flag to plot spin up band structure separately.",
        )

        parser_bands.add_argument(
            "-spindown",
            action="store_true",
            help="Flag to plot spin down band structure separately.",
        )
        parser_bands.add_argument(
            "-compare",
            action="store_true",
            help="Compare with DFT band structure (requires KPOINTS and EIGENVAL).",
        )

        parser_bands.add_argument(
            "-plotpartial",
            action="store_true",
            help="Flag to plot projected band structure",
        )
        parser_bands.add_argument(
            "-wo",
            "--wanorbs",
            default=[1, 4],
            type=int,
            nargs="+",
            help="List of Wannier orbitals to project. Ordering follows atom order in structure and the Wannier orbital order. Starts from 1.",
        )
        parser_bands.add_argument(
            "-vlim",
            type=float,
            nargs=2,
            help="Spectral intensity range. If -normalize flag is set,\
            this will correspond to the min and max values of the normalization range.",
        )
        parser_bands.add_argument(
            "-normalize",
            action="store_true",
            help="Normalize spectral intensity range. -vlim sets range.",
        )
        parser_bands.add_argument(
            "-show", action="store_true", help="Display the bands"
        )
        parser_bands.set_defaults(func=PostProcess().bands)

        args = parser.parse_args()
        args.func(args)

    else:
        print("Usage: postDMFT.py -h")
