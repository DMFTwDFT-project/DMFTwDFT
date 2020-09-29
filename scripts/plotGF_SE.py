#!/usr/bin/env python

""" Plotting DMFTwDFT Green's functions and Self Energies.

This script plots the real and imaginary parts of the Green's function (G_loc.out files)
and the self-energies (sig.inp files) calculated from the DMFTwDFT code.

The order of the columns are:

|---------------------|-----------------------|----------------------------|
| Matsubara Frequency | Real part of GF or SE | Imaginary part of GF or SE |  x Repeats for each group
|---------------------|-----------------------|----------------------------|

of "cor_orb".

Run this script inside the DMFT folder.

Currently, this script is for correlated d-orbitals considering eg and t2g i.e.

"cor_orb":   [[['d_z2','d_x2y2'],['d_xz','d_yz','d_xy']]]

as set in INPUT.py. Please modify according to your requirements.


"""

import matplotlib.pyplot as plt
import glob


# plotting Green's functions
gf_file = sorted(
    sorted(glob.glob("G_loc.out.*")), key=lambda x: (len(x), float(x[10:]))
)[-1]

with open(gf_file, "r") as f:
    lines = f.readlines()
    x = [float(line.split()[0]) for line in lines]
    y1_real = [float(line.split()[1]) for line in lines]  # eg
    y1_imag = [float(line.split()[2]) for line in lines]
    y2_real = [float(line.split()[3]) for line in lines]  # t2g
    y2_imag = [float(line.split()[4]) for line in lines]

plt.figure(1)
plt.plot(x, y1_real, label="eg")
plt.plot(x, y2_real, label="t2g")
plt.title("Real Green function")
plt.xlabel("$i{\omega_n}$")
plt.ylabel("Re $G(i{\omega_n})$")
plt.legend()
plt.savefig("Gf_real.png")
#plt.show()

plt.figure(2)
plt.plot(x, y1_imag, label="eg")
plt.plot(x, y2_imag, label="t2g")
plt.title("Imaginary Green function")
plt.legend()
plt.xlabel("$i{\omega_n}$")
plt.ylabel("Im $G(i{\omega_n})$")
plt.savefig("Gf_imag.png")
#plt.show()
f.close()

# plotting Self-energies
se_file = sorted(
    sorted(glob.glob("sig.inp.*"))[1:], key=lambda x: (len(x), float(x[8:]))
)[-1]

with open(se_file, "r") as f:
    for i in range(5):
        f.readline()
    lines = f.readlines()
    x = [float(line.split()[0]) for line in lines]
    y1_real = [float(line.split()[1]) for line in lines]  # eg
    y1_imag = [float(line.split()[2]) for line in lines]
    y2_real = [float(line.split()[3]) for line in lines]  # t2g
    y2_imag = [float(line.split()[4]) for line in lines]

plt.figure(3)
plt.plot(x, y1_real, label="eg")
plt.plot(x, y2_real, label="t2g")
plt.title("Real Self-Energy")
plt.xlabel("$i{\omega_n}$")
plt.ylabel("Re $\Sigma(i{\omega_n})$")
plt.legend()
plt.savefig("Selfenergy_real.png")
#plt.show()

plt.figure(4)
plt.plot(x, y1_imag, label="eg")
plt.plot(x, y2_imag, label="t2g")
plt.title("Imaginary Self-Energy")
plt.legend()
plt.xlabel("$i{\omega_n}$")
plt.ylabel("Im $\Sigma(i{\omega_n})$")
plt.savefig("Selfenergy_imaginary.png")
#plt.show()
f.close()

# plotting analytically continued Self-energies
with open("./ac/Sig.out", "r") as f:
    for i in range(2):
        f.readline()
    lines = f.readlines()
    x = [float(line.split()[0]) for line in lines]
    y1_real = [float(line.split()[1]) for line in lines]  # eg
    y1_imag = [float(line.split()[2]) for line in lines]
    y2_real = [float(line.split()[3]) for line in lines]  # t2g
    y2_imag = [float(line.split()[4]) for line in lines]

plt.figure(5)
plt.plot(x, y1_real, label="eg")
plt.plot(x, y2_real, label="t2g")
plt.title("Real Self-Energy Analytically continued")
plt.legend()
# plt.ylim(-0.1,0.1)
plt.xlabel("$\omega$")
plt.ylabel("Re $\Sigma(\omega})$")
plt.savefig("RealSelf-EnergyAnalyticallycontinued.png")
#plt.show()

plt.figure(6)
plt.plot(x, y1_imag, label="eg")
plt.plot(x, y2_imag, label="t2g")
plt.title("Imaginary Self-Energy Analytically continued")
plt.legend()
# plt.ylim(-0.1,0.1)
plt.xlabel("$\omega$")
plt.ylabel("Im $\Sigma(\omega})$")
plt.savefig("Imaginary Self-Energy Analytically continued.png")
#plt.show()
f.close()
