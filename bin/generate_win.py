#!/usr/bin/env python2
import VASP
import Struct

if __name__ == '__main__':

  p = {"Niter":     1,               # Number of DFT+DMFT iterations
     "atomnames": ['Ni','O'],       # The name of atoms"],
     "orbs"     : ['d','p'],       # The name  of orbitals
     "L_rot"     : [1,0],       # Whether rotate local axis or not
     "ewin":      [-8,3.0],           # Energy Window
     "cor_at":    [['Ni1','Ni2']],      # Correlated atoms, put degenerate atoms in the same list
     "cor_orb":   [[['d_z2','d_x2y2'],['d_xz','d_yz','d_xy']]], # DMFT orbitals, other orbitals are treated by HF"],
  }

  TB=Struct.TBstructure('POSCAR',p['atomnames'],p['orbs'])
  TB.Compute_cor_idx(p['cor_at'],p['cor_orb'])
  print TB.TB_orbs
  DFT=VASP.VASP_class()
  DFT.NBANDS=72
  DFT.Create_win(TB,p['atomnames'],p['orbs'],p['L_rot'],DFT.NBANDS,DFT.EFERMI+p['ewin'][0],DFT.EFERMI+p['ewin'][1])