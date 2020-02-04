#!/usr/bin/env python2

from scipy import *
import copy, Fileio, re
from scipy import interpolate

if __name__=='__main__':

   emin=-7.0;emax=3.0;rom=3000;broaden=0.03
   headerline=2
   om,Sig=Fileio.Read_complex_multilines('Sig.out',headerline)
   s_oo = None
   Vdc = None
   fi=open('Sig.out','r')
   for i in range(headerline):
      line=fi.readline()
      m=re.search('#(.*)',line)
      exec(m.group(1).strip())
   s_oo_Vdc=array(s_oo)-array(Vdc)

   ommesh=linspace(emin,emax,rom)
   Sig_tot=zeros((len(Sig),rom),dtype=complex)
   for i in range(len(Sig)):
      SigSpline = interpolate.splrep(om, Sig[i].real, k=1, s=0)
      Sig_tot[i,:] += interpolate.splev(ommesh, SigSpline)
      SigSpline = interpolate.splrep(om, Sig[i].imag, k=1, s=0)
      Sig_tot[i,:] += 1j*interpolate.splev(ommesh, SigSpline)

   header1='# nom,ncor_orb= '+str(len(ommesh))+' '+str(len(Sig_tot))
   #header2='# T= %18.15f'%(1.0/pC['beta'][0])#+str(self.T)
   header2='# T= %18.15f'%(broaden)#+str(self.T)
   header3='# s_oo-Vdc= '
   for i in range(len(s_oo_Vdc)):
      header3+='%18.15f '%(s_oo_Vdc[i])
   header4='# s_oo= '+str(s_oo)
   header5='# Vdc= '+str(Vdc)

   Fileio.Print_complex_multilines(Sig_tot,ommesh,'sig.inp_real',[header1,header2,header3,header4,header5])

