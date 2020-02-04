#!/usr/bin/env python2
from scipy import *
import os,subprocess
from mysub import *
from scipy import linalg
import time
def now(): return time.strftime("at %H:%M:%S, %D")

class WANNIER:
   def __init__(self,seedname):
      self.name=seedname
      
   def load_eig(self):
      if os.path.exists(self.name+'.eig'): print "eig file exists!"
      else: print "eig file does not exist! We are exiting"; exit()
      fi_eig=open(self.name+'.eig','r')
      eigvals=[float(line.split()[2]) for line in fi_eig.readlines()] 
      self.num_kpts=int(line.split()[1]);self.num_bands=int(line.split()[0])
      self.eigvals=array(eigvals).reshape(self.num_kpts,self.num_bands)
      fi_eig.close()
   def load_chk(self,path):
      """ Load seedname.chk file and store """
      def list_to_complex(a): return complex(float(a[0]),float(a[1]))
      if os.path.exists(self.name+'.chk.fmt'): print "formmated chk file exists!"
      elif os.path.exists(self.name+'.chk'):
         print "unformmated chk file exists! we will change it to formmted file."
         cmd = path+"w90chk2chk.x -export "+self.name
         out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
         print out
      else: print "chk file does not exist! We are exiting"; exit() 
      fi_chk=open(self.name+'.chk.fmt','r')
      header=fi_chk.readline()
      num_bands=int(fi_chk.readline())
      num_exclude_bands=int(fi_chk.readline())
      if num_exclude_bands>0: excl_bands=map(int,fi_chk.readline().split())
      real_latt=map(float,fi_chk.readline().split())
      real_latt=array(real_latt).reshape((3,3),order='F')
      recip_latt=map(float,fi_chk.readline().split())
      self.recip_latt=array(recip_latt).reshape((3,3),order='F')
      num_kpts=int(fi_chk.readline())
      self.num_kpts=num_kpts
      mp_grid=map(int,fi_chk.readline().split())
      self.mp_grid=array(mp_grid)
      kpt_latt=[map(float,fi_chk.readline().split()) for i in range(num_kpts)]
      self.kpt_latt=array(kpt_latt)
      self.ikpt=array([map(int,map(round,array(mp_grid)*array(kp)))%array(mp_grid) for kp in kpt_latt])
      nn_tot=int(fi_chk.readline())
      num_wann=int(fi_chk.readline())
      self.num_wann=num_wann
      checkpoint=fi_chk.readline()
      have_disentangled=int(fi_chk.readline())
      if have_disentangled: 
         omega_invariant=float(fi_chk.readline())
         lwindow=[int(fi_chk.readline()) for i in range(num_kpts) for j in range(num_bands)]
         lwindow=array(lwindow).reshape(num_kpts,num_bands)
         ndimwin=[int(fi_chk.readline()) for i in range(num_kpts)]
         u_matrix_opt=[list_to_complex(fi_chk.readline().split()) for i in range(num_kpts) for j in range(num_wann) for k in range(num_bands)]
         u_matrix_opt=array(u_matrix_opt).reshape(num_kpts,num_wann,num_bands)
      u_matrix=[list_to_complex(fi_chk.readline().split()) for i in range(num_kpts) for j in range(num_wann) for k in range(num_wann)] 
      u_matrix=array(u_matrix).reshape(num_kpts,num_wann,num_wann)
      fi_chk.close()
      #m_matrix=[list_to_complex(fi_chk.readline().split()) for i in range(num_kpts) for j in range(nn_tot) for k in range(num_wann) for l in range(num_wann)] 
      #m_matrix=array(m_matrix).reshape(num_kpts,nn_tot,num_wann,num_wann)
      #wannier_centres=[map(float,fi_chk.readline().split()) for i in range(num_wann)]
      #wannier_centres=array(wannier_centres).reshape(num_wann,3)
      #wannier_spreads=[float(fi_chk.readline()) for i in range(num_wann)]

      num_tot_bands=num_exclude_bands+num_bands
      lexclude_band=zeros(num_tot_bands,dtype=int)
      for i in range(num_exclude_bands): lexclude_band[excl_bands[i]]=1
      num_band_max=1; band_win=[]
      for nkp in range(num_kpts):
         IP=0;nbmin=num_tot_bands-1;nbmax=0
         for nb in range(num_tot_bands):
            if lexclude_band[nb]: continue
            IP+=1
            if not lwindow[nkp,IP-1]: continue
            if nb+1<nbmin: nbmin=nb
            if nb+1>nbmax: nbmax=nb
         band_win.append([nbmin,nbmax])
         if num_band_max<nbmax-nbmin+1: num_band_max=nbmax-nbmin+1
      self.num_band_max=num_band_max
      self.band_win=array(band_win)
      self.WANU=zeros((num_kpts,num_wann,num_band_max),dtype=complex)
      for i in range(num_kpts): self.WANU[i]=u_matrix[i].dot(u_matrix_opt[i][:,:num_band_max])
   
   def Compute_dU(self):
      ########## READ AMN file ##########
      def list_to_complex(a): return complex(float(a[0]),float(a[1]))
      if os.path.exists(self.name+'.amn'): print "amn file exists!"
      else: print "amn file does not exist! We are exiting"; exit()
      fi_amn=open(self.name+'.amn','r')
      fi_amn.readline()
      line=fi_amn.readline().split()
      num_bands=int(line[0]); num_kpts=int(line[1]);  num_proj=int(line[2])
      AMN=[list_to_complex(fi_amn.readline().split()[3:5]) for i in range(num_kpts) for j in range(num_proj) for k in range(num_bands)]
      fi_amn.close()
      AMN=array(AMN).reshape(num_kpts,num_proj,num_bands)
      ######## READ dAMN file #############
      if os.path.exists(self.name+'.damn'): print "damn file exists!"
      else: print "damn file does not exist! We are exiting"; exit()
      fi_amn=open(self.name+'.damn','r')
      fi_amn.readline()
      line=fi_amn.readline().split()
      num_bands=int(line[0]); num_kpts=int(line[1]);  num_proj=int(line[2]); num_dir=int(line[3])
      Atom=[[0,5],[5,8],[8,11],[11,14]]
      self.num_dir=num_dir; self.num_fdir=num_dir*len(Atom)
      dAMN=[list_to_complex(fi_amn.readline().split()[4:6]) for i in range(num_kpts) for j in range(num_proj) for k in range(num_bands) for l in range(num_dir)]
      dAMN=array(dAMN).reshape(num_kpts,num_proj,num_bands,num_dir)
      fi_amn.close()

      self.dWANU=zeros((num_kpts,num_proj,self.num_band_max,len(Atom)*num_dir),dtype=complex)
      for i in range(num_kpts):
         nbmin=self.band_win[i][0]; nbmax=self.band_win[i][1]
         num_band_max=nbmax-nbmin+1
         FAC=array(matrix(self.WANU[i][:,:num_band_max])*matrix(AMN[i][:,nbmin:nbmax+1]).I)
         for ii in range(len(Atom)):
            for k in range(num_dir):
               dU=zeros((num_proj,num_band_max),dtype=complex)
               dU[Atom[ii][0]:Atom[ii][1],:]=dAMN[i][Atom[ii][0]:Atom[ii][1],nbmin:nbmax+1,k]
               self.dWANU[i,:,:num_band_max,ii*num_dir+k]=FAC.dot(dU)
                  



   def Compute_HamR0(self):
      self.HamR0=zeros((self.num_wann,self.num_wann),dtype=complex)
      for i in range(self.num_kpts):
         nbmin=self.band_win[i][0]; nbmax=self.band_win[i][1]
         num_band_max=nbmax-nbmin+1
         Uwanband=matrix(self.WANU[i][:self.num_wann,:num_band_max])
         self.HamR0+=Uwanband*diag(self.eigvals[i][nbmin:nbmax+1])*Uwanband.H
      self.HamR0=self.HamR0/self.num_kpts
         
   def Check_Unitarity(self):
      for i in range(self.num_kpts):
         nbmin=self.band_win[i][0]; nbmax=self.band_win[i][1]
         num_band_max=nbmax-nbmin+1
         Uwanband=matrix(self.WANU[i][:5,:num_band_max])
         print sum(diag(Uwanband*Uwanband.H))

if __name__=='__main__':
   WAN=WANNIER('wannier90')
   WAN.load_chk('/home/uthpala/Documents/Research/projects/DMFTwDFT/bin/')
   WAN.load_eig()
   WAN.Check_Unitarity() 
   print WAN.ikpt
   print WAN.recip_latt
   exit()

   print "Start", now()
   WAN.Compute_HamR0()
      #Ham0=zeros((WAN.num_wann,WAN.num_wann),dtype=complex)
      #for i in range(WAN.num_kpts): 
      #   #WANU=matrix(WAN.u_matrix[i])*matrix(WAN.u_matrix_opt[i][:,:WAN.num_band_max])
      ###   #WANU=matrix(WAN.u_matrix[i].dot(WAN.u_matrix_opt[i]))
      #   nbmin=WAN.band_win[i][0]; nbmax=WAN.band_win[i][1]
      #   num_band_max=nbmax-nbmin+1
      #   Sigma=linalg.diagsvd(WAN.eigvals[i][nbmin:nbmax+1],num_band_max,num_band_max)
      #   #Sigma=ones((num_band_max,num_band_max))
      #   Ham0+=matrix(WAN.WANU[i][:num_band_max,:num_band_max])*matrix(Sigma)*matrix(WAN.WANU[i][:num_band_max,:num_band_max]).H
      #   #mul=[sum(WAN.WANU[i][ii][:num_band_max]*WAN.eigvals[i][nbmin:nbmax+1]*(WAN.WANU[i][jj][:num_band_max]).conjugate()) for ii in range(WAN.num_wann) for jj in range(WAN.num_wann)]
      #   #Ham0+=array(mul).reshape(WAN.num_wann,WAN.num_wann)
   print "end", now()
   #print linalg.eigvalsh(Ham0)
   pmatrix(array(WAN.HamR0)[:5,:5].real)
   ##pmatrix(WAN.u_matrix[0])
   #pmatrix(WAN.u_matrix_opt[0])
   ##WAN.print_H0()
   execfile('rham.py')
   #num_wann=len(Hopping[(0,0,0)])
   #Hamk=zeros((num_wann,num_wann),dtype=complex)
   #for rvec in Hopping.keys():
   #   Hamk+=Hopping[rvec]
   ##Hamk=Hamk/len(Hopping.keys())
   #print linalg.eigvalsh(Hamk)
   pmatrix(array(Hopping[(0,0,0)])[:5,:5].real)
