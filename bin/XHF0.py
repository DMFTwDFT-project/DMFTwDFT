#!/usr/bin/env python

import sys, subprocess
import re, os, glob, shutil, socket, time
from os.path import getsize
from scipy import *
import Struct,WAN90,Fileio
from scipy import linalg
from mysub import *
from mpi4py import MPI
import MATRIX_util
import fort_kpt_tools as fkpt
import copy

###########################################################################
#### This program executes VASP+DMFT using CTQMC impurity solver ##########
#### This part controls the overall self-consistent steps        ##########
############################################  Dr. Hyowon Park    ##########

         
def HF_eigsystem(KPTS,TOT_NKPTS,NBANDMAX,EIGVAL,BAND_WIN,WANU):
   """ This function perform the k-sum of HF Self-consistent calculation
   """
   evalHF=zeros((len(KPTS),NBANDMAX),dtype=float)
#   ncor_orb=len(SigMdc)

   for ikp,kp in enumerate(KPTS):
      if kp<TOT_NKPTS+1: 
         nbmin=BAND_WIN[ikp][0]; nbmax=BAND_WIN[ikp][1]
         num_band_max=nbmax-nbmin+1
         evalHF[ikp,:num_band_max]=array(EIGVAL[ikp][nbmin:nbmax+1])
   return evalHF
   

def Ksum_HF(KPTS,TOT_NKPTS,NBANDMAX,NWANN,EIGVAL,BAND_WIN,WANU,evalHF,KWEIGHT):
   """Compute dm and EKIN"""
   DM=zeros(NWANN,dtype=float)
   EKIN=0.0
   for ikp,kp in enumerate(KPTS):
      if kp<TOT_NKPTS+1:
         nbmin=BAND_WIN[ikp][0]; nbmax=BAND_WIN[ikp][1]
         num_band_max=nbmax-nbmin+1
         Uwanband=matrix(WANU[ikp][:NWANN,:num_band_max])
         EKIN+=2*dot(EIGVAL[ikp][nbmin:nbmax+1],KWEIGHT[ikp][:num_band_max]).real
         #EKIN+=dot(EIGVAL[ikp][nbmin:nbmax+1],diag(Uwanband2*DMFTW_perk*Uwanband2.H)).real
         DM+=2*diag(Uwanband*diag(KWEIGHT[ikp][:num_band_max])*Uwanband.H).real
   return DM,EKIN       
   
def Compute_totN_HF2(KPTS,TOT_NKPTS,BAND_WIN,KWEIGHT):
   totN=0.
   for ikp,kp in enumerate(KPTS):
      if kp<TOT_NKPTS+1:
         nbmin=BAND_WIN[ikp][0]; nbmax=BAND_WIN[ikp][1]
         num_band_max=nbmax-nbmin+1
         for iband in range(num_band_max):
            totN+=KWEIGHT[ikp,iband]
   return totN


def Get_Tetra(recip_latt,nk,ikpt):
   vec=transpose(recip_latt.copy())
   verts=((0,0,0),(1,0,0),(0,1,0),(1,1,0),(0,0,1),(1,0,1),(0,1,1),(1,1,1))
   diags=((0,7),(1,6),(4,3))
   d=sum(dot(vec,diags[0][0]-diags[0][1])**2)
   # find the smallest diagonal...
   mdiag=sorted([(sum(dot(vec,array(verts[ss[0]])-array(verts[ss[1]]))**2),s) for s,ss in enumerate(diags)])[0][-1]
   mdiag=(verts[diags[mdiag][0]],verts[diags[mdiag][1]])
   # get the primitive tetrahedra...
   temp=[s for s in verts if s!=mdiag[0] and s!=mdiag[1]]
   ptet=[]
   for i,ii in enumerate(temp):
     for j,jj in enumerate(temp[i:]):
       if sum((array(ii)-array(jj))**2)==1:ptet.append(sorted([mdiag[0],mdiag[1],ii,jj]))
   nk=tuple(nk)
   ntet=nk[0]*nk[1]*nk[2]*6
   sym=[array([[1, 0, 0],[0, 1, 0],[0, 0, 1]])]
   ki=[(i,j,k) for k in range(nk[2]+1) for j in range(nk[1]+1) for i in range(nk[0]+1)]
   kptr,gptr=fkpt.get_kptr(nk,sym,ki)
   kptr=[tuple(s) for s in kptr]
   kibz=sorted(set(kptr))
   tet_idx=[]
   for kk in kibz:
      for ik,ikk in enumerate(ikpt):
         if kk[0]==ikk[0] and kk[1]==ikk[1] and kk[2]==ikk[2]:
            tet_idx.append(ik)
            break
   tetkptr,tet=fkpt.get_itet(nk,ptet,kibz,kptr,ntet)

   return ntet,kibz,tetkptr,tet,gptr,tet_idx

def Compute_TETRA(ef,TOT_NKPTS,NBANDMAX,tetkptr,tet_idx,evalHF_reciv,TET_KPTS,ntet):
   WEIGHT=zeros((TOT_NKPTS,NBANDMAX),dtype=float)
   for ikp,kp in enumerate(TET_KPTS):
      if kp<ntet+1:
         tetk=tetkptr[ikp]
         for b in range(NBANDMAX):
            E1=evalHF_reciv[tet_idx[tetk[0]-1],b]
            E2=evalHF_reciv[tet_idx[tetk[1]-1],b]
            E3=evalHF_reciv[tet_idx[tetk[2]-1],b]
            E4=evalHF_reciv[tet_idx[tetk[3]-1],b]
            Elist=[[E1,0],[E2,1],[E3,2],[E4,3]]
            Eorder=array(sorted(Elist))[:,0]
            order=map(int,array(sorted(Elist))[:,1])
            i1=tet_idx[tetk[order[0]]-1]
            i2=tet_idx[tetk[order[1]]-1]
            i3=tet_idx[tetk[order[2]]-1]
            i4=tet_idx[tetk[order[3]]-1]
            e1=Eorder[0];e2=Eorder[1];e3=Eorder[2];e4=Eorder[3]
            cw1=0;cw2=0;cw3=0;cw4=0
            if ef>e1 and ef<e2:
               c=0.25*(ef-e1)**3/((e2-e1)*(e3-e1)*(e4-e1))
               cw1=c*(4-(ef-e1)*(1/(e2-e1)+1/(e3-e1)+1/(e4-e1)))
               cw2=c*(ef-e1)/(e2-e1)
               cw3=c*(ef-e1)/(e3-e1)
               cw4=c*(ef-e1)/(e4-e1)
            elif ef>e2 and ef<e3:
               c1=0.25*(ef-e1)**2/((e4-e1)*(e3-e1))
               c2=0.25*((ef-e1)*(ef-e2)*(e3-ef))/((e4-e1)*(e3-e2)*(e3-e1))
               c3=0.25*(ef-e2)**2*(e4-ef)/((e4-e2)*(e3-e2)*(e4-e1))
               cw1=c1+(c1+c2)*(e3-ef)/(e3-e1)+(c1+c2+c3)*(e4-ef)/(e4-e1)
               cw2=c1+c2+c3+(c2+c3)*(e3-ef)/(e3-e2)+c3*(e4-ef)/(e4-e2)
               cw3=(c1+c2)*(ef-e1)/(e3-e1)+(c2+c3)*(ef-e2)/(e3-e2)
               cw4=(c1+c2+c3)*(ef-e1)/(e4-e1)+c3*(ef-e2)/(e4-e2)
            elif ef>e3 and ef<e4:
               c=0.25*(e4-ef)**3/((e4-e1)*(e4-e2)*(e4-e3))
               cw1=0.25-c*(e4-ef)/(e4-e1)
               cw2=0.25-c*(e4-ef)/(e4-e2)
               cw3=0.25-c*(e4-ef)/(e4-e3)
               cw4=0.25-c*(4-(1/(e4-e1)+1/(e4-e2)+1/(e4-e3))*(e4-ef))
            elif ef>e4:
               cw1=0.25;cw2=0.25;cw3=0.25;cw4=0.25
            WEIGHT[i1,b]+=cw1/6
            WEIGHT[i2,b]+=cw2/6
            WEIGHT[i3,b]+=cw3/6
            WEIGHT[i4,b]+=cw4/6
   return WEIGHT   

def Compute_EPOT(DM,MOM,ncor_orb,norb,atm_idx,U,Up,Jh):
   EPOT=0.0
   Nd=zeros((len(atm_idx)))
   for i,ats in enumerate(atm_idx):
      fi=open('UC'+str(ats)+'.dat','r')
      UC=[]
      for line in fi.readlines():
         UC.append(map(float,line.split()))
      UC=array(UC)+U[i]-diag(ones(len(UC))*U[i])
      OCC=array(list((DM[norb*i:norb*(i+1)]+MOM[norb*i:norb*(i+1)])/2)+list((DM[norb*i:norb*(i+1)]-MOM[norb*i:norb*(i+1)])/2))
      Nd[i]=sum(DM[norb*i:norb*(i+1)])
      #print shape(UC), shape(OCC)
      EPOT+=0.5*dot(OCC,dot(UC,OCC))
      EPOT-=(Up[i]*Nd[i]*(Nd[i]-1)/2.0-Jh[i]*Nd[i]*(Nd[i]-2)/4.0)
   return Nd,EPOT

if __name__ == '__main__':
   
   comm = MPI.COMM_WORLD
   size = comm.Get_size()
   rank = comm.Get_rank()
   name = MPI.Get_processor_name()

   print "Hello, World! I am process %d of %d on %s." % (rank, size, name)

   HF_iter=open('INFO_KSUM','a')
   #DM_iter=open('DMINFO','a')

   fi=open('DFT_mu.out','r')
   mu=float(fi.readline().split()[0])
   fi.close()
#   if (os.path.exists('SigMdc.out')):
#      SigMdc=Fileio.Read_float('SigMdc.out')
#   else: print "SigMdc.out is missing"; exit()

   if (os.path.exists('dmft_params.dat')):
      fi=open('dmft_params.dat','r')
      fi.readline()
      fi.readline()
      fi.readline()
      n_tot=float(fi.readline().strip())
      fi.readline()
      fi.readline()
      fi.readline()
      mu_iter=int(fi.readline().strip())
      fi.readline()
      fi.readline()
      fi.readline()
      fi.readline()
      fi.readline()
      n_orbs=int(fi.readline().strip())
   else: print "dmft_params.dat is missing"; exit()

   #if nspin!=2: print "nspin must be 2 for HF calculations"; exit()
#   if ncor_orb!=len(SigMdc): print "ncor_orb is not consistent with SigMdc"; exit()
#   if ncor_orb/norb!=len(U): print "ncor_orb is not consistent with U"; exit()

   execfile('INPUT.py')
   TB=Struct.TBstructure('POSCAR',p['atomnames'],p['orbs'])
   cor_at=p['cor_at']; cor_orb=p['cor_orb']
   Natom=0
   for ats in cor_at:
      Natom+=len(ats)
   #TB.Compute_cor_idx(cor_at,cor_orb)

   if rank==0:
      WAN=WAN90.WANNIER('wannier90')
      WAN.load_chk(p['path_bin'])
      WAN.load_eig()
      WAN.Compute_HamR0()
      ntet,kibz,tetkptr,tet,gptr,tet_idx=Get_Tetra(WAN.recip_latt,WAN.mp_grid,WAN.ikpt)
      ed=[]
      for i,ats in enumerate(cor_at):
         ed.append([])
         for j,orbs in enumerate(cor_orb[i]):
            ed_loc=0
            for at in ats:
               for orb in orbs:
                  ed_loc+=WAN.HamR0[TB.idx[at][orb],TB.idx[at][orb]].real
            ed[i].append(ed_loc/len(ats)/len(orbs))
      #savetxt('ED.out',ed)

      NK_proc=WAN.num_kpts//size
      if WAN.num_kpts%size>0: NK_proc+=1
      KPTS_send=arange(NK_proc*size,dtype=int)+1
      TOT_NKPTS=WAN.num_kpts
      assert len(KPTS_send)>=TOT_NKPTS
      TET_NK_proc=ntet//size
      if ntet%size>0: TET_NK_proc+=1
      TET_KPTS_send=arange(TET_NK_proc*size,dtype=int)+1
      assert len(TET_KPTS_send)>=ntet
      tetkptr_send=zeros((len(TET_KPTS_send),4),dtype=int)
      tetkptr_send[:ntet,:]=array(tetkptr)[:,:]
      KVEC_send=zeros((len(KPTS_send),3),dtype=float)
      KVEC_send[:TOT_NKPTS,:]=WAN.kpt_latt[:,:]
      NBANDMAX=WAN.num_band_max
      NBANDS=WAN.num_bands
      NWANN=WAN.num_wann
      BAND_WIN_send=zeros((len(KPTS_send),2),dtype=int)
      BAND_WIN_send[:TOT_NKPTS,:]=array(WAN.band_win[:,:])
      KWEIGHT_send=zeros((len(KPTS_send),NBANDMAX),dtype=float)
      KWEIGHT_dn_send=zeros((len(KPTS_send),NBANDMAX),dtype=float)
      #EIGVALTET_send=zeros((len(KPTS_send),4,NBANDS),dtype=float)
      EIGVAL_send=zeros((len(KPTS_send),NBANDS),dtype=float)
      EIGVAL_send[:TOT_NKPTS,:]=array(WAN.eigvals[:,:])
      WANU_send=zeros((len(KPTS_send),NWANN,NBANDMAX),dtype=complex)
      WANU_send[:TOT_NKPTS,:,:]=array(WAN.WANU[:,:,:])
      
   else: 
      ntet=None; tet_idx=None; TET_NK_proc=None; TET_KPTS_send=None; tetkptr_send=None; KVEC_send=None; NK_proc=None; NKTET_proc=None; KPTS_send=None; NWANN=None; NBANDS=None; TOT_NKPTS=None; NBANDMAX=None; BAND_WIN_send=None; NBANDMAX_send=None; #EIGVALTET_send=None;  
      EIGVAL_send=None; WANU_send=None; KWEIGHT_send=None;KWEIGHT_dn_send=None; 
      
   ######## Broadcast data for mpi #####
   NK_proc = comm.bcast(NK_proc, root=0)   
   NBANDS = comm.bcast(NBANDS, root=0)   
   NWANN = comm.bcast(NWANN, root=0)   
   TOT_NKPTS = comm.bcast(TOT_NKPTS, root=0)   
   NBANDMAX = comm.bcast(NBANDMAX, root=0)
   TET_NK_proc = comm.bcast(TET_NK_proc, root=0)   
   tet_idx = comm.bcast(tet_idx, root=0)
   ntet = comm.bcast(ntet, root=0)
   evalHF_reciv=zeros((NK_proc*size,NBANDMAX),dtype=float)
   evalHF_dn_reciv=zeros((NK_proc*size,NBANDMAX),dtype=float)
   ######## Scatter data for mpi #####
   KPTS = zeros(NK_proc,dtype=int)
   BAND_WIN = zeros((NK_proc,2),dtype=int)
   EIGVAL = zeros((NK_proc,NBANDS),dtype=float)
   KWEIGHT = zeros((NK_proc,NBANDMAX),dtype=float)
   KWEIGHT_dn = zeros((NK_proc,NBANDMAX),dtype=float)
   WANU = zeros((NK_proc,NWANN,NBANDMAX),dtype=complex)
   KVEC = zeros((NK_proc,3),dtype=float)
   comm.Scatter(KPTS_send, KPTS, root=0)
   comm.Scatter(BAND_WIN_send, BAND_WIN, root=0)
   comm.Scatter(EIGVAL_send, EIGVAL, root=0)
   comm.Scatter(WANU_send, WANU, root=0)
   comm.Scatter(KVEC_send, KVEC, root=0)
   tetkptr=zeros((TET_NK_proc,4),dtype=int)
   comm.Scatter(tetkptr_send, tetkptr, root=0)
   TET_KPTS=zeros(TET_NK_proc,dtype=int)
   comm.Scatter(TET_KPTS_send, TET_KPTS, root=0)
         

   evalHF = HF_eigsystem(KPTS,TOT_NKPTS,NBANDMAX,EIGVAL,BAND_WIN,WANU)

   comm.Allgather(evalHF, evalHF_reciv)

   ### Compute mu ###
   low_mu=mu;high_mu=mu
   LOWB=False;HIGHB=False
   for i in range(mu_iter):
      mu=(low_mu+high_mu)/2.0

      TET_KWEIGHT_sum=zeros((TOT_NKPTS,NBANDMAX),dtype=float)
      TET_KWEIGHT=zeros((TOT_NKPTS,NBANDMAX),dtype=float)
      TET_KWEIGHT=Compute_TETRA(mu,TOT_NKPTS,NBANDMAX,tetkptr,tet_idx,evalHF_reciv,TET_KPTS,ntet)
      comm.Reduce(TET_KWEIGHT,TET_KWEIGHT_sum,op=MPI.SUM)
      if rank==0: KWEIGHT_send[:TOT_NKPTS,:]=TET_KWEIGHT_sum[:,:]

      comm.Scatter(KWEIGHT_send, KWEIGHT, root=0)
      totN=Compute_totN_HF2(KPTS,TOT_NKPTS,BAND_WIN,KWEIGHT)
      totN=comm.allreduce(totN,op=MPI.SUM)
      totN=2*totN/TOT_NKPTS
      
      if abs(totN-n_tot)<1e-10: break
      if totN<n_tot:
         if HIGHB==False:
            high_mu=mu+(n_tot-totN)/(0.1*n_tot)
         low_mu=mu
         LOWB=True
      else:
         if LOWB==False:
            low_mu=mu+(n_tot-totN)/(0.1*n_tot)
         high_mu=mu
         HIGHB=True
      if rank==0: print mu, totN
#   savetxt('DMFT_mu.out',[mu])

   DM,EKIN = Ksum_HF(KPTS,TOT_NKPTS,NBANDMAX,NWANN,EIGVAL,BAND_WIN,WANU,evalHF,KWEIGHT)
   DM_sum=zeros(shape(DM))
   EKIN=comm.allreduce(EKIN,op=MPI.SUM)
   comm.Allreduce(DM,DM_sum,op=MPI.SUM)
   
   if rank==0:
      EKIN=EKIN/TOT_NKPTS
      DM_sum=DM_sum/TOT_NKPTS
#      MOM_sum=zeros(shape(DM_sum))
      #print DM_sum, MOM_sum; exit()
#      Nd,EPOT=Compute_EPOT(DM_sum,MOM_sum,ncor_orb,norb,atm_idx,U,Up,Jh)
      
      HF_iter.write( '%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f' % (mu, totN, sum(DM_sum[:n_orbs]), sum(DM_sum[n_orbs*(Natom-1):n_orbs*Natom]), EKIN, 0.0) )
      HF_iter.write('\n')
      HF_iter.flush()

      #for i in range(ncor_orb):
      #   DM_iter.write( '%18.12f' % (DM_sum[i]) )
      #for i in range(ncor_orb):
      #   DM_iter.write( '%18.12f' % (MOM_sum[i]) )
      #DM_iter.write('\n')
      #DM_iter.flush()
