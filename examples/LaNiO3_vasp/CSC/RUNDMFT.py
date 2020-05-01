#!/usr/bin/env python

import sys, subprocess
import re, os, glob, shutil, socket, time
from os.path import getsize
import Struct,Fileio,VASP, DMFT_MOD, IMP_SOLVER
from scipy import *
import shutil,copy,scipy.interpolate

###########################################################################
#### This program executes VASP+DMFT using CTQMC impurity solver ##########
#### This part controls the overall self-consistent steps        ##########
############################################  Dr. Hyowon Park    ##########

def now(): return time.strftime(' at %c') #%Y:%H HH:%M MM:%S SS')


def CreateINCAR(params_vasp):
   " Creates input file (INCAR) for vasp"
   f = open('INCAR', 'w')
   for p in params_vasp:
       print >> f, p, '  ', params_vasp[p][0], '\t', params_vasp[p][1]
   f.close()   



if __name__ == '__main__':
   
   execfile('INPUT.py') # Read input file

   main_out=open('INFO_TIME','w')
   main_iter=open('INFO_ITER','w')
   DMFT_iter=open('INFO_KSUM','w')
   DM_iter=open('INFO_DM','w')
   E_iter=open('INFO_ENERGY','w')
   DFT_iter=open('INFO_DFT_loop','w')


   if os.path.exists("para_com.dat"):
      fipa=open('para_com.dat','r')
      para_com=str(fipa.readline())[:-1]
      fipa.close()
   else:
      para_com=""
   if os.path.exists("para_com2.dat"):
      fipa=open('para_com2.dat','r')
      para_com2=str(fipa.readline())[:-1]
      fipa.close()
   else:
      para_com2=""

 ############ Initial Preparation ########################

   if p['Nit']>0 and p['Niter']>1: main_out.write( '-----------Charge Self-consistent DFT+DMFT calculations-----------' )
   if p['Nit']>0 and p['Niter']==1: main_out.write( '----------Non-charge Self-consistent DFT+DMFT calculations-----------' )
   main_out.write('\n')
   main_out.flush()

   main_out.write( 'Caculation Starts'+now() )
   main_out.write('\n')
   main_out.flush()

   DMFT_iter.write( '%12s %12s %12s %12s %12s %12s %12s' % ('# mu','TOT_ELEC','Nd[0]','Nd[-1]','EKIN','Sigoo-Vdc','<d_epsk>') )
   DMFT_iter.write('\n')
   DMFT_iter.flush()

   DM_iter.write( '%12s' % ('# Occupancy') )
   DM_iter.write('\n')
   DM_iter.flush()

   main_iter.write( '%3s %6s %8s %8s %16s %16s %16s %16s' % ('#','# DMFT','Nd_latt','Nd_imp','(Sigoo-Vdc)_latt','(Sigoo-Vdc)_imp','TOT_E(Tr(SigG))','TOT_E(EPOT_imp)') )
   main_iter.write('\n')
   main_iter.flush()

   E_iter.write( '%12s %18s %18s %12s' % ('# DFT_E','<epsk>-<epsk>_o','1/2*Tr(Sigma*G)','EPOT_imp'))
   E_iter.write('\n')
   E_iter.flush()

   DFT_iter.write( '%10s %10s %10s %10s' % ('# mu','TOT_ELEC','Nd[0]','EKIN') )
   DFT_iter.write('\n')
   DFT_iter.flush()

   CHGDIFF=0.0 
   #### Read POSCAR ########
   TB=Struct.TBstructure('POSCAR',p['atomnames'],p['orbs'])
   cor_at=p['cor_at']; cor_orb=p['cor_orb']
   TB.Compute_cor_idx(cor_at,cor_orb)
   print TB.TB_orbs
   if TB.LHF=='.TRUE.': p['Nd_qmc']=0
   U=p['U'];J=p['J']
   T=1.0/pC['beta'][0]
   noms=p['noms']
#   T_high=p['T_high']
#   if T_high<1.0/pC['beta'][0]: T_high=1.0/pC['beta'][0]
#   noms_high=int(p['noms']/(T_high*pC['beta'][0]))
   ########### Create symmetry index ################
   N_atoms=0
   for i,ats in enumerate(cor_at):
      N_atoms+=len(ats)
   atm_idx=zeros(N_atoms,dtype=int)
   loc_idx=-1
   for i,ats in enumerate(cor_at):
      loc_idx+=1
      for at in ats:
         at_idx=int(at[-1])-1
         atm_idx[at_idx]=loc_idx
   #print atm_idx
   sym_idx=[]
   loc_idx=0
   for i,ats in enumerate(cor_at):
      d_orb=TB.TB_orbs[ats[0]]
      idx=zeros(len(d_orb),dtype=int)
      for j,orbs in enumerate(cor_orb[i]):
         loc_idx+=1
         for orb in orbs:
            idx[d_orb.index(orb)]=loc_idx
      sym_idx.append(idx.tolist())
   #print sym_idx
#      for at in ats:
#   for i,ats in enumerate(cor_at):
#      for ii,at in enumerate(ats):
#         at_idx=int(at[-1])-1
#         d_orb=TB.TB_orbs[at]
#         for j,orbs in enumerate(cor_orb[i]):
#            if ii==0: loc_idx+=1
#            for orb in orbs:
#               idx[d_orb.index(orb)]=loc_idx
#         sym_idx[at_idx].append(idx)

   DMFT=DMFT_MOD.DMFT_class(p,pC,TB)
   DFT=VASP.VASP_class()

   ETOT_old=0.0;ETOT2_old=0.0;ETOT=0.0;ETOT2=0.0
   CHGDIFF=0.;CHGDIFF2=0.
   shutil.copy2('sig.inp','sig.inp.0')

   for itt in range(p['Niter']):
      main_out.write( '--- Starting Charge loop '+str(itt+1)+now()+'---' )
      main_out.write('\n')
      main_out.flush()

      #if itt<p['Niter']-p['Nforce']: Fileio.Create_INPUT(p,pC,TB,T,noms,'.FALSE.')
      #else: Fileio.Create_INPUT(p,pC,TB,T,noms,'.TRUE.')
      Fileio.Create_dmft_params(p,pC,N_atoms,atm_idx,sym_idx)

      for it in range(p['Nit']):
         main_out.write( '--- Starting DMFT loop '+str(it+1)+now()+'---' )
         main_out.write('\n')
         main_out.flush()
         
         if itt==0 and it==0:
            for i,ats in enumerate(cor_at):
               if p['nspin']==2: pD['para=']=0 
               else: pD['para=']=1
               if p['orbs'][0]=='f': pD['l=']=3
               elif p['orbs'][0]=='d': pD['l=']=2
               pD['J=']=float(J[i])
               pD['Eimp=']=zeros(p['nspin']*len(cor_orb[i]))#array(ed[i])-ed[i][0]+array(DMFT.sig_st[i])-DMFT.sig_st[i][0]
               IMP_SOLVER.Create_atomd(pD)
               IMP_SOLVER.Create_Trans(TB.ncor_orb_list[i],p['nspin'],ats[0],cor_orb[i],TB)
               cmd = 'cp Trans.dat Trans'+str(i+1)+'.dat'
               print os.popen(cmd).read()
               cmd = p['path_bin']+"atom_d.py 'atom_d.inp' 2> atom_d.error || { echo 'ATOMD run failed!'; exit 1; }"
               #cmd = "./atom_d.py 'atom_d.inp' 2> atom_d.error || { echo 'ATOMD run failed!'; exit 1; }"
               out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
               print out #, err
               cmd = 'cp UC.dat UC'+str(i+1)+'.dat'
               print os.popen(cmd).read()
#            DMFT.Read_Sig(TB,p['nspin'])
#            #print DMFT.Ekin,DMFT.Eimp,DMFT.imp_mu
#            DMFT.Nd_latt=copy.deepcopy(DMFT.Nd_imp)
#            DMFT.Compute_HF(1,p,TB)
#            DMFT.Cmp_Sig_highT(T_high,noms_high,p['nspin'])
#            #DMFT.Mix_Sig(p['mix_sig'])
#            if TB.LHF=='.FALSE.':
#               Fileio.Print_complex_multilines(DMFT.Sig,DMFT.ommesh,'SigMoo.out')
#               if p['nspin']==2: Fileio.Print_complex_multilines(DMFT.Sig_dn,DMFT.ommesh,'SigMoo_dn.out')
#               Fileio.Print_complex_multilines(DMFT.Sig_highT,DMFT.ommesh_highT,'SigMoo_highT.out')
#               if p['nspin']==2: Fileio.Print_complex_multilines(DMFT.Sig_dn_highT,DMFT.ommesh_highT,'SigMoo_dn_highT.out')
#
#            if (os.path.exists('SigMdc.out')):
#               pass
#            else:
#               fi=open('SigMdc.out','w')
#               for data in DMFT.SigMdc:
#                  print >>fi, data,
#               print >>fi, ''
#               fi.close()
#            if p['nspin']==2: 
#               if (os.path.exists('SigMdc_dn.out')):
#                  pass
#               else:
#                  fi=open('SigMdc_dn.out','w')
#                  for data in DMFT.SigMdc_dn:
#                     print >>fi, data,
#                  print >>fi, ''
#                  fi.close()
#   
   
         if it==0:
            DMFT.EKIN0=0
            cmd = para_com+" "+p['path_bin']+"XHF0.py > ksum_output 2> ksum_error"
            out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
            fiEKIN=open('INFO_KSUM','r')
            lines=fiEKIN.readlines()
            DMFT.EKIN0=float(lines[-1].split()[4])
            fiEKIN.close()

         if TB.LHF=='.FALSE.':
            #cmd = para_com+" "+p['path_bin']+"dmft_ksum_sp > ksum_output 2> ksum_error"
            cmd = para_com+" "+p['path_bin']+"dmft.x > ksum_output 2> ksum_error"
            out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
         else: 
            cmd = para_com+" "+p['path_bin']+"XHF.py > ksum_output 2> ksum_error"
            out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
         print out #, err

         fi=open('DMFT_mu.out','r')
         DMFT.mu=float(fi.readline())
         fi.close()

         DMFT.Update_Sigoo_and_Vdc(TB,'sig.inp')
         DMFT.Update_Nlatt(TB,p['nspin'])
   
#         if TB.LHF=='.FALSE.': DMFT_mu,ed,TrDeltaG=DMFT.Compute_Delta(DMFT.T,p['nspin'],cor_at,cor_orb,TB,p['nom'])
#         
#         if itt==0 and it==0:
#            fi=open('VASP_mu.out','w')
#            print >>fi, DMFT_mu      
#            fi.close()
   
         main_out.write( '--- DMFT Ksum is finished '+now()+'---' )
         main_out.write('\n')
         main_out.flush()
   
   
         #for i in range(len(cor_at)):
         #   DMFT.sig_st[i][:]+=0.5*(DMFT.Nd_imp[i]-DMFT.Nd_latt[i])
         #DMFT.Impurity_Solver
         #if itt%3==0:

         om_loc,Delta=Fileio.Read_complex_multilines('Delta.out')
         loc_idx=-1
         for i,ats in enumerate(cor_at):
            Delta_loc=zeros((len(cor_orb[i]),len(om_loc)),dtype=complex)
            for j,orbs in enumerate(cor_orb[i]):
               loc_idx+=1
               Delta_loc[j,:]=copy.deepcopy(Delta[loc_idx,:])
               Fileio.Print_complex_multilines(Delta_loc,om_loc,'Delta'+str(i+1)+'.inp')
       # shutil.copy2('Delta.out','Delta1.inp')
         Ed=array([loadtxt('Ed.out')])
         if TB.LHF=='.FALSE.': IMP_SOLVER.RUN_CTQMC(p,pC,pD,it,itt,para_com,DMFT.mu,Ed,DMFT.Vdc)
         main_out.write( '--- Impurity solver is finished '+now()+'---' )
         main_out.write('\n')
         main_out.flush()

################# Mix sig.inp ###############################3

         DMFT.Read_Sig(TB,p['nspin'])
         DMFT.Compute_Energy(DFT,TB,Ed)
         DMFT.Compute_Sigoo_and_Vdc(p,TB)
         DMFT.Mix_Sig_and_Print_sig_inp(TB,p['Nd_qmc'],p['mix_sig'],'sig.inp')
         shutil.copy2('sig.inp','sig.inp.'+str(itt+1)+'.'+str(it+1))
         shutil.copy2('G_loc.out','G_loc.out.'+str(itt+1)+'.'+str(it+1)) 

         DMFT.ETOT=DFT.E-DMFT.EKIN0+DMFT.EKIN+DMFT.EPOT-DMFT.Edc
         DMFT.ETOT2=DFT.E-DMFT.EKIN0+DMFT.EKIN+DMFT.EPOT2-DMFT.Edc_imp
         # Print out results
         main_iter.write( '%3d %3d %12.6f %12.6f %14.6f %14.6f %16.6f %16.6f %16.6f' %(itt+1,it+1,DMFT.Nd_latt[0],DMFT.Nd_imp[0],DMFT.Sigoo[0][0]-DMFT.Vdc[0],DMFT.Sigoo_imp[0][0]-DMFT.Vdc_imp[0],DMFT.ETOT, DMFT.ETOT2, CHGDIFF) )
         main_iter.write('\n')
         main_iter.flush()

         E_iter.write( '%14.6f %14.6f %14.6f %14.6f' %(DFT.E, DMFT.EKIN-DMFT.EKIN0, DMFT.EPOT-DMFT.Edc, DMFT.EPOT2-DMFT.Edc_imp) )
         E_iter.write('\n')
         E_iter.flush()
 
#      n_avg=int(p['Nit'])
#      fiDMFT=open('ENERGY','r')
#      Eline=fiDMFT.readlines()[-n_avg:]
#      EAVG=0.; EAVG2=0.; 
#      for i in range(n_avg):
#         EAVG+=float(Eline[i].split()[2])
#         EAVG2+=float(Eline[i].split()[4])
#      fiDMFT.close()
#      main_iter.write( '--------------------------\n' )
#      main_iter.write( '%14.6f %14.6f %10.6f ' %(EAVG/n_avg, EAVG2/n_avg, CHGDIFF) )
#      main_iter.write('\n')
#      main_iter.write( '--------------------------\n' )
#      main_iter.flush()
#      E_iter.write( '%14.6f %14.6f' %(EAVG/n_avg, EAVG2/n_avg) )
#      E_iter.write('\n')
#      E_iter.flush()


##################### CHARGE UPDATE ################
#      if p['Nrelax']>p['Nforce']: print "Nrelax should not be larger than Nforce"; exit()
#      if itt>=p['Niter']-p['Nrelax']:
#         TB.Update_POSCAR('POSCAR')
#         pV['NSW= ']=[5,'# NSW']
#         pV['IBRION= ']=[1,'# IBRION']
#         pV['ISIF= ']=[2,'# ISIF']
#         pV['EDIFFG= ']=[-0.01,'# EDIFFG']
#
#
      if itt<p['Niter']-1:
         main_out.write( '--- Running vaspCHG '+now()+'---' )
         main_out.write('\n')
         main_out.flush()
   
   #      #if itt>0:
   #      #   cmd = para_com+" "+p['path_bin']+"XHF_fix_Nd.py > ksum_output 2> ksum_error"
   #      #   out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
   #
   #         if pV.keys().count('NBANDS='): DFT.NBANDS=pV['NBANDS='][0]
   #      DFT.Create_win(TB,p['atomnames'],p['orbs'],p['L_rot'],DFT.NBANDS,DFT.EFERMI+p['ewin'][0],DFT.EFERMI+p['ewin'][1])
   #      if os.path.exists('CHGCAR') and itt==0:
   #         CreateINCAR(pV)
         if itt==0:
            f = open('INCAR', 'a')
   #         #print >> f, "LDMFT= .TRUE."
            print >> f, "NELM= "+str(p['Ndft'])
            f.close()

         cmd = para_com+" "+p['path_bin']+"vaspDMFT > vasp.out 2> vasp.error || { echo 'Parallel run failed!'; exit 1; }"
         out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
         shutil.copy2('CHGCAR','CHGCAR.'+str(itt))
         #shutil.copy2('OSZICAR','OSZICAR.'+str(itt))
         if itt>0:
            DFT.Read_NELECT()
            CHGDIFF=DFT.Diff_CHGCAR('CHGCAR.'+str(itt-1),'CHGCAR.'+str(itt))
         DFT.Read_NBANDS()
         DFT.Read_EFERMI()
         DFT.Update_win(DFT.NBANDS,DFT.EFERMI+p['ewin'][0],DFT.EFERMI+p['ewin'][1])
   #      if itt==0:
   #
   #      if itt>0:# or (not os.path.exists('CHGCAR')):
         print os.popen("rm wannier90.chk").read()
         print os.popen("rm wannier90.chk.fmt").read()
   ##         print os.popen("rm wannier90.eig").read()
   #
   #      #if itt<p['Niter']-1:
         main_out.write( '-------------- Running wannier 90 '+str(itt+1)+'----------------' )
         main_out.write('\n')
         main_out.flush()
         cmd = para_com+" "+p['path_bin']+"wannier90.x wannier90"
         out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
         print out #, err

   main_out.write( 'Caculation Ends'+now() )
   main_out.write('\n')
   main_out.flush()


