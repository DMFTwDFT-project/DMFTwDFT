from scipy import *
import os,sys,copy

class TBstructure:
   def __init__(self,poscar,atomnames,orbs):
      self.orb_idx={}
      self.orb_idx['s']=['s']
      self.orb_idx['p']=['p_z','p_x','p_y']
      self.orb_idx['d']=['d_z2','d_xz','d_yz','d_x2y2','d_xy']
#     self.orb_idx['f']=['f_z3','f_xz2','f_yz2','f_z(x2-y2)','f_xyz','f_x(x2-3y2)']
      self.orb_idx['f']=['f1','f2','f3','f4','f5','f6','f7']
      self.L={}
      self.L['s']=0;self.L['p']=1;self.L['d']=2;self.L['f']=3
      fi=open(poscar,'r')
      self.title=fi.readline().strip()
      self.scale=float(fi.readline().strip())
      self.vec=zeros((3,3),float)
      for i in range(3):self.vec[i]=array([float(s) for s in fi.readline().split()])
      self.vec=self.vec*self.scale
      linenum=0;lines=[]
      for i,line in enumerate(fi.readlines()):
         if line=='\n': break
         if i==0: line0=line
         elif i==1: line1=line
         elif i==2: line2=line
         else: linenum+=1
         lines.append(line)
      self.atomname=[]
      if str(line1.split()[0])[0]=='D' or str(line1.split()[0])[0]=='C' or str(line1.split()[0])[0]=='d' or str(line1.split()[0])[0]=='c': print "INPUT atom name info in POSCAR!!"; exit() #self.specpoint=array([int(s) for s in line0.split()]); linenum+=1
      elif str(line2.split()[0])[0]=='D' or str(line2.split()[0])[0]=='C' or str(line2.split()[0])[0]=='d' or str(line2.split()[0])[0]=='c': self.specpoint=array([int(s) for s in line1.split()]);self.atomname=[s for s in line0.split()]
      if linenum<sum(self.specpoint): print "Number of atoms for poscar is not right!"; exit()
      self.nspecies=len(self.specpoint)
      #if self.nspecies!=len(atomnames): print "Number of atoms for poscar and for INPUT.py is not consistent!"; exit()
      self.coord=line2.strip().lower()[0]
      self.at=zeros((sum(self.specpoint),3),float)
      for i in range(sum(self.specpoint)): self.at[i]=array([float(s) for s in lines[i+3].split()[:3]])
      if self.coord=='c':self.atc=self.at*self.scale; self.atd=dot(self.atc,inv(self.vec)) ; self.at=self.atc
      elif self.coord=='d':self.atd=self.at; self.atc=dot(self.atd,self.vec)
      else: print "ERROR: coordinates must be direct or cartesian";sys.exit()
      self.scale=1.0 # the vec and atc were rescaled above... set scale to 1...
      #print self.atc;exit()

      self.TB_orbs={}
      self.at_idx={}
      self.num_atoms={}
      self.names=[]
      self.at_list=[]
      for ia,atname in enumerate(atomnames):
         chkpoint=0
         idxidx=0
         for ispec,nspec in enumerate(self.specpoint):
            if self.atomname[ispec]==atname: 
               self.num_atoms[atname]=nspec
               for i in range(nspec):
                  at_name=atname+str(i+1)
                  self.names.append(at_name)
                  self.TB_orbs[at_name]=self.orb_idx[orbs[ia]]
                  self.at_idx[at_name]=idxidx
                  idxidx+=1
               chkpoint=1
            else: idxidx+=nspec
         if chkpoint==0: print "Couldn't find the atom name", atname, "in POSCAR"; exit()
      self.natoms=len(self.TB_orbs.keys())
      #print self.at_idx; exit()
      #print "atomic names=", self.names
      
      idx={}
      ix=0
      for at in self.names:
         orb=self.TB_orbs[at]
         if not at in idx : idx[at]={}
         for iorb in orb:
            idx[at][iorb]=ix
            ix=ix+1
      self.idx=idx
      self.Msize=ix
      #print "atomic indices=", self.idx

   def Update_POSCAR(self,poscar):
      fi=open(poscar,'r')
      self.title=fi.readline().strip()
      self.scale=float(fi.readline().strip())
      self.vec=zeros((3,3),float)
      for i in range(3):self.vec[i]=array([float(s) for s in fi.readline().split()])
      self.vec=self.vec*self.scale
      linenum=0;lines=[]
      for i,line in enumerate(fi.readlines()):
         if line=='\n': break
         if i==0: line0=line
         elif i==1: line1=line
         elif i==2: line2=line
         else: linenum+=1
         lines.append(line)
      self.atomname=[]
      if str(line1.split()[0])[0]=='D' or str(line1.split()[0])[0]=='C' or str(line1.split()[0])[0]=='d' or str(line1.split()[0])[0]=='c': print "INPUT atom name info in POSCAR!!"; exit() #self.specpoint=array([int(s) for s in line0.split()]); linenum+=1
      elif str(line2.split()[0])[0]=='D' or str(line2.split()[0])[0]=='C' or str(line2.split()[0])[0]=='d' or str(line2.split()[0])[0]=='c': self.specpoint=array([int(s) for s in line1.split()]);self.atomname=[s for s in line0.split()]
      if linenum<sum(self.specpoint): print "Number of atoms for poscar is not right!"; exit()
      self.nspecies=len(self.specpoint)
      self.coord=line2.strip().lower()[0]
      self.at=zeros((sum(self.specpoint),3),float)
      for i in range(sum(self.specpoint)): self.at[i]=array([float(s) for s in lines[i+3].split()[:3]])
      if self.coord=='c':self.atc=self.at*self.scale; self.atd=dot(self.atc,inv(self.vec)) ; self.at=self.atc
      elif self.coord=='d':self.atd=self.at; self.atc=dot(self.atd,self.vec)
      else: print "ERROR: coordinates must be direct or cartesian";sys.exit()
      self.scale=1.0 # the vec and atc were rescaled above... set scale to 1...


   def Compute_cor_idx(self,cor_at,cor_orb):
      if len(cor_at)!=len(cor_orb): print "The number of cor_at should be same as the number of cor_orb"; exit()
      self.ncor_orb=0;self.ncor_orb_list=[];self.max_cor_orb=0
      cor_at_flat=[item for sublist in cor_at for item in sublist]
      for at in cor_at_flat: self.ncor_orb+=len(self.TB_orbs[at]); self.max_cor_orb=max(self.max_cor_orb,len(self.TB_orbs[at]))
      for at in cor_at: self.ncor_orb_list.append(len(self.TB_orbs[at[0]]))
      cor_orb_flat=[]
      for i in range(len(cor_orb)):
         cor_orb_flat.append([item for sublist in cor_orb[i] for item in sublist])
   
      #ncor_orb=0 # total number of orbitals treated within correlated atoms
      #for at in cor_at_flat: ncor_orb+=len(TB.TB_orbs[at])
      self.cor_idx=zeros(self.ncor_orb,dtype=int)
   
      ######## Give the index for correlated orbitals ##########
      # 0: s,p, 1: d(HF), 2,3,4,...: d(DMFT for non-eqival. orbitals=> reads from Sig.out file)
      #Nd_atom=[];
   
      idx1=0
      self.LHF='.FALSE.'
      for i,ats in enumerate(cor_at):
         #Nd_at=0
         for at in ats:
            for orb in self.TB_orbs[at]:
               idx=self.idx[at][orb]
               if self.orb_idx['s'].count(orb)==1 or self.orb_idx['p'].count(orb)==1 :
                  print at," atom ",orb," orbital is treated by LDA!"
                  self.cor_idx[idx]=0
               elif cor_orb_flat[i].count(orb)==0 and (self.orb_idx['d'].count(orb)==1 or self.orb_idx['f'].count(orb)==1):
                  if rank==0: print at," atom ",orb," orbital is treated by HF!"
                  self.cor_idx[idx]=1
                  self.LHF='.TRUE.'
                  #Nd_at+=dm[idx]
               elif cor_orb_flat[i].count(orb)==1 and (self.orb_idx['d'].count(orb)==1 or self.orb_idx['f'].count(orb)==1):
                  if rank==0: print at," atom ",orb," orbital is treated by DMFT!"
                  for j,cor in enumerate(cor_orb[i]):
                    if cor.count(orb):
                       self.cor_idx[idx]=2+idx1+j; break
                  #Nd_at+=dm[idx]
               else: print at," atom ",orb," orbital is not supported!"; exit()
         idx1+=len(cor_orb[i])
         #Nd_atom.append(Nd_at/len(ats))

   def Rot_axis(self,rot_at):
      rot_mat={}
      cut=3.0;transx=2;transy=2;transz=2
      natoms=self.natoms; at=self.atc; v=self.vec; 
      # loop over atoms in the unit cell...
      comp=[];
      at1=self.at_idx[rot_at]
      for tx in range(-transx,transx+1):
        for ty in range(-transy,transy+1):
          for tz in range(-transz,transz+1):
            for at2 in range(sum(self.specpoint)):
              dv=at[at1]-(at[at2]+tx*v[0]+ty*v[1]+tz*v[2])
              dist=sqrt(dot(dv,dv)) ; dv=list(dv)
              # need to format dv for printing
              if 0.01<dist<cut:
                comp=comp+[[dist,dv[0],dv[1],dv[2]]]
      #comp.sort(cmp_list)
      comp.sort()
      n3=array(comp[-1][1:])/comp[-1][0]
      del comp[-1]
      comp.sort(cmp_list2)
      for cp in comp:
         n1=array(cp[1:])/cp[0]
         if abs(dot(n3,n1))<0.1: break
      n_2=cross(n3,n1)
      n_1=cross(n_2,n3)
      #print 'Atom=',rot_at
      return '%.8f,%.8f,%.8f'%(self.atd[at1][0],self.atd[at1][1],self.atd[at1][2]),'z=%.8f,%.8f,%.8f:x=%.8f,%.8f,%.8f' %(n3[0],n3[1],n3[2],n_1[0],n_1[1],n_1[2]) 

def off_matrix(angle,H0):
   dim=len(H0)
   U=Rot.Cmp_U((dim-1)/2,angle[0],angle[1],angle[2])
   HU=matrix(U)*matrix(H0)*matrix(U).T
   off_term=0.0
   for i in range(dim):
      for j in range(dim):
         if i!=j: off_term+=HU[i,j]**2
   return off_term

def cmp_list(x,y,idx=3):
   a=abs(x[idx]);b=abs(y[idx])
   return cmp(abs(x[idx]),abs(y[idx]))

def cmp_list2(x,y,idx=1):
   return cmp(abs(x[idx]),abs(y[idx]))
