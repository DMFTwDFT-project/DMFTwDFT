#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
# 

pycxx_available=True
symeig_available=False

from scipy import *
import sys, re, os
import copy
import getopt
import pickle
import glob
import scipy.weave as weave
from numpy import linalg


ROOT = os.environ.get('WIEN_DMFT_ROOT')
if ROOT is not None:
    sys.path.append( ROOT )
else:
    print >> sys.stderr, "Environment variable WIEN_DMFT_ROOT must be set!"
    print "Environment variable WIEN_DMFT_ROOT must be set!"
    sys.exit(1)

import gaunt                
if pycxx_available: import gutils
if symeig_available: import symeig


import numpy
if numpy.__version__ == '1.0.1':
    loadtxt = io.read_array
    savetxt = io.write_array
             

def union(data1, data2):
    " Takes a union of two lists"
    res = data1
    for d in data2:
        if d not in res: res.append(d)
    return res
def overlap(data1, data2):
    " Checks if there is any overlap between data1 and data2"
    for d in data2:
        if d in data1:
            return True
    return False

def compres(groups):
    loopsDone = True
    while (loopsDone):
        loopsDone = False
        for i in range(len(groups)):
            if loopsDone: break
            for j in range(i+1,len(groups)):
                if loopsDone: break
                if overlap(groups[i],groups[j]):
                    groups[i] = union(groups[i],groups[j])
                    del groups[j]
                    loopsDone = True
    for g in groups: g.sort()
    groups.sort(cmp=lambda x,y: x[0]-y[0])
    return groups

if pycxx_available:
    compress = gutils.compres
else:
    compress = compres

def cprint(fh, U):
    for i in range(shape(U)[0]):
        for j in range(shape(U)[1]):
            f = U[i,j]
            if abs(f)<1e-10: f = 0j
            print >> fh, "%7.4f %7.4f*i " % (f.real, f.imag),
        print >> fh
    print >> fh

def mprint(fh, U):
    for i in range(shape(U)[0]):
        for j in range(shape(U)[1]):
            f = U[i,j]
            if abs(f)<1e-10: f = 0
            print >> fh, "%7.4f " % f,
        print >> fh
    print >> fh


def TS2C_p():
    s2 = 1./sqrt(2.)
    T2 = [[],[],[]]
    ##   m        x      y       z   
    T2[-1+1]= [    s2, -1j*s2,     0.0]
    T2[ 0+1]= [   0.0,    0.0,     1.0]
    T2[ 1+1]= [    s2,  1j*s2,     0.0]
    return array(T2)
   
def TS2C_d():
  """ Generates transformation matrix from complex        
  spherical harmonic to cubic harmonic for l=2.

  A^{cub} = T^+ * A^{spher} * T
  
  Complex spherical harmonics are (after Varshalovich)

  order of cubics = [yz, zx, xy, x^2-y^2, 3z^2-r^2]

 Spheric harmonics are:
  
  Y   (r)=a/sqrt(2)*(x^2-y^2-2ixy)/r^2      Y   (r)=a/sqrt(2)*(2zx-2izy)/r^2
   2-2                                       2-1
   
  Y   (r)=a/sqrt(6)*(3z^2/r^2-1)
   2 0                                                 
   
  Y   (r)=-a/sqrt(2)*(2zx+2izy)/r^2          Y   (r)=a/sqrt(2)*(x^2-y^2+2ixy)/r^2
   2 1                                        2 2
          
  Cubic harmonics are compatible with Wien2K code:

  Y   (r)= a*(3z^2/r^2-1)/sqrt(6) 
   2 1
   
  Y   (r)= a*(x^2-y^2)/r^2     
   2 2

  Y   (r)= 2*a*yz/r^2     
   2 3                    

  Y   (r)= 2*a*zx/r^2   
   2 4                  

  Y   (r)= a*2xy/r^2    
   2 5                   

  where a=sqrt(3*5/pi)/4 is a normalization constant.  
                                                       
  Transformation matrix T(mu,i) is given by            
       _   --                                       
  Y   (r)= >   T(m,i)*Y  (r)                          
   2m      --          2i                             
          i=xyz                                        
  """
  s2 = 1./sqrt(2.)
  T2 = [[],[],[],[],[]]
  ##   m       z^2 x^2-y^2     yz    zx       xy   
  T2[-2+2]= [  0.0,    s2,    0.0,  0.0,  -1j*s2]   
  T2[-1+2]= [  0.0,   0.0, -1j*s2,   s2,     0.0]   
  T2[ 0+2]= [  1.0,   0.0,    0.0,  0.0,     0.0]   
  T2[ 1+2]= [  0.0,   0.0, -1j*s2,  -s2,     0.0]   
  T2[ 2+2]= [  0.0,    s2,    0.0,  0.0,   1j*s2]
  return array(T2)

def TS2C(l):
    if l==0: return array([[1]])
    if l==1: return TS2C_p()
    if l==2: return TS2C_d()
    #if l==3: return Ts2C_f()
    

def CoulUs(T2C, l):
    # Precomputes gaunt coefficients for speed
    # shape(gck)=(4,7,7,4). It contains gck(l, m4, m1, k)
    gck = gaunt.cmp_all_gaunt()
    nw = 2*l+1    
    T2Cp = conj(T2C.transpose())
    
    UC = zeros((l+1, nw, nw, nw, nw), dtype=complex)
    shft = 3-l
    for i4 in range(nw):
        for i3 in range(nw):
            for i2 in range(nw):
                for i1 in range(nw):
                    for k in range(l+1):
                        dsum = 0
                        for m4 in range(nw):
                            for m3 in range(nw):
                                for m2 in range(nw):
                                    for m1 in range(nw):
                                        if (m1+m2!=m3+m4): continue
                                        dsum += T2Cp[i4,m4] * gck[l,shft+m4,shft+m1,k] * T2C[m1,i1] * T2Cp[i3,m3] * gck[l,shft+m2,shft+m3,k] * T2C[m2,i2]
                        UC[k,i4,i3,i2,i1] = dsum

    # 3-l,3+l
    # 3-l,
    #print
    #for i4 in range(nw):
    #    for i1 in range(nw):
    #        for i3 in range(nw):
    #            for i2 in range(nw):
    #                f = UC[2,i4,i3,i2,i1]
    #                print "%6.3f " % f.real,
    #            print
    #        print
    #    print
        
    return UC

def CoulUsC1(l, T2C):
    # Precomputes gaunt coefficients for speed
    # shape(gck)=(4,7,7,4). It contains gck(l, m4, m1, k)
    gck = gaunt.cmp_all_gaunt()

    #print 'Gaunt coefficients precomputed - shape(gck)', shape(gck)
    nw = 2*l+1
    T2Cp = conj(T2C.transpose())
    
    UC = zeros((l+1, nw, nw, nw, nw), dtype=complex)
    shft = 3-l
    
    source="""
    using namespace std;
    for (int i4=0; i4<nw; i4++){
        for (int i3=0; i3<nw; i3++){
            for (int i2=0; i2<nw; i2++){
                for (int i1=0; i1<nw; i1++){
                    for (int k=0; k<l+1; k++){
                        complex<double> dsum = 0;
                        for (int m4=0; m4<nw; m4++){
                            for (int m3=0; m3<nw; m3++){
                                for (int m2=0; m2<nw; m2++){
                                    for (int m1=0; m1<nw; m1++){
                                        if (m1+m2!=m3+m4) continue;
                                        dsum += T2Cp(i4,m4)*gck(l,shft+m4,shft+m1,k)*T2C(m1,i1) * T2Cp(i3,m3)*gck(l,shft+m2,shft+m3,k)*T2C(m2,i2);
                                    }
                                }
                            }
                        }
                        UC(k,i4,i3,i2,i1) = dsum;
                    }
                }
            }
        }
    }
    """

    weave.inline(source, ['UC', 'gck', 'l', 'T2C', 'T2Cp', 'shft', 'nw'],
                 type_converters=weave.converters.blitz, compiler = 'gcc')

    return UC


def CoulUsC2(l, T2C):
    # Precomputes gaunt coefficients for speed
    # shape(gck)=(4,7,7,4). It contains gck(l, m4, m1, k)
    gck = gaunt.cmp_all_gaunt()

    #print 'Gaunt coefficients precomputed - shape(gck)', shape(gck)
    mw = 2*l+1
    if len(T2C) == mw:
        nw = mw
        ns = 1
    elif len(T2C) == 2*mw:
        nw = 2*(2*l+1)
        ns = 2
    else:
        print "ERROR in atom_d.py: T2C has wrong shape"
        sys.exit(0)
    
    T2Cp = conj(T2C.transpose())
    
    UC = zeros((l+1, nw, nw, nw, nw), dtype=complex)
    shft = 3-l
    shft2 = 2*l+1
    Sum1 = zeros((nw,nw,shft2*2),dtype=complex)
    Sum2 = zeros((nw,nw,shft2*2),dtype=complex)
    
    source="""
    using namespace std;
    for (int k=0; k<l+1; k++){
        Sum1=0;
        for (int i4=0; i4<nw; i4++){
            for (int i1=0; i1<nw; i1++){
                for (int m4=0; m4<mw; m4++){
                    for (int m1=0; m1<mw; m1++){
                        for (int s=0; s<ns; s++) Sum1(i4,i1,m1-m4+shft2) += T2Cp(i4,m4+s*mw)*gck(l,shft+m4,shft+m1,k)*T2C(m1+s*mw,i1);
                    }
                }
            }
        }
        Sum2=0;
        for (int i3=0; i3<nw; i3++){
            for (int i2=0; i2<nw; i2++){
                for (int m3=0; m3<mw; m3++){
                    for (int m2=0; m2<mw; m2++){
                        for (int s=0; s<ns; s++) Sum2(i3,i2,m3-m2+shft2) += T2Cp(i3,m3+s*mw)*gck(l,shft+m2,shft+m3,k)*T2C(m2+s*mw,i2);
                    }
                }
            }
        }
        for (int i4=0; i4<nw; i4++){
            for (int i3=0; i3<nw; i3++){
                for (int i2=0; i2<nw; i2++){
                    for (int i1=0; i1<nw; i1++){
                        complex<double> csum=0.0;
                        for (int dm=0; dm<shft2*2; dm++) csum += Sum1(i4,i1,dm)*Sum2(i3,i2,dm);
                        UC(k,i4,i3,i2,i1) = csum;
                    }
                }
            }
        }
    }
    """

    weave.inline(source, ['UC', 'gck', 'l', 'T2C', 'T2Cp', 'shft', 'nw', 'mw', 'ns', 'shft2', 'Sum1', 'Sum2'],
                 type_converters=weave.converters.blitz, compiler = 'gcc')

    return UC

def CoulUsC(l, T2C, op):
    # Precomputes gaunt coefficients for speed
    # shape(gck)=(4,7,7,4). It contains gck(l, m4, m1, k)
    gck = gaunt.cmp_all_gaunt()

    #print 'shape(gck)=', shape(gck)
    #print 'Gaunt precomputed'
    
    nw = 2*l+1
    # creating large T2C base
    if ( len(T2C) < 2*nw ):
        if (len(T2C)!=nw): print 'ERROR: Something wrong with size of T2C'
        T2Cl = zeros((2*nw,2*nw),dtype=complex)
        T2Cl[:nw,:nw] = T2C
        T2Cl[nw:,nw:] = T2C
    else:
        T2Cl = T2C
    
    T2Cp = conj(T2Cl.transpose())
    
    UC = zeros((l+1, 2*nw, 2*nw, 2*nw, 2*nw), dtype=complex)
    shft = 3-l

    bi = array(op.bi)
    sz = array(op.sz)
    #print 'bi=', bi
    #print 'sz=', sz
    #print 'nw=', nw
    
    source="""
    using namespace std;
    for (int i4=0; i4<2*nw; i4++){
        for (int i3=0; i3<2*nw; i3++){
            for (int i2=0; i2<2*nw; i2++){
                for (int i1=0; i1<2*nw; i1++){
                    for (int k=0; k<l+1; k++){
                        complex<double> dsum = 0;
                        for (int ms4=0; ms4<2*nw; ms4++){
                            int m4 = bi(ms4);
                            int s4 = sz(ms4);
                            //cout<<"ms4="<<ms4<<" "<<m4<<" "<<s4<<endl;
                            for (int ms3=0; ms3<2*nw; ms3++){
                                int m3 = bi(ms3);
                                int s3 = sz(ms3);
                                //cout<<"ms3="<<ms3<<" "<<m3<<" "<<s3<<endl;
                                for (int ms2=0; ms2<2*nw; ms2++){
                                    int m2 = bi(ms2);
                                    int s2 = sz(ms2);
                                    //cout<<"ms2="<<ms2<<" "<<m2<<" "<<s2<<endl;
                                    for (int ms1=0; ms1<2*nw; ms1++){
                                        int m1 = bi(ms1);
                                        int s1 = sz(ms1);
                                        //cout<<"ms1="<<ms1<<" "<<m1<<" "<<s1<<endl;
                                        if (m1+m2!=m3+m4) continue;
                                        if (s1!=s4 || s2!=s3) continue;
                                        dsum += T2Cp(i4,ms4)*gck(l,shft+m4,shft+m1,k)*T2C(ms1,i1) * T2Cp(i3,ms3)*gck(l,shft+m2,shft+m3,k)*T2C(ms2,i2);
                                    }
                                }
                            }
                        }
                        UC(k,i4,i3,i2,i1) = dsum;
                    }
                }
            }
        }
    }
    """

    
    
    weave.inline(source, ['UC', 'gck', 'l', 'T2C', 'T2Cp', 'shft', 'nw', 'bi', 'sz'],
                 type_converters=weave.converters.blitz, compiler = 'gcc')
    
    return UC



class operateLS(object):
    def __init__ (self, baths, T2C, Q3d):
        self.baths = baths
        self.Nband = self.baths/2
        self.N = self.baths
        
        self.mask=[]
        for i in range(self.N): self.mask.append(1<<i);

        self.T2C = T2C
        self.Q3d = Q3d

        if not self.Q3d:
        #############################################
        # Here for 5d's where spin-orbit is kept    #
        #############################################
            self.Q3d=False
            M2=[]
            l=(self.Nband-1)/2
            print 'L is here ', l
            for s in [0.5,-0.5]:
                for m in range(-l,l+1):
                    M2.append( (m+2*s)/2.)
            #print 'M2=',M2
            self.M2a=zeros((len(M2),len(M2)),dtype=float)
            for a in range(len(M2)):
                for b in range(len(M2)):
                    for ms in range(len(M2)):
                        self.M2a[a,b] += real(conj(T2C[ms,a])*T2C[ms,b]*M2[ms])
            #print 'M2a=', self.M2a        
        else:
        ####################################################
        # Here only for 3d's where spin-orbit is neglected #
        ####################################################
            self.Q3d=True
            self.bi=[] # band index
            self.sz=[] # sz
            for i in range(self.Nband):
                self.sz.append(1);
                self.bi.append(i);
            for i in range(self.Nband):
                self.bi.append(i)
                self.sz.append(-1)
                
            self.mask_u = []
            self.mask_d = []
            for i in range(self.Nband):
                self.mask_u.append(self.mask[i])
            for i in range(self.Nband):
                self.mask_d.append(self.mask[self.Nband+i])
        
    def Nel(self, state):
        n=0
        for k in self.mask:
            if (k&state): n+=1
        return n
    def occup(self, state):
        """ gives a list of occupancies per band [n_{band1},n_{band2},...]
        """
        oc=[]
        for i in range(self.N):
            if state & self.mask[i]: oc.append(1)
            else: oc.append(0)
        return oc
        
    def sign(self, state, mask_min, mask_max):
        """ Sign when electron hops from mask_min to mask_max
        Counts number of electrons between the two spaces
        """
        # mask will run between mask_min to mask_max
        mask = mask_min<<1 
        n=0           # number of electrons between mask_min and mask_max
        while (mask<mask_max):    # loop to mask_max
            if (mask&state): n+=1 # found electron between the two places
            mask = mask<<1        # increment the mask
        return 1-2*(n%2)          # (-1)^n
    
    def sign_(self, state, mask_max):
        """ Sign when electron is added to the state (from the left)
        """
        # mask will run between mask_min to mask_max
        mask = 1
        n=0           # number of electrons between mask_min and mask_max
        while (mask<mask_max):    # loop to mask_max
            if (mask&state): n+=1 # found electron between the two places
            mask = mask<<1        # increment the mask
        return 1-2*(n%2)          # (-1)^n
    
    def N_el_before(self, state, i):
        n=0
        for q in range(i):
            if self.mask[q]&state: n+=1
        return n

    def DM(self, state):
        "Density matrix"
        DenM = [[[] for i in range(self.N)] for j in range(self.N)]
        for j in range(self.N):
            if not self.mask[j]&state: continue # c_j operator
            jsig = self.sign_(state, self.mask[j])
            nst = state^self.mask[j]
            for i in range(self.N):
                if self.mask[i]&nst: continue   # c_i^\dagger operator
                nstate = nst^self.mask[i]
                isig = self.sign_(nst, self.mask[i])
                DenM[i][j].append( (nstate, jsig*isig) )
        return DenM
    
    def Fp(self, state, ib):
        """ This implements psi^dagger_{ib} operator acting on state
        indexes are:
          ib - band+spin index
        """
        if state&self.mask[ib]: return (0,1)  # This state is already occupied
        newstate = state^self.mask[ib]
        sig = self.sign_(state, self.mask[ib])
        return (newstate, sig)

        
    def CoulombU(self, state, UC, FkoJ, Ising=False):
        sts=[]
        ni=-1
        maxk=l+1
        if (self.Q3d):
            ### will evaluate  again <sts| U(b,a,j,i) psi^+_b psi^+_a psi_j psi_i | state>, but
            ### this time U is defined in (m4,m3,m2,m1) basis only, and is missing the spin component.
            ### Need to make sure that s_1==s_4 and s_2==s_3
            ### hence  <sts| U(m4,m3,m2,m1) psi^+_{m4,s} psi^+_{m3,s'} psi_{m2,s'} psi_{m1,s} | state>
            for i in range(self.N):
                if not(self.mask[i]&state) : continue # (i,m1) does not exists
                ni+=1
                state1 = state^self.mask[i]
                m1 = self.bi[i]
                s1 = self.sz[i]
                nj=-1
                for j in range(self.N):
                    if not(self.mask[j]&state1) : continue  # (j,m2) does not exists
                    nj+=1
                    # here we have: mask[i]&state && mask[j]&state
                    state2 = state1^self.mask[j]
                    m2 = self.bi[j]
                    s2 = self.sz[j]
                    for a in range(self.N): # (a,m3) exists
                        if self.mask[a]&state2: continue
                        if self.sz[a]!=s2 : continue # s3 == s2
                        na = self.N_el_before(state2,a)
                        state3 = state2^self.mask[a]
                        m3 = self.bi[a]
                        s3 = self.sz[a]
                        for b in range(self.N): # (b,m4) exists
                            if self.mask[b]&state3 : continue
                            if self.sz[b]!=s1: continue # s4 == s1
                            nb = self.N_el_before(state3,b)
                            state4 = state3^self.mask[b]
                            m4 = self.bi[b]
                            s4 = self.sz[b]
            
                            if Ising and state4!=state: continue
                            
                            sign = 1-2*((ni+nj+na+nb)%2)
                            U0 = sign*UC[0,m4,m3,m2,m1]*FkoJ[0]
                            
                            dsum=0
                            for k in range(1,maxk):
                                dsum += UC[k,m4,m3,m2,m1]*FkoJ[k]                            
                            U1 = sign*dsum
            
                            if (abs(U0)>1e-6 or abs(U1)>1e-6): sts.append([state4, [U0,U1]])
        
        else :  # This is used for 5d, but not for 3d
            ### will evaluate  <sts| U(b,a,j,i) psi^+_b psi^+_a psi_j psi_i | state>
            for i in range(self.N):
                if not(self.mask[i]&state) : continue # i should exist, otherwise continue
                ni+=1
                state1 = state^self.mask[i]
                nj=-1
                for j in range(self.N):
                    if not(self.mask[j]&state1) : continue  # j should exist, otherwise continue
                    nj+=1
                    # here we have: mask[i]&state && mask[j]&state
                    state2 = state1^self.mask[j]
                    for a in range(self.N): 
                        if self.mask[a]&state2: continue # a should not exist exist
                        na = self.N_el_before(state2,a)
                        state3 = state2^self.mask[a]
                        for b in range(self.N): 
                            if self.mask[b]&state3: continue # b should not exist
                            nb = self.N_el_before(state3,b)
                            state4 = state3^self.mask[b]
                            
                            if Ising and state4!=state: continue
                            
                            sign = 1-2*((ni+nj+na+nb)%2)
                            U0 = sign*UC[0,b,a,j,i]*FkoJ[0]
                            
                            dsum=0
                            for k in range(1,maxk):
                                dsum += UC[k,b,a,j,i]*FkoJ[k]
                            U1 = sign*dsum
            
                            if (abs(U0)>1e-6 or abs(U1)>1e-6): sts.append([state4, [U0,U1]])
        return sts

    def printn(self, state):
        sstate=''
        if self.Q3d:
            for i in range(self.Nband):
                if (state & self.mask_u[i]) and (state & self.mask_d[i]) : sstate += '2'
                elif (state & self.mask_u[i]) : sstate += 'u'
                elif (state & self.mask_d[i]) : sstate += 'd'
                else : sstate += '0'
        else:
            for i in range(self.N):
                if state & self.mask[i]: sstate += '1'
                else : sstate += '0'
        return sstate
    
    def Mz(self,state):
        m2=0.0
        for i in range(self.N):
            if state&self.mask[i]: m2 += self.M2a[i,i]
        return m2

    def OrbDiff(self, state, orb1=2, orb2=3):
        n1 = 0
        if (state&self.mask_u[orb1]): n1+=1
        if (state&self.mask_d[orb1]): n1+=1
        n2 = 0
        if (state&self.mask_u[orb2]): n2+=1
        if (state&self.mask_d[orb2]): n2+=1
        return n1-n2
    #########################
    
    def Sz(self, state):
        if not self.Q3d: return 0
        nu = 0
        nd = 0
        for i in range(self.Nband):
            if state&self.mask_u[i] : nu += 1
            if state&self.mask_d[i] : nd += 1
        return nu-nd
    
    def S2(self, state):
        l2p1 = self.Nband
        sts=[]
        # diagonal part
        dd=0;
        for ilz in range(l2p1):
            up=0; dn=0
            if self.mask_u[ilz] & state: up = 1
            if self.mask_d[ilz] & state: dn = 1
            # if only up or only down in certain lz
            if up+dn==1: dd += 0.5
        # Sz^2
        fct = (0.5*self.Sz(state))**2 + dd
        # store diagonal
        sts.append([state,fct])
        # off diagonal
        for ilz in range(l2p1):
            im1 = self.mask_u[ilz]
            im2 = self.mask_d[ilz]
            ib1 = bool(state & im1)
            ib2 = bool(state & im2)
            if ib1 and not ib2: # S^-_i gives nonzero
                isig = self.sign(state, min(im1,im2), max(im1,im2))
                istate = state^im1^im2
                for jlz in range(l2p1):
                    if (ilz==jlz): continue
                    jm1 = self.mask_d[jlz]
                    jm2 = self.mask_u[jlz]
                    jb1 = bool(state & jm1)
                    jb2 = bool(state & jm2)
                    if jb1 and not jb2: # S^+_j gives nonzero
                        jsig = self.sign(istate, min(jm1,jm2), max(jm1,jm2))
                        jstate = istate^jm1^jm2
                        sts.append([jstate, isig*jsig])
        return sts

    
    def PairHop(self, state):
        """ Computes the pair-hopping term:  D_a^\dagger D_a , where D_a creates or
        anhilates a double occupied site. There is no minus sign in this term!
        """
        doubles=[]
        empty=[]
        for i in range(self.Nband):
            if state & self.mask_u[i] and state & self.mask_d[i]:
                doubles.append(i)
            elif not(state & self.mask_u[i]) and not(state & self.mask_d[i]):
                empty.append(i)

        rst=[]
        for id in doubles:
            nst1 = state^self.mask_u[id]^self.mask_d[id]
            for ie in empty:
                nst2 = nst1^self.mask_u[ie]^self.mask_d[ie]
                rst.append(nst2)
        return rst
                
    def NDouble(self, state):
        ne = 0
        for i in range(self.Nband):
            if state & self.mask_u[i] and state & self.mask_d[i]: ne += 1
        return ne
            
    def OneBodyNab(self, state, Sab):
        """ computing the term Sab[a,i] f^+_a f_i
        returns all matrix elements generated by the above one-body term
        when acting on state
        """
        sts=[]
        ni=-1
        for i in range(self.baths):
            if not(self.mask[i]&state) : continue
            ni+=1
            state1 = state^self.mask[i]
            m1 = self.bi[i]
            s1 = self.sz[i]
            # here we have: mask[i]&state
            for a in range(self.baths):
                if self.mask[a]&state1 : continue # sz_a == sz_j
                # here we have: state&mask[i] and not(state1&mask[a])
                na = self.N_el_before(state1,a)
                state2 = state1^self.mask[a]
                m2 = self.bi[a]
                s2 = self.sz[a]
                        
                sign = 1-2*((ni+na)%2)
                
                nab = sign*Sab[a,i]
                        
                if (abs(nab)>1e-6): sts.append([state2, nab])
        return sts

#class operateLS(object):
#    def __init__ (self, Nband):
#        self.Nband = Nband
#        self.baths = 2*Nband
#
#        self.N = self.baths
#        
#        self.mask=[]
#        for i in range(self.N): self.mask.append(1<<i);
#        
#        self.bi=[] # band index
#        self.sz=[] # sz
#        for i in range(self.Nband):
#            self.sz.append(1);
#            self.bi.append(i);
#        for i in range(self.Nband):
#            self.bi.append(i)
#            self.sz.append(-1)
#            
#        self.mask_u = []
#        self.mask_d = []
#        for i in range(self.Nband):
#            self.mask_u.append(self.mask[i])
#        for i in range(self.Nband):
#            self.mask_d.append(self.mask[self.Nband+i])
#        
#        
#    def printn(self, state):
#        sstate=''
#        for i in range(self.Nband):
#            if (state & self.mask_u[i]) and (state & self.mask_d[i]) : sstate += '2'
#            elif (state & self.mask_u[i]) : sstate += 'u'
#            elif (state & self.mask_d[i]) : sstate += 'd'
#            else : sstate += '0'
#            #sstate += ' '
#        return sstate
#    
#    def Nel(self, state):
#        n=0
#        for k in self.mask:
#            if (k&state): n+=1
#        return n
#    
#    def Sz(self, state):
#        nu = 0
#        nd = 0
#        for i in range(self.Nband):
#            if state&self.mask_u[i] : nu += 1
#            if state&self.mask_d[i] : nd += 1
#        return nu-nd
#    
#    def occup(self, state):
#        """ gives a list of occupancies per band [n_{band1},n_{band2},...]
#        """
#        oc=[]
#        for i in range(self.Nband):
#            ne = 0
#            if state & self.mask_u[i] : ne += 1
#            oc.append(ne)
#        for i in range(self.Nband):
#            ne = 0
#            if state & self.mask_d[i] : ne += 1
#            oc.append(ne)
#        return oc
#        
#    def sign(self, state, mask_min, mask_max):
#        """ Sign when electron hops from mask_min to mask_max
#        Counts number of electrons between the two spaces
#        """
#        # mask will run between mask_min to mask_max
#        mask = mask_min<<1 
#        n=0           # number of electrons between mask_min and mask_max
#        while (mask<mask_max):    # loop to mask_max
#            if (mask&state): n+=1 # found electron between the two places
#            mask = mask<<1        # increment the mask
#        return 1-2*(n%2)          # (-1)^n
#    
#    def sign_(self, state, mask_max):
#        """ Sign when electron is added to the state (from the left)
#        """
#        # mask will run between mask_min to mask_max
#        mask = 1
#        n=0           # number of electrons between mask_min and mask_max
#        while (mask<mask_max):    # loop to mask_max
#            if (mask&state): n+=1 # found electron between the two places
#            mask = mask<<1        # increment the mask
#        return 1-2*(n%2)          # (-1)^n
#    
#    def N_el_before(self, state, i):
#        n=0
#        for q in range(i):
#            if self.mask[q]&state: n+=1
#        return n
#
#    def DM(self, state):
#        "Density matrix"
#        DenM = [[[] for i in range(self.baths)] for j in range(self.baths)]
#        for j in range(self.baths):
#            if not self.mask[j]&state: continue # c_j operator
#            jsig = self.sign_(state, self.mask[j])
#            nst = state^self.mask[j]
#            for i in range(self.baths):
#                if self.mask[i]&nst: continue   # c_i^\dagger operator
#                nstate = nst^self.mask[i]
#                isig = self.sign_(nst, self.mask[i])
#                DenM[i][j].append( (nstate, jsig*isig) )
#        return DenM
#    
#        
#    def S2(self, state):
#        l2p1 = self.Nband
#        sts=[]
#        # diagonal part
#        dd=0;
#        for ilz in range(l2p1):
#            up=0; dn=0
#            if self.mask_u[ilz] & state: up = 1
#            if self.mask_d[ilz] & state: dn = 1
#            # if only up or only down in certain lz
#            if up+dn==1: dd += 0.5
#        # Sz^2
#        fct = (0.5*self.Sz(state))**2 + dd
#        # store diagonal
#        sts.append([state,fct])
#        # off diagonal
#        for ilz in range(l2p1):
#            im1 = self.mask_u[ilz]
#            im2 = self.mask_d[ilz]
#            ib1 = bool(state & im1)
#            ib2 = bool(state & im2)
#            if ib1 and not ib2: # S^-_i gives nonzero
#                isig = self.sign(state, min(im1,im2), max(im1,im2))
#                istate = state^im1^im2
#                for jlz in range(l2p1):
#                    if (ilz==jlz): continue
#                    jm1 = self.mask_d[jlz]
#                    jm2 = self.mask_u[jlz]
#                    jb1 = bool(state & jm1)
#                    jb2 = bool(state & jm2)
#                    if jb1 and not jb2: # S^+_j gives nonzero
#                        jsig = self.sign(istate, min(jm1,jm2), max(jm1,jm2))
#                        jstate = istate^jm1^jm2
#                        sts.append([jstate, isig*jsig])
#        return sts
#
#    
#    def Fp(self, state, ib):
#        """ This implements psi^dagger_{ib} operator acting on state
#        indexes are:
#          ib - band+spin index
#        """
#        if state&self.mask[ib]: return (0,1)  # This state is already occupied
#        newstate = state^self.mask[ib]
#        sig = self.sign_(state, self.mask[ib])
#        return (newstate, sig)
#
#    def PairHop(self, state):
#        """ Computes the pair-hopping term:  D_a^\dagger D_a , where D_a creates or
#        anhilates a double occupied site. There is no minus sign in this term!
#        """
#        doubles=[]
#        empty=[]
#        for i in range(self.Nband):
#            if state & self.mask_u[i] and state & self.mask_d[i]:
#                doubles.append(i)
#            elif not(state & self.mask_u[i]) and not(state & self.mask_d[i]):
#                empty.append(i)
#
#        rst=[]
#        for id in doubles:
#            nst1 = state^self.mask_u[id]^self.mask_d[id]
#            for ie in empty:
#                nst2 = nst1^self.mask_u[ie]^self.mask_d[ie]
#                rst.append(nst2)
#        return rst
#                
#    def NDouble(self, state):
#        ne = 0
#        for i in range(self.Nband):
#            if state & self.mask_u[i] and state & self.mask_d[i]: ne += 1
#        return ne
#            
#    def CoulombU(self, state, UC, FkoJ, Ising=False):
#        sts=[]
#        ni=-1
#        maxk=l+1
#        #maxk=2
#        for i in range(self.baths):
#            if not(self.mask[i]&state) : continue # (i,m1) does not exists
#            ni+=1
#            state1 = state^self.mask[i]
#            m1 = self.bi[i]
#            s1 = self.sz[i]
#            nj=-1
#            for j in range(self.baths):
#                if not(self.mask[j]&state1) : continue  # (j,m2) does not exists
#                nj+=1
#                # here we have: mask[i]&state && mask[j]&state
#                state2 = state1^self.mask[j]
#                m2 = self.bi[j]
#                s2 = self.sz[j]
#                for a in range(self.baths): # (a,m3) exists
#                    if self.mask[a]&state2 or self.sz[a]!=s2 : continue # s3 == s2
#                    na = self.N_el_before(state2,a)
#                    state3 = state2^self.mask[a]
#                    m3 = self.bi[a]
#                    s3 = self.sz[a]
#                    for b in range(self.baths): # (b,m4) exists
#                        if self.mask[b]&state3 or self.sz[b]!=s1: continue # s4 == s1
#                        nb = self.N_el_before(state3,b)
#                        state4 = state3^self.mask[b]
#                        m4 = self.bi[b]
#                        s4 = self.sz[b]
#
#
#                        if Ising and state4!=state: continue
#                        
#                        sign = 1-2*((ni+nj+na+nb)%2)
#                        U0 = sign*UC[0,m4,m3,m2,m1]*FkoJ[0]
#                        
#                        dsum=0
#                        for k in range(1,maxk):
#                            dsum += UC[k,m4,m3,m2,m1]*FkoJ[k]                            
#                        U1 = sign*dsum
#
#                        if (abs(U0)>1e-6 or abs(U1)>1e-6): sts.append([state4, [U0,U1]])
#        return sts
#
#    def OneBodyNab(self, state, Sab):
#        """ computing the term Sab[a,i] f^+_a f_i
#        returns all matrix elements generated by the above one-body term
#        when acting on state
#        """
#        sts=[]
#        ni=-1
#        for i in range(self.baths):
#            if not(self.mask[i]&state) : continue
#            ni+=1
#            state1 = state^self.mask[i]
#            m1 = self.bi[i]
#            s1 = self.sz[i]
#            # here we have: mask[i]&state
#            for a in range(self.baths):
#                if self.mask[a]&state1 : continue # sz_a == sz_j
#                # here we have: state&mask[i] and not(state1&mask[a])
#                na = self.N_el_before(state1,a)
#                state2 = state1^self.mask[a]
#                m2 = self.bi[a]
#                s2 = self.sz[a]
#                        
#                sign = 1-2*((ni+na)%2)
#                
#                nab = sign*Sab[a,i]
#                        
#                if (abs(nab)>1e-6): sts.append([state2, nab])
#        return sts
        

def baseN(Nband, prop,Q3d):
    Ntot = len(prop)
    wstates=[]

    if Q3d:
        for n1 in range(Nband*2+1):
            for sz1 in range(-n1,n1+1,2):
                states=[]
                for i in range(Ntot):
                    if prop[i][0]==n1 and prop[i][1]==sz1:
                        states.append(i)
                        
                if (len(states)>0): wstates.append([n1, sz1, states])
    else:
        for n1 in range(Nband*2+1):
            states=[]
            for i in range(Ntot):
                if prop[i][0]==n1:
                    states.append(i)
            if (len(states)>0): wstates.append([n1, 0, states])
        
    return wstates

def list_to_string(x):
    return str(array(x).flatten().tolist())



def analizeGroups(A, small = 1e-4):
    groups=[]
    for i in range(shape(A)[0]):    
        nonz=[]
        for j in range(shape(A)[1]):
            if abs(A[i,j])>small : nonz.append(j)
        if (len(nonz)>0): groups.append(nonz)
    
    groups0 = compress(groups)
    
    groups=[]
    for i in range(shape(A)[1]):
        nonz=[]
        for j in range(shape(A)[0]):
            if abs(A[j,i])>small : nonz.append(j)
        if (len(nonz)>0): groups.append(nonz)
    
    groups1 = compress(groups)
    return (groups1,groups0)

def coupled(A, groups0, groups1, small = 1e-4):
    #ng0 = len(array(groups0).flatten().tolist())
    #ng1 = len(array(groups1).flatten().tolist())

    fpair = [-1]*len(groups0)
    #pairs=[]
    for ii,ig0 in enumerate(groups0):
        nonz=[]
        for ir0 in ig0:
            for q in range(shape(A)[1]):
                if abs(A[ir0,q])>small : nonz.append(q)
        for jj,jg1 in enumerate(groups1):
            if overlap(nonz,jg1):
                #pairs.append([ii,jj])
                fpair[ii] = jj
    return fpair

def comp(x, y):
    if x[2]!=y[2]: return int(x[2]-y[2])
    else:
        if abs(x[3]-y[3])<1e-5: return 0
        elif (x[3]<y[3]): return -1
        else: return 1



def SpinOrbitM(l,T2C):
    # one electron |l,m,s> base
    ms_base=[]
    for s in [1/2.,-1/2.]:
        for m in range(-l,l+1):
            ms_base.append([m,s])
    #print 'ms_base=', ms_base
    
    # one electron |j,mj> base
    pj = [l-1/2.,l+1/2.]
    if l==0 : pj = [0.5]
    jj_base=[]
    for j in pj:
        for mj in arange(-j,j+1):
            jj_base.append([j, mj])
    #print 'jj_base=', jj_base
    
    # transforms between |lms> and |jmj> base
    Tjls = zeros((len(ms_base),len(jj_base)))
    for ijj,jj in enumerate(jj_base):
        for ims,ms in enumerate(ms_base):
            Tjls[ijj,ims] = gaunt.clebschg(jj[0], jj[1], l, ms[0], 1/2., ms[1])

    # the one-body operator l*s in matrix form
    # in the j-j base
    jSO = zeros((len(jj_base),len(jj_base)))
    for ijj,jj in enumerate(jj_base):
        jSO[ijj,ijj] = 0.5*(jj[0]*(jj[0]+1) - l*(l+1) - 3/4.)


    #print 'jSO=', jSO
    #mprint(jSO)
    
    # changing to lms base
    mSO = matrix(Tjls.transpose())*matrix(jSO)*matrix(Tjls)

    # creating large T2C base
    if ( len(T2C) < len(jj_base) ):
        T2Cl = zeros(tuple(array(shape(T2C))*2),dtype=complex)
        T2Cl[:len(T2C),:len(T2C)] = T2C
        T2Cl[len(T2C):,len(T2C):] = T2C
    else:
        T2Cl = T2C
    
    # changing to cubic harmonics base
    cSO = matrix(conj(T2Cl.transpose())) * mSO * matrix(T2Cl)

    print 'spin-orbit='
    mprint(sys.stdout,real(cSO))
    
    return cSO


def Diagonalize(Ham, small=1e-4, fh=sys.stdout):
    """ Diagonalization is done in blocks. This is not because of efficiency but because
    the resulting eigenvectors must not mix states of direct base if not absolutely necessary.
    If brute force diagonalization is used in large scale problems, eigenvectors can be seriously
    mix direct states with different symmetry.
    """
    diff = sum(Ham-transpose(conj(Ham)))
    if abs(diff)>1e-6:
        print 'H NOT HERMITIAN!'
    
    # Check block structure of Hamiltonian
    # States which are mixed in Ham will have the same blck[i]
    ndim = len(Ham)
    blck=range(ndim)
    for i in range(ndim):
        for j in range(i+1,ndim):
            if (abs(Ham[i][j])>small):
                commonb = min(blck[i],blck[j])
                for k in range(ndim):
                    if blck[k] in [blck[i],blck[j]]: blck[k]=commonb
        #print ('%2d'%i), 'current blck=', '%2d,'*len(blck) % tuple(blck)
        
    # Having blck[i] a new array block[:][:] is created, which contains indexes to all blocks
    # for example [[1,2,3],[4,5,6]] for Hamiltonian containing two blocks
    block=[]
    for i in range(ndim):
        bb=[]
        for j in range(ndim):
            if blck[j]==i: bb.append(j)
        if len(bb)>0:
            block.append(bb)
    #print 'block=', block
    
    # Here we go over all blocks and diagonalize each one.
    eigv=[] # contains all eigenvalues
    eigx=[] # contains full eigenvectors
    for ibl,bl in enumerate(block):
        hs = zeros((len(bl),len(bl)), dtype=complex)
        for i,ib in enumerate(bl):
            for j,jb in enumerate(bl):
                hs[i,j] = Ham[ib,jb]

        if symeig_available:
            eigy = symeig.symeig(hs) # diagonalization of small block
        else:
            eigy = linalg.eigh(hs)

        print >> fh, 'Eigenvalues[',bl,']=',eigy[0]
        
        # Checking if eigenvectors are complex!
        for l in range(len(eigy[1])):
            imax=0
            #print 'shape(eigy[1])', shape(eigy[1])
            for iu in range(len(eigy[0])):
                #print iu, imax
                if abs(eigy[1][iu,l])>abs(eigy[1][imax,l]): imax=iu
            z=eigy[1][imax,l]
            phi=math.atan2(z.imag,z.real)
            eigy[1][:,l] *= exp(-phi*1j)
            
            ime = sum([abs(x.imag) for x in eigy[1][:,l]])
            if (abs(ime))<1e-10: ime=0
            print >> fh, 'im=%2d %2d %f' % (ibl, l, ime)
            #ime = sum([abs(eigy[1][u,l].imag) for u in range(len(eigy[1]))])
        #    if ime>1e-7: print 'TROUBLES!!! Complex eigenvector! You sould improve that!'
        
        # Creating a big eigenvector with all components
        for l in range(len(eigy[1])):
            large_eig=zeros(ndim, dtype=complex)
            small_eig = eigy[1][:,l]
            for m,mb in enumerate(bl):  large_eig[mb] = small_eig[m]
            eigx.append(large_eig)
        eigv += eigy[0].tolist()

    # Now we need to sort eigenvectors and eigenvalues
    # index is created for sorting
    indx=range(ndim)    
    indx.sort(lambda a,b: cmp(eigv[a],eigv[b]))
    # and actual sorting is performed
    seigv=[]
    seigx=[]
    for i in range(ndim):
        seigv.append(eigv[indx[i]])
        seigx.append(eigx[indx[i]])

    # Eigenvectors should be in the form Ham*v[:,i] = w[i]*v[:,i]
    # which means that we need to transpose the list of eigenvectors
    seigx = array(seigx).transpose()

    # We also do a brute force diagonalization just to check if something goes wrong with block diagonalization
    # Note that the two resulting eigensystems are not necessary the same due to freedom in choosing eigenvectors
    eig = linalg.eigh(Ham)

    # If eigenvalues from block diagonalization and full diagonalization are different, something is wrong
    if sum(map(abs,eig[0]-array(seigv)))>small:
        print '!!!!!TEZAVE!'
        print 'The right eigenvalues are:', eig[0]
    
    return [seigv,seigx]

def EquivalentBaths(Eimp):
    """ Finds which baths are equivalent from impurity levels"""
    wE = [(i,Eimp[i]) for i in range(len(Eimp))]
    #print 'wE=', wE
    kbths=[]
    while len(wE)>0:
        En = wE[0][1]
        j=0
        rr=[]
        while j < len(wE):
            if abs(En-wE[j][1])<1e-10:
                rr.append(wE[j][0])
                del wE[j]
            else: j+=1
            #print 'w', j, rr, En, wE[j][1]
        kbths.append(rr)
    
    bathis=range(len(Eimp))
    for ik,k in enumerate(kbths):
        for ij in k: bathis[ij]=ik

    #print 'kbths=', kbths
    Ed=[]
    for ik,k in enumerate(kbths):
        Ed.append(Eimp[k[0]])
        
    return (bathis,kbths,Ed)

def thesame(mx,my,small=1e-3):
    if mx.keys() != my.keys(): return False
    for k in mx.keys():
        if abs(mx[k]-my[k])>small: return False
    return True
    
def VEquivalentStates(mps,ind):
    """ Finds which states have the same bubbles """
    wx = [(i,mps[i]) for i in range(len(mps))]
    iequiv=[]
    while len(wx)>0:
        mx = wx[0][1]
        j=0
        rr=[]
        while j < len(wx):
            if thesame(mx,wx[j][1]):
                rr.append(wx[j][0])
                del wx[j]
            else: j+=1
        iequiv.append(rr)
    
    for ik in range(len(iequiv)):
        for ij in range(len(iequiv[ik])):
            iequiv[ik][ij] = ind[iequiv[ik][ij]]
            
    return iequiv

def AverageBubbles(tmps):
    """ Compute average over 'almost' equivalent states """
    trmp=[]
    for mps in tmps:
        all_keys=[]
        for mp in mps:
            all_keys = union(all_keys,mp.keys())
        
        rmp={}
        for k in all_keys:
            sm=0.0
            for mp in mps:
                if mp.has_key(k):
                    sm += mp[k]
            #sm/=len(mps)
            rmp[k]=sm
        trmp.append(rmp)    
    return trmp


def EquivalentStates(ipE, ipN):
    iequiv=[]
    equiv = range(len(ipE))
    leq=0
    ju=0
    Nmax = ipN[-1]
    for Ni in range(Nmax+1):
        # all states of the same N are in the interval [ju,je]
        je=ju
        while je<len(ipN) and ipN[je]==Ni: je+=1

        ind = range(ju,je)
        ind.sort(lambda x,y: cmp(ipE[x],ipE[y]) )

        #print Ni
        i0=0
        while (i0<len(ind)):
            Ec = ipE[ind[i0]]
            ieq=[]
            #print 'Ec=', Ec
            while i0<len(ind) and abs(ipE[ind[i0]]-Ec)<1e-10:
                #print ind[i0], ipE[ind[i0]], leq
                equiv[ind[i0]] = leq
                ieq.append(ind[i0])
                i0+=1
            leq += 1
            iequiv.append(ieq)
            #print
        #print
        ju=je
    return (equiv, iequiv)

def RenumberStates(pseudostates, Enes, wstates, S2ws):
    # renumbers states such that each of 1024 states has unique index
    # also remembers energy and N for each state
    ij=0
    puniq={}
    ipuniq=[]
    ipE=[]
    ipN=[]
    ipS=[]
    for ii,iwp in enumerate(pseudostates):
        wdim = len(Enes[ii])
        for j in range(wdim):
            puniq[(ii,j)]=ij
            ipuniq.append((ii,j))
            ipE.append(Enes[ii][j])
            ipS.append(S2ws[ii][j])
            wN = sum(wstates[iwp[0]][0])
            ipN.append(wN)
            ij+=1
    return (puniq, ipE, ipN, ipS)


def CreateEmpty3D_Dict(n0,n1,n2):
    return [[[{} for i2 in range(n2)] for i1 in range(n1)] for i0 in range(n0)]
def CreateEmpty2D_Dict(n0,n1):
    return [[{} for i1 in range(n1)] for i0 in range(n0)]


def ReadTrans(filename, fh_info):
    """Read the self-energy index file Sigind and the local transformation matrix CF from a file"""
    fh = open(filename, 'r')
    data = fh.readlines()
    
    (n1,n2) = map(int, data[0].split()[:2])

    Sigind=[]
    for i in range(n1):
        Sigind.append( map(int, data[i+2].split()[:n2]) )
    Sigind = array(Sigind)

    print >> fh_info, 'len(data)', len(data)
    print >> fh_info, 'n1=', n1
    if len(data) >= n1+n1+3:
        n2 = n1
        CF=[]
        for i in range(n2):
            cl = array(map(float, data[n1+3+i].split()))
            CF.append( cl[0::2]+cl[1::2]*1j )
        CF = array(CF)
    elif len(data)>=n1+n1/2+3:
        n2 = n1/2
        CF=[]
        for i in range(n2):
            cl = array(map(float, data[n1+3+i].split()))
            CF.append( cl[0::2]+cl[1::2]*1j )
        CF = array(CF)
        CFN = zeros((2*n2,2*n2), dtype=complex)
        CFN[:n2,:n2] = CF
        CFN[n2:,n2:] = CF
        CF = CFN
    else:
        CF = identify(n1)
    print >> fh_info, 'CF=', CF
    
    return (Sigind, CF)

def SlaterF(U, J, l):
    Fk = zeros((4,4), dtype=float)
    if l==0:
        # F0 for s-electrons
        Fk[0,0] = U
    elif l==1:
        # F2 for p-electrons
        Fk[0,1] = U
        if type(J) is list:
            Fk[1,1] = 5*J[0]
        else:
            Fk[1,1] = 5*J
    elif l==2:
        # F2 and F4 for d-electrons
        Fk[0,2] = U
        if type(J) is list:
            Fk[1,2] = 14./1.625 * J[0]
            Fk[2,2] = 14.*0.625/1.625 * J[1]
        else:
            Fk[1,2] = 14./1.625 * J
            Fk[2,2] = 14.*0.625/1.625 * J
    elif l==3:
        # F2, F4 and F6 for f-electrons
        Fk[0,3] = U
        if type(J) is list:
            Fk[1,3] = 6435./(286+195*0.668+250*0.494) * J[0]
            Fk[2,3] = 0.668*6435./539.76 * J[1]
            Fk[3,3] = 0.494*6435./539.76 * J[2]
        else:
            Fk[1,3] = 6435./(286+195*0.668+250*0.494) * J
            Fk[2,3] = 0.668*6435./539.76 * J
            Fk[3,3] = 0.494*6435./539.76 * J
    return Fk



def Check_T2C_Real(T2C, l, fh_info, small):
    """
     Here we added a routine which checks that cubic harmonics are real.
     Only in this case operators F^+ and F used in ctqmc will be real.
     Otherwise these matrix elements might be complex.
     The condition for cubic harmonics to be real is:
    
      Imag( T2C[m,i] + (-1)**m * T2C[-m,i] ) == 0
        and
      Real( T2C[m,i] - (-1)**m * T2C[-m,i] ) ==0
    
     which follows from the requirement: \sum_m T2C[m,i]*exp(i*m*phi) is real for any phi
    
     We are free to add any phase to cubic harmonics, hence T2C[m,i] -> T2C[m,i]*exp(i*phi_i)
       with phi_i being arbitrary
    
     This leads to the following 2x2 system of equations:
      ( Rp[m,i], Qp[m,i] ) ( sin(phi_i) )
      ( Qm[m,i], Rm[m,i] ) ( cos(phi_i) ) = 0
    
     where
          Qp[m,i] = Imag( T2C[m,i] + (-1)**m * T2C[-m,i] )
          Rm[m,i] = Real( T2C[m,i] - (-1)**m * T2C[-m,i] )
          Rp[m,i] = Real( T2C[m,i] + (-1)**m * T2C[-m,i] )
          Qm[m,i] = Imag(-T2C[m,i] + (-1)**m * T2C[-m,i] )
    """
    
    for i in range(2*l+1):
        ctg=None
        for m in range(0,l+1):
            Qp = T2C[m+l,i].imag + (-1)**m * T2C[-m+l,i].imag
            Rm = T2C[m+l,i].real - (-1)**m * T2C[-m+l,i].real
            Rp = T2C[m+l,i].real + (-1)**m * T2C[-m+l,i].real
            Qm =-T2C[m+l,i].imag + (-1)**m * T2C[-m+l,i].imag

            if abs(Qp) > small or abs(Rm) > small:
                if abs(Qp) > small :
                    ctg = -Rp/Qp
                    xb = -Rp
                    xa = Qp
                if abs(Rm) > small :
                    ctg = -Qm/Rm
                    xb = -Qm
                    xa = Rm
                
        if ctg is not None:
            for m in range(0,l+1):
                Qp = T2C[m+l,i].imag + (-1)**m * T2C[-m+l,i].imag
                Rm = T2C[m+l,i].real - (-1)**m * T2C[-m+l,i].real
                Rp = T2C[m+l,i].real + (-1)**m * T2C[-m+l,i].real
                Qm =-T2C[m+l,i].imag + (-1)**m * T2C[-m+l,i].imag

                if abs(Rp + Qp * ctg)>small or abs(Qm + Rm * ctg)>small:
                    print 'ERROR: Could not find an angle to make all cubic harmonics real'

            phi = arctan2(xa, xb)                        
            #print i, ctg, exp(phi*1j)
            print >> fh_info, 'Correcting T2C because original cubic harmonics were not real'
            print >> fh_info, 'T2C before correction:'
            cprint(fh_info, T2C)
            T2C[:,i] = T2C[:,i] * exp(phi*1j)
            print >> fh_info, 'T2C after correction:'
            cprint(fh_info, T2C)
            

if __name__ == '__main__':
    """ Help here"""
    
    n=[1,2,3]  # occupanices used for OCA
    l=1        # angular momentum
    J = 0.3    # Hunds coupling
    qOCA=1     # OCA diagrams are computed 
    Eoca=10.   # Energy window for OCA diagrams
    mOCA=1e-3  # matrix element for OCA should be greater than that
    Ncentral=[5] # OCA diagrams are selected such that central occupancy is in Ncentral
    Ewindow = [-1000,1000]
    max_M_size=500
    add_occupancy=True
    CoulombF = 'Full' #'Bulla_Jarrel' #'Full' # 'Bulla_Jarrel' # 'Full' 'Oles'
    OCA_G=True
    PrintReal=True
    HB2 = False
    Eimp = [0.0]
    Nmax = 1024

    small = 1e-6
    
    args = sys.argv[1:]
    if ('-h' in args) or ('--help' in args):
        print """Code for generating impurity cix file for a d-type of material
                 The output cix-file can be used for OCA or CTQMC solvers
                 The outputs:
                    out.cix       -- the oca cix file
                    impurity.cix  -- the ctqmc cix file
                 The input is:
                    filename      -- another python script, which contains some definitions of parameters
                    Eimp          -- list of impurity levels (i.e., [0,0,0...] )
                    Sigind        -- The symmetry of the self-energy and impurity levels given as a 2D-list
                    CF            -- The local rotation matrix given as a 2D list or array
                    n             -- a list of occupancies use in oca output [1,2,3]
                    l             -- angular momentum (l=0, l=1, l=2 supported)
                    J             -- Hund's coupling
                    qOCA          -- OCA diagrams included if qOCA=1
                    Eoca          -- OCA diagrams cutof: if any of the atomic states has for Eoca higher energy than the ground state for particular occupancy, the diagram is dropped
                    mOCA          -- If matrix element is smaller than mOCA, the diagram is dropped
                    Ewindow       -- Energy window for the states kept (used in ctqmc only)
                    max_M_size    -- maximum matrix size kept for ctqmc
                    Ncentral      -- a list of central occupancies for OCA [1]
                    OCA_G         -- bool - comput input for OCA as well
                 """
        sys.exit(0)
    Q3d=None
    for arg in args:
        if os.path.isfile(arg):
            execfile(arg)
            print 'Executed file', arg
        else:
            exec(arg)

    #if CoulombF=='Ising': HB2 = True
    if CoulombF=='Georges': CoulombF='Ising'
            
    fh_info = open('info_atom_d.dat','w')
    print >> fh_info, ' '.join(sys.argv)
    print >> fh_info, 'Eimp=', '[','%f, '*len(Eimp) % tuple(Eimp),']'
    print >> fh_info, 'n=', n
    print >> fh_info, 'l=', l
    print >> fh_info, 'J=', J
    print >> fh_info, 'qOCA=', qOCA
    print >> fh_info, 'Eoca=', Eoca
    print >> fh_info, 'mOCA=', mOCA
    print >> fh_info, 'Ewindow=', Ewindow
    print >> fh_info, 'max_M_size=', max_M_size


    ftrans='Trans.dat'
    if (len(glob.glob(ftrans))!=0):
        print >> fh_info, 'Reading file', ftrans
        (Sigind, CF) = ReadTrans(ftrans, fh_info)
        
        if len(Sigind)==(2*l+1):
            dn = 2*l+1
            SigindN = zeros((2*dn, 2*dn), dtype=int)
            SigindN[:dn,:dn] = Sigind
            SigindN[dn:,dn:] = Sigind
            Sigind=SigindN
        
        if (len(Sigind)!= 2*(2*l+1)):
            print 'ERROR: Sigind from file', ftrans, 'does not have correct dimension!'
            sys.exit(1)

        if len(CF)==2*l+1:
            if Q3d==None: Q3d=True # this must be 3d orbital
        elif len(CF)==2*(2*l+1):
            dn=2*l+1
            off_diag = CF[dn:,:dn]
            off = sum(sum(abs(off_diag)))
            print off
            if Q3d==None:
                if off>1e-5:
                    Q3d=False
                else:
                    Q3d=True
        else:
            print 'ERROR: Transformation CF=T2C does not have correct dimension'
            sys.exit(1)

        # If any diagonal entry in Sigind[] is equal to zero, we want to project it out.
        # This is achieved by making impurity levels Eimp very large for this orbital
        Simax = max(map(max,Sigind))
        for i in range(len(Sigind)):
            if Sigind[i,i]==0:
                Sigind[i,i]=Simax+1
                if len(Eimp)<=Simax : Eimp += [2000.]
                
        
        # Matrix contained in Trans.dat should be T(i,m):
        #   - i runs over real harmonics (z^2, x^2-y^2, yz, xz, xy)
        #   - m runs over complex harmonics (-2, -1, 0, 1, 2)
        # The matrix in Trans.dat, read into CF, is the transpose of
        # what we want in T2C, hence the T2C = transpose(CF) statement
        if Q3d:
            ##### change 2013
            T2C = transpose(CF[:(2*l+1),:(2*l+1)])
            Check_T2C_Real(T2C, l, fh_info, small)
        else:
            #CFN = zeros(shape(CF),dtype=complex)
            #CFN[:,:(2*l+1)] = CF[:,(2*l+1):]
            #CFN[:,(2*l+1):] = CF[:,:(2*l+1)]
            CFN = CF
            T2C = transpose(CFN)
        ##### change 2013, you will need to generalize this
        
    else:
        """
        Cubic harmonics have the order compatible with the Wien2K package:
        z^2, x^2-y^2, yz, zx, xy
        """
        print ftrans, 'file not found; generating Sigind & complex-to-real spherical harmonics transformation inside atom_d.'
        T2C = TS2C(l)
        if Q3d==None: Q3d=True

        Sigind = zeros((2*(2*l+1),2*(2*l+1)), dtype=int)
        if Q3d:
            for i in range(2*l+1):
                Sigind[i,i] = i+1
                Sigind[i+(2*l+1),i+(2*l+1)] = i+1
        else:
            zr = zeros(2*l+1)
            CF = transpose(T2C)
            CFn=[]
            for i in range(2*l+1):
                CFn.append( CF[i].tolist()+zr.tolist()  )
                CFn.append( zr.tolist()+CF[i].tolist()  )
                Sigind[2*i,2*i] = i+1
                Sigind[2*i+1,2*i+1] = i+1
            CFn = array(CFn)
            T2C = transpose(CFn)
            
    if Q3d:
        global_flip = range(2*l+1) + range(2*l+1)
    else:
        global_flip=[]
        for i in range(2*l+1):
            global_flip += [i,i]
            
    print 'Sigind=', Sigind
    
    print 'T2C='
    for i in range(len(T2C)):
        for j in range(len(T2C)):
            print "%6.3f %6.3f   " % (T2C[i,j].real, T2C[i,j].imag),
        print
    
    print 'global_flip=', global_flip
    
    print >> fh_info, 'Impurity level structure Sigind is:'
    print >> fh_info, Sigind

    print >> fh_info, 'T2C follows:'
    print >> fh_info, '\n'.join('   '.join('%10f %10f' % (x.real, x.imag) for x in row) for row in T2C)
    print >> fh_info, 'shape(T2C)=', shape(T2C)
    print >> fh_info, 'T2C is Unitary=', sum(abs(matrix(T2C) * matrix(T2C).H - identity(len(T2C))))
    
    Nitt=1  # To figure out the symmetry, we itterate Nitt times
    
    Jc = J
    cx = 0. # No spin orbit at the moment

    # Ratio between F2,F4,F6 and J! At the end of the day, we want to use U and J only!
    #Fk = gaunt.slaterf(1., Jc)
    U0=1.
    Fk = SlaterF(U0, Jc, l)
    print >> fh_info, 'Slater integrals F^k are ', Fk[:,l]


    # one electron base
    baths=[]
    for s in [1,-1]:
        for b in range(2*l+1):
            baths.append([b,s])
    bathi=[Sigind[b,b]-1 for b in range(len(Sigind))]

    dkbth={}
    for i in range(len(bathi)):
        if dkbth.has_key(bathi[i]):
            dkbth[bathi[i]].append(i)
        else:
            dkbth[bathi[i]]=[i]
    kbth=[]
    for k in sort(dkbth.keys()):
        kbth.append(dkbth[k])
        
    kbth0=[]
    for i in range(len(baths)): kbth0.append([i])


    bkeep=[]
    for b in range(len(bathi)):
        if Eimp[bathi[b]]<1000:  bkeep.append(b)
    tEimp = filter(lambda x: x<1000,Eimp)
    tkbth=[]
    for k in kbth:
        if k[0] in bkeep: tkbth.append(k)
        
    print >> fh_info, 'Some other info in ED:'
    print >> fh_info, 'bathi=', bathi
    print >> fh_info, 'kbth=', kbth
    print >> fh_info, 'tkbth=', tkbth
    print >> fh_info, 'Eimp=', Eimp
    print >> fh_info, 'tEimp=', tEimp
    print >> fh_info, 'bkeep=', bkeep
    
    Ntot = 2**(len(baths)) # size of the direct base
    
    op = operateLS(2*(2*l+1), T2C, Q3d) # operators class
    
    
    if op.Q3d:
        print >> fh_info, 'baths bi=', op.bi
        print >> fh_info, 'spin Sz=', op.sz
        print >> fh_info, 'mask-down=', op.mask_d
        print >> fh_info, 'mask-up  =', op.mask_u
    
    # some properties of integers which will serve a direct base - partial occupancy and Sz
    prop=[]
    for i in range(Ntot):
        #### 2013 ### all these are wrong for 5d
        occ = op.occup(i)
        prop.append([sum(occ), op.Sz(i),occ])
    # creates direct base from integers having correct properties
    # wstates contains [N, Sz, [all direct states with this N and Sz]]
    wstates = baseN(2*l+1,prop,op.Q3d)

    SO = SpinOrbitM(l,T2C) # Spin orbit matrix
    
    UC = CoulUsC2(l,T2C)     # Coulomb repulsion matrix

    if os.path.isfile('../Uc.dat') and os.path.getsize('../Uc.dat')>0:
        Uc = loadtxt('../Uc.dat')
        for m1 in range(5):
            for m2 in range(5):
                UC[0,m1,m2,m2,m1] = Uc[m1,m2]
                #if abs(UC[0,m1,m2,m2,m1])>1e-3:
                #    print "%2d %2d %2d %2d  %5.2f " % (m1, m2, m2, m1, UC[0,m1,m2,m2,m1])
    else:
        UC[0,:,:,:,:]=0.0
        
    
    # New for self-energy sampling
    UHa=zeros((2*l+1,2*l+1,2*l+1))
    UFo=zeros((2*l+1,2*l+1,2*l+1))
    for m1 in range(2*l+1):
        for m2 in range(2*l+1):
            for m3 in range(2*l+1):
                for k in range(0,l+1):
                    UHa[m1,m2,m3] += real(UC[k,m1,m2,m3,m1])*Fk[k,l]
                    UFo[m1,m2,m3] += real(UC[k,m1,m2,m1,m3])*Fk[k,l]
    #print 'Together='
    for bs1 in baths:
        for bs2 in baths:
            for bs3 in baths:
                m1 = bs1[0]
                s1 = bs1[1]
                m2 = bs2[0]
                s2 = bs2[1]
                m3 = bs3[0]
                s3 = bs3[1]
                Uc = 0.0
                if s2==s3:
                    if s1==s2: Uc = UHa[m1,m2,m3]-UFo[m1,m2,m3] # Equal spins: Hartree and Fock
                    else: Uc = UHa[m1,m2,m3], # Opposite spins: Hartree Only
                #print "%10.6f" % Uc,
            #print
        
    indx={}
    for ni,ns in enumerate(wstates):
        indx[(ns[0],ns[1])] = ni  # index according to N and Sz of the state

    print >> fh_info, 'indx='
    print >> fh_info, indx
    kindx = indx.keys()
    
    print >> fh_info, 'Stage0: Exact diagonalization of the atom'

    mxs = max(map(max,Sigind))
    if len(Eimp)<mxs:
        print 'ERROR: The dimension of the Eimp should be equal to the maximum index of Sigind->',mxs
        sys.exit(1)
    
        
    Eimpc = zeros((2*(2*l+1), 2*(2*l+1)), dtype=complex)
    for ic in range(len(Sigind)):
        Eimpc[ic,ic] = Eimp[Sigind[ic,ic]-1]

    print >> fh_info, 'impurity levels Eimpc0='
    mprint(fh_info, real(Eimpc))
    #Eimpc = matrix(T2C) * matrix(Eimpc) * matrix(T2C).H
    #print 'Eimpc1='
    #mprint(Eimpc)


    Ene=[] # Energy
    Te=[]  # eigenvectors
    S2w=[] # Spin
    for ni,ns in enumerate(wstates):

        #print 'Br:', 'n=', ns[0], 'sz=', ns[1]/2.
        
        print >> fh_info, '----------------------------------------------------'
        print >> fh_info, 'n=', ns[0], 'sz=', ns[1]/2.
        states = ns[2]
        # printing all states in this sector
        print >> fh_info, 'states=',
        for ist,st in enumerate(states): print >> fh_info, ('%d:'%ist),op.printn(st),
        print >> fh_info
        
        S2 = zeros((len(states),len(states)),dtype=complex)
        if CoulombF != 'Ising' and op.Q3d:
            # Computes matrix of S^2
            for js,st in enumerate(states):
                #print st, op.printn(st), "    ",
                stn = op.S2(st)
                #print stn
                for ps in stn:
                    ii = ps[0]
                    iu = states.index(ii)
                    S2[js,iu] += ps[1]
        
        Ham = zeros((len(states),len(states)),dtype=complex)

        for js,st in enumerate(states):
            # on-site Energies in base of cubic harmonics
            # contain crystal-field splittings
            DM = op.DM(st)
            for i in range(len(Eimpc)):
                for j in range(len(Eimpc)):
                    if abs(Eimpc[i,j])<1e-5:continue
                    for p in DM[i][j]:
                        iu = states.index(p[0])
                        Ham[js,iu] += p[1]*Eimpc[i,j]
            
            if CoulombF=='Full':
                ## Coulomb interaction including F2 and F4
                cst = op.CoulombU(st, UC, Fk[:,l])
                for cs in cst:
                    ii = cs[0]
                    U0 = cs[1][0]
                    U1 = cs[1][1]
                    iu = states.index(ii)
                    Ham[js,iu] +=  0.5*U1 # adding only F2,F4,... but not F0
                    Ham[js,iu] +=  0.5*U0 # adding only F2,F4,... but not F0
                    
            elif CoulombF=='Ising':
                ## Coulomb interaction including F2 and F4
                cst = op.CoulombU(st, UC, Fk[:,l],Ising=True)
                for cs in cst:
                    ii = cs[0]
                    U0 = cs[1][0]
                    U1 = cs[1][1]
                    iu = states.index(ii)
                    Ham[js,iu] +=  0.5*U1 # adding only F2,F4,... but not F0
                    Ham[js,iu] +=  0.5*U0 # adding only F2,F4,... but not F0
                    
            elif CoulombF in ['Bulla_Jarrel', 'Oles']:
                occ = op.occup(st)
                # Model Coulomb interaction -J*S^2
                cst = op.S2(st)
                #print 'S2=', cst
                for cs in cst:
                    ii = cs[0]
                    iu = states.index(ii)
                    Ham[js,iu] += -Jc*cs[1]
                nd = sum(occ)
                
                #print 'nd=', nd
                #print 'h1='
                #mprint(Ham)
                
                if CoulombF == 'Bulla_Jarrel':
                    #### Model Coulomb interaction is -J*S^2 + J*N - 1/4*J*N^2
                    Ham[js,js] += Jc*nd*(1-0.25*nd)
                else :
                    ## Model Coulomb interaction is -J*S^2 - J[5/4*N^2-3/2*N-ND] + J D_a^+ D_b
                    ND = op.NDouble(st)
                    # See e.g.: PRB 72, 214431 (2005), Eq.2.5.
                    Ham[js,js] += Jc*(-5/4.*nd**2 + 2.*nd + ND)
                    if ND>0:
                        ph = op.PairHop(st)
                        for p in ph:
                            iu = states.index(p)
                            Ham[js,iu] += Jc
            else:
                print 'Not yet implemented!'
                sys.exit(1)
                
            # Spin-orbit interaction
            if cx>1e-5:
                cst = op.OneBodyNab(st, SO)
                for cs in cst:
                    iu = states.index(cs[0])
                    #if (js>=len(states) or iu>=len(states)): print 'Tezave!'
                    Ham[js,iu] += cs[1]*cx
            
        #if (cx>1e-5): cprint(Ham)
        #else:
        #print >> fh_info, 'H='
        #mprint(fh_info, Ham)
        
        if CoulombF != 'Ising':
            eig = Diagonalize(Ham, small, fh_info)  # Block diagonalization is better!!
            Ex = eig[0]
            Tx = eig[1]
            Ene.append( Ex )
            Te.append( Tx )
        else:
            Ex = [real(Ham[i,i]) for i in range(len(Ham))]
            Tx = eye(len(Ham),len(Ham))
            Ene.append( Ex )
            Te.append( Tx )
        
        if CoulombF != 'Ising' and op.Q3d:
            # Here we compute matrix of S^2 in eigenbase. Should be diagonal if no spin-orbit coupling
            S2e = matrix(conj(Tx.transpose())) * S2 * matrix(Tx)

            printS=False
            trouble=[]
            for i0 in range(shape(S2e)[0]):
                for i1 in range(shape(S2e)[1]):
                    if i0!=i1 and abs(S2e[i0,i1])>1e-6 :
                        print 'WARNING: Troubles->', i0, i1, S2e[i0,i1]
                        printS=True
                        trouble.append(i0)
            
            printS=False # BRISISSS
            if printS:
                print >> fh_info, 'S2='
                mprint(fh_info, S2e)
                print >> fh_info, 'H='
                cprint(fh_info, Ham)
                for it,t in enumerate(trouble):
                    print >> fh_info, 'A[%d]=' % t
                    print >> fh_info, Tx[t]
            
            S2w.append([0.5*int(round(-1+sqrt(1+4*S2e[i,i].real))) for i in range(len(S2e))]) # Spin is computed using formula s(s+1)
        else:
            S2w.append( [0 for i in range(len(S2))] )
        

        print >> fh_info, 'E=', '%f '*len(Ex) % tuple(Ex)
        


    #print 'kindx=', kindx
    #print 'wstates=', wstates

    print 'Br:', 'kindx=', kindx
    
    # Here we create index for psi^dagger
    iFi = zeros((len(wstates),len(baths)),dtype=int)
    for ni,ns in enumerate(wstates):
        for ib,be in enumerate(baths):
            if op.Q3d:
                st = (ns[0]+1, ns[1] + be[1])  # (n+1,sz+s)
                if st in kindx:
                    iFi[ni,ib] = indx[st]
            else:
                st = (ns[0]+1, 0)  # (n+1,sz+s)
                if st in kindx:
                    iFi[ni,ib] = indx[st]
    
    wgr=[]
    for iw in range(len(wstates)): wgr.append([])

    print 'Br:', 'Stage1: Computing F^ in direct base'
    print >> fh_info, 'Stage1: Computing F^ in direct base'
    
    # Below we compute matrix elements of F^ (FKP)
    kindx = indx.keys()
    FKP = []
    for ni,ns in enumerate(wstates):
        states = ns[2]
        bFp=[]

        if CoulombF == 'Ising':
            wgr[ni] += [[ist] for ist in range(len(states))]
            
        for ib,wib in enumerate(baths):
            inew = iFi[ni,ib]
            if inew==0:
                bFp.append([])
                continue
            
            newstates = wstates[inew][2]

            Fp = zeros((len(states), len(newstates)), dtype=complex)

            for js,st in enumerate(states):
                (newst, sig) = op.Fp(st, ib)
                if newst>0:
                    ii = newstates.index(newst)
                    Fp[js,ii] += sig
                #print 'state=', st, newst, ii

            if CoulombF == 'Ising':
                bFp.append(Fp)
            else:
                
                Fpn = matrix(conj(Te[ni].transpose())) * matrix(Fp) * matrix(Te[inew])
                
                # Set to zero small values
                for i0 in range(shape(Fpn)[0]):
                    for i1 in range(shape(Fpn)[1]):
                        if abs(Fpn[i0,i1])<small: Fpn[i0,i1]=0.0
                
                gr = analizeGroups(Fpn, small)
                
                # |inew> = F^+ |ni>
                wgr[ni] += gr[0]   # which states are coupled by F^ in |ni>
                wgr[inew] += gr[1] # which states are coupled by F^ in |inew>
            
                bFp.append(Fpn)
            
        FKP.append(bFp)

    #FKP created!
    
    print 'Br:', 'Stage2: Compressing F^+ according to its block diagonal form'
    print >> fh_info, 'Stage2: Compressing F^+ according to its block diagonal form'
    
    for i in range(len(wstates)):
        wgr[i] = compress(wgr[i])
        print >> fh_info, i+1, wgr[i]
    
    print >> fh_info, 'Stage3: Renumbers states -- creates superstates for ctqmc'
    print 'Br:', 'Stage3: Renumbers states -- creates superstates for ctqmc'

    # Here we store ground state and N for each superstates to be used later for sorting
    tstates=[]
    for iw in range(len(wstates)):
        Nt = sum(wstates[iw][0])
        for ip in range(len(wgr[iw])):
            Eg = Ene[iw][wgr[iw][ip][0]]
            tstates.append([iw,ip,Nt,Eg])
    tstates.sort(comp)
    # tstates contains [index-to-wstates, index-to-state-inside-wstates, N, E]
    
    # superstates == pseudostates are defined
    pseudostates=[]
    indpseudo={}
    jj=0
    for st in tstates:
        iw = st[0]
        ip = st[1]
        pseudostates.append([iw,ip])
        indpseudo[(iw,ip)] = jj
        jj+=1


    iFi_inside=[]
    for iw in range(len(wstates)):
        biFi=[]
        for ib in range(len(baths)):
            ifi = iFi[iw,ib]
            if ifi>0:
                fpair = coupled(FKP[iw][ib], wgr[iw], wgr[ifi], small)
                biFi.append(fpair)
            else:
                biFi.append([])
        iFi_inside.append(biFi)
        

    # creates arrays containing Energy, occupancy and index table for all superstates
    iFinal = zeros((len(pseudostates),len(baths)),dtype=int)
    Enes = []
    S2ws = []
    Occ=[]
    for ii,iwp in enumerate(pseudostates):
        iw = iwp[0]
        ip = iwp[1]
        wstate = wstates[iw]
        group = wgr[iw][ip]

        for ib in range(len(baths)):
            ifi = iFi[iw,ib]
            ifinal = -1
            if (ifi>0):
                ifi_ins = iFi_inside[iw][ib][ip]
                if ifi_ins>=0:
                    ifinal = indpseudo[(ifi,ifi_ins)]
            iFinal[ii,ib] = ifinal

        Ens=[]
        occ=[]
        S2s=[]
        for iq,q in enumerate(group):
            Ens.append(Ene[iw][q])
            occ.append(wstate[0])
            S2s.append(S2w[iw][q])
            
        Enes.append(Ens)
        Occ.append(occ)
        S2ws.append(S2s)

    #print 'pseu=', pseudostates
    #print 'Enes=', Enes
    #print 'Occ=', Occ
    #print 'S2=', S2ws

    print >> fh_info, 'Stage4: F^dagger matrices between superstates evaluated'
    print 'Br:', 'Stage4: F^dagger matrices between superstates evaluated'
    
    # creates small F^dagger matrices between superstates
    maxs = 0
    rFKP = []
    rNn=[]  # This is the matrix F*F^ == 1-N

    for ii,iwp in enumerate(pseudostates):
        iw = iwp[0]
        ip = iwp[1]
        wstate = wstates[iw]
        group = wgr[iw][ip]

        bNn=[]
        bFKP=[]
        for ib in range(len(baths)):
            ifi = iFi[iw,ib]
            ifinal = iFinal[ii,ib]
            if (ifi>0): ifi_ins = iFi_inside[iw][ib][ip]

            if ifinal>=0:
                M = zeros((len(group),len(wgr[ifi][ifi_ins])),dtype=complex)
                #Nn =zeros(len(group),dtype=float) 
                for ii0,i0 in enumerate(group):
                    for jj0,j0 in enumerate(wgr[ifi][ifi_ins]):
                        M[ii0,jj0] = FKP[iw][ib][i0,j0]
                        
                    #Nn[ii0] = sum(map(lambda x: x**2, M[ii0]))
                if max(shape(M)) > maxs : maxs = max(shape(M))

                Nn = zeros((len(group),len(group)))
                for ii0,i0 in enumerate(group):
                    for ii1,i1 in enumerate(group):
                        Nn[ii0,ii1]=sum(M[ii0]*M[ii1]).real
                        
                #print 'ii=', ii, 'ib=', ib, 'ifinal=', ifinal, 'M=', M, 'Nn=', Nn, 'Nt=', Nt
                
            else:
                M = array([])
                Nn = array([])
            bFKP.append(M)
            bNn.append(Nn)
        rFKP.append(bFKP)
        rNn.append(bNn)


    ############################################################################################
    ######## The part of the code between this symbols,  generates input for OCA solver  #######
    ############################################################################################
    if (OCA_G):
        print >> fh_info, 'Stage5: Renumber states for oca - each atomic states has unique number'
        # renumbers states such that each of 1024 states has unique index
        # also remembers energy and N for each state
        (puniq, ipE, ipN, ipS) = RenumberStates(pseudostates, Enes, wstates, S2ws)
        
        # bubbles will be accessed by bubbles[ib][ii][ij]
        # where ib is integer running over all baths, ii is integer running over all possible states
        # and ij is integer which can be accessed by F^{+,ib}|ii>
        bubbles = CreateEmpty2D_Dict(len(kbth), len(ipE))
        FpF = zeros((len(ipE),len(kbth)))
        smallb=1e-4
        for ib,bs in enumerate(kbth): # over all different baths
            bubl=bubbles[ib]
            for ii,iwp in enumerate(pseudostates): # over all pseudostates
                for iib in bs: # over equivalent baths
                    ifinal = iFinal[ii,iib]
                    if ifinal>=0: 
                        dims = shape(rFKP[ii][iib])
                        for i0 in range(dims[0]):
                            istr = puniq[(ii,i0)]
                            for j0 in range(dims[1]):
                                iend = puniq[(ifinal,j0)]
                                FpF[istr][ib] += abs(rFKP[ii][iib][i0,j0])**2
                                if (abs(rFKP[ii][iib][i0,j0])>smallb):
                                    if bubl[istr].has_key(iend):
                                        bubl[istr][iend] += abs(rFKP[ii][iib][i0,j0])**2
                                    else:
                                        bubl[istr][iend] = abs(rFKP[ii][iib][i0,j0])**2
                                    
        
        (equiv, iequiv) = EquivalentStates(ipE, ipN)
        # Now we have all the bubbles.
        # We need to find which states are equivalent
        for tt in range(Nitt):
        
            # We want to merge bubbles that we believe are equivalent
            ebubbles=[]
            for ib,bs in enumerate(kbth):
                bubl = bubbles[ib]
            
                nbubl=[]
                for i in range(len(bubl)): nbubl.append({})
                
                for i0 in range(len(bubl)):
                    for i1 in bubl[i0].keys():
                        if nbubl[i0].has_key(equiv[i1]):
                            nbubl[i0][equiv[i1]] += bubl[i0][i1]
                        else:
                            nbubl[i0][equiv[i1]] = bubl[i0][i1]
                ebubbles.append(nbubl)
            
            # Here we check if the states, which we idenfified above as equivalent,
            # really are equivalent, i.e., have the same type of bubble
            new_iequiv=[]
            back_bubs=[]
            for ii in iequiv:
                cbubs=[]
                for ij in ii:
                    cbub={}
                    for ib in range(len(kbth)):
                        cbub.update(ebubbles[ib][ij])
                    cbubs.append(cbub)
        
                abubs = AverageBubbles([[ebubbles[ib][ij] for ij in ii] for ib in range(len(kbth))])
                
                ieqs = VEquivalentStates(cbubs,ii)
                
                back_bubs.append(abubs)
                
                new_iequiv += ieqs
                
            
            new_equiv=range(len(equiv))
            for i,ii in enumerate(new_iequiv):
                for j in ii: new_equiv[j]=i
                
        
            if len(iequiv)==len(new_iequiv) or tt+1==Nitt: break
            equiv = new_equiv
            iequiv = new_iequiv
        
        print >> fh_info, 'before qOCA'
        
        if qOCA:  # Here we add second order OCA diagrams to NCA bubbles
            # Second order diagramms
            # OCAdiag will be accessed by OCAdiag[ib1][ib2][ii][ij]
            # where ib1, ib2 is integer running over all baths, ii is integer running over all possible states
            # and ij is integer which can be accessed by F^{+,ib}|ii>
            OCAdiag = CreateEmpty3D_Dict(len(kbth), len(kbth), len(iequiv))
                
            for ib1,bs1 in enumerate(kbth):     # over all baths ones
                for ib2,bs2 in enumerate(kbth): # over all baths twice
                    OCAs=OCAdiag[ib1][ib2]
                    for ii,iwp in enumerate(pseudostates): # over all pseudostates
                        for iib1 in bs1:        # over equivalent baths ones
                            for iib2 in bs2:    # over equivalent baths twice
                                ifinal_j = iFinal[ii,iib1]
                                ifinal_l = iFinal[ii,iib2]
                                if ifinal_j>=0 and ifinal_l>=0:
                                    ifinal_k = iFinal[ifinal_j,iib2]
                                    ifinal_k2 = iFinal[ifinal_l,iib1]
                                    if ifinal_k>=0 and ifinal_k2>=0 and ifinal_k==ifinal_k2:
                                        Fij = rFKP[ii][iib1]
                                        Fil = rFKP[ii][iib2]
                                        Fjk = rFKP[ifinal_j][iib2]
                                        Flk = rFKP[ifinal_l][iib1]
                                        (dims_i, dims_j)  = shape(Fij)
                                        (dims_i2,dims_l)  = shape(Fil)
                                        (dims_j2, dims_k) = shape(Fjk)
                                        (dims_l2,dims_k2) = shape(Flk)
                                        if dims_i != dims_i2: print 'Troubles i'
                                        if dims_j != dims_j2: print 'Troubles j'
                                        if dims_l != dims_l2: print 'Troubles l'
                                        if dims_k != dims_k2: print 'Troubles k'
                                        
                                        for i0 in range(dims_i):
                                            iu = equiv[puniq[(ii,i0)]]
                                            for j0 in range(dims_j):
                                                ju = equiv[puniq[(ifinal_j,j0)]]
                                                if (abs(Fij[i0,j0])<smallb): continue
                                                for k0 in range(dims_k):
                                                    ku = equiv[puniq[(ifinal_k,k0)]]
                                                    if (abs(Fjk[j0,k0])<smallb): continue
                                                    for l0 in range(dims_l):
                                                        lu = equiv[puniq[(ifinal_l,l0)]]
                                                        if (abs(Fil[i0,l0])<smallb): continue
                                                        if (abs(Flk[l0,k0])<smallb): continue
                                                        contr = -Fij[i0,j0]*Fjk[j0,k0]*Flk[l0,k0]*Fil[i0,l0]
                                                        akey = (ju,ku,lu)
                                                        if OCAs[iu].has_key(akey):
                                                            OCAs[iu][akey] += contr
                                                        else:
                                                            OCAs[iu][akey] = contr
            # OCA diagramms are renumbered to (i0,i1,i2,i3,ib1,ib2) where
            # i0,i1,i2,i3 are atomic states involved in the diagram and ib1,ib2 are the two bath
            # propagators in the diagram.
            OCAf={}
            for i in range(len(iequiv)):
                for ib1 in range(len(kbth)):     
                    for ib2 in range(len(kbth)): 
                        for ks in OCAdiag[ib1][ib2][i].keys():
                            if abs(OCAdiag[ib1][ib2][i][ks])<1e-10: continue
                            if (ib2<=ib1):
                                new_key = (i,) + ks + (ib1,ib2)
                                OCAf[new_key] = OCAdiag[ib1][ib2][i][ks]
                            else: # Due to time invariance, some diagrams are equivalent. For example
                                  #  0 (b1) 1 (b2) 4 (b1) 2 (b2) 0    and      0 (b2) 2 (b1) 4 (b2) 1 (b1) 0
                                new_key = (i,) + (ks[2],ks[1],ks[0]) + (ib2,ib1)
                                OCAf[new_key] = OCAdiag[ib1][ib2][i][ks]
        
        # Here we regroup N_a =F^+ F  in more convenient way    
        rFpF = zeros((len(iequiv),len(kbth)))
        for i in range(len(equiv)):
            df = FpF[0]-FpF[i]
            for j in range(len(df)):  # if we don't do that, the occupancy can be annoyingly negative
                if abs(df[j])<1e-10: df[j]=0
            rFpF[equiv[i]] += df
        for i,ii in enumerate(iequiv):
            rFpF[i]/=len(iequiv[i]) # it has to be average, not the sum        
        
        # Bubble contains diagrams named b (back). We need to construct from these also the other diagramms
        # which are called f (forward).
        forw_bubs = CreateEmpty2D_Dict(len(iequiv), len(kbth))
        for i in range(len(iequiv)):
            for b in range(len(kbth)):
                for ks in back_bubs[i][b]:
                    forw_bubs[ks][b][i]= back_bubs[i][b][ks]
        
        # We want to have energy of each OCA-pseudoparticle ready
        Eq=zeros(len(iequiv))
        Egs=zeros(len(iequiv))
        Nq=zeros(len(iequiv),dtype=int)
        Nc0=0; Eg0=ipE[0]
        for i,ii in enumerate(iequiv):
            Nq[i] = ipN[ii[0]]
            Eq[i] = ipE[ii[0]]
            Egs[i] = Eg0   # ground state in this sector
            if (Nq[i]>Nc0):
                Nc0=Nq[i]
                Eg0=Eq[i]
                Egs[i]=Eg0
        
        # last renumbering for printing!
        # Only occupancies in n=[....] need to be kept. The rest of the diagrams is ignored.
        pu=-ones(len(iequiv),dtype=int)
        pl=0
        for i,ii in enumerate(iequiv):
            if Nq[i] in n:
                pu[i]=pl
                pl+=1
        
        # Printing output for OCA
        foca=open('out.cix', 'w')
        print >> foca, '# Input file for OCA impurity solver.', 'l=', l, 'J=', Jc, 'Eimp=', Eimp, 'c=', cx,  'mOCA=', mOCA,  'Eoca=', Eoca
        print >> foca, len(kbth), (("%d "*len(kbth)) % tuple(map(len,kbth))), pl, 0,
        print >> foca, '# Number of baths it\'s degeneracy and number of local valence and local core states'
        
        print >> foca, '%3s' % '#',
        print >> foca, ("%6s")*len(kbth) % tuple(map(lambda x: 'N'+str(x), range(len(kbth)))),
        print >> foca, "%4s" % 'Mtot', '%4s' % 'deg', '%10s' % 'Eatom',
        print >> foca, ("%3s"%'#b')*len(back_bubs[i]),
        print >> foca, ("%3s"%'#f')*len(forw_bubs[i])
        
        for i,ii in enumerate(iequiv):
            if pu[i] <0 : continue  # state not used
            
            print >> foca, "%3d" % pu[i], (("%6.2f")*len(kbth) % tuple(rFpF[i])), "%4d" % Nq[i], #ipN[ii[0]],
            print >> foca, "%4d" % len(ii),
        
            Eatom = ipE[ii[0]]  # This is the atomic energy
            for ib in range(len(kbth)): Eatom -= rFpF[i][ib]*Eimp[ib]  # This part will be added back inside the impurity solver, therefore the energy should be subtracted
            
            print >> foca, "%10.4f" % Eatom,
        
            for b in range(len(kbth)): # delete diagrams which include states that were removed
                for ks in back_bubs[i][b].keys():
                    if pu[ks]<0: # this diagram involves state not considered
                        del back_bubs[i][b][ks]
                for ks in forw_bubs[i][b].keys():
                    if pu[ks]<0: # this diagram involves state not considered
                        del forw_bubs[i][b][ks]
                
            
            print >> foca, ("%3d"*len(back_bubs[i])) % tuple(map(len,back_bubs[i])),
            print >> foca, ("%3d"*len(forw_bubs[i])) % tuple(map(len,forw_bubs[i])),
            print >> foca, '  ',
            
            for b in range(len(kbth)):
                for ks in back_bubs[i][b]:
                    print >> foca, "%6.2f x %-3d  " % (back_bubs[i][b][ks], pu[ks]),
                    
            for b in range(len(kbth)):
                for ks in forw_bubs[i][b]:
                    print >> foca, "%6.2f x %-3d  " % (forw_bubs[i][b][ks], pu[ks]),
        
            print >> foca, '# S ', ipS[ii[0]], ' Eatom=', ipE[ii[0]]
        
        if qOCA:
            print >> foca, '# OCA diagrams, information is (pp0,pp1,pp2,pp3) (b1,b2) fact , where pp is pseudoparticle and b is bath'
            OCAF = OCAf.items()
            OCAF.sort(lambda x,y: cmp(y[1],x[1]))
            for i in range(len(OCAF)):
                excitedE = [Eq[j]-Egs[j] for j in OCAF[i][0]]
                states_involved = [pu[l] for l in OCAF[i][0][:4]]
                #print states_involved
                if (-1 in states_involved): continue  # One of the states is not considered
                if max(excitedE)>Eoca:  continue      # We take it into account only if all states that are involved, have energy close to the ground state energy for this occupancy
                if abs(OCAF[i][1])<mOCA: continue     # Matrix element negligible
                if not (Nq[OCAF[i][0][1]] in Ncentral): continue
                
                print >> foca, "%3d %3d %3d %3d   " % tuple(states_involved), #tuple([pu[l] for l in OCAF[i][0][:4]]),
                print >> foca, "%2d %2d" % tuple(OCAF[i][0][4:]),
                print >> foca, real(OCAF[i][1]), 
                print >> foca, '   #', [Eq[j]-Egs[j] for j in OCAF[i][0]]
                
        ############################################################################################
        ######## End of the part which generates input for OCA solver                        #######
        ############################################################################################


    # Extract low energy states
    lowE=[]
    low_maxsize=0
    for ii in range(len(Enes)):
        size=0
        plowE=[]
        for iq in range(min(len(Enes[ii]),max_M_size)):
            if Enes[ii][iq]>=Ewindow[0] and Enes[ii][iq]<Ewindow[1] and iq<Nmax:
                plowE.append(iq)
                size += 1
        if size>low_maxsize: low_maxsize = size
        if len(plowE)>0:
            lowE.append((ii,plowE))

    # Creates index array between all states and low energy ones
    inv_lowE1={-1:-1}
    for i in range(len(pseudostates)): inv_lowE1[i]=-1
    for i in range(len(lowE)):
        ii = lowE[i][0]
        inv_lowE1[ii]=i


    fcix = open('actqmc.cix', 'w')
    # ---------------------------------------------------------------------------------------
    # -------------- Below is printing for ctqmc  solver ------------------------------------
    # ---------------------------------------------------------------------------------------
    print >> fcix, '# CIX file for ctqmc! '
    print >> fcix, '# cluster_size, number of states, number of baths, maximum_matrix_size'
    print >> fcix, 1, len(lowE), len(bkeep), low_maxsize
    print >> fcix, '# baths, dimension, symmetry'
    
    for ib in range(len(bkeep)):
        print >> fcix, ib, '  ', 1, Sigind[bkeep[ib],bkeep[ib]]-1, '  ', global_flip[bkeep[ib]]
    
    print >> fcix, '# cluster energies for non-equivalent baths, eps[k]'
    for E in tEimp: print >> fcix, E,
    print >> fcix
    print >> fcix, '#   N   K   Sz size'

 
    for i in range(len(lowE)):
        ii = lowE[i][0]
        iwp = pseudostates[ii]
        iw = iwp[0]
        ip = iwp[1]
        wstate = wstates[iw]

        if op.Q3d:
            Mz = sum(wstate[1])/2.
            gs=wstate[2][ip] 
            Mz = op.OrbDiff(gs)    #### Here we compute N_{xz}-N_{yz}
        else:
            gs=wstate[2][ip]
            Mz = op.Mz(gs)
        
        print >> fcix, "%3d  %2d %2d %6.3f %2d " % (i+1, sum(wstate[0]), 0, Mz, len(lowE[i][1])),#len(Enes[ii])),
        for ib in bkeep:
            ifinal = iFinal[ii,ib]
            print >> fcix, "%3d" % (inv_lowE1[ifinal]+1),
        print >> fcix, "  ",
        for iq in lowE[i][1]:
            print >> fcix, "%10.6f" % (Enes[ii][iq],),
        print >> fcix, "  ",
        for iq in lowE[i][1]:
            print >> fcix, S2ws[ii][iq],
        if CoulombF == 'Ising':
            print >> fcix, "  # ", op.printn(gs),
        print >> fcix
        
    print >> fcix, '# matrix elements'

    for i in range(len(lowE)):
        ii = lowE[i][0]
        for ib in bkeep:
            ifinal = iFinal[ii,ib]
            low_ifinal = inv_lowE1[ifinal]
            print >> fcix, "%3d %3d " % (i+1, low_ifinal+1),
            if low_ifinal>=0: 
                ind0 = lowE[i][1]
                ind1 = lowE[low_ifinal][1]
                print >> fcix, "%2d %2d" % (len(ind0), len(ind1)), 
                for i0 in ind0:
                    for j0 in ind1:
                        x = rFKP[ii][ib][i0,j0]
                        if abs(x.imag)<1e-4 or PrintReal:
                            print >> fcix, x.real,
                        else:
                            print >> fcix, x,
            else:
                print >> fcix, "%2d %2d" % (0, 0),
            print >> fcix

    if HB2 : print >> fcix, 'HB2'
    else: print >> fcix, 'HB1'

    if (HB2):
        print >> fcix, "# Uc = U[m1,m2,m3,m1]-U[m1,m2,m1,m3] ; loops [m1,m2,m3]"
        for bs1 in baths:
            for bs2 in baths:
                for bs3 in baths:
                    m1 = bs1[0]
                    s1 = bs1[1]
                    m2 = bs2[0]
                    s2 = bs2[1]
                    m3 = bs3[0]
                    s3 = bs3[1]
                    Uc = 0.0
                    if s2==s3:
                        if s1==s2: Uc = UHa[m1,m2,m3]-UFo[m1,m2,m3] # Equal spins: Hartree and Fock
                        else: Uc = UHa[m1,m2,m3], # Opposite spins: Hartree Only
                    print >> fcix, "%10.6f" % Uc,
                print >> fcix

    print >> fcix, '# number of operators needed'

    if not add_occupancy:
        print >> fcix, '0'
    else:
        print >> fcix, '1'
        print >> fcix, '# Occupancy '
        
        for i in range(len(lowE)):
            ii = lowE[i][0]
            ind0 = lowE[i][1]

            #tkbth = kbth
            #if HB2: tkbth=kbth0
            
            for ikb,bt in enumerate(tkbth):
                Oub = zeros((len(ind0),len(ind0)),dtype=float)
                for ib in bt:
                    Nm = zeros((len(ind0),len(ind0)),dtype=float)
                    if len(rNn[ii][ib])>0:  
                        #Nm = rNn[ii][ib]
                        for j in range(len(ind0)):
                            for k in range(len(ind0)):
                                Nm[j,k] = rNn[ii][ib][ind0[j],ind0[k]]
                                
                    Oub += identity(len(ind0))-Nm
                print >> fcix, ("%3d " % (i+1)),
                print >> fcix, "%2d %2d" % (len(ind0), len(ind0)), 
                for iw,i0 in enumerate(ind0):
                    for iz,j0 in enumerate(ind0):
                        ff = Oub[iw,iz]
                        if abs(ff)<small: ff=0.0
                        print >> fcix, ff,
                print >> fcix
    
    
    print >> fcix, '# Data for HB1'
    
    print >> fcix, 1, len(pseudostates), len(bkeep), maxs
    print >> fcix, '# ind   N   K   Jz size'
    
    for ii,iwp in enumerate(pseudostates):
        iw = iwp[0]
        ip = iwp[1]
        wstate = wstates[iw]
        print >> fcix, "%3d  %3d  %2d %2d %4.1f %2d " % (ii+1, inv_lowE1[ii]+1, sum(wstate[0]), 0, sum(wstate[1])/2., len(Enes[ii])),
        for ib in bkeep:
            print >> fcix, "%3d" % (iFinal[ii,ib]+1),
        print >> fcix, "  ",
        for iq in range(len(Enes[ii])):
            print >> fcix, Enes[ii][iq],
        print >> fcix, "  ",
        for iq in range(len(Enes[ii])):
            print >> fcix, 0,
        print >> fcix, "  ",
        print >> fcix
        
    print >> fcix, '# matrix elements'

    for ii in range(len(pseudostates)):
        for ib in bkeep:

                print >> fcix, "%3d %3d " % (ii+1, iFinal[ii,ib]+1), 

                ffp = zeros(len(Enes[ii]),dtype=float)

                
                if iFinal[ii,ib]>=0:
                    (dim0, dim1) = shape(rFKP[ii][ib])
                    print >> fcix, "%2d %2d" % (dim0,dim1), 
                    for i0 in range(dim0):
                        for j0 in range(dim1):
                            x = rFKP[ii][ib][i0,j0]
                            if abs(x.imag)<1e-4 or PrintReal:
                                print >> fcix, x.real,
                            else:
                                print >> fcix, x,
                            #print >> fcix, rFKP[ii][ib][i0,j0],



                    
                    for i0 in range(dim0):
                        dsum=0
                        for j0 in range(dim1):
                            dsum += abs(rFKP[ii][ib][i0][j0])**2
                        ffp[i0] += dsum

                        
                    
                else:
                    print >> fcix, "%2d %2d" % (0, 0),
                print >> fcix

            
        
