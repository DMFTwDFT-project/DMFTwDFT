#!/usr/bin/env python2
from scipy import *
from scipy import linalg
from mysub import *

def Make_orthogonal(Egval,EgvecL,EgvecR,smalleps=1e-5):
   for i in range(len(Egval)):
      for j in range(i):
         if abs(Egval[i]-Egval[j])<smalleps and dot(conjugate(EgvecL[:,j]),EgvecR[:,i])>smalleps:
            print "EigenvectorR are not orthogonal!"
            EgvecR[:,i]=EgvecR[:,i]-dot(conjugate(EgvecL[:,j]),EgvecR[:,i])/dot(conjugate(EgvecL[:,j]),EgvecR[:,j])*EgvecR[:,j]
      for j in range(i):
         if abs(Egval[i]-Egval[j])<smalleps and dot(conjugate(EgvecL[:,i]),EgvecR[:,j])>smalleps:
            print "EigenvectorL are not orthogonal!"
            EgvecL[:,i]=EgvecL[:,i]-dot(EgvecL[:,i],conjugate(EgvecR[:,j]))/dot(EgvecL[:,j],conjugate(EgvecR[:,j]))*EgvecL[:,j]
   ctmp=sqrt(diag(matrix(EgvecL).H*matrix(EgvecR)))
   for i in range(len(Egval)):
      EgvecL[:,i]=EgvecL[:,i]/conjugate(ctmp[i])
      EgvecR[:,i]=EgvecR[:,i]/ctmp[i]
   return EgvecL,EgvecR 

  
if __name__=='__main__':
   a=(arange(64)+1j*arange(64)).reshape(8,8)
   a[0,0]+=10
   eval,evecL,evecR=linalg.eig(a,left=True)
   evecL,evecR=Make_orthogonal(eval,evecL,evecR) 
   pmatrix(array(matrix(evecR)*matrix(evecL).H).real)
   pmatrix(array(matrix(evecR)*diag(eval)*matrix(evecL).H).imag)
