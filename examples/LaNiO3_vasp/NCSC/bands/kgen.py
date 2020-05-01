#!/usr/bin/env python
from scipy import *

def Create_kpath(KPoints,nk_band):
   def dist(a,b): return sqrt(sum((array(a)-array(b))**2))
   # returns the distance of given a and b points
   #print KKPoints
   path_len=[]
   for i in range(len(KPoints)-1):
      path_len.append(dist(KPoints[i+1],KPoints[i]))
   #path_len=[1.0,1.0,sqrt(2.0)]
   path_nk=map(int,nk_band*array(path_len)/sum(path_len))
   print path_nk
   klist=[]; dist_K=[0.]; dist_SK=[0.]
   for i,nkk in enumerate(path_nk):
      for n in range(nkk):
         klist.append(KPoints[i]+(KPoints[i+1]-KPoints[i])*n/nkk)
         if len(klist)>1: dist_K.append(dist_K[-1]+dist(klist[-1],klist[-2]))
      dist_SK.append(dist_SK[-1]+path_len[i])
   # Add the ending point
   klist.append(KPoints[-1])
   dist_K.append(dist_K[-1]+dist(klist[-1],klist[-2]))
   return array(klist), array(dist_K), array(dist_SK)

if __name__=='__main__':
   KPoints=[[0,0,0],[0.5,0.5,0],[0.5,0,0],[0,0,0]]
   KPoints=[matrix([[0.66972822, -1.16000331, 0.48613733],[0.66972822, 1.16000331, 0.48613733],[-1.33945645, 0.00000000, 0.48613733]]).dot(array(KPoints[i])) for i in range(len(KPoints))]
   KPoints=array(KPoints)
   SKPoints=['$\Gamma$','X','M','$\Gamma$']
   nk_band=500
   klist, dist_K, dist_SK = Create_kpath(KPoints,nk_band)
   klist=[matrix([[0.66972822, -1.16000331, 0.48613733],[0.66972822, 1.16000331, 0.48613733],[-1.33945645, 0.00000000, 0.48613733]]).I.dot(array(klist[i][0,:])) for i in range(len(klist))]
   #print klist,dist_K
   fi=open('klist.dat','w')
   for i in range(nk_band):
      kcheck=0
      for j,d in enumerate(dist_SK):
         if abs(dist_K[i]-d)<1e-4: 
            print >>fi, '%.14f  %.14f  %.14f  %.14f  %s' %(dist_K[i],klist[i][0,0],klist[i][0,1],klist[i][0,2],SKPoints[j])
            kcheck=1
            break
      if kcheck==0:
         print >>fi, '%.14f  %.14f  %.14f  %.14f' %(dist_K[i],klist[i][0,0],klist[i][0,1],klist[i][0,2])

   #print klist,dist_K,dist_SK
