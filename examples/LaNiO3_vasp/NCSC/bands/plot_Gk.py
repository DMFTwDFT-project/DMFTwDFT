#!/usr/bin/env python

from scipy import *
import matplotlib
matplotlib.use('ps')
from matplotlib.font_manager import fontManager, FontProperties
from pylab import *
import scipy.interpolate
import re

vmm = [0,15.0]

nk=0
SKP=[]
SKPoints=[]
distk=[]
kpts=[]
fi=open('klist.dat','r')
for line in fi.readlines():
   line=line.split()
   distk.append(float(line[0]))
   kpts.append([float(line[1]),float(line[2]),float(line[3])])
   if len(line)==5:
      SKP.append(float(line[0]))
      SKPoints.append(line[4])
#print distk
#print kpts
print SKP
print SKPoints
fi.close() 

fi=open('EIGENVAL','r')
klist=[]
kxais=[]
for i in range(6):
   fi.readline()
    
jj=0
eigval0=[]
for j in range(48):
   fi.readline()
   line=fi.readline()
   klist.append(map(float,line.split()[0:3]))
   if j>=0 and j<=15:
      kxais.append((SKP[1]-SKP[0])*jj/15.0)
      jj+=1
      eigval0.append([])
      for i in range(80):
         lines=fi.readline().split() 
         eigval0[jj-1].append(float(lines[1]))   
   elif j>=17 and j<=31:
      kxais.append(SKP[1]+(SKP[2]-SKP[1])*(jj-15.0)/15.0)
      jj+=1
      eigval0.append([])
      for i in range(80):
         lines=fi.readline().split() 
         eigval0[jj-1].append(float(lines[1]))   
   elif j>=33 and j<=47:
      kxais.append(SKP[2]+(SKP[3]-SKP[2])*(jj-30.0)/15.0)
      jj+=1
      eigval0.append([])
      for i in range(80):
         lines=fi.readline().split() 
         eigval0[jj-1].append(float(lines[1]))   
   else: 
      for i in range(80):
         lines=fi.readline().split()
      
fi.close()
eigval0=array(eigval0)
#print kxais
#print eigval0
#exit()
fi=open('ksum.input','r')
numk=int(fi.readline())
nom=int(fi.readline())
fi.close()

A_k=[]
dist_k=[]
om=[]
kpts=[]
fi=open('Gk.out','r')
for i in range(numk):
   kpts.append(map(float,fi.readline().split()[1:]))
   A_k.append([]) 
   om.append([]) 
   for j in range(nom):
      line=map(float,fi.readline().split())
      A_k[i].append(-1*line[2]/3.14159265)
      om[i].append(line[0])
   A_k[i]=A_k[i][::-1]
fi.close()
A_k=array(matrix(array(A_k)[::-1]).T)
#print A_k
#exit()

#om=[]
#for line in file.readlines()[::-1]:
#   dat = map(float,line.split())
#   om.append(dat[0])
#   t_A_k.append(dat[1])
#A_k.append(t_A_k)
#
#for i in range(1,num_k+1):
#   file=open('G_k.out.%03d'%(i),'r')
#   t_A_k=[]
#   om=[]
#   for line in file.readlines()[::-1]:
#      dat = map(float,line.split())
#      om.append(dat[0])
#      t_A_k.append(dat[1])
#   A_k.append(t_A_k)
#A_k=array(matrix(array(A_k)).T)
# 
#k_eig=[]
#eigval0=[]
#file=open('eigval0.dat','r')
#for i,line in enumerate(file.readlines()):
#   dat = map(float,line.split())
#   eigval0.append([])
#   k_eig.append(dat[0])
#   for j in range(len(dat)-1):
#      eigval0[i].append(dat[j+1])
#eigval0=array(eigval0)
#k_eig=2*array(k_eig)/(float(k_eig[-1])+0.0)

axhline(y=0,color='white',ls='--')
axvline(x=SKP[1],color='white',ls='--')
axvline(x=SKP[2],color='white',ls='--')

(ymin,ymax) = (om[0][0],om[0][-1])
#print ymin,ymax
(xmin,xmax) = (distk[0],distk[-1])#(0, (shape(A_k)[1]-1)/(numk+0.0)*3.14159265)
#print shape(A_k),xmin,xmax
#vmm = [0,max(map(max,A_k))*itensity]

im=imshow(A_k,cmap=cm.hot, vmin=vmm[0], vmax=vmm[1], extent=[xmin,xmax,ymin,ymax],aspect='auto')
#im=imshow(A_k,cmap=cm.hot, interpolation='bilinear', vmin=vmm[0], vmax=vmm[1], extent=[0.0,10.0,-1.0,1.0])
colorbar(im,orientation='vertical',pad=0.05,shrink=1.0,ticks=arange(0,15.0,1.0))
ax1=subplot(111)
for i in range(15,shape(eigval0)[1]):
   ax1.plot(kxais,eigval0[:,i]-7.5434        ,color='green')
##ax1.plot([k_eig[319],k_eig[319]], [ymin,ymax], 'w-')
#xticks([0,1,2],['$\Gamma$','X','M'])
xticks(SKP,SKPoints)

xlabel('k-path',fontsize='xx-large')
ylabel('Energy',fontsize='xx-large')
ax1.set_ylim(ymin,ymax)
ax1.set_xlim(xmin,xmax)


show()
savefig('A_k.eps')
exit()


