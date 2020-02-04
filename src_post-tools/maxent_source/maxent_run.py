#!/usr/bin/env python
""" This script takes self-energy of DFT+DMFT code,
and performs max-entropy on auxiliary Green's function of the form
Gc(iom) = 1/(iom-Sigma+s_oo)
Then it performs Kramars-Kronig on resulting Dos, to obtain auxiliary G(omega).
From G(omega), one can obtain Sigma(omega) by
Sigma(omega) = omega+s_oo-1/G(omega)
which is written to 'Sig.out'
"""
import sys, shutil, os
from scipy import *
import shutil
from maxentropy import *

if len(sys.argv)<2:
    print 'give input file Sigma(iom)'
    sys.exit(0)

Sfile = sys.argv[1]
fi = open(Sfile, 'r')
firstlines = [fi.next(),fi.next()]
fi.close()

Sdata = loadtxt(Sfile).transpose()
# number of baths
# number of baths
nb2=0; nz2=0;
for i in range(1,len(Sdata)):
    if sum(abs(Sdata[i]))>0:
        nb2 +=1
    else:
        nz2 +=1
nb = nb2/2
nz = nz2/2
print 'nb=',nb, 'nz=',nz


# path to skrams executable
#exepath = os.environ.get('WIEN_DMFT_ROOT')
# need 'maxent_params.dat'
execfile('maxent_params.dat')


iom = Sdata[0]
beta = pi/iom[0]
tau = linspace(0,beta,params['Ntau']+1)

print 'beta=', beta
Sigt=[]
for b in range(nb):
    Gm = 1/(iom*1j-Sdata[1+2*b]-Sdata[2+2*b]*1j)
    
    Gt = InverseFourier(Gm, iom, tau, beta, params['Nf'])
    savetxt('gt0.'+str(b), vstack((tau,Gt)).transpose())
    
    (Aw, omega) = MaximumEntropy(params, tau, Gt)
    savetxt('dos.out.'+str(b), vstack((omega,Aw)).transpose())

    shutil.copy2('gtn', 'gtn.'+str(b))
    
    # removes zero from mesh, because Kramars-Kronig can be be done with zero
    #izero = omega.tolist().index(0)
    izero = argmin(abs(omega))
    if abs(omega[izero])<1e-6:
        omega_n = hstack([omega[:izero],omega[izero+1:]])
        Aw_n = hstack([Aw[:izero],Aw[izero+1:]])
    else:
        omega_b = omega
        Aw_n = Aw
    savetxt('dosn', vstack((omega_n,Aw_n)).transpose())
    
    # Performs Kramars-Kronig
    cmd = 'skrams -cn 2 -s -pi dosn > Gc'
    print cmd
    print os.popen(cmd).read()
    Gdata = loadtxt('Gc').transpose()
    om = Gdata[0]
    Sc = om-1/(Gdata[1]+Gdata[2]*1j)
    savetxt('sig.'+str(b), array([om,real(Sc),imag(Sc)]).transpose())
    if b==0: Sigt.append(om)
    Sigt.append(real(Sc))
    Sigt.append(imag(Sc))

for z in range(nz):
    Sigt.append( zeros(len(Sigt[0])) )
    Sigt.append( zeros(len(Sigt[0])) )

Sigt = array(Sigt).transpose()
fo = open('Sig.out', 'w')
print >> fo, firstlines[0].strip()
print >> fo, firstlines[1].strip()
for sg in Sigt:
    for b in sg:
        print >> fo, b,
    print >> fo
