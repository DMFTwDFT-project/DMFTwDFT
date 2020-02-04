#!/usr/bin/env python2
# @Copyright 2007 Kristjan Haule
# 

from sys import *
from scipy import *
from pylab import *
from scipy import random
from scipy import interpolate
from scipy import integrate
from scipy import optimize
import os
import subprocess

def GiveTanMesh(x0,L,Nw):
    def fun(x,x0,L,Nw):
        "x[0]=d, x[1]=w"
        d=x[0]
        w=x[1]
        #print 'd=', d, 'w=', w
        return array([L-w/tan(d), x0-w*tan(pi/(2*Nw)-d/Nw) ])
    
    xi=x0/L
    d0 = Nw/2.*(tan(pi/(2*Nw))-sqrt(tan(pi/(2*Nw))**2 - 4*xi/Nw))
    w0 = L*d0

    sol=optimize.root(fun, [d0,w0], args=(x0,L,Nw) )
    (d,w) = sol.x
    om = w*tan(linspace(0,1,2*Nw+1)*(pi-2*d) -pi/2+d)
    return om

def FindHighFrequency(Gm,om,Nf):
    S=0.; Sx=0.; Sy=0.; Sxx=0.; Sxy=0;
    for j in range(len(om)-Nf,len(om)):
        x = om[j]
        y = Gm[j].imag * x
        x2= x**2
        Sy += y
        Sx += 1/x2
        Sxx += 1/x2**2
        Sxy += y/x2
        S += 1

    dd = S*Sxx-Sx**2
    a = (Sxx*Sy-Sx*Sxy)/dd
    bx = (S*Sxy-Sx*Sy)/dd
    ah = -a;
    if abs(ah-1.0)<1e-3: ah=1.0
    return ah


def GiveMatsubaraKernel(om,Niom,omega):
    """ omega     -- real axis mesh
        om[:Niom] -- imaginary axis Matsubara points to be considered
    """
    # Kernel on imaginary axis
    Ker = zeros((Niom,len(omega)),dtype=complex)
    for iw in range(Niom):
        Ker[iw,:]=-(omega+om[iw]*1j)*dom/(omega**2+om[iw]**2)
    return Ker

icounter=0
def MaxentFunctional(A,alpha,dom,Ker,Gm,mod,error):
    global icounter
    
    A *= 1./sum(A*dom)
    gm = dot(Ker,A)
    dg = Gm-gm
    
    Niom = len(gm)
    chi2 = sum( abs(dg)**2/error**2 )/Niom
    Se = sum(dom*A*log(abs(A/mod)))
    L = chi2 + alpha*Se
    
    if (icounter%200==0):
        print icounter, 'chi2=', chi2, 'S=', Se, 'L=', L
    icounter = icounter + 1
    return L

def dMaxentFunctional(A,alpha,dom,Ker,Gm,mod,error):
    A *= 1./sum(A*dom)
    gm = dot(Ker,A)
    dg = Gm-gm

    Niom = len(gm)
    dchi2 = -2*dot(conj(dg)/error**2,Ker).real/Niom
    dSe = dom*(log(abs(A/mod))+1.)
    dL = dchi2 + alpha*dSe

    return dL - dom*sum(A*dL)


def MaxEntropyMatsubara(alpha, om, Gm0, error, mod, dom, Ker, Nf=20):
    
    # The tail, which determines normalization
    ah = FindHighFrequency(Gm0,om,Nf)

    Niom=len(Ker)
    # Normalizing such that G \propto 1/omega
    Gm = Gm0[:Niom]/ah
    
    # We will start with default model
    A = copy(mod)
    # Values of real axis spectra need to be positive
    abounds=[(1e-15,100.) for i in range(len(A))]
    
    # non-linear optimization
    (Af, fmin, dinf) = optimize.fmin_l_bfgs_b(MaxentFunctional, A, fprime=dMaxentFunctional,args=(alpha,dom,Ker,Gm,mod,error),bounds=abounds, factr=1e7)
    icounter=0
    # Checking what are all values at the end.
    MaxentFunctional(Af,alpha,dom,Ker,Gm,mod,error)
    
    # This is our best approximation for G
    gm = dot(Ker, Af)*ah
    # Data normalized
    Gm *= ah
    
    return (Af*ah, gm, ah)


      
if __name__=='__main__':
    # default value of some parameters
    p={'statistics': 'fermi', # fermi/bose
        'x0'        : 0.02,    # low energy cut-off
        'L'         : 40.0,    # cutoff frequency on real axis
        'Nw'        : 250,     # number of frequency points on real axis
        'Niom'      : 250,     # Number of frequency points on imaginary axix
        'gwidth'    : 20.0,    # width of gaussian
        'idg'       : 1,       # error scheme: idg=1 -> sigma=deltag ; idg=0 read from the file
        'deltag'    : 0.001,   # error
        'alpha'     : 0.1,     # How much entropy to add
        'iflat'     : 1,       # iflat=1 : gaussian of width gwidth, iflat=2 : input using file model.dat
        'Nf'        : 20,      # to perform inverse Fourier, high frequency limit is computed from the last Nf points
        'PLOT'      : True,    # whether to plot data with matplotlib
    }

    
    # need 'maxent_paramsi.dat'
    fprs='maxent_paramsi.dat'
    if os.path.exists(fprs):
        execfile(fprs)
        p.update(params)
    
    if len(sys.argv)<2:
        print 'give input file Sigma(iom)'
        sys.exit(0)

    imAxisData = sys.argv[1]

    
    fi = open(imAxisData, 'r')
    firstlines = [fi.next(),fi.next()]
    fi.close()
        
    alpha=p['alpha']
    Niom = p['Niom']

    # Real axis mesh mesh
    omega = GiveTanMesh(p['x0'],p['L'],p['Nw'])
    dom = array([(omega[1]-omega[0])]+[0.5*(omega[i+1]-omega[i-1]) for i in range(1,len(omega)-1)]+[(omega[-1]-omega[-2])])

    # Set the model
    if p['iflat']==0:
        mod = ones(len(omega))
    elif p['iflat']==1:
        mod = exp(-omega**2/p['gwidth']**2)
    else:
        dmod=loadtxt('model.dat').transpose()
        # Bringing model to the current real axis mesh
        fmod = interpolate.UnivariateSpline(dmod[0],dmod[1],s=0)
        mod = fmod(omega)
    
    # normalization of the model
    mod *= 1./sum(mod*dom)

            
    # Reading Im Axis data    
    dat=loadtxt(imAxisData).transpose()
    om=dat[0]
    beta = pi/om[0]
    print 'beta=', beta
    Gmdata=dat[1::2]+dat[2::2]*1j

    # Gets Kernel for maxent
    Ker = GiveMatsubaraKernel(om,Niom,omega)
    
    if p['idg']==1:
        norb = len(Gmdata)
    else:
        norb = len(Gmdata)/2

    Sg_imag=[]
    for iorb in range(norb):
        Gm0 = Gmdata[iorb]
        
        # Setting the error
        if p['idg']==1:
            error = ones(Niom)*p['deltag']
        else:
            error = abs(Gmdata[norb+iorb])[:Niom]
        print 'average-error=', sum(error)/Niom

        # Main call to maximum entropy
        (Af, gm, ah) = MaxEntropyMatsubara(alpha, om, Gm0, error, mod, dom, Ker, p['Nf'])
        
        Sg_imag.append( -Af*pi )
        
        Gm=Gm0[:Niom]
        # Save the result on imaginar axis
        savetxt('siwn.'+str(iorb), vstack((om[:Niom], gm.real, gm.imag)).transpose() )
        # This is the origonal data from the file
        savetxt('siw0.'+str(iorb), vstack((om[:Niom], Gm.real, Gm.imag)).transpose() )
        # The best real axis data
        savetxt('dos.out.'+str(iorb), vstack((omega,Af*pi)).transpose())
        
        if p['PLOT']:
            subplot(2,1,1)
            plot(omega, Af*pi,label='minimization')
            plot(omega, mod*ah*pi,'k:',label='model')
            legend(loc='best')

            subplot(2,1,2)
            plot(om[:Niom], Gm.real,'or',label='data')
            plot(om[:Niom], Gm.imag,'or')
            plot(om[:Niom], gm.real, 'g-', label='fit')
            plot(om[:Niom], gm.imag, 'g-')
            legend(loc='best')
            show()

    # Below we perform Kramars-Kronig transformation
    fh_info = open('log.out','w')
    #exepath = os.environ.get('WIEN_DMFT_ROOT')
    Sigt=[]
    for iorb in range(norb):
        Sgi = Sg_imag[iorb]

        # Removing zero from the mesh, so that Kramars-Kronig does not fail
        izero = argmin(abs(omega))
        if abs(omega[izero])<1e-16:
            omega_n = hstack([omega[:izero],omega[izero+1:]])
            Sgin = hstack([Sgi[:izero],Sgi[izero+1:]])
        else:
            omega_n = omega
            Sgin = Sgi
        
        savetxt('dosn', vstack((omega_n,Sgin)).transpose())
        
        # Performs Kramars-Kronig
        cmd = 'skrams -cn 2 dosn > Sc'
        print cmd
        subprocess.call(cmd,shell=True,stdout=fh_info,stderr=fh_info)
        
        Sdata = loadtxt('Sc').transpose()
        om = Sdata[0]
        Sdata[1]+Sdata[2]*1j
        
        if iorb==0: Sigt.append(om)
        Sigt.append(Sdata[1])
        Sigt.append(Sdata[2])

Sigt = array(Sigt).transpose()
fo = open('Sig.out', 'w')
print >> fo, firstlines[0].strip()
print >> fo, firstlines[1].strip()
for sg in Sigt:
    for b in sg:
        print >> fo, b,
    print >> fo
