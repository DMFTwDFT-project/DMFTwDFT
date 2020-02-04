#!/usr/bin/env python
from sys import *
from scipy import *
from pylab import *
from scipy import random
from scipy import interpolate
from scipy import integrate
import maxent_routines as me
import shutil

def Broad(width, om, fw):
    " Broadens the data with gaussian of width=width"
    def MakeTanMesh(N, tanc, tanw, b0, b1):
        if not(b0<b1): print "Relation must hold: b0<b1!"
        if not(b0<tanw and tanw<b1): print "Relation mesu hold: b0<tanw<b1!"
        if not(b0>0): print "b0 must be positive!"
        du = arctan(((tanc-b0)/tanw))
        b1n = arctan((b1-tanc)/tanw)+du
        m0 = [tanc + tanw * tan(b1n*(i-1)/(N-2)-du) for i in range(1,N)]
        return hstack( (-array(m0[::-1]), array([0]+m0) ) )
    
    w=width
    x = MakeTanMesh(200,0.0,w,w/50,w*20)
    fwi = interpolate.interp1d(om, fw)
    fwn=[]
    for im in range(len(om)):
        x1 = filter(lambda t: t>=om[im]-x[-1] and t<=om[im]-x[0], om)
        x2 = filter(lambda t: t>=om[0] and t<=om[-1], x+om[im])
        eps = sorted(hstack((x1, x2)))
        x3 = om[im]-eps
        gs = exp(-x3**2/(2*w**2))/(sqrt(2*pi)*w)
        yn = integrate.trapz(fwi(eps) * gs, x=eps)
        fwn.append(yn)
    return array(fwn)

def MaximumEntropy(p, tau, Gt):

    beta = tau[-1]

    random.seed( 1 ) # seed for random numbers
    
    omega = linspace(-p['L'],p['L'],2*p['Nw']+1)
    f0,f1,f2 = me.initf0(omega)
    fsg=1
    if p['statistics']=='fermi':
        Gt = -Gt
        fsg=-1
        normalization = Gt[0]+Gt[-1]
        Ker = me.initker_fermion(omega,beta,tau)        
    elif p['statistics']=='bose':
        normalization = integrate.trapz(Gt,x=tau)
        Ker = me.initker_boson(omega,beta,tau)
    
    print 'beta=', beta
    print 'normalization=', normalization

    # Set error
    if p['idg']:
        sxt = ones(len(tau))/(p['deltag']**2)
    else:
        sxt = Gt*p['deltag']
        for i in range(len(sxt)):
            if sxt[i]<1e-5: sxt[i]=1e-5
        sxt = 1./sxt**2
    
    # Set model
    if p['iflat']==0:
        model = normalization*ones(len(omega))/sum(f0)
    elif p['iflat']==1:
        model = exp(-omega**2/p['gwidth'])
        model *= normalization/dot(model,f0)
    else:
        dat=loadtxt('model.dat').transpose()
        fm=interpolate.interp1d(dat[0],dat[1])
        model = fm(omega)
        model *= normalization/dot(model,f0)
        #savetxt('brisi_test', vstack((tau, fsg*dot(model,Ker))).transpose())

        
    print 'Model normalization=', dot(model,f0)

    # Set starting Aw(omega)
    Aw = random.rand(len(omega))
    Aw = Aw * (normalization/dot(Aw,f0))
    print 'Aw normalization=', dot(Aw,f0)

    dlda = me.initdlda(omega,Ker,sxt)
    
    temp=10.
    rfac=1.
    alpha=p['alpha0']
    
    for itt in range(p['Nitt']):
        print itt, 'Restarting maxent with rfac=', rfac, 'alpha=', alpha
        iseed = random.randint(0,maxint)
    
        me.maxent(Aw,rfac,alpha,temp,Ker,sxt,Gt,model,f0,p['Asteps'],iseed)
        S = me.entropy(Aw,model,f0)
        Trc = me.lambdac(alpha,Aw,omega,dlda)
        
        ratio = -2*S*alpha/Trc
        print 'Finished maxent with alpha=',alpha,'-2*alpha*S=',-2*alpha*S,'Trace=',Trc
        print '   ratio=', ratio

        savetxt('dos_'+str(itt), vstack((omega,Aw)).transpose())
        temp=0.001
        rfac=0.05
    
        if abs(ratio-1)<p['min_ratio']: break
    
        if (abs(ratio)<0.05):
            alpha *= 0.5
        else:
            alpha *= (1.+0.001*(random.rand()-0.5))/ratio
        
    for itt in range(p['Nr']):
        print 'Smoothing itt ', itt
        Aw = Broad(p['bwdth'],omega,Aw)
        Aw *= (normalization/dot(Aw,f0)) # Normalizing Aw
        
        savetxt('dos_'+str(p['Nitt']), vstack((omega,Aw)).transpose())
        
        temp=0.005
        rfac=0.005
        iseed = random.randint(0,maxint)
        me.maxent(Aw,rfac,alpha,temp,Ker,sxt,Gt,model,f0,p['Asteps'],iseed)
        
        S = me.entropy(Aw,model,f0)
        Trc = me.lambdac(alpha,Aw,omega,dlda)
        ratio = -2*S*alpha/Trc
        print 'Finished smoothing run with alpha=',alpha,'-2*alpha*S=',-2*alpha*S,'Trace=',Trc
        print '   ratio=', ratio
        
    savetxt('gtn', vstack((tau, fsg*dot(Aw,Ker))).transpose())
    Aw = Broad(p['bwdth'],omega,Aw)
    savetxt('dos.out', vstack((omega,Aw)).transpose())
    return (Aw, omega)

def Pade(om, Gm, x, gamma, Norder):
    zn = om[:Norder]*1j
    gn = Gm[:Norder]
    Pt = me.padecof(gn, zn)
    Gx = array([me.padeg(w+gamma*1j, zn, Pt) for w in x])
    return Gx

def InverseFourier(Gm, om, tau, beta, Nf=40):
    """Inverse Fourier transform which
       computes G(tau) from G(iom)
    """
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
    
    ah = FindHighFrequency(Gm,om,Nf)
    Gtau = zeros(len(tau),dtype=float)
    for it,t in enumerate(tau):
        Gtau[it] = me.fourpart(t,Gm,om,ah,beta)
    # Correction 1/omega^2 (sometimes usefull)
    df = Gm[-1].real*om[-1]/pi
    Gtau[0] += df
    Gtau[-1] -= df
    return Gtau

####################################
# Limited use -- only for rear cases
####################################
def MaximumEntropyTest(p, tau, Gt):
    
    def MEStep(alpha,rfac,Aw,temp,Ker,sxt,Gt,model,f0,Asteps,itt,reset=True):
        if (reset):
            temp=0.001
            rfac=0.05
        print 'Restarting maxent with rfac=', rfac, 'alpha=', alpha
        iseed = random.randint(0,maxint)
        me.maxent(Aw,rfac,alpha,temp,Ker,sxt,Gt,model,f0,Asteps,iseed)
        S = me.entropy(Aw,model,f0)
        Trc = me.lambdac(alpha,Aw,omega,dlda)
        ratio = -2*S*alpha/Trc
        print 'Finished maxent with alpha=',alpha,'-2*alpha*S=',-2*alpha*S,'Trace=',Trc
        print '   ratio=', ratio
        savetxt('dos_'+str(itt), vstack((omega,Aw)).transpose())
        return ratio
    
    beta = tau[-1]

    random.seed( 1 ) # seed for random numbers
    
    omega = linspace(-p['L'],p['L'],2*p['Nw']+1)
    f0,f1,f2 = me.initf0(omega)
    fsg=1
    if p['statistics']=='fermi':
        Gt = -Gt
        fsg=-1
        normalization = Gt[0]+Gt[-1]
        Ker = me.initker_fermion(omega,beta,tau)        
    elif p['statistics']=='bose':
        normalization = integrate.trapz(Gt,x=tau)
        Ker = me.initker_boson(omega,beta,tau)
    
    print 'beta=', beta
    print 'normalization=', normalization

    # Set error
    if p['idg']:
        sxt = ones(len(tau))/(p['deltag']**2)
    else:
        sxt = Gt*p['deltag']
        for i in range(len(sxt)):
            if sxt[i]<1e-5: sxt[i]=1e-5
        sxt = 1./sxt**2
    
    # Set model
    if p['iflat']==0:
        model = normalization*ones(len(omega))/sum(f0)
    elif p['iflat']==1:
        model = exp(-omega**2/p['gwidth'])
        model *= normalization/dot(model,f0)
    else:
        dat=loadtxt('model.dat').transpose()
        fm=interpolate.interp1d(dat[0],dat[1])
        model = fm(omega)
        model *= normalization/dot(model,f0)
        #savetxt('brisi_test', vstack((tau, fsg*dot(model,Ker))).transpose())

        
    print 'Model normalization=', dot(model,f0)

    # Set starting Aw(omega)
    Aw = random.rand(len(omega))
    Aw = Aw * (normalization/dot(Aw,f0))
    print 'Aw normalization=', dot(Aw,f0)

    dlda = me.initdlda(omega,Ker,sxt)
    
    temp=10.
    rfac=1.
    alpha=p['alpha0']

    for itt in range(10):
        ratio = MEStep(alpha,rfac,Aw,temp,Ker,sxt,Gt,model,f0,p['Asteps'],itt,itt!=0)
        if abs(ratio-1)<p['min_ratio']: break
        if (ratio<0.05):
            if ratio>0:
                alpha *= 1.1
            else:
                alpha /= 10.
        else:
            alpha *= (1.+0.001*(random.rand()-0.5))/ratio
            
    if abs(ratio-1)>2*p['min_ratio']:
        alpha=1.
        p['Asteps'] *= 1.5
        for itt in range(3,6):
            ratio = MEStep(alpha,rfac,Aw,temp,Ker,sxt,Gt,model,f0,p['Asteps'],itt)
    
    for itt in range(p['Nr']):
        print 'Smoothing itt ', itt
        Aw = Broad(p['bwdth'],omega,Aw)
        Aw *= (normalization/dot(Aw,f0)) # Normalizing Aw
        ratio = MEStep(alpha,rfac,Aw,temp,Ker,sxt,Gt,model,f0,p['Asteps'],p['Nitt']+itt)
        
        
    savetxt('gtn', vstack((tau, fsg*dot(Aw,Ker))).transpose())
    Aw = Broad(p['bwdth'],omega,Aw)
    savetxt('dos.out', vstack((omega,Aw)).transpose())
    return (Aw, omega)

def PadeTest(om, Gm, x, gamma, Norder):
    from scipy import poly1d
    
    zn = om[:Norder]*1j
    gn = Gm[:Norder]
    an = me.padecof(gn, zn)

    print 'zn=', zn
    print 'an=', an
    
    Aq_0 = poly1d([0])
    Aq_1 = poly1d([an[0]])
    Bq_0 = poly1d([1.])
    Bq_1 = poly1d([1.])
    for i in range(Norder-1):
        Aq_2 = Aq_1 +poly1d([1,-zn[i]])*an[i+1]*Aq_0
        Aq_0 = Aq_1
        Aq_1 = Aq_2
        Bq_2 = Bq_1 +poly1d([1,-zn[i]])*an[i+1]*Bq_0
        Bq_0 = Bq_1
        Bq_1 = Bq_2
        
    poles = sorted( roots(Bq_2), key=real )
    ezeros = sorted( roots(Aq_2), key=real )
    Bqp = poly1d(poles, r=True)
    Cbq = Bq_2(0.0)/Bqp(0.0)
    print 'ratio=', Cbq
    
    wgh=[]
    print 'poles='
    for i in range(len(poles)):
        Rq = poly1d(poles[:i] + poles[i+1:], r=True)
        wi = Aq_2(poles[i])/(Cbq*Rq(poles[i]))
        wgh.append(wi)
        if (poles[i].imag<0): used='Yes'
        else: used='No'
        print "%2d %12.4g %12.4g   %12.4g %12.4g  used=%s" % (i+1,poles[i].real, poles[i].imag, wi.real, wi.imag, used)
    wgh=array(wgh)
    #print 'zeros='
    #for i in range(len(ezeros)):
    #    print i+1, ezeros[i]
    #print 'weights=', wgh

    yt = zeros(len(om), dtype=complex)
    for i in range(len(poles)):
        if (poles[i].imag<0):
            yt += real(wgh[i])/(om*1j-poles[i])
    
    normy = sum(real(yt))
    normg = sum(real(Gm))
    if normy>1e-6:
        print 'Warning: norm mismatch: ratio=', normg/normy
        print 'Renormalizing'
        wgh *= (normg/normy)
    else:
        print 'Warning: Not enough poles. Bailing out'
    
    yt = zeros(len(om), dtype=complex)
    for i in range(len(poles)):
        if (poles[i].imag<0):
            yt += real(wgh[i])/(om*1j-poles[i])
    yr = zeros(len(x), dtype=complex)
    for i in range(len(poles)):
        if (poles[i].imag<0):
            yr += real(wgh[i])/(x-poles[i]+gamma*1j)
    
    savetxt('pade.corrected', vstack((x,real(yr),imag(yr))).transpose())
    G0 = Aq_2(zn)/Bq_2(zn)
    print 'G0=', G0

    subplot(2,1,1)
    plot(om, real(Gm), 'o')
    plot(om, real(yt), ':')
    #show()
    subplot(2,1,2)
    plot(x, imag(yr))
    #xlim([0,2.4])
    #ylim([0,0.6])
    show()
    
    
    Gx = array([me.padeg(w+gamma*1j, zn, an) for w in x])
    return Gx
    
def test1(filename,stat,deltag=0.006,L=4):
    params={'statistics': stat,    # fermi/bose
            'L'         : L,       # cutoff frequency on real axis
            'gwidth'    : 2*L,     # width of gaussian
            'Nw'        : 201,     # number of frequencies
            'idg'       : 1,       # error scheme: idg=1 -> sigma=deltag ; idg=0 -> sigma=deltag*G(tau)
            'deltag'    : deltag,  # error
            'Asteps'    : 4000,    # anealing steps
            'alpha0'    : 1000,    # starting alpha
            'min_ratio' : 0.001,    # condition to finish, what should be the ratio
            'iflat'     : 1,       # iflat=0 : constant model, iflat=1 : gaussian of width gwidth, iflat=2 : input using file model.dat
            'Nitt'      : 1000,    # maximum number of outside iterations
            'Nr'        : 0,       # number of smoothing runs
            'bwdth'     : 0.03,    # smoothing width
            }
    
    Gtdat = loadtxt(filename).transpose()
    tau = Gtdat[0]
    Gt  = Gtdat[1]
    (Aw, omega) = MaximumEntropy(params, tau, Gt)

def test3(filename='Giom.dat',stat='fermi',Ntau=200,L=20.):
    params={'statistics': stat,    # fermi/bose
            'L'         : L,       # cutoff frequency on real axis
            'gwidth'    : 2*L,     # width of gaussian
            'Nw'        : 200,     # number of frequencies
            'idg'       : 1,       # error scheme: idg=1 -> sigma=deltag ; idg=0 -> sigma=deltag*G(tau)
            'deltag'    : 0.005,   # error
            'Asteps'    : 4000,    # anealing steps
            'alpha0'    : 1000,    # starting alpha
            'min_ratio' : 0.01,    # condition to finish, what should be the ratio
            'iflat'     : 1,       # iflat=0 : constant model, iflat=1 : gaussian of width gwidth, iflat=2 : input using file model.dat
            'Nitt'      : 1000,    # maximum number of outside iterations
            'Nr'        : 0,       # number of smoothing runs
            'bwdth'     : 0.03,    # smoothing width
            }

    Gdat = loadtxt(filename).transpose()
    om = Gdat[0]
    Gm = Gdat[1]+Gdat[2]*1j
    beta = pi/om[0]
    tau = linspace(0,beta,Ntau+1)
    
    Gt = InverseFourier(Gm, om, tau, beta, Nf=40)
    (Aw, omega) = MaximumEntropy(params, tau, Gt)

    
def test5(filename, Norder, L=6.):
    data = loadtxt(filename).transpose()
    tau = data[0]
    Gt = data[1]
    beta = tau[-1]
    x = linspace(-L,L,801)
    gamma = 0.005
    
    fGt = interpolate.UnivariateSpline(tau, Gt, s=0)
    tau2 = linspace(0,beta,len(tau)*10)
    fGt2 = fGt(tau2)
    savetxt('Gt_interpolated', vstack((tau2,fGt2)).transpose())
    Gm=[]
    om=[]
    for im in range(Norder):
        omi = 2*im*pi/beta
        om.append(omi)
        Gm.append( integrate.trapz( fGt2*cos(tau2*omi), x=tau2 ) + 0j )
    Gm = array(Gm)
    om = array(om)

    savetxt('Gm', vstack((om,real(Gm))).transpose())

    Gx = Pade(om, Gm, x, gamma, Norder)
    savetxt('dos.pade', vstack((x,real(Gx),imag(Gx))).transpose())

    Ax = zeros(len(x),dtype=float)
    for ix,xx in enumerate(x):
        if (abs(xx)>1e-12):
            Ax[ix] = imag(Gx[ix])/(pi*xx)
        else:
            Ax[ix] = 0.5*imag(Gx[ix-1])/(pi*x[ix-1]) + 0.5*imag(Gx[ix+1])/(pi*x[ix+1]) 
    savetxt('Ax', vstack((x,Ax)).transpose())
    Ker = me.initker_boson(x,beta,tau)
    savetxt('gtn', vstack((tau, dot(Ax,Ker))).transpose())

def test5b(filename, Norder, L=6.):
    data = loadtxt(filename).transpose()
    tau = data[0]
    Gt = data[1]
    beta = tau[-1]
    x = linspace(-L,L,801)
    gamma = 0.005
    
    fGt = interpolate.UnivariateSpline(tau, Gt, s=0)
    tau2 = linspace(0,beta,len(tau)*10)
    fGt2 = fGt(tau2)
    savetxt('Gt_interpolated', vstack((tau2,fGt2)).transpose())
    Gm=[]
    om=[]
    for im in range(Norder*3):
        omi = 2*im*pi/beta
        om.append(omi)
        Gm.append( integrate.trapz( fGt2*cos(tau2*omi), x=tau2 ) + 0j )
    Gm = array(Gm)
    om = array(om)

    savetxt('Gm', vstack((om,real(Gm))).transpose())

    Gx = PadeTest(om, Gm, x, gamma, Norder)
    savetxt('dos.pade', vstack((x,real(Gx),imag(Gx))).transpose())

    Ax = zeros(len(x),dtype=float)
    for ix,xx in enumerate(x):
        if (abs(xx)>1e-12):
            Ax[ix] = imag(Gx[ix])/(pi*xx)
        else:
            Ax[ix] = 0.5*imag(Gx[ix-1])/(pi*x[ix-1]) + 0.5*imag(Gx[ix+1])/(pi*x[ix+1]) 
    savetxt('Ax', vstack((x,Ax)).transpose())
    Ker = me.initker_boson(x,beta,tau)
    savetxt('gtn', vstack((tau, dot(Ax,Ker))).transpose())

def test6(filename,Norder=100,L=10):
    
    Gdat = loadtxt(filename).transpose()
    om = Gdat[0]
    Gm = Gdat[1]+Gdat[2]*1j
    beta = pi/om[0]
    
    x = linspace(-L,L,501)
    gamma = 0.001
    
    Gx = Pade(om, Gm, x, gamma, Norder)
    savetxt('dos.pade', vstack((x,real(Gx),imag(Gx))).transpose())

def test6b(filename,Norder=100,L=10):
    
    Gdat = loadtxt(filename).transpose()
    om = Gdat[0]
    Gm = Gdat[1]+Gdat[2]*1j
    beta = pi/om[0]
    
    x = linspace(-L,L,501)
    gamma = 0.001
    
    Gx = PadeTest(om, Gm, x, gamma, Norder)
    savetxt('dos.pade', vstack((x,real(Gx),imag(Gx))).transpose())


def test2(fin, fout, Norder):
    # for the one band model on Bethe lattice
    def Gz(z):
        return 2*(z - sign(z.real+1e-16)*sqrt(z**2-1.))
    
    x = linspace(-6,6,1001)
    gamma = 0.001
    data = loadtxt(fin).transpose()
    om=data[0]
    Sm = data[2]*1j
    beta = pi/om[0]
    print 'beta=',beta

    Sx = Pade(om, Sm, x, gamma, Norder)
    Ax = -imag(Gz(x-Sx))/pi

    savetxt('sig.pade', vstack((x,real(Sx),imag(Sx))).transpose())
    savetxt(fout, vstack((x,Ax)).transpose())

    
if __name__ == '__main__':

    #test1('Augtau_nonint.dat','bose')
    
    #test1('Gtau.dat','fermi',0.006,4)  # test2

    #test6('Gf.out_T_0.005', 54, 20.)

    #fs=['U_1.0', 'U_2.0', 'U_2.4i', 'U_2.4m', 'U_3.0']
    #for fi in fs:
    #    test2('Sig.out_'+fi, 'DOS_'+fi, 72)
    
    test5b('Augtau_U_0.dat', 38)
    shutil.copy('pade.corrected', 'Pade_U_0.0')
    test5b('Augtau_U_1.0', 36) 
    shutil.copy('pade.corrected', 'Pade_U_1.0')
    test5b('Augtau_U_2.0', 36) 
    shutil.copy('pade.corrected', 'Pade_U_2.0')
    test5b('Augtau_U_2.4m', 36) 
    shutil.copy('pade.corrected', 'Pade_U_2.4m')
    test5b('Augtau_U_3.0', 24) 
    shutil.copy('pade.corrected', 'Pade_U_3.0')
    
    test1('FeAsAugtau_U_5.0','bose',0.01,15)
    
