#!/usr/bin/env python2
# @Copyright 2007 Kristjan Haule
import sys,re,os,shutil
import optparse, subprocess
from scipy import *
from scipy import optimize
import numpy

nv = map(int,numpy.__version__.split('.'))
if (nv[0],nv[1]) < (1,6):
    loadtxt = io.read_array
    #
    def savetxt(filename, data):
        io.write_array(filename, data, precision=16)

def union(data):
    " Takes a union of array or list"
    c = []
    for d in data:
        if d not in c:
            c.append(d)
    return c

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
    
def default_realaxis_mesh():
    # Probably should allow user to control the mesh generation
    # to some extent view program options...
    
    # use gaumesh.py to generate standard mesh
    gaumesh = os.path.join( utils.DmftEnvironment().ROOT, 'gaumesh.py' )
    args = [
        'x0=[0,0]',
        'dx0=[0.25,0.01]',
        'fwhm=[30, 3]',
        'xmin=-20',
        'xmax=20'
        ]
    stdoutdata, stderrdata = subprocess.Popen([gaumesh] + args, stdout=subprocess.PIPE).communicate()
    return [float(x) for x in stdoutdata.split('\n') if x.strip() and not x.strip().startswith('#')]


if __name__=='__main__':
    """ Create a zero input self-energy for dmft0, dmft1 and dmft0
    """
    usage = """usage: %prog [ options ]

    Create a zero input self-energy for dmft0, dmft1 and dmft0
    """
    # n==(200/T-1)/2.
    if os.path.exists('INPUT.py'): 
       execfile('INPUT.py')
       nc=0
       for i,ats in enumerate(p['cor_at']):
          for j,orbs in enumerate(p['cor_orb'][i]):
             if len(orbs)>0: nc+=1
       nc*=p['nspin']
       Vdc = p['U'][0]*(p['nf']-0.5)-0.5*p['J'][0]*(p['nf']-1.)
       T = 1.0/pC['beta'][0]
       parser = optparse.OptionParser(usage)
       parser.add_option("-c", "--Nc",  dest="Nc",    type="int", default=nc, help="Number of correlated orbitals")
       parser.add_option("-e", "--Edc",  dest="Edc",    type="float", default=Vdc, help="Starting double counting and Hartree value of the self-energy. Should be close to U(n-1/2)")
       parser.add_option("-i", "--sinp", dest="insig",  default='sig.inp', help="the mesh will be used from this file.", metavar="FILE")
       parser.add_option("-o", "--sout", dest="outsig", default='sig.inp', help="the result will be written to this file")
       parser.add_option("-T", "--Temp",  dest="T",     type="float", default=T, help="Temperature")
       parser.add_option("-n", "--nom",  dest="nom",    type="int", default=None, help="Number of frequency points")
       parser.add_option("-L", "--range", dest="L",     type="float", default=20., help="energy range on real axis")
       parser.add_option("-x", "--x0",    dest="x",     type="float", default=0.05, help="energy range on real axis")
       parser.add_option("-N", "--Nom",  dest="Nom",    type="int", default=400, help="Number of frequency points on real axis")

    else: 
       parser = optparse.OptionParser(usage)
       parser.add_option("-c", "--Nc",  dest="Nc",    type="int", default=1, help="Number of correlated orbitals")
       parser.add_option("-e", "--Edc",  dest="Edc",    type="float", default=0.0, help="Starting double counting and Hartree value of the self-energy. Should be close to U(n-1/2)")
       parser.add_option("-i", "--sinp", dest="insig",  default='sig.inp', help="the mesh will be used from this file.", metavar="FILE")
       parser.add_option("-o", "--sout", dest="outsig", default='sig.inp', help="the result will be written to this file")
       parser.add_option("-T", "--Temp",  dest="T",     type="float", default=0.0, help="Temperature")
       parser.add_option("-n", "--nom",  dest="nom",    type="int", default=None, help="Number of frequency points")
       parser.add_option("-L", "--range", dest="L",     type="float", default=20., help="energy range on real axis")
       parser.add_option("-x", "--x0",    dest="x",     type="float", default=0.05, help="energy range on real axis")
       parser.add_option("-N", "--Nom",  dest="Nom",    type="int", default=400, help="Number of frequency points on real axis")

    # Next, parse the arguments
    (options, args) = parser.parse_args()
    
    if options.nom==None and options.T>0:
        options.nom = int((300./options.T-1)/(2.*pi))+1
        
#    env = utils.W2kEnvironment()
#    case = env.case
#
#    print 'Edc=%s case=%s, insig=%s, outsig=%s, nom=%s T=%s' %  (options.Edc, case, options.insig, options.outsig, options.nom, options.T)
#
#    inl = indmffile.Indmfl(case)
#    inl.read()
#    m_extn = 'dn' if os.path.exists(case+'.indmfl'+'dn') else ''
#    if m_extn:
#        print 'INFO: case.indmfldn present => magnetic calculation with two dmft2 steps'
#        inldn = indmffile.Indmfl(case, 'indmfl'+m_extn)
#        inldn.read()
#
#    if options.T>0 and inl.matsubara:
#        print '..... creating new matsubara mesh of size '+str(options.nom)+' of T='+str(options.T)
#        omega = (arange(1,options.nom,1)*2-1)*pi*options.T
#    elif os.path.isfile(options.insig) and os.path.getsize(options.insig)>0:
#        # Read the input file
#        sigdata = loadtxt(options.insig)  # self-energy from 'sig.inp' on long mesh
#        if len(shape(sigdata))==1:
#            omega = sigdata
#        else:
#            omega = (sigdata.transpose())[0]
#    else:
#        if not inl.matsubara:
#            print '..... Could not find '+options.insig+'. Generating default real axis mesh.'
#            omega = GiveTanMesh(options.x, options.L,options.Nom/2)
#        else:
#            Found=False
#            if os.path.isfile('params.dat'):
#                execfile('params.dat')
#                if iparams0.has_key('beta'): 
#                    Found=True
#                    beta = iparams0['beta'][0]
#                    if options.nom==None: options.nom = (300*beta-1)/(2.*pi)
#                    print '..... creating new matsubara mesh of size '+str(options.nom)+' of T='+str(1./beta)
#                    omega = (arange(1,options.nom,1)*2-1)*pi/beta
#            if not Found:
#                print '..... Could not find '+options.insig+'. Do not know the temperature. Can not create self-energy!'
#                print '..... Boiling out.....'
#                sys.exit(1)
#        
#    print '..... Going over all correlated blocks'
#    cols=[]
#    for icix in inl.siginds.keys():    # over all imp.problems, even those that are equivalent
#        Sigind = inl.siginds[icix]
#        cols = sort(union(array(Sigind).flatten().tolist()+cols)).tolist()
#    if m_extn:
#        for icix in inldn.siginds.keys():    # over all imp.problems, even those that are equivalent
#            Sigind = inldn.siginds[icix]
#            cols = sort(union(array(Sigind).flatten().tolist()+cols)).tolist()
#    if 0 in cols: cols.remove(0)
#
#    print 'cols=', cols
#    Nc = cols[-1]
#
#    # saving the original self-energy if necessary
#    if options.insig==options.outsig and os.path.isfile(options.insig) and os.path.getsize(options.insig)>0:
#        shutil.copy2(options.insig, options.insig+'.bak')
#
#
#    print 'om=', omega
#
#    if options.Edc==0 and os.path.isfile('params.dat'):
#        execfile('params.dat')
#        if iparams0.has_key('U') and iparams0.has_key('J') and iparams0.has_key('nf0'):
#            U = iparams0['U'][0]
#            J = iparams0['J'][0]
#            nf = iparams0['nf0'][0]
#            options.Edc = U*(nf-0.5)-0.5*J*(nf-1.)
            
        
    omega = (arange(1,options.nom,1)*2-1)*pi*options.T
    # writing to the file
    fo = open(options.outsig, 'w')
    print >> fo, '# nom,ncor_orb=', options.nom-1, options.Nc
    print >> fo, '# T=', (options.T)
    print >> fo, '# s_oo-Vdc=', ("0.0 "*(options.Nc))
    print >> fo, '# s_oo=', (ones(options.Nc)*options.Edc).tolist()
    print >> fo, '# Vdc=', (ones(options.Nc)*options.Edc).tolist()
    for iom,om in enumerate(omega):
        print >> fo, ("%20.15f "%om), ("0.0 "*(2*options.Nc))
        
    print(options.outsig,'written to the disc.')
