#!/usr/bin/env python2
# @Copyright 2007 Kristjan Haule
import sys,re,os
import optparse
from scipy import *

import numpy
nv = map(int,numpy.__version__.split('.'))
if (nv[0],nv[1]) < (1,6):
    loadtxt = io.read_array
    def savetxt(filename, data):
        io.write_array(filename, data, precision=16)  

if __name__=='__main__':
    """ Takes several self-energy files and produces an average over these self-energy files
    """
    usage = """usage: %prog [ options ] argumens

    The script takes several self-energy files and produces an average self-energy

    arguments  -- all input self-energies
    option -o  -- output self-energy file
    """

    parser = optparse.OptionParser(usage)
    parser.add_option("-o", "--osig", dest="osig", default='sig.inpx', help="filename of the output self-energy file. Default: 'sig.inp'")
    parser.add_option("-l", "--lext", dest="m_extn", default='', help="For magnetic calculation, it can be 'dn'.")
    parser.add_option("-n", "--nexecute", dest="nexecute", action="store_true", default=False,  help="Do not execute the comments")
    parser.add_option("-s", dest="stdev", action='store_true', default=False, help="Computes also the standard deviation - the error of self-energy")

    # Next, parse the arguments
    (options, args) = parser.parse_args()

    print 'files to average over:', args
    print 'output: ', options.osig
    ws_oo=[]
    wVdc=[]
    wdata=[]
    for f in args:
        fi=open(f,'r')
        header=[]
        dat = fi.readlines()
        s_oo = None
        Vdc = None
        data=[]
        for line in dat:
            m = re.search('#(.*)',line)
            if m is not None: 
                #print m.group(1).strip()[:5]
                if m.group(1).strip()[:5]=="s_oo=" or m.group(1).strip()[:5]=="Vdc= ":
                   if not options.nexecute: exec(m.group(1).strip())
                else:
                   header.append(line)
            else:
                data.append( map(float, line.split() ) )

        if s_oo is not None: ws_oo.append(s_oo)
        if Vdc is not None: wVdc.append(Vdc)
        wdata.append(data)    
        fi.close()
    
    fout = open(options.osig, 'w')
    
    #if len(header):
    #    for i in range(len(header)):
    #       print >> fout, header[i],

    if len(ws_oo):
        ws_oo = array(ws_oo)
        as_oo=[]
        for i in range(shape(ws_oo)[1]):  as_oo.append( sum(ws_oo[:,i])/len(ws_oo) )
        print 's_oo=', as_oo
        print >> fout, '# s_oo=', as_oo

    if len(wVdc):
        wVdc = array(wVdc)
        aEdc=[]
        for i in range(shape(wVdc)[1]):  aEdc.append( sum(wVdc[:,i])/len(wVdc) )
        print 'Vdc=', aEdc
        print >> fout, '# Vdc=', aEdc

    wdata = array(wdata)
    wres = zeros(shape(wdata)[1:], dtype=float)
    for i in range(len(wdata)): wres[:,:] += wdata[i,:,:]
    wres *= 1./len(wdata)
    
    if options.stdev:
        sw = shape(wdata)
        wstd = zeros((sw[1],sw[2]-1), dtype=float) # no frequency
        for i in range(len(wdata)): wstd[:,:] += wdata[i,:,1:]**2
        wstd *= 1./len(wdata)
        wstd[:,:] = sqrt(wstd[:,:] - wres[:,1:]**2)
        
        wres = hstack( (wres, wstd) )
        
    savetxt(fout,wres)
    
    
    
