#!  /usr/bin/env python

"""
Calculate distribution of link SNRs for a given topology.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-10-19 17:04:02 -0500 (Wed, 19 Oct 2011) $
* $LastChangedRevision: 5220 $

:author: Ketan Mandke <kmandke@mail.utexas.edu>

:copyright:
    Copyright 2009 The University of Texas at Austin

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
"""
__docformat__ = "restructuredtext en"

from SimPy.Simulation import *
from scapy.all import *
from wins import *
from wins.ieee80211 import *
from copy import copy, deepcopy

import sys
from optparse import OptionParser

import numpy as np
import struct
import gc
import time

from pylab import *

def read_topo(options, topofile):
    f = file(topofile, 'r')
    s = f.readline()
    topo = {'border':None, 'layout':None}
    done = not (s)
    while not done:
        # covert s to dict (check for border and layout)
        try:
            d = eval(s)
            assert isinstance(d, dict)
            assert ('border' in d) and ('layout' in d)
        except:
            d = None
        # add dict to datain
        if d: topo = d
        # get next input
        s = f.readline()
        done = not(s)
    return topo

def topo_to_snr(options, topofile):
    # set up channel
    modeltype  = options.tgnmodel       # default -> LOS Channel
    alpha      = options.alpha
    # set topology
    topo = read_topo(options, topofile)
    border = topo['border']
    layout = topo['layout']
    xmin, xmax, ymin, ymax = border[0:4]
    numnodes = len(layout)
    # set up pathloss parameters for calculating SNR
    Ptx = DOT11N_MAXPOWER
    Lrx, Ltx = Dot11NRadio.Lrx, Dot11NRadio.Ltx
    Grx, Gtx = Dot11NRadio.Grx, Dot11NRadio.Gtx
    No = Dot11NRadio.thermalnoise(DOT11N_BANDWIDTH) + DOT11N_NOISEFIGURE
    Psignal = Ptx + Gtx + Grx - Ltx - Lrx
    # use breakpoint distance as refdist
    refdist = Dot11NChannel(modeltype=modeltype, n=alpha).bpdist
    refloss = Propagation.freespace(refdist, n=2.0, fc=Dot11NChannel.fc)
    # use topology to calculate corresponding link SNRs
    linksnr, allsnr = [], []
    nlinks = int(numnodes/2)  # first half are sender, second half are receivers
    for n in range(numnodes):
        for m in range(numnodes):
            if n==m: continue
            pa, pb = layout[n], layout[m]
            ds = Propagation.distance(pa,pb)
            PL = refloss + alpha*linear2db(ds/refdist)
            snr = Psignal - PL - No
            if (n<nlinks): allsnr.append(snr)   # all links
            if (n<nlinks) and (m<nlinks): linksnr.append(snr)   # direct links
    # output link snrs and all snrs to stdout
    linkdata = {'data':linksnr, 'parameters':{'label':'Link SNRs'}}
    alldata = {'data':allsnr,  'parameters':{'label':'All SNRs'}}
    if options.verbose:
        f = sys.stdout
        f.write("%s\n"%linkdata)
        f.write("%s\n"%alldata)

    # plot histogram
    plines, labels = [], []
    bins, fmt = 40, ['c ', 'g/']
    normed = options.normed
    # plot all
    if options.plotall:
        x = np.array(allsnr)
        f = fmt.pop(0)
        c, h = f[0], f[1]
        fmt.append(f)
        n, bins, p = hist(x, range=(0,x.max()), bins=40, facecolor=c, hatch=h, normed=normed)
        plines.append(p[0]), labels.append('All Links')
    # plot direct
    if options.plotlink:
        x = np.array(linksnr)
        if plines: hold(True)
        f = fmt.pop(0)
        c, h = f[0], f[1]
        fmt.append(f)
        n, bins, p = hist(x, range=(0,x.max()), bins=40, facecolor=c, hatch=h, normed=normed)
        plines.append(p[0]), labels.append('Direct Links')
    # show plot
    if (options.plotlink or options.plotall):
        grid(True)
        if (len(plines)>1): legend(plines, labels)
        xlabel('SNR (dB)')
        ylabel('Number of Links')
        # save to figure and show plot
        if options.savefig:
            savefig(options.savefig,transparent=True)
        show()
    
def main():
    usage = "%prog [OPTIONS] TOPOLOGY"
    parser = OptionParser(usage=usage)
    # channel parameters
    parser.add_option("", "-v", dest="verbose", action="store_true",  \
            default=False, help="Print SNR data to standard output.")
    parser.add_option("", "--tgn-model", dest="tgnmodel", \
            default=None, help="Specify TGn model.")
    parser.add_option("", "--alpha", dest="alpha", type="float", \
            default=2.0, help="Specify pathloss exponent [default=%default]")
    parser.add_option("", "--plot-link", dest="plotlink", action="store_true",  \
            default=False, help="Plot SNRs of directs links.")
    parser.add_option("", "--plot-all", dest="plotall", action="store_true",  \
            default=False, help="Plot SNR of all links.")
    parser.add_option("", "--normed", dest="normed", action="store_true",  \
            default=False, help="Normalize plots.")
    parser.add_option("", "--save-fig", dest="savefig", \
                      help="Save figure to file.")
    parser.add_option("-s", "--fscale", dest="fscale", type="float", \
            default=1.0, help="Scale figure text/markers/lines when saving figure [default=%default]")
    (options, args) = parser.parse_args()

    if len(args)<1:
        print "Invalid number of arguments."
        parser.print_help()
        raise SystemExit

    # get topology file
    tfile = args[0]

    # set parameters for saving figure to file
    if options.savefig:
        if (options.savefig.lower().find('svg')>-1):
            rcParams.update({'backend': 'svg'})
        elif (options.savefig.lower().find('pdf')>-1):
            rcParams.update({'backend': 'pdf'})
        else:
            print "Cannot create non-pdf or non-svg file %s"%(options.savefig)
            raise SystemExit
        fscale = options.fscale
        rcParams.update({'font.size': 12.0*fscale, \
                         'lines.linewidth': 1.0*fscale, \
                         'lines.markeredgewidth': 0.5*fscale, \
                         'lines.markersize': 6*fscale, \
                         'axes.linewidth': 1.0*fscale, \
                         'patch.linewidth': 1.0*fscale, \
                         'grid.linewidth': 0.5*fscale})

    # plot snr data
    topo_to_snr(options, tfile)

if __name__ == '__main__':
    main()
