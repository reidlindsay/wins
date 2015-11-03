#!  /usr/bin/env python

"""
Check topology for sufficient connectivity.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-12-12 00:22:16 -0600 (Mon, 12 Dec 2011) $
* $LastChangedRevision: 5370 $

:author: Ketan Mandke <kmandke@mail.utexas.edu>

:copyright:
    Copyright 2009-2011 The University of Texas at Austin

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

from optparse import OptionParser
import numpy as np

def read_topo(options, topofile):
    """Read topology layout from file."""
    f = file(topofile, 'r')
    s = f.readline()
    topo = {'border':None, 'layout':None}
    done = not (s)
    while not done:
        # convert s to dict (check for border and layout)
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
    f.close()
    return topo

def get_topology(options, topofile):
    """Get topology."""
    # load topology from file
    topo = read_topo(options, topofile)
    border = topo['border']
    layout = topo['layout']
    xmin, xmax, ymin, ymax = border[:4]
    return topo

def main():
    usage = "%prog [OPTIONS] TOPOFILE"
    parser = OptionParser(usage=usage)
    parser.add_option("-n", "--neighbors", dest="neighbors", type="int", \
            default=1, help="Minimum number of neighbors for each node [default=%default]")
    parser.add_option("-d", "--dist", dest="dist", type="float", \
            default=1.0, help="Cutoff distance to be considered a neighbor [default=%default]")
    parser.add_option("-c", "--connect", dest="nconnect", type="int", \
            default=5, help="Number of sources used to make connections [default=%default]")
    (options, args) = parser.parse_args()

    if len(args)<>1:
        print "Invalid number of arguments."
        parser.print_help()
        raise SystemExit

    # get parameters
    topofile = args[0]
    topo = get_topology(options, topofile)
    layout = topo['layout']
    border = topo['border']

    # determine number of neighbors for all nodes
    neighbors = 0*np.arange(len(layout))
    for n in range(len(layout)):
        neighbors[n] = 0
        p = np.array(layout[n])
        for m in range(len(layout)):
            x = np.array(layout[m])
            if (m==n): continue
            if (np.linalg.norm(p-x)<options.dist):
                neighbors[n] += 1
    pconnected  = []
    pdisconnect = []
    psources    = []
    for n in range(len(neighbors)):
        checkdist = (neighbors>n)
        checksrc  = (neighbors[:options.nconnect]>n)
        pconnected.append(checkdist.mean())
        psources.append(checksrc.mean())

    NN = len(neighbors)
    NB = options.neighbors
    print "num sources = %d"%(options.nconnect)
    print "num nodes   = %d"%(NN)
    print "              %s"%(["%4s"%k for k in range(1,NB+1)])
    print "Pconnected  = %s"%(["%.2f"%(v) for v in pconnected[:NB]])
    print "Psources    = %s"%(["%.2f"%(v) for v in psources[:NB]])
    #print "              %s"%(["%4s"%(k+1) for k in range(NB-1, 2*NB)])
    #print "Pdisconnect = %s"%(["%.2f"%(1-v) for v in pconnected[NB-1:2*NB]])
    print "avg neighbors = %.2f"%(neighbors.mean())

if __name__ == '__main__':
    main()
