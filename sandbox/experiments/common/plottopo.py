#!  /usr/bin/env python

"""
Plot data from topology file.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-11-11 19:39:03 -0600 (Fri, 11 Nov 2011) $
* $LastChangedRevision: 5332 $

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
from pylab import *

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
    parser.add_option("-c", "--num-connections", dest="nconnect", type="int", \
            default=1, help="Set number of active connections [default=%default]")
    parser.add_option("", "--xmin", dest="xmin", type="float", \
            default=None, help="Set minimum value on x-axis [default=%default]")
    parser.add_option("", "--xmax", dest="xmax", type="float", \
            default=None, help="Set maximum value on x-axis [default=%default]")
    parser.add_option("", "--ymin", dest="ymin", type="float", \
            default=None, help="Set minimum value on y-axis [default=%default]")
    parser.add_option("", "--ymax", dest="ymax", type="float", \
            default=None, help="Set maximum value on y-axis [default=%default]")
    parser.add_option("-s", "--fscale", dest="fscale", type="float", \
            default=1.0, help="Scale figure text/markers/lines when saving figure [default=%default]")
    parser.add_option("", "--save-fig", dest="savefig", \
                      help="Save figure to file.")
    (options, args) = parser.parse_args()

    if len(args)<>1:
        print "Invalid number of arguments."
        parser.print_help()
        raise SystemExit

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
                         'legend.fontsize': 'medium', \
                         'lines.linewidth': 1.0*fscale, \
                         'lines.markeredgewidth': 0.5*fscale, \
                         'lines.markersize': 6*fscale, \
                         'axes.linewidth': 1.0*fscale, \
                         'patch.linewidth': 1.0*fscale, \
                         'grid.linewidth': 0.5*fscale, \
                         'xtick.major.size' : 4*fscale, \
                         'xtick.minor.size' : 2*fscale, \
                         'xtick.major.pad'  : 4*fscale, \
                         'xtick.minor.pad'  : 4*fscale, \
                         'ytick.major.size' : 4*fscale, \
                         'ytick.minor.size' : 2*fscale, \
                         'ytick.major.pad'  : 4*fscale, \
                         'ytick.minor.pad'  : 4*fscale, \
                         })

    # get parameters
    nconnect = options.nconnect
    topofile = args[0]
    topo = get_topology(options, topofile)
    layout = topo['layout']
    border = topo['border']
    options.xmin = border[0]
    options.xmax = border[1]
    options.ymin = border[2]
    options.ymax = border[3]
    fig = gcf()

    # plot senders and receivers
    numnodes = len(layout)
    sx = [layout[k][0]    for k in range(numnodes) if (k<numnodes/2)]
    sy = [layout[k][1]    for k in range(numnodes) if (k<numnodes/2)]
    rx = [layout[-k-1][0] for k in range(numnodes) if (k<numnodes/2)]
    ry = [layout[-k-1][1] for k in range(numnodes) if (k<numnodes/2)]
    plot(sx, sy, 'bo')
    hold(True), grid(True)
    plot(rx, ry, 'r^', markerfacecolor="None")

    # plot number of connections
    for k in range(nconnect):
        cx = [sx[k], rx[-k-1]]
        cy = [sy[k], ry[-k-1]]
        plot(cx, cy, 'g--')
        mx = (cx[1] - cx[0])/2 + cx[0]
        my = (cy[1] - cy[0])/2 + cy[0]
        text(mx, my, "C%d"%(k))

    # adjust axis if xmin/xmax/ymin/ymax is set
    xmin, xmax = options.xmin, options.xmax
    ymin, ymax = options.ymin, options.ymax
    z = list(axis())
    if xmin is not None: z[0] = xmin
    if xmax is not None: z[1] = xmax
    if ymin is not None: z[2] = ymin
    if ymax is not None: z[3] = ymax
    axis(z)

    # adjust aspect ratio
    dx = border[1] - border[0]
    dy = border[3] - border[2]
    h = fig.get_figheight()
    w = fig.get_figwidth()
    fig.set_figwidth(w*sqrt(dx/dy))
    if options.savefig:
        savefig(options.savefig, transparent=True)

    show()

if __name__ == '__main__':
    main()
