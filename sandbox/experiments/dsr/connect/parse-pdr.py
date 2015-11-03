#!  /usr/bin/env python

"""
Parse PDR (vs. number of nodes) for DSR simulation.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-12-18 16:57:02 -0600 (Sun, 18 Dec 2011) $
* $LastChangedRevision: 5381 $

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

from wins import *
from wins.ieee80211 import *

from optparse import OptionParser
import sys

import numpy as np

def read_trace(options, tracefile):
    # load trace from file
    tr = Trace()
    tr.read(tracefile)
    # return trace
    return tr

def hasfield(event, *args):
    """Check if event has all the specified arguments."""
    ok = all([(a in event) for a in args])
    return ok

def parse_info(options, trace):
    # (data, param) variables
    data = []
    param = {'plen': None}
    rreqrate = None
    envspeed = None
    for e in trace.events:
        # check timestamp - [tmin, tmax]
        ts = float(e['ts'])
        if (ts<options.tmin): continue
        elif (ts>options.tmax): break
        # get event/object/packet parameters
        obj, evt = e['obj'], e['event']
        # check for MC PDR
        if (obj=="MC") and (evt=="PDR"):
            assert hasfield(e, 'pdr', 'numnodes')
            pdr = float(e['pdr'])
            nn  = int(e['numnodes']) - 2
            dp = {'x':nn, 'y':pdr, 'ndata': 1}
            data.append(dp)
        # check for other parameters
        if (obj=="MON") and (evt=="RATE"):
            assert hasfield(e, 'rreqrate')
            rate = int(e['rreqrate'])
    # set remaining parameters
    param['use-yerror'] = True
    param['rreqrate'] = True
    return data, param

def parse_pdr():
    usage = "%prog [OPTIONS] TRACEFILE1 [TRACEFILE2 ...]\n" + \
            "      Writes parsed data to standard output."
    parser = OptionParser(usage=usage)
    parser.add_option("", "--tmin", dest="tmin", \
            type="float", default=-np.inf, \
            help="Specify simulation time to start from [default=%default]")
    parser.add_option("", "--tmax", dest="tmax", \
            type="float", default=np.inf, \
            help="Specify simulation time to stop at [default=%default]")
    parser.add_option("-f", "--fmt", dest="fmt", \
            default="bo", help="Format of data parsed from tracefiles [default=%default]")
    parser.add_option("", "--label", dest="label", \
            default="Hydra PHY", help="Specify label for newly parsed data [default=%default]")
    (options, args) = parser.parse_args()

    if len(args)<1:
        print "Insufficient number of arguments."
        parser.print_help()
        raise SystemExit

    # get trace files
    tracefile = args[0:]
    numtraces = len(tracefile)

    # set default parameters
    defaultparam = {'xlabel': "Number of Nodes", \
                    'ylabel': "% connected", \
                    'title': "Connected Topologies vs. Number of Nodes", \
                    'label': None, \
                    'source': None, \
                    'format': None}

    # define formats
    formats = [(options.fmt, options.fmt[0]+"x:")]     # use same format for all
    #formats = [('ro','rx:'), ('bo', 'bx:'), ('go', 'gx:')]

    # parse data from each data set
    for k in range(numtraces):
        # get tracefile
        tfile = tracefile[k]
        trace = read_trace(options, tfile) # treat as normal wins trace file
        if not trace: continue
        # get other parameters
        fmt = formats[k%len(formats)]
        # parse trace to get (data, param)
        sys.stderr.write("Parsing trace from %s ... \n"%(tfile))
        data, dparam = parse_info(options, trace)
        if data:
            param = defaultparam.copy()
            param.update(dparam)
            param['format'] = fmt[0]
            param['source'] = tfile
            param['label']  = options.label
            pdata = {'data': data, 'parameters':param}
            sys.stdout.write("%s\n"%(pdata))

if __name__ == '__main__':
    parse_pdr()
