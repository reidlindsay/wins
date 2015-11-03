#!  /usr/bin/env python

"""
Parse PDR of DSR simulation.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-12-08 12:26:08 -0600 (Thu, 08 Dec 2011) $
* $LastChangedRevision: 5357 $

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

FILTER_DISCONNECTED = 1

def parse_info(options, trace):
    # (data, param) variables
    data = []
    param = {'plen': None}
    envspeed = None
    # event parameters
    macargs     = ('mac-txts', )
    agtargs     = ('net-root', )
    # paramaters for logging data
    sendbuffer = {}     # log sent packets using AGT name for index
    sendnodes  = {}     # nodes that sent data : number delivered
    # default arguments for buffer entry
    defaultpkt = {'txts':None, 'rxts':None, 'sender':None}
    for e in trace.events:
        # check timestamp - [tmin, tmax]
        ts = float(e['ts'])
        if (ts<options.tmin): continue
        elif (ts>options.tmax): break
        # get event/object/packet parameters
        obj, evt = e['obj'], e['event']
        pname, pid = e['packet'], e['pid']
        rname, rid = e['root'], e['rid']
        nname, nid = e['node'], e['nid']
        pkt = "%s(%s)"%(pname, pid)
        root = "%s(%s)"%(rname, rid)
        node = "%s(%s)"%(nname, nid)
        # check for AGT SND
        if (obj=="AGT") and (evt=="SND"):
            assert (pkt not in sendbuffer)
            sendbuffer[pkt] = defaultpkt.copy()
            sendbuffer[pkt]['txts'] = float(ts)
            sendbuffer[pkt]['sender'] = node
            if node not in sendnodes: sendnodes[node] = 0
        # check for AGT RCV
        if (obj=="AGT") and (evt=="RCV"):
            assert (pkt in sendbuffer)
            sendbuffer[pkt]['rxts'] = float(ts)
            sender = sendbuffer[pkt]['sender']
            assert (sender in sendnodes)
            sendnodes[sender] += 1
        # check for other parameters
        if (obj=="AGT") and (evt=="SND"):
            assert hasfield(e, 'plen')
            plen = int(e['plen'])
            if param['plen'] is None: param['plen'] = plen
            else:                     assert (plen==param['plen'])
        if (obj=="MON") and (evt=="MODEL"):
            assert hasfield(e, 'environmentspeed')
            envspeed = float(e['environmentspeed'])
    # skip data from disconnected topologies
    if FILTER_DISCONNECTED:
        #sys.stderr.write("WARNING: FILTERING DISCONNECTED!\n")
        isdisconnected = False
        for s, v in sendnodes.items():
            if (v==0): isdisconnected = True
            #sys.stderr.write("  %10s delivered %d\n"%(s,v))
        if isdisconnected:
            sys.stderr.write("WARNING: Skipping data for disconnected topology!\n")
            return data, param
    # calculate PDR from sendbuffer
    assert (envspeed is not None)
    envspeed = max(envspeed, const.EPSILON)
    for p,d in sendbuffer.items():
        plen = param['plen']
        assert (plen is not None)
        if (d['rxts'] is None):
            # packet was not delivered
            x, y = envspeed, 0
        else:
            # packet was delivered
            x, latency = envspeed, d['rxts']-d['txts']
            L = 8*plen               # length in bits
            tput = (1e-6*L/latency)     # throughput in Mbps
            y = 1.0
        dp = {'x': x, 'y': y, 'ndata': 1}
        data.append(dp)
    # set remaining parameters
    param['use-yerror'] = True
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
    defaultparam = {'xlabel': "$v_0$ (km/hr)", \
                    'ylabel': "PDR", \
                    'title': "Packet-Delivery Ratio vs. Environment Speed", \
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
