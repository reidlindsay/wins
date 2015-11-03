#!  /usr/bin/env python

"""
Parse distribution of hop-count in DSR simulation.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-12-10 17:09:21 -0600 (Sat, 10 Dec 2011) $
* $LastChangedRevision: 5368 $

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

from config_parse import FILTER_DISCONNECTED

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
    rate = None
    envspeed = None
    # event parameters
    macargs     = ('mac-txts', )
    agtargs     = ('net-root', )
    # paramaters for logging data
    sendbuffer = {}     # log sent packets using AGT name for index
    sendnodes  = {}     # nodes that sent data : number delivered
    hopcount = {0: 0}   # log number of hopcounts in simulation
    # default arguments for buffer entry
    defaultpkt = {'txts':None, 'rxts':None, 'hopcount':None, 'sender':None}
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
        # check for DSR SND
        hasagt = hasfield(e, 'agt-root')
        if (obj=="DSR") and (evt=="SND") and hasagt:
            agt = e['agt-root']
            assert (agt in sendbuffer)
            nhops = None
            if   (rname=="DSR|SRCROUTE"):
                assert hasfield(e, 'srcroute.addresses')
                srcroute = eval(e['srcroute.addresses'])
                nhops = len(srcroute) + 1
            elif (rname=="DSR"):
                nhops = 1
            assert (nhops is not None)
            sendbuffer[agt]['hopcount'] = nhops
        # check for other parameters
        if (obj=="AGT") and (evt=="SND"):
            assert hasfield(e, 'plen')
            plen = int(e['plen'])
            if param['plen'] is None: param['plen'] = plen
            else:                     assert (plen==param['plen'])
        if (obj=="MON") and (evt=="RATE"):
            assert hasfield(e, 'rate')
            rate = float(e['rate'])
        if (obj=="MON") and (evt=="MODEL"):
            assert hasfield(e, 'environmentspeed')
            envspeed = float(e['environmentspeed'])
    # skip data from disconnected topologies
    if FILTER_DISCONNECTED:
        isdisconnected = False
        for s, v in sendnodes.items():
            if (v==0): isdisconnected = True
        if isdisconnected:
            sys.stderr.write("WARNING: Skipping data for disconnected topology!\n")
            return data, param
    # calculate hop histogram from sendbuffer
    assert (rate is not None)
    nsent = 0
    for p,d in sendbuffer.items():
        plen = param['plen']
        assert (plen is not None)
        nhops = d['hopcount']
        if nhops is None: continue
        if nhops not in hopcount: hopcount[nhops] = 0
        # update hopcount histogram
        hopcount[nhops] += 1
        nsent += 1
    # get hop-count data
    for h, n in hopcount.items():
        y = 1.0*n/nsent
        dp = {'x': h, 'y': y, 'ndata': 1}
        data.append(dp)
    # set remaining parameters
    param['use-yerror'] = True
    param['environmentspeed'] = envspeed
    param['rate'] = rate
    return data, param

def parse_delay():
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
    defaultparam = {'xlabel': "Hop Count", \
                    'ylabel': "PDF", \
                    'title': "Hop Count Distribution", \
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
    parse_delay()
