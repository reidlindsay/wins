#!  /usr/bin/env python

"""
Parse throughput of RBAR from trace files.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-10-23 19:18:52 -0500 (Sun, 23 Oct 2011) $
* $LastChangedRevision: 5287 $

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
    rmsdelay = []
    # event parameters
    phyargs     = ('dot11n-sinr', 'dot11n-channel-fading')
    macargs     = ('mac-txts', )
    agtargs     = ('net-root', )
    # paramaters for logging data
    avgsnr = None       # average SNR reported by MC
    currpkt = None      # current IP packet being sent
    sendbuffer = {}     # log sent packets using IP name for index
    # default arguments for buffer entry
    defaultpkt = {'agt':None, 'mac': None}
    defaultpkt.update({'txts':None, 'rxts':None, 'snr':None})
    for e in trace.events:
        # check timestamp - [tmin, tmax]
        ts = float(e['ts'])
        if (ts<options.tmin): continue
        elif (ts>options.tmax): break
        # get event/object/packet parameters
        obj, evt = e['obj'], e['event']
        pname, pid = e['packet'], e['pid']
        rname, rid = e['root'], e['rid']
        pkt = "%s(%s)"%(pname, pid)
        root = "%s(%s)"%(rname, rid)
        # check for MC AVGSNR
        if (obj=="MC") and (evt=="AVGSNR"):
            avgsnr = float(e['snr'].lower().replace("db","") )
        # check for RTG SND
        if (obj=="RTG") and (evt=="SND"):
            assert (avgsnr is not None)
            assert (pkt not in sendbuffer)
            # new packet being sent from AGT
            sendbuffer[pkt] = defaultpkt.copy()
            sendbuffer[pkt]['agt'] = root
            sendbuffer[pkt]['snr'] = float(avgsnr)
        # check for RBAR BACKOFF
        if (obj=="RBAR") and (evt=="BACKOFF"):
            assert hasfield(e, 'net-root')
            ip = e['net-root']
            if (sendbuffer[ip]['mac'] is None):
                # new packet being sent by MAC
                sendbuffer[ip]['txts'] = float(ts)
                sendbuffer[ip]['mac'] = pkt
                sendbuffer[ip]['snr'] = float(avgsnr)
                currpkt = ip
        # check for AGT RCV
        if (obj=="AGT") and (evt=="RCV"):
            assert hasfield(e, 'net-root')
            ip = e['net-root']
            assert (ip==currpkt)
            assert (ip in sendbuffer)
            # set rxts
            assert hasfield(e, *agtargs)
            assert (sendbuffer[ip]['rxts'] is None)
            sendbuffer[ip]['rxts'] = float(ts)
        # check for other parameters
        rcv = (evt=="RCV")
        if (obj=="80211N") and rcv:
            assert hasfield(e, 'dot11n-rmsdelay')
            tau = float(e['dot11n-rmsdelay'])
            rmsdelay.append(tau)
        if (obj=="AGT") and rcv:
            assert hasfield(e, 'plen')
            plen = int(e['plen'])
            if param['plen'] is None: param['plen'] = plen
            else:                     assert (plen==param['plen'])
    # calculate throughput from sendbuffer
    for p,d in sendbuffer.items():
        plen = param['plen']
        assert (plen is not None)
        if (d['rxts'] is None):
            x, y = d['snr'], 0       # log with zero throughput
        else:
            x, latency = d['snr'], d['rxts']-d['txts']
            L = 8*plen               # length in bits
            y = (1e-6*L/latency)     # throughput in Mbps
        dp = {'x': x, 'y': y, 'ndata': 1}
        data.append(dp)
    # set remaining parameters
    if rmsdelay:
        d = np.array(rmsdelay).mean()
        param['rmsdelay'] = d
    param['use-yerror'] = True
    return data, param

def parse_per():
    usage = "%prog [OPTIONS] TRACEFILE1 [TRACEFILE2 ...]\n" + \
            "      Writes parsed data to standard output."
    parser = OptionParser(usage=usage)
    parser.add_option("", "--tmin", dest="tmin", \
            type="float", default=-np.inf, \
            help="Specify simulation time to start from [default=%default]")
    parser.add_option("", "--tmax", dest="tmax", \
            type="float", default=np.inf, \
            help="Specify simulation time to stop at [default=%default]")
    parser.add_option("", "--use-average", dest="use_average", \
            action="store_true", default=False, \
            help="Bin data to average (not instantaneous) SNR [default=%default]")
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
    defaultparam = {'xlabel': "SNR (dB)", \
                    'ylabel': "${\\rm T}_{put}$ (Mbps)", \
                    'title': "Throughput vs. SNR", \
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
    parse_per()
