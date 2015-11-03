#!  /usr/bin/env python

"""
Parse throughput of ALOHA network from trace files.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-10-20 21:39:01 -0500 (Thu, 20 Oct 2011) $
* $LastChangedRevision: 5236 $

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
    param = {'mcs': None, 'plen': None, 'T': None}
    rmsdelay = []
    # event parameters
    phyargs     = ('dot11n-sinr', 'dot11n-channel-fading')
    hydraerror  = ('dot11n-error',)
    modelerror  = ('dot11n-model-per', 'dot11n-model-error')
    loadargs    = ('G', 'T', )
    # paramaters for logging data
    G = None
    nrcvd = 0
    stoptime = None
    for e in trace.events:
        # check timestamp - [tmin, tmax]
        ts = float(e['ts'])
        if (ts<options.tmin): continue
        elif (ts>options.tmax): break
        # get event/object/packet parameters
        obj, evt = e['obj'], e['event']
        pkt = "%s.%s"%(e['packet'], e['pid'])
        root = "%s.%s"%(e['root'], e['rid'])
        # parse tput
        rcv = (evt=="RCV")
        # check for ALOHA RCV
        if (obj=="ALOHA") and rcv: nrcvd += 1
        # check for other parameters
        if (evt=="WORKLOAD"):
            assert hasfield(e, 'G', 'T')
            g = float(e['G'])
            tstring = e['T']
            if 'usec' in tstring:
                t = float(tstring.replace("usec",""))*1e-6
            elif 'msec' in tstring:
                t = float(tstring.replace("msec",""))*1e-3
            else:
                t = float(tstring)
            if G is None: G = g
            else:         assert (abs(g-G)<const.EPSILON)
            if param['T'] is None: param['T'] = t
            else:                  assert (abs(t-param['T'])<const.EPSILON)
        if (evt=="STOPTIME"):
            assert hasfield(e, 'stoptime')
            t = float(e['stoptime'])
            if stoptime is None: stoptime = t
            else:                assert (abs(t-stoptime)<const.EPSILON)
        if (obj=="80211N") and rcv:
            assert hasfield(e, 'phy-rate', 'dot11n-rmsdelay')
            mcs = int(e['phy-rate'])
            if param['mcs'] is None: param['mcs'] = mcs
            else:                    assert (mcs==param['mcs'])
            tau = float(e['dot11n-rmsdelay'])
            rmsdelay.append(tau)
        if (obj=="AGT") and rcv:
            assert hasfield(e, 'plen')
            plen = int(e['plen'])
            if param['plen'] is None: param['plen'] = plen
            else:                     assert (plen==param['plen'])
    # create data point
    if stoptime>0:
        assert (G>0)
        T = param['T']
        assert (T>0)
        nslots = (1.0*stoptime/T)
        tput = nrcvd/nslots
        dp = {'x':G, 'y':tput, 'ndata': nslots}
        data.append(dp)
        param['use-yerror'] = True
    # set remaining parameters
    if rmsdelay:
        d = np.array(rmsdelay).mean()
        param['rmsdelay'] = d
    return data, param

def parse_tput():
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
            default="hydra", help="Name used in label subscript [default=%default]")
    (options, args) = parser.parse_args()

    if len(args)<1:
        print "Insufficient number of arguments."
        parser.print_help()
        raise SystemExit

    # get trace files
    tracefile = args[0:]
    numtraces = len(tracefile)

    # set default parameters
    defaultparam = {'xlabel': "G", \
                    'ylabel': "${\\rm T}_{put}$", \
                    'title': "Throughput vs. Workload", \
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
        sys.stderr.write("Parsing tput data from %s ... \n"%(tfile))
        data, dparam = parse_info(options, trace)
        if data:
            param = defaultparam.copy()
            param.update(dparam)
            param['format'] = fmt[0]
            param['source'] = tfile
            param['label']  = "${\\rm T}_{put}$"
            if 'mcs' in dparam: param['label'] += " (MCS %d)"%dparam['mcs']
            pdata = {'data': data, 'parameters':param}
            sys.stdout.write("%s\n"%(pdata))

if __name__ == '__main__':
    parse_tput()
