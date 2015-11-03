#!  /usr/bin/env python

"""
Parse PER data for ALOHA network from trace files.

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
    cdata = {}
    param = {'mcs': None, 'plen': None}
    cdata['all'] = {'data': []}
    rmsdelay = []
    # event parameters
    phyargs     = ('dot11n-sinr', 'dot11n-channel-fading')
    hydraerror  = ('dot11n-error', 'cif-collision')
    modelerror  = ('dot11n-model-per', 'dot11n-model-error')
    # paramaters for logging data
    maxcoll = options.maxcoll
    for e in trace.events:
        # check timestamp - [tmin, tmax]
        ts = float(e['ts'])
        if (ts<options.tmin): continue
        elif (ts>options.tmax): break
        # get event/object/packet parameters
        obj, evt = e['obj'], e['event']
        pkt = "%s.%s"%(e['packet'], e['pid'])
        # check for 80211N RCV
        if (obj=="80211N") and (evt=="RCV"):
            assert hasfield(e, *phyargs)
            # get SINR
            sinr = float(e['dot11n-sinr'].lower().replace("db","") )
            if options.use_average:
                sinr -= float(e['dot11n-channel-fading'].lower().replace("db","") )
            # get PER
            per = 0.0
            assert hasfield(e, *hydraerror)
            perror = eval(e['dot11n-error'])
            if perror: per = 1.0
            # create new data point
            dp = {'x':sinr, 'y':per, 'ndata': 1}
            # check for number of collisions
            ncoll = len(eval(e['cif-collision']))
            if (ncoll<=maxcoll):
                # add new data point to collision-specifc entry
                if ncoll not in cdata:
                    cdata[ncoll] = {'data':[]}
                cdata[ncoll]['data'].append(dp.copy())
            # add new data to "all" entry
            cdata['all']['data'].append(dp.copy())
        # check for other parameters
        rcv = (evt=="RCV")
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
    # set remaining parameters
    if rmsdelay:
        d = np.array(rmsdelay).mean()
        param['rmsdelay'] = d
    param['use-yerror'] = True
    # update combined data with parameters
    for key, pdata in cdata.items():
        assert (param['mcs'] is not None)
        pdata['parameters'] = param.copy()
        labelname = options.label.lower()
        label  = "${\\rm PER}_{%s}$"%(labelname)
        label += " (MCS %d, $N_{coll}=%s$)"%(param['mcs'], key)
        pdata['parameters']['label'] = label
        pdata['parameters']['ncollision'] = key
    return cdata

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
    parser.add_option("", "--max-coll", dest="maxcoll", type="int", \
            default=3, help="Maximum number of collisions to consider [default=%default]")
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
    defaultparam = {'xlabel': "SINR (dB)", \
                    'ylabel': "PER", \
                    'title': "Packet-Error Rate vs. SINR", \
                    'label': None, \
                    'source': [], \
                    'format': None}

    # define formats
    formats = [(options.fmt, options.fmt[0]+"x:")]  # use same format for all
    #formats = [('ro','rx:'), ('bo', 'bx:'), ('go', 'gx:')]

    # parse data from each data set
    for k in range(numtraces):
        # get tracefile
        tfile = tracefile[k]
        trace = read_trace(options, tfile) # treat as normal wins trace file
        if not trace: continue
        # get other parameters
        fmt = formats[0]
        # parse trace to get combined data
        sys.stderr.write("Parsing per data from %s ...\n"%(tfile))
        cdata = parse_info(options, trace)
        for key, pdata in cdata.items():
            assert ('data' in pdata)
            assert ('parameters' in pdata)
            data = pdata['data']
            param = defaultparam.copy()
            param.update(pdata['parameters'])
            param['format'] = fmt[0]
            param['source'] = tfile
            param['ncollision'] = key
            tdata = {'data': data, 'parameters': param}
            sys.stdout.write("%s\n"%(tdata))
            if isinstance(key, int):
                assert (key<=options.maxcoll), "Invalid Key (%d)!"%(key)

if __name__ == '__main__':
    parse_per()
