#!  /usr/bin/env python

"""
Parse FDR data from trace files.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-10-19 15:19:51 -0500 (Wed, 19 Oct 2011) $
* $LastChangedRevision: 5216 $

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

def parse_info(options, trace, usemodel=False):
    # (data, param) variables
    data = []
    param = {'mcs': None, 'plen': None}
    rmsdelay = []
    # event parameters
    phyargs     = ('dot11n-sinr', 'dot11n-channel-fading')
    hydradetect = ('dot11n-detect', )
    modeldetect = ('dot11n-model-fdr', 'dot11n-model-detect')
    # paramaters for logging data
    currpkt = None  # current packet being considered by parsing loop
    for e in trace.events:
        # check timestamp - [tmin, tmax]
        ts = float(e['ts'])
        if (ts<options.tmin): continue
        elif (ts>options.tmax): break
        # get event/object/packet parameters
        obj, evt = e['obj'], e['event']
        pkt = "%s.%s"%(e['packet'], e['pid'])
        # check for 80211N SND
        if (obj=="80211N") and (evt=="SND"):
            currpkt = pkt
        # check for 80211N DETECT or IGNORE
        detect, ignore = (evt=="DETECT"), (evt=="IGNORE")
        if (obj=="80211N") and (detect or ignore):
            assert (pkt==currpkt)
            assert hasfield(e, *phyargs)
            # get SINR
            sinr = float(e['dot11n-sinr'].lower().replace("db","") )
            if options.use_average:
                sinr -= float(e['dot11n-channel-fading'].lower().replace("db","") )
            # get FDR
            fdr = 0.0
            if usemodel:
                assert hasfield(e, *modeldetect)
                pdetect = eval(e['dot11n-model-detect'])
                if pdetect: fdr = 1.0
            else:
                assert hasfield(e, *hydradetect)
                pdetect = eval(e['dot11n-detect'])
                if pdetect: fdr = 1.0
            # create new data point
            dp = {'x':sinr, 'y':fdr, 'ndata': 1}
            data.append(dp)
        # check for other parameters
        rcv, drp = (evt=="RCV"), (evt=="DRP")
        if (obj=="80211N") and (rcv or drp):
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
            else:                    assert (plen==param['plen'])
    # set remaining parameters
    if rmsdelay:
        d = np.array(rmsdelay).mean()
        param['rmsdelay'] = d
    if not usemodel:
        param['use-yerror'] = True
    return data, param

def parse_fdr():
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
                    'ylabel': "FDR", \
                    'title': "Frame-Detection Rate vs. SNR", \
                    'label': None, \
                    'source': None, \
                    'format': None}

    # define formats
    formats = [('ro','r:'), ('bo', 'b:'), ('go', 'g:')]

    # parse data from each data set
    for k in range(numtraces):
        # get tracefile
        tfile = tracefile[k]
        trace = read_trace(options, tfile) # treat as normal wins trace file
        if not trace: continue
        # get other parameters
        fmt = formats[k%len(formats)]
        # parse trace to get (data, param) for HYDRA
        sys.stderr.write("Parsing trace from %s: \n"%(tfile))
        sys.stderr.write("      hydra ...\n")
        data, dparam = parse_info(options, trace, usemodel=False)
        if data:
            param = defaultparam.copy()
            param.update(dparam)
            param['format'] = fmt[0]
            param['source'] = tfile
            param['label']  = "Hydra PHY"
            pdata = {'data': data, 'parameters':param}
            sys.stdout.write("%s\n"%(pdata))
        # parse trace to get (data, param) for MODEL
        sys.stderr.write("      model ...\n")
        data, dparam = parse_info(options, trace, usemodel=True)
        if data:
            param = defaultparam.copy()
            param.update(dparam)
            param['format'] = fmt[1]
            param['source'] = tfile
            param['label']  = "WiNS Model"
            pdata = {'data': data, 'parameters':param}
            sys.stdout.write("%s\n"%(pdata))

if __name__ == '__main__':
    parse_fdr()
