#!  /usr/bin/env python

"""
Parse FDR data from trace files.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-11-01 19:59:54 -0500 (Tue, 01 Nov 2011) $
* $LastChangedRevision: 5321 $

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
    cdata['all'] = {'data': [], 'parameters': {}}
    rmsdelay = []
    # event parameters
    phyargs     = ('dot11n-sinr', 'dot11n-channel-fading')
    cifargs     = ('cif-collision', )
    # paramaters for logging data
    stoptime = None
    G, T = None, None
    senderid = None
    sendbuffer = {}
    curreth = None
    for e in trace.events:
        # check timestamp - [tmin, tmax]
        ts = float(e['ts'])
        if (ts<options.tmin): continue
        elif (ts>options.tmax): break
        # get event/object/packet parameters
        obj, evt = e['obj'], e['event']
        pkt = "%s.%s"%(e['packet'], e['pid'])
        root = "%s.%s"%(e['root'], e['rid'])
        try:    nid = int(e['nid'])
        except: nid = None
        # get for other parameters
        if (evt=="STOPTIME"):
            assert hasfield(e, 'stoptime')
            stoptime = float(e['stoptime'])
        if (obj=="80211N") and (evt=="RCV"):
            assert hasfield(e, 'phy-rate', 'dot11n-rmsdelay')
            mcs = int(e['phy-rate'])
            if param['mcs'] is None: param['mcs'] = mcs
            else:                    assert (mcs==param['mcs'])
            tau = float(e['dot11n-rmsdelay'])
            rmsdelay.append(tau)
        if (obj=="AGT") and (evt=="SND"):
            assert hasfield(e, 'plen')
            plen = int(e['plen'])
            if param['plen'] is None: param['plen'] = plen
            else:                     assert (plen==param['plen'])
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
            G, T = g, t
        # parse FDR -> check for sender ID
        if (evt=="SENDER"):
            assert hasfield(e, 'nodeid')
            senderid = int(e['nodeid'])
        if (senderid is None): continue
        if (G is None) or (T is None): continue
        # check for ALOHA SND
        if (obj=="ALOHA") and (evt=="SND") and (senderid==nid):
            curreth = pkt
        # check for 80211N SND
        if (obj=="80211N") and (evt=="SND"):
            if (root==curreth):
                assert (pkt not in sendbuffer)
                # add new packet to sendbuffer
                sendbuffer[pkt] = {'eth':root, 'detect':None, 'ncoll':None}
        # check for 80211N IGNORE
        if (obj=="80211N") and (evt=="IGNORE"):
            assert hasfield(e, *cifargs)
            ncoll = len(eval(e['cif-collision']))
            if (pkt in sendbuffer):
                assert (sendbuffer[pkt]['detect'] is  None)
                sendbuffer[pkt]['detect'] = False
                sendbuffer[pkt]['ncoll'] = ncoll
        # check for 80211N DETECT
        if (obj=="80211N") and (evt=="DETECT"):
            assert hasfield(e, 'dot11n-detect', *cifargs)
            ncoll = len(eval(e['cif-collision']))
            if (pkt in sendbuffer):
                assert (sendbuffer[pkt]['detect'] is None)
                sendbuffer[pkt]['detect'] = eval(e['dot11n-detect'])
                sendbuffer[pkt]['ncoll'] = ncoll
        # check for 80211N DRP
        if (obj=="80211N") and (evt=="DRP"):
            assert hasfield(e, 'drop', *cifargs)
            ncoll = len(eval(e['cif-collision']))
            if (pkt in sendbuffer) and (sendbuffer[pkt]['detect'] is None):
                # assume packet was not detected
                sendbuffer[pkt]['detect'] = False
                sendbuffer[pkt]['ncoll'] = ncoll
            # otherwise -> assume header decoding failed
            pass
    # calculate FDR data from sendbuffer
    assert (G is not None)
    nodetect = 0    # number of undetected packets
    for p, d in sendbuffer.items():
        assert hasfield(d, 'eth', 'detect', 'ncoll')
        eth, detect, ncoll = d['eth'], d['detect'], d['ncoll']
        maxcoll = options.maxcoll
        if (detect is None): nodetect += 1
        if (eth!=curreth) and (detect is not None):
            # create new data point
            assert (ncoll is not None)
            if detect: fdr = 1.0
            else:      fdr = 0.0
            dp = {'x': G, 'y': fdr, 'ndata': 1}
            # add to 'all', 'max', and/or ncoll
            cdata['all']['data'].append(dp)
            if (maxcoll is not None):
                if (ncoll<=maxcoll): key = ncoll
                else:                key = 'max'
                if (key not in cdata):
                    cdata[key] = {'data': [], 'parameters': {}}
                cdata[key]['data'].append(dp)
    # set remaining parameters
    if rmsdelay:
        d = np.array(rmsdelay).mean()
        param['rmsdelay'] = d
    param['use-yerror'] = True
    param['nodetect'] = nodetect
    assert (param['mcs'] is not None)
    assert (param['plen']>0)
    # update combined data with parameters
    for key, pdata in cdata.items():
        assert (param['mcs'] is not None)
        pdata['parameters'].update(param)
        label = options.label
        clabel = "= %s"%(key)
        if (key=='max'): clabel = "> %s"%(options.maxcoll)
        label += " (MCS %d, $N_{coll}$ %s)"%(param['mcs'], clabel)
        pdata['parameters']['label'] = label
        pdata['parameters']['ncollision'] = key
    return cdata

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
    parser.add_option("-f", "--fmt", dest="fmt", \
            default="bo", help="Format of data parsed from tracefiles [default=%default]")
    parser.add_option("", "--label", dest="label", \
            default="Hydra PHY", help="Name used in label subscript [default=%default]")
    parser.add_option("", "--max-coll", dest="maxcoll", type="int", \
            default=10, help="Maximum number of colliders to consider [default=%default]")
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
                    'ylabel': "FDR", \
                    'title': "Frame-Detection Ratio vs. Workload", \
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
        fmt = formats[0]
        # parse trace to get combined data
        sys.stderr.write("Parsing FDR data from %s ...\n"%(tfile))
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

if __name__ == '__main__':
    parse_fdr()
