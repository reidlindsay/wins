#!  /usr/bin/env python

"""
Parse throughput data from trace files.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-10-19 17:04:02 -0500 (Wed, 19 Oct 2011) $
* $LastChangedRevision: 5220 $

:author: Ketan Mandke <kmandke@mail.utexas.edu>

:copyright:
    Copyright 2009 The University of Texas at Austin

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
from copy import copy
import numpy as np

def read_trace(options, tracefile):
    # load trace from file
    tr = Trace()
    tr.read(tracefile)
    # return trace
    return tr

def parse_tput_info(trace, options, **kwargs):
    # initialize parameters
    ninfo = {}
    mcs, load, duration, stoptime = None, None, None, None
    for p in ['mcs', 'load', 'duration', 'stoptime']:
        exec("if '%s' in kwargs: %s = kwargs['%s']"%(p,p,p))
    ncollision = options.ncollision
    # parse trace
    for e in trace.events:
        # get event parameters
        ts = float(e['ts'])
        n, obj, evt = e['nid'], e['obj'], e['event']
        uid, pkt = e['uid'], e['packet']
        try:    plen = e['plen']
        except: plen = None
        # check timestamp
        if (ts<options.tmin):
            continue    # skip if ts < Tmin
        elif (ts>options.tmax):
            stoptime = ts
            break       # stop if ts > Tmax
        # add new node (if needed)
        if n not in ninfo:
            ninfo[n] = {'nsent':0, 'nrcvd':0, 'nignore':0, 'sinr': None}
        node = ninfo[n]     # node info to modify
        # update sender info
        G, T, R = None, None, None
        if 'G' in e: G = float(e['G'])
        if 'duration' in e: T = 1e-6*float(e['duration'].replace("usec", "") )
        if 'phy-rate' in e: R = int(e['phy-rate'])
        if (evt=="LOAD"):
            if load is None: load = G
            else: assert (abs(load-G)<const.EPSILON)
        if (obj=="AGT") and (evt=="WAIT"):
            if duration is None: duration = T
            else: assert (abs(duration-T)<const.EPSILON)
        if (obj=="80211N") and (evt=="SND"):
            if mcs is None: mcs = R
            else: assert (mcs==R)
        # skip while load or duration are not set
        if (load is None): continue
        if (duration is None): continue
        # update ninfo with AGT info
        if (obj=="AGT") and (evt=="SND"):
            node['nsent'] += 1
        if (obj=="AGT") and (evt=="RCV"):
            rcvd = True
            if ncollision is not None:
                if 'cif-collision' in e:
                    coll = eval(e['cif-collision'])
                    assert isinstance(coll, list)
                    if (len(coll)>ncollision): rcvd = False
            if rcvd: node['nrcvd'] += 1
            else:    node['nignore'] += 1
            if 'sinr' in e: node['sinr'] = e['sinr']
        # update stoptime
        if 'stoptime' in e: stoptime = float(e['stoptime'])
    # return node info and parameters
    for p in ['mcs', 'load', 'duration', 'stoptime']:
        exec("assert (%s is not None), '%s is None, ts = %.5f'"%(p, p, ts))
    param = {'mcs':mcs, 'load':load, 'duration':duration, 'stoptime':stoptime}
    return ninfo, param

def calc_network_tput(ninfo, param):
    # calculate actual network throughput -> return datapoint {x,y,ndata}
    load = param['load']
    duration = param['duration']
    stoptime = param['stoptime']
    # calculate Tput = Nrcvd/Nslots
    nslots = (stoptime/duration)
    nrcvd, nignore  = 0, 0
    for nid, data in ninfo.items():
        nrcvd += data['nrcvd']
        nignore += data['nignore']
    tput = nrcvd/nslots
    dp = {}
    if ninfo: dp = {'x':load, 'y':tput, 'ndata': nslots, 'nignore':nignore}
    return dp

def parse_tput():
    usage = "%prog [OPTIONS] TRACEFILE1 [TRACEFILE2 ...]\n" + \
            "      Writes parsed data to standard output."
    parser = OptionParser(usage=usage)
    parser.add_option("", "--tmin", dest="tmin", \
            type="float", default=0, \
            help="Minimum simulation timestamp to parse [default=%default]")
    parser.add_option("", "--tmax", dest="tmax", \
            type="float", default=np.inf, \
            help="Maximum simulation timestamp to parse [default=%default]")
    parser.add_option("-c", "--ncollision", dest="ncollision", type="int", \
            default=None, help="Filter results based on maximum number of collisions allowed. [default=%default]")
    (options, args) = parser.parse_args()

    if len(args)<1:
        print "Insufficient number of arguments."
        parser.print_help()
        raise SystemExit

    tracefile = args[0:]
    numtraces = len(tracefile)

    # set parameters
    default_parameters = {'xlabel': "G", \
                          'ylabel': "Throughput", \
                          'title': "Throughput vs. Traffic Load", \
                          'label': None, \
                          'source': None, \
                          'format': None}

    data, mcs, duration = [], None, None
    kwargs = {'stoptime':1.25*0.1}
    def check_param(p, pname, val):
        if pname in param:
            if p is None: p = val
            else: assert (abs(p-val)<const.EPSILON)
        return p

    for k in range(numtraces):
        tfile = tracefile[k]
        # treat as normal wins trace file
        trace = read_trace(options, tfile)
        if not trace: continue
        sys.stderr.write("Parsing trace from %s ...\n"%(tfile))
        # parse Tput info from trace
        ninfo, param = parse_tput_info(trace, options, **kwargs)
        # calculate network tput -> output data point
        dp = calc_network_tput(ninfo, param)
        if dp: data.append(dp)
        duration = check_param(duration, 'duration', param['duration'])
        mcs = check_param(mcs, 'mcs', param['mcs'])
    # set parameters and log data
    if data:
        param = {'format':'bo'}
        param['label'] = "${\\rm T}_{put}$  ${\\rm MCS = %d}$"%(mcs)
        parameters = copy(default_parameters)
        parameters.update(param)
        parameters['source'] = tfile
        parameters['title'] += " (${\\rm T}_{packet}$ = %.2f usec)"%(duration*1e6)
        assert (param['label'] is not None)
        parsed_data = {'parameters': parameters, 'data': data}
        sys.stdout.write("%s\n"%(parsed_data) )

if __name__ == '__main__':
    parse_tput()
