#!  /usr/bin/env python

"""
Parse throughput data for RBAR from trace files.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-10-21 15:51:37 -0500 (Fri, 21 Oct 2011) $
* $LastChangedRevision: 5246 $

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
    packetlength = None
    for p in ['packetlength']:
        exec("if '%s' in kwargs: %s = kwargs['%s']"%(p,p,p))
    pinfo, pktname, avgsnr = None, None, None
    snrdata, avgsnr = {}, None
    agtsnd = []
    # parse trace
    for e in trace.events:
        # get event parameters
        ts = float(e['ts'])
        n, obj, evt = e['nid'], e['obj'], e['event']
        uid, pkt, pid = e['uid'], e['packet'], e['pid']
        try:    plen = int(e['plen'])
        except: plen = None
        try:    pid = int(e['pid'])
        except: pid = None
        try:    sinr = float(e['sinr'].replace("dB", ""))
        except: sinr = None
        # check timestamp
        if (ts<options.tmin):
            continue    # skip if ts < Tmin
        elif (ts>options.tmax):
            stoptime = ts
            break       # stop if ts > Tmax
        # check if motion control information is found
        if (obj=="MC") and (evt=="AVGSNR"):
            avgsnr = e['snr']
            if (pinfo is None): pinfo = {}
            snrdata[avgsnr] = {'ts': [], 'snr': []}
        if (pinfo is None): continue
        # check for new packet being sent
        if (obj=="AGT") and (evt=="SND") and (pid>=0) and (plen>0):
            assert (avgsnr is not None)
            pktname = "%s.%s"%(obj, pid)
            assert (pktname not in pinfo)
            agtsnd.append(pktname)  # buffer pktname
            pinfo[pktname] = {}
            pinfo[pktname]['tstart'] = float(ts)
            pinfo[pktname]['tstop'] = None
            pinfo[pktname]['latency'] = None
            pinfo[pktname]['plen'] = plen
            pinfo[pktname]['avgsnr'] = float(avgsnr.replace("dB", ""))
            pinfo[pktname]['snr'] = None
            packetlength = plen
        # check for packet being delivered to receiver
        if (obj=="AGT") and (evt=="RCV") and (pid>=0) and (plen>0):
            pktname = "%s.%s"%(obj, pid)
            repeat = (pktname in pinfo) and (pinfo[pktname]['latency'] is not None)
            assert agtsnd or repeat
            assert (pktname in pinfo)
            assert (sinr is not None)
            if not repeat:
                pname = agtsnd.pop(0)
                while (pname!=pktname): pname = agtsnd.pop(0)
            pinfo[pktname]['tstop'] = float(ts)
            pinfo[pktname]['latency'] = float(ts)-pinfo[pktname]['tstart']
            pinfo[pktname]['snr'] = float(sinr)
        # update snr data
        if ('phy-sinr' in e):
            snrdata[avgsnr]['ts'].append(float(ts))
            snrdata[avgsnr]['snr'].append(float(e['phy-sinr'].replace("dB","")))
        elif ('dot11n-sinr' in e):
            snrdata[avgsnr]['ts'].append(float(ts))
            snrdata[avgsnr]['snr'].append(float(e['dot11n-sinr'].replace("dB","")))
        # update stoptime
        if 'stoptime' in e: stoptime = float(e['stoptime'])
    # return node info and parameters
    for p in ['packetlength']:
        exec("assert (%s is not None), '%s is None, ts = %.5f'"%(p, p, ts))
    param = {'packetlength':packetlength}
    return pinfo, snrdata, param

def calc_avgsnr(snrdata):
    # calculate actual average SNR for every expected SNR
    snr = {}
    for a, d in snrdata.items():
        x = np.array(d['ts'])
        y = np.array(d['snr'])
        snr[a] = a
        assert (len(x)==len(y))
        if (len(x)>1):
            dx = x[1:] - x[:-1]             # dt
            T = 1.0*(x[-1] - x[0])          # T
            avg = (0.5*(y[:-1]+y[1:])*dx).sum()/T       # integrate
            snr[a] = "%.4f dB"%(avg)
        elif (len(x)>0):
            snr[a] = "%.4f dB"%(d['snr'][0])
    return snr

def calc_tput(pinfo, options):
    # calculate tput for each average SNR value
    data, plen = [], None
    for p,d in pinfo.items():
        plen = d['plen']
        avgsnr = d['avgsnr']
        latency = d['latency']
        # get x-value: SNR
        x = d['snr']
        if options.use_average or (x is None): x = avgsnr
        # get y-value: Tput
        y = 0
        if latency is not None: y = 1e-6*(8.0*plen/latency)    # record Mbps
        # record data point
        data.append({'x':x, 'y':y, 'ndata':1})
    return data

def parse_tput():
    usage = "%prog [OPTIONS] TRACEFILE1 [TRACEFILE2 ...]\n" + \
            "      Writes parsed data to standard output."
    parser = OptionParser(usage=usage)
    parser.add_option("", "--tmin", dest="tmin", \
            type="float", default=-np.inf, \
            help="Minimum simulation timestamp to parse [default=%default]")
    parser.add_option("", "--tmax", dest="tmax", \
            type="float", default=np.inf, \
            help="Maximum simulation timestamp to parse [default=%default]")
    parser.add_option("", "--use-average", dest="use_average", \
            action="store_true", default=False, \
            help="Bin data to average or expected SNR [default=%default]")
    (options, args) = parser.parse_args()

    if len(args)<1:
        print "Insufficient number of arguments."
        parser.print_help()
        raise SystemExit

    tracefile = args[0:]
    numtraces = len(tracefile)

    # set parameters
    default_parameters = {'xlabel': "SNR", \
                          'ylabel': "Throughput (Mbps)", \
                          'title': "Throughput vs. SNR", \
                          'label': None, \
                          'source': None, \
                          'format': None}

    formats = ['bo', 'ro', 'go']
    for k in range(numtraces):
        tfile = tracefile[k]
        # treat as normal wins trace file
        trace = read_trace(options, tfile)
        if not trace: continue
        sys.stderr.write("Parsing trace from %s ...\n"%(tfile))
        # parse Tput info from trace
        pinfo, snrdata, param = parse_tput_info(trace, options)
        # calculate throughput data -> output new data and parameters
        data = calc_tput(pinfo, options)
        if data:
            plen = param['packetlength']
            param.update(default_parameters.copy())
            param['format'] = formats[k%len(formats)]
            param['source'] = tfile
            param['label'] = " ${\\rm T}_{put}$, $L = %d$ bytes"%(plen)
            assert (param['label'] is not None)
            parsed_data = {'parameters': param, 'data': data}
            sys.stdout.write("%s\n"%(parsed_data) )
        # calculate average SNR -> write to stderr
        actual = calc_avgsnr(snrdata)
        sys.stderr.write("Expected/Actual SNR:\n  %s\n"%(actual))

if __name__ == '__main__':
    parse_tput()
