#!  /usr/bin/env python

"""
Parse PER vs. SINR data from trace files.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-10-19 17:04:02 -0500 (Wed, 19 Oct 2011) $
* $LastChangedRevision: 5220 $

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
from copy import copy

from numpy import array

def read_trace(options, tracefile):
    # load trace from file
    tr = Trace()
    tr.read(tracefile)
    # return trace
    return tr

DETECTFAIL1 = "not detected in LISTEN"
HEADERFAIL1 = "header parameters failed"
HEADERFAIL2 = "header decoding failed"
IGNOREFAIL1 = "ignore rxdata in DECODE"
IGNOREFAIL2 = "ignore detect in DECODE"

def parse_per_info(options, trace, fmt='bo', usemodel=False):
    # initialize parameters
    param, data = {}, []
    mcs, rmsdelay = None, []
    ncollision = options.ncollision
    # parse trace
    for e in trace.events:
        obj, evt = e['obj'], e['event']
        # check for MCS parameter
        if ('phy-rate' in e):
            rate = int(e['phy-rate'])
            hparamfail = ('drop' in e) and (e['drop']==HEADERFAIL1)
            if not hparamfail:
                if mcs is None: mcs = rate
                else:           assert (mcs == rate)
        # check for 802.11n RCV & DRP events
        if (obj=="80211N"):
            rcv, drp = (evt=="RCV"), (evt=="DRP")
            x, y = None, None
            if drp:
                drop = e['drop']
                notdetected = (drop==DETECTFAIL1)
                hparamfail  = (drop==HEADERFAIL1)
                headerfail  = (drop==HEADERFAIL2)
                ignorefail  = (drop==IGNOREFAIL1) or (drop==IGNOREFAIL2)
                assert (notdetected or hparamfail or headerfail or ignorefail), "%s"%(e)
                #sinr = float(e['dot11n-sinr'].lower().replace("db","") )
                #x, y = sinr, 1.0       # log header drop as a packet error also
            elif rcv:
                sinr = float(e['dot11n-sinr'].lower().replace("db","") )
                err  = e['crc']
                haserror = (err=="FAIL")
                noerror  = (err=="OK")
                assert (haserror or noerror)
                if usemodel:
                    per = float(e['dot11n-model-per'])
                else:
                    if haserror: per = 1.0
                    else:        per = 0.0
                # check if ncollision matches
                keepdata = True
                if (ncollision is not None):
                    keepdata = False
                    if 'cif-collision' in e:
                        coll = eval(e['cif-collision'])
                        assert isinstance(coll, list)
                        keepdata = (len(coll) == ncollision)
                if keepdata:
                    x, y = sinr, per
            # log data point
            if (x is not None) and (y is not None):
                dp = {'x':x, 'y':y, 'ndata': 1}
                data.append(dp)
            # check for RMS delay
            if (rcv or drp):
                tau = float(e['dot11n-rmsdelay'])
                rmsdelay.append(tau)
    # check parameters
    assert (rmsdelay)
    assert (mcs is not None)
    avgdelay = array(rmsdelay).mean()
    pertype = "actual"
    if usemodel: pertype = "model"
    # return param and data
    param['mcs']      = mcs
    param['rmsdelay'] = avgdelay
    param['format']   = fmt
    label             = "${\\rm PER}_{%s}$ ${\\rm (MCS = %d}$, "%(pertype,mcs)
    if ncollision is not None: label +="$N_{coll} = %d$, "%(ncollision)
    label            += "$\\sigma_{rms} = %.3g ns)$"%(avgdelay*1e9)
    param['label']    = label
    return param, data

def parse_per():
    usage = "%prog [OPTIONS] TRACEFILE1 [TRACEFILE2 ...]\n" + \
            "      Writes parsed data to standard output."
    parser = OptionParser(usage=usage)
    parser.add_option("-c", "--ncollision", dest="ncollision", type="int", \
            default=None, help="Filter results using number of collisions. [default=%default]")
    (options, args) = parser.parse_args()

    if len(args)<1:
        print "Insufficient number of arguments."
        parser.print_help()
        raise SystemExit

    tracefile = args[0:]
    numtraces = len(tracefile)

    # set parameters
    default_parameters = {'xlabel': "SINR (dB)", \
                  'ylabel': "PER", \
                  'title': "PER vs. SINR", \
                  'label': None, \
                  'source': None, \
                  'format': None}

    lgd, formats = [], [('ro','r:'), ('bo', 'b:'), ('go', 'g:')]
    for k in range(numtraces):
        tfile = tracefile[k]
        # treat as normal wins trace file
        trace = read_trace(options, tfile)
        fmt = formats[k%len(formats)]
        if not trace: continue
        sys.stderr.write("Parsing trace from %s ...\n"%(tfile))
        # parse actual PER from trace
        param, data = parse_per_info(options, trace)
        if data:
            parameters = copy(default_parameters)
            parameters.update(param)
            parameters['source'] = tfile
            parameters['format'] = fmt[0]
            assert (param['label'] is not None)
            parsed_data = {'parameters': parameters, 'data': data}
            sys.stdout.write("%s\n"%(parsed_data) )
        # parse model PER from trace
        param, data = parse_per_info(options, trace, usemodel=True)
        if data:
            parameters = copy(default_parameters)
            parameters.update(param)
            parameters['source'] = tfile
            parameters['format'] = fmt[1]
            assert (param['label'] is not None)
            parsed_data = {'parameters': parameters, 'data': data}
            sys.stdout.write("%s\n"%(parsed_data) )

if __name__ == '__main__':
    parse_per()
