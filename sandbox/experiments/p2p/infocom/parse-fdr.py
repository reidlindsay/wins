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

def parse_per_info(trace, fmt='bo', usemodel=False):
    # initialize parameters
    param, data = {}, []
    mcs, rmsdelay = None, []
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
                sinr = float(e['dot11n-sinr'].lower().replace("db","") )
                drop = e['drop']
                notdetected = (drop==DETECTFAIL1)
                hparamfail  = (drop==HEADERFAIL1)
                headerfail  = (drop==HEADERFAIL2)
                ignorefail  = (drop==IGNOREFAIL1) or (drop==IGNOREFAIL2)
                assert (notdetected or hparamfail or headerfail)
                if usemodel:
                    fdr = float(e['dot11n-fdr'])
                else:
                    fdr = 0.0
                    if (hparamfail or headerfail): fdr = 1.0
                # log outcome of frame detection
                x, y = sinr, fdr
            elif rcv:
                sinr = float(e['dot11n-sinr'].lower().replace("db","") )
                if usemodel:
                    fdr = float(e['dot11n-fdr'])
                else:
                    fdr = 1.0
                # log outcome of frame detection
                x, y = sinr, fdr
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
    param['label']    = "${\\rm FDR}_{%s}$ ${\\rm (MCS = %d}$, "%(pertype,mcs) + \
                        "$\\sigma_{rms} = %.3g ns)$"%(avgdelay*1e9)
    return param, data

def parse_per():
    usage = "%prog [OPTIONS] TRACEFILE1 [TRACEFILE2 ...]\n" + \
            "      Writes parsed data to standard output."
    parser = OptionParser(usage=usage)
    (options, args) = parser.parse_args()

    if len(args)<1:
        print "Insufficient number of arguments."
        parser.print_help()
        raise SystemExit

    tracefile = args[0:]
    numtraces = len(tracefile)

    # set parameters
    default_parameters = {'xlabel': "SNR (dB)", \
                  'ylabel': "FDR", \
                  'title': "Frame Detection Rate vs. SNR", \
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
        param, data = parse_per_info(trace)
        if data:
            parameters = copy(default_parameters)
            parameters.update(param)
            parameters['source'] = tfile
            parameters['format'] = fmt[0]
            assert (param['label'] is not None)
            parsed_data = {'parameters': parameters, 'data': data}
            sys.stdout.write("%s\n"%(parsed_data) )
        # parse model PER from trace
        param, data = parse_per_info(trace, usemodel=True)
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
