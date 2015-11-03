#!  /usr/bin/env python

"""
Runs tests for `digital` package.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-09-28 21:43:47 -0500 (Wed, 28 Sep 2011) $
* $LastChangedRevision: 5169 $

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

import unittest

from SimPy.Simulation import *
from scapy.all import *
from wins import *
from copy import copy, deepcopy

class TestDigital(unittest.TestCase):
    """Test elements in `digital` package."""

    def setUp(self):
        Trace.Global.reset()

    def test_mqam(self):
        """Test for `MQAM`."""
        verbose = 100
        mqam = MQAM(verbose=verbose)
        snr = 5
        #snr = numpy.array([0, 4, 8, 12])
        for mod in ["BPSK", "QPSK", "4QAM", "16qam", "64qam"]:
            mqam.mtype = mod
            uber = mqam.ber(snr)
            mqam.log(mod, ber="%.5g"%(uber) )
            #mqam.log(mod, ber=["%.5g"%(u) for u in uber] )
        #mqam.trace.output()
        # check output
        qpsk, qam4 = 0, 1
        bpsk, qam16, qam64 = None, None, None
        for e in mqam.trace.events:
            if (e['event']=="BPSK"):  bpsk = float(e['ber'])
            if (e['event']=="QPSK"):  qpsk = float(e['ber'])
            if (e['event']=="4QAM"):  qam4 = float(e['ber'])
            if (e['event']=="16QAM"): qam16 = float(e['ber'])
            if (e['event']=="64QAM"): qam64 = float(e['ber'])
        self.assertAlmostEquals(qpsk, qam4, 5)      # (QPSK==4QAM) error!
        self.assertTrue(bpsk<=qpsk<=qam16<=qam64)   # M-QAM ber() failed!

    def test_rcpc(self):
        """Test for `RCPC`."""
        verbose = 100
        rcpc = RCPC(verbose=verbose)
        mqam = MQAM(verbose=verbose)
        mod, snr = "BPSK", 7
        mqam.mtype = mod
        uber = mqam.ber(snr)
        #uber = numpy.array([0, 1e-3, 0.5, 1.0])
        for rate in ["1/2", "2/3", "3/4", "4/5"]:
            rcpc.rate = rate
            cber = rcpc.ber(uber)
            rcpc.log(rate, ber="%.5g"%(cber) )
            #rcpc.log(rate, ber=["%.5g"%(c) for c in cber] )
        #rcpc.trace.output()
        # check output
        r12, r23, r34, r45 = None, None, None, None
        for e in rcpc.trace.events:
            if (e['event']=="1/2"):  r12 = float(e['ber'])
            if (e['event']=="2/3"):  r23 = float(e['ber'])
            if (e['event']=="3/4"):  r34 = float(e['ber'])
            if (e['event']=="4/5"):  r45 = float(e['ber'])
        errmsg = "r12 = %.5g, r23 = %.5g, r34 = %.5g, r45 = %.5g, uber = %.5g"%(r12,r23,r34,r45,uber)
        self.assertTrue(r12<=r23<=r34<=r45<uber, errmsg)    # RCPC ber() failed!

    def test_rs(self):
        """Test for `ReedSolomon`."""
        verbose = 100
        n, k = 255, 239
        plen = 1024
        rs = ReedSolomon(n, k, verbose=verbose)
        self.assertAlmostEqual(rs.rate, (1.0*k/n) )
        mqam = MQAM(verbose=verbose)
        mod, snr = "BPSK", 5
        mqam.mtype = mod
        uber = mqam.ber(snr)
        uper = mqam.per(plen, snr)
        rates = [(255,239)]
        for (n,k) in rates:
            rs.rate = (n,k)
            self.assertAlmostEqual(rs.rate, (1.0*k/n) )
            per = rs.per(plen, uber)
            self.assertTrue(per<uper)   # RS per() failed!
