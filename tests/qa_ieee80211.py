#!  /usr/bin/env python

"""
Runs tests for `ieee80211` package.

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
from wins.ieee80211 import Dot11APHY, Dot11ARadio
from wins.ieee80211 import Dot11NPHY, Dot11NRadio
from copy import copy, deepcopy

from wins.ieee80211 import Dot11NChannel

import numpy

DATARATE  = 4
PACKETLEN = 1000

class Agent(Element):
    """Simple agent to act as source and sink."""
    name="agent"
    tracename="AGT"
    def __init__(self, **kwargs):
        self.nsent, self.nrcvd = 0, 0
        self.dest = None
        Element.__init__(self, **kwargs)
    def configure(self, dest=None, rate=None, **kwargs):
        if rate is None: rate=1.0
        if dest: self.dest = dest
        self.addport("TX")
        self.addport("RX")
        f = self.newchild("txfsm", FSM, tracename=self.tracename+".TX")
        g = self.newchild("rxfsm", FSM, tracename=self.tracename+".RX")
        f.goto(self.SEND, rate)
        g.goto(self.RECV)
    def connect(self, net):
        self.TX.connect(net.RXU)
        net.TXU.connect(self.RX)
    def SEND(self, fsm, rate):
        global DATARATE, PACKETLEN
        delay = 0
        if rate>0: delay = 1.0/rate
        while fsm.active() and (delay>0):
            self.log("wait", delay=delay)
            yield hold, fsm, delay
            s = "helloworld"
            p = Packet()/(s*(PACKETLEN/len(s) ) )
            p.setanno('phy-rate', DATARATE)
            crcupdate(p)
            if self.dest: p.setanno('net-dst', self.dest)
            self.log("snd", p)
            yield self.TX.send(fsm, [p])
            self.nsent += 1
        yield fsm.stop()
    def RECV(self, fsm):
        while fsm.active():
            yield self.RX.recv(fsm, "all")
            for p in fsm.got:
                err = CRC32.haserror(p)
                if err: self.log("drp", p, crcerror=err)
                else:   self.log("rcv", p, crcerror=err)
                if not err: self.nrcvd += 1
        yield fsm.stop()

class TestIEEE80211(unittest.TestCase):
    """Test elements in `ieee80211` package."""

    def setUp(self):
        Trace.Global.reset()

    def test_dot11n(self):
        """Test `dot11n` module."""
        global DATARATE
        DATARATE = 2

        Dot11NPHY.usewaveform = False

        class Node(Element):
            name = "node"
            tracename = "NODE"
            def __init__(self, **kwargs):
                Element.__init__(self, **kwargs)
            def configure(self, rate=None, pos=None, useshared=False, **kwargs):
                cif = self.newchild('cif', Dot11NRadio)
                phy = self.newchild('phy', Dot11NPHY, radio=cif)
                mac = self.newchild('mac', Aloha, phy=phy)
                net = self.newchild('net', Routing)
                arp = self.newchild('arp', ARP, useshared=useshared)
                agt = self.newchild('agent', Agent, rate=rate)
                mobi = self.newchild('motion', Motion, pos=pos)
                # connect ports
                agt.connect(net)
                arp.connect(net, mac)
                mac.connect(phy)
                phy.connect(cif)

        # set up parameters
        initialize()
        stoptime = 2.1
        verbose = 80
        useshared = True
        kwargs = {'verbose':verbose, 'useshared': useshared}
        ch = Channel(model=Dot11NChannel, **kwargs)
        n0 = Node(rate=2.0, pos=(0,  0), **kwargs)
        n1 = Node(rate=0.0, pos=(0, 39), **kwargs)
        n2 = Node(rate=0.0, pos=(0, 40), **kwargs)
        #n3 = Node(rate=0.0, pos=(0,200), **kwargs)
        nlist = [n0, n1, n2]
        # connect interfaces
        for n in nlist:
            for m in nlist:
                if (n is not m): ch.add_edge(n.cif, m.cif, n=2.5)
        # set agent destination
        d0 = n0.net.address
        d1 = n1.net.address
        d2 = n2.net.address
        n0.agent.dest = n0.net.broadcast #d2
        n0.net.addroute(d2, d1)     # n0->n1
        n1.net.addroute(d2)         # n1->n2

        # run simulation
        print ""
        simulate(until=stoptime)
        """
        if verbose>40:
            ch.trace.output()
            
            # report errors
            for n in [n0, n1, n2]:
                ch.stderr("%s: nsent = %s\n"%(n, n.agent.nsent) )
                ch.stderr("%s  nrcvd = %s\n"%(" "*len(str(n)), n.agent.nrcvd) )
        """

        # check output
        peranno = 'dot11n-per'
        prob = [None, None, None]
        nid0, nid1, nid2 = n0.uid, n1.uid, n2.uid
        aid0, aid1, aid2 = n0.phy.uid, n1.phy.uid, n2.phy.uid
        for e in ch.trace.events:
            if (e['event']=="RCV") and (e['obj']=="80211N"):
                for (id,aid,k) in [(nid0,aid0,0), (nid1,aid1,1), (nid2,aid2,2)]:
                    if (int(e['nid'])==id) and (int(e['uid'])==aid):
                        if (peranno in e): prob[k] = float(e[peranno])
        # check PER
        nsent0 = n0.agent.nsent
        nrcvd1 = n1.agent.nrcvd
        nrcvd2 = n2.agent.nrcvd
        std = [None]*len(prob)
        for k in range(len(prob)):
            p = prob[k]
            if p is None: continue
            std[k] = (p*(1-p))**0.5+const.EPSILON
        per1 = 1.0 - (1.0*nrcvd1/nsent0)
        per2 = 1.0 - (1.0*nrcvd2/nsent0)
        beta = 0.50
        """
        ch.stderr("node %d: prob = %.2g, per = %.2g, std = %.2g\n"%(1, prob[1], per1, std[1]))
        ch.stderr("node %d: prob = %.2g, per = %.2g, std = %.2g\n"%(2, prob[2], per2, std[2]))
        """
        self.assertTrue(abs(prob[1]-per1)<beta*std[1])  # Dot11NPHY PER failed
        self.assertTrue(abs(prob[2]-per2)<beta*std[2])  # Dot11NPHY PER failed

    def not_test_dot11a(self):
        """Test `dot11a` module."""
        class Node(Element):
            name = "node"
            tracename = "NODE"
            def __init__(self, **kwargs):
                Element.__init__(self, **kwargs)
            def configure(self, rate=None, pos=None, useshared=False, **kwargs):
                cif = self.newchild('cif', Dot11ARadio)
                phy = self.newchild('phy', Dot11APHY, radio=cif)
                mac = self.newchild('mac', Aloha, phy=phy)
                net = self.newchild('net', Routing)
                arp = self.newchild('arp', ARP, useshared=useshared)
                agt = self.newchild('agent', Agent, rate=rate)
                mobi = self.newchild('motion', Motion, pos=pos)
                # connect ports
                agt.connect(net)
                arp.connect(net, mac)
                mac.connect(phy)
                phy.connect(cif)

        # set up parameters
        initialize()
        stoptime = 20.1
        verbose = 80
        useshared = True
        kwargs = {'verbose':verbose, 'useshared': useshared}
        ch = Channel(model=Propagation, **kwargs)
        n0 = Node(rate=2.0, pos=(0,  0), **kwargs)
        n1 = Node(rate=0.0, pos=(0, 39), **kwargs)
        n2 = Node(rate=0.0, pos=(0, 40), **kwargs)
        #n3 = Node(rate=0.0, pos=(0,200), **kwargs)
        nlist = [n0, n1, n2]
        # connect interfaces
        for n in nlist:
            for m in nlist:
                if (n is not m): ch.add_edge(n.cif, m.cif, n=2.5)
        # set agent destination
        d0 = n0.net.address
        d1 = n1.net.address
        d2 = n2.net.address
        n0.agent.dest = n0.net.broadcast #d2
        n0.net.addroute(d2, d1)     # n0->n1
        n1.net.addroute(d2)         # n1->n2

        # run simulation
        simulate(until=stoptime)
        #ch.trace.output()

        # check for errors
        """
        for n in [n0, n1, n2]:
            ch.stderr("%s: nsent = %s\n"%(n, n.agent.nsent) )
            ch.stderr("%s  nrcvd = %s\n"%(" "*len(str(n)), n.agent.nrcvd) )
        """

        # check output
        peranno = 'dot11a-per'
        prob = [None, None, None]
        nid0, nid1, nid2 = n0.uid, n1.uid, n2.uid
        aid0, aid1, aid2 = n0.phy.uid, n1.phy.uid, n2.phy.uid
        for e in ch.trace.events:
            if (e['event']=="RCV") and (e['obj']=="80211A"):
                for (id,aid,k) in [(nid0,aid0,0), (nid1,aid1,1), (nid2,aid2,2)]:
                    if (int(e['nid'])==id) and (int(e['uid'])==aid):
                        if (peranno in e): prob[k] = float(e[peranno])
        # check PER
        nsent0 = n0.agent.nsent
        nrcvd1 = n1.agent.nrcvd
        nrcvd2 = n2.agent.nrcvd
        std = [None]*len(prob)
        for k in range(len(prob)):
            p = prob[k]
            if p is None: continue
            std[k] = (p*(1-p))**0.5+const.EPSILON
        per1 = 1.0 - (1.0*nrcvd1/nsent0)
        per2 = 1.0 - (1.0*nrcvd2/nsent0)
        beta = 0.50
        """
        ch.stderr("node %d: prob = %.2g, per = %.2g, std = %.2g\n"%(1, prob[1], per1, std[1]))
        ch.stderr("node %d: prob = %.2g, per = %.2g, std = %.2g\n"%(2, prob[2], per2, std[2]))
        """
        self.assertTrue(abs(prob[1]-per1)<beta*std[1])  # Dot11APHY PER failed
        self.assertTrue(abs(prob[2]-per2)<beta*std[2])  # Dot11APHY PER failed
