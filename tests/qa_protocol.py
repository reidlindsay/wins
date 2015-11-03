#!  /usr/bin/env python

"""
Runs tests for `protocol` package.

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

class Agent(Element):
    """Simple agent to act as source and sink."""
    name="agent"
    tracename="AGT"
    def __init__(self, **kwargs):
        self.nsent, self.nrcvd = 0, 0
        self.dest = None
        Element.__init__(self, **kwargs)
    def configure(self, dest=None, rate=None, duration=None, **kwargs):
        if rate is None: rate=1.0
        if dest: self.dest = dest
        self.addport("TX")
        self.addport("RX")
        f = self.newchild("txfsm", FSM, tracename=self.tracename+".TX")
        g = self.newchild("rxfsm", FSM, tracename=self.tracename+".RX")
        f.goto(self.SEND, rate, duration)
        g.goto(self.RECV)
    def connect(self, net):
        self.TX.connect(net.RXU)
        net.TXU.connect(self.RX)
    def SEND(self, fsm, rate, duration):
        delay = 0
        if rate>0: delay = 1.0/rate
        if duration is None: duration = 1.0e-3
        while fsm.active() and (delay>0):
            self.log("wait", delay=delay)
            yield hold, fsm, delay
            p = Packet()/"helloworld"
            crcupdate(p)
            p.setanno('cif-duration', duration)
            if self.dest: p.setanno('net-dst', self.dest)
            self.log("snd", p, duration=time2usec(duration) )
            yield self.TX.send(fsm, [p])
            self.nsent += 1
            yield hold, fsm, duration
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

class MyARP(Element):
    name = "arp"
    tracename = "ARP"
    _share = {NET.broadcast:MAC.broadcast} # translate net->mac address
    def __init__(self, **kwargs):
        self.__netaddr = None
        self.__macaddr = None
        Element.__init__(self, **kwargs)

    netaddr = property(fget=lambda self: self.__netaddr)
    macaddr = property(fget=lambda self: self.__macaddr)

    def configure(self, netaddr=None, macaddr=None, **kwargs):
        if netaddr and macaddr: self.setaddr(netaddr, macaddr)
        # add ports and FSM
        self.addport("RXU"), self.addport("TXU")
        self.addport("TXD"), self.addport("RXD")
        txfsm = self.newchild("txfsm",FSM,tracename=self.tracename+".TX")
        rxfsm = self.newchild("rxfsm",FSM,tracename=self.tracename+".RX")
        # init states
        txfsm.goto(self.SEND)
        rxfsm.goto(self.RECV)

    def setaddr(self, netaddr, macaddr):
        if self.netaddr or self.macaddr: del self._share[self.netaddr]
        self.__netaddr = netaddr
        self.__macaddr = macaddr
        # register with _share
        self._share[netaddr] = macaddr

    def SEND(self, fsm):
        while fsm.active():
            yield self.RXU.recv(fsm, 1)
            p = fsm.got[0]
            pdst = p.dst #ANNO.supported(p) and p.getanno('net-dst')
            if p.hasanno('net-dst'): pdst = p.getanno('net-dst')
            assert (pdst in self._share)   # p.dst in _share?
            src, dst, etype = self.macaddr, self._share[pdst], const.ETHERTYPE_IP
            eth = Ether(src=src,dst=dst,type=etype)/p
            self.log("snd", eth, src=src, dst=dst, type=etype)
            yield self.TXD.send(fsm, [eth])
        yield fsm.stop()

    def RECV(self, fsm):
        while fsm.active():
            yield self.RXD.recv(fsm, "all")
            for p in fsm.got:
                self.log("rcv", p)
                pay = p.payload
                crcremove(pay)
                yield self.TXU.send(fsm, [pay])
        yield fsm.stop()

class TestProtocol(unittest.TestCase):
    """Test elements in `protocol` package."""

    def setUp(self):
        Trace.Global.reset()

    def test_routing(self):
        """Test `Routing`."""
        class Node(Element):
            name = "node"
            tracename = "NODE"
            def __init__(self, **kwargs):
                Element.__init__(self, **kwargs)
            def configure(self, rate=None, duration=None, pos=None, **kwargs):
                agt = self.newchild('agent', Agent, rate=rate, duration=duration)
                net = self.newchild('net', Routing)
                cif = self.newchild('cif', Radio)
                mobi = self.newchild('motion', Motion, pos=pos)
                # connect ports
                agt.connect(net)
                net.TXD.connect(cif.RXU)
                cif.TXU.connect(net.RXD)

        # set up parameters
        initialize()
        stoptime = 10.1
        verbose = 80
        kwargs = {'verbose':verbose}
        ch = Channel(model=Propagation, **kwargs)
        n0 = Node(rate=0.5, pos=(0,  0), **kwargs)
        n1 = Node(rate=0.0, pos=(0, 50), **kwargs)
        n2 = Node(rate=0.0, pos=(0,100), **kwargs)
        # connect to channel
        ch.add_edge(n0.cif, n1.cif)
        ch.add_edge(n1.cif, n2.cif)
        # set agent destination
        d0 = n0.net.address
        d1 = n1.net.address
        d2 = n2.net.address
        n0.agent.dest = d2
        n0.net.addroute(d2, d1)
        n1.net.addroute(d2)

        # run simulation
        simulate(until=stoptime)
        #ch.trace.output()

        """
        # check for errors
        for n in [n0, n1, n2]:
            ch.stderr("%s: nsent = %s\n"%(n, n.agent.nsent) )
            ch.stderr("%s  nrcvd = %s\n"%(" "*len(str(n)), n.agent.nrcvd) )
        """
        # check output
        nsent = [0, 0, 0]
        nrcvd = [0, 0, 0]
        nfwrd = [0, 0, 0]
        ndrop = [0, 0, 0]
        nid0, nid1, nid2 = n0.uid, n1.uid, n2.uid
        rtgid0 = n0.net.uid
        for e in ch.trace.events:
            if (e['event']=="SND") and (e['obj']=="RTG"):
                if (int(e['nid'])==nid0): nsent[0] += 1
                if (int(e['nid'])==nid1): nsent[1] += 1
                if (int(e['nid'])==nid2): nsent[2] += 1
            if (e['event']=="RCV") and (e['obj']=="RTG") and (e['packet']=='Packet'):
                if (int(e['nid'])==nid0): nrcvd[0] += 1
                if (int(e['nid'])==nid1): nrcvd[1] += 1
                if (int(e['nid'])==nid2): nrcvd[2] += 1
            if (e['event']=="FWD") and (e['obj']=="RTG"):
                if (int(e['nid'])==nid0): nfwrd[0] += 1
                if (int(e['nid'])==nid1): nfwrd[1] += 1
                if (int(e['nid'])==nid2): nfwrd[2] += 1
            if (e['event']=="DRP") and (e['obj']=="RTG"):
                if (int(e['uid'])==rtgid0) and (e['drop']=="not for me"): ndrop+=1
        self.assertEqual(n0.agent.nsent, n2.agent.nrcvd) # routing failed
        self.assertEqual(nsent[0], nrcvd[2])    # send error
        self.assertEqual(nsent[0], nfwrd[1])    # forward error
        self.assertEqual(nsent[2], 0)           # send error
        self.assertEqual(nrcvd[2], n2.agent.nrcvd) # did not send upstream

    def test_aloha(self):
        """Test `Aloha`."""
        class Node(Element):
            name = "node"
            tracename = "NODE"
            def __init__(self, **kwargs):
                Element.__init__(self, **kwargs)
            def configure(self, rate=None, duration=None, pos=None, **kwargs):
                agt = self.newchild('agent', Agent, rate=rate, duration=duration)
                net = self.newchild('net', Routing)
                mac = self.newchild('mac', Aloha)
                arp = self.newchild('arp', MyARP, \
                                    netaddr=net.addr, macaddr=mac.addr)
                cif = self.newchild('cif', Radio)
                mobi = self.newchild('motion', Motion, pos=pos)
                # connect ports
                agt.connect(net)
                net.TXD.connect(arp.RXU), arp.TXU.connect(net.RXD)
                arp.TXD.connect(mac.RXU), mac.TXU.connect(arp.RXD)
                mac.connect(cif)
                #arp.TXD.connect(cif.RXU), cif.TXU.connect(arp.RXD)

        # set up parameters
        initialize()
        stoptime = 10.1
        verbose = 80
        kwargs = {'verbose':verbose}
        ch = Channel(model=Propagation, **kwargs)
        n0 = Node(rate=0.5, pos=(0,  0), **kwargs)
        n1 = Node(rate=0.0, pos=(0, 50), **kwargs)
        n2 = Node(rate=0.0, pos=(0,100), **kwargs)
        #n3 = Node(rate=0.0, pos=(0,200), **kwargs)
        ch.add_edge(n0.cif, n1.cif)
        ch.add_edge(n1.cif, n2.cif)
        ch.add_edge(n1.cif, n0.cif)
        ch.add_edge(n0.cif, n2.cif)
        #ch.add_edge(n3.cif, n0.cif)
        # set agent destination
        d0 = n0.net.address
        d1 = n1.net.address
        d2 = n2.net.address
        n0.agent.dest = n2.net.address
        n0.net.addroute(d2, d1)
        n1.net.addroute(d2)

        # run simulation
        simulate(until=stoptime)
        #ch.trace.output()

        # check output
        nsent = [0, 0, 0]
        nrcvd = [0, 0, 0]
        ndrop = [0, 0, 0]
        nid0, nid1, nid2 = n0.uid, n1.uid, n2.uid
        aid0, aid1, aid2 = n0.mac.uid, n1.mac.uid, n2.mac.uid
        for e in ch.trace.events:
            if (e['event']=="SND") and (e['obj']=="ALOHA"):
                if (int(e['nid'])==nid0): nsent[0] += 1
                if (int(e['nid'])==nid1): nsent[1] += 1
                if (int(e['nid'])==nid2): nsent[2] += 1
            if (e['event']=="RCV") and (e['obj']=="ALOHA"):
                if (int(e['nid'])==nid0): nrcvd[0] += 1
                if (int(e['nid'])==nid1): nrcvd[1] += 1
                if (int(e['nid'])==nid2): nrcvd[2] += 1
            if (e['event']=="DRP") and (e['obj']=="ALOHA"):
                if (int(e['uid'])==aid0) and (e['drop']=="not for me"): ndrop[0]+=1
                if (int(e['uid'])==aid1) and (e['drop']=="not for me"): ndrop[1]+=1
                if (int(e['uid'])==aid2) and (e['drop']=="not for me"): ndrop[2]+=1
        self.assertEqual(n0.agent.nsent, n2.agent.nrcvd) # pass-upstream failed
        self.assertEqual(nsent[0], nrcvd[1])    # send error
        self.assertEqual(nsent[1], nrcvd[2])    # forward error
        self.assertEqual(ndrop[0], nsent[1])    # address match error
        self.assertEqual(nsent[0], ndrop[2])    # address match error

    def test_arp(self):
        """Test `ARP` module."""
        class Node(Element):
            name = "node"
            tracename = "NODE"
            def __init__(self, **kwargs):
                Element.__init__(self, **kwargs)
            def configure(self, rate=None, duration=None, pos=None, \
                                useshared=False, **kwargs):
                agt = self.newchild('agent', Agent, rate=rate, duration=duration)
                net = self.newchild('net', Routing)
                mac = self.newchild('mac', Aloha)
                arp = self.newchild('arp', ARP, useshared=useshared)
                cif = self.newchild('cif', Radio)
                mobi = self.newchild('motion', Motion, pos=pos)
                # connect ports
                agt.connect(net)
                arp.connect(net, mac)
                mac.connect(cif)
                # add kludge
                mac.duration = lambda *args,**kwargs: self.duration(*args,**kwargs)
            def duration(self, p):
                d = 1e-4 #const.EPSILON
                if p.hasanno('cif-duration'): d = p.getanno('cif-duration')
                p.setanno('cif-duration', d)
                return d

        # set up parameters
        initialize()
        stoptime = 20.1
        verbose = 80
        useshared = True
        kwargs = {'verbose':verbose, 'useshared': useshared}
        ch = Channel(model=Propagation, **kwargs)
        n0 = Node(rate=0.5, pos=(0,  0), **kwargs)
        n1 = Node(rate=0.0, pos=(0, 50), **kwargs)
        n2 = Node(rate=0.0, pos=(0,100), **kwargs)
        #n3 = Node(rate=0.0, pos=(0,200), **kwargs)
        nlist = [n0, n1, n2]
        # connect interfaces
        for n in nlist:
            for m in nlist:
                if (n is not m): ch.add_edge(n.cif, m.cif)
        # set agent destination
        d0 = n0.net.address
        d1 = n1.net.address
        d2 = n2.net.address
        n0.agent.dest = d2
        n0.net.addroute(d2, d1)     # n0->n1
        n1.net.addroute(d2)         # n1->n2

        # run simulation
        simulate(until=stoptime)
        #ch.trace.output()
        """
        # check for errors
        for n in [n0, n1, n2]:
            ch.stderr("%s: nsent = %s\n"%(n, n.agent.nsent) )
            ch.stderr("%s  nrcvd = %s\n"%(" "*len(str(n)), n.agent.nrcvd) )
        """

        # check output
        nsent = [0,0,0]
        nrcvd = [0,0,0]
        ndrop = [0,0,0]
        nsreq, nrreq = [0,0,0], [0,0,0]     # requests sent/rcvd
        nsrep, nrrep = [0,0,0], [0,0,0]     # replies sent/rcvd
        ndreq, ndrep = [0,0,0], [0,0,0]     # dropped requests/replies
        nid0, nid1, nid2 = n0.uid, n1.uid, n2.uid
        aid0, aid1, aid2 = n0.arp.uid, n1.arp.uid, n2.arp.uid
        reason1 = "dst not found in IPSEND"
        reason2 = "ARP request not for me in ARPRECV"
        for e in ch.trace.events:
            if (e['event']=="SND") and (e['obj']=="ARP"):
                for (id, k) in [(nid0,0), (nid1,1), (nid2,2)]:
                    if (int(e['nid'])==id):
                        if (e['root']=="ARP Request"): nsreq[k] += 1
                        elif (e['root']=="ARP Reply"): nsrep[k] += 1
                        else:                          nsent[k] += 1
            if (e['event']=="RCV") and (e['obj']=="ARP"):
                for (id, k) in [(nid0,0), (nid1,1), (nid2,2)]:
                    if (int(e['nid'])==id):
                        if (e['packet']=="ARP Request"): nrreq[k] += 1
                        elif (e['packet']=="ARP Reply"): nrrep[k] += 1
                        else:                            nrcvd[k] += 1
            if (e['event']=="DRP") and (e['obj']=="ARP"):
                for (id, k) in [(aid0,0), (aid1,1), (aid2,2)]:
                    if (int(e['uid'])==id):
                        if (e['packet']=="ARP Request"): ndreq[k] += 1
                        elif (e['packet']=="ARP Reply"): ndrep[k] += 1
                        else:                            ndrop[k] += 1
        if useshared:
            self.assertEqual(ndrop[0], 0)   # ARP lookup failed in IPSEND
            self.assertEqual(ndrop[1], 0)   # ARP lookup failed in IPSEND
            for k in [0,1,2]:
                self.assertEqual(nsreq[k], 0)   # Shared ARP fail -> requests!
                self.assertEqual(nsrep[k], 0)   # Shared ARP fail -> replies!
        else:
            self.assertEqual(ndrop[0], 1)   # ARP lookup failed in IPSEND
            self.assertEqual(ndrop[1], 1)   # ARP lookup failed in IPSEND
            self.assertEqual(nsreq[0], 1)   # failed to send ARP request
            self.assertEqual(nsrep[1], 1)   # failed to send ARP reply
            self.assertEqual(nrrep[0], 2)   # failed to recv ARP reply
            self.assertEqual(nsreq[1], 1)   # failed to send ARP request
            self.assertEqual(nsrep[2], 1)   # failed to send ARP reply
            self.assertEqual(nrrep[1], 1)   # failed to recv ARP reply
            self.assertEqual(ndreq[2], 1)   # drop 'not for me' request failed
            self.assertEqual(ndreq[0], 1)   # drop 'not for me' request failed
            for k in [0,1,2]:
                self.assertEqual(ndrep[k], 0)   # dropped ARP Reply!
        tdrop = 0
        for x in ndrop: tdrop += x
        self.assertEqual(n0.agent.nsent, n2.agent.nrcvd+tdrop) # end-to-end fail
