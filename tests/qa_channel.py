#!  /usr/bin/env python

"""
Runs tests for `Channel`.

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
from wins import *
from scapy.all import *
from copy import copy, deepcopy

class Agent(Element):
    """Simple agent to act as source and sink."""
    name="agent"
    tracename="AGT"
    def __init__(self, **kwargs):
        self.nsent, self.nrcvd = 0, 0
        Element.__init__(self, **kwargs)
    def configure(self, rate=None, duration=None, **kwargs):
        if rate is None: rate=1.0
        self.addport("TX")
        self.addport("RX")
        f = self.newchild("txfsm", FSM, tracename=self.tracename+".TX")
        g = self.newchild("rxfsm", FSM, tracename=self.tracename+".RX")
        f.goto(self.SEND, rate, duration)
        g.goto(self.RECV)
    def SEND(self, fsm, rate, duration):
        delay = 0
        if rate>0: delay = 1.0/rate
        if duration is None: duration = 1.0
        while fsm.active() and (delay>0):
            yield hold, fsm, delay
            p = Packet()/"helloworld"
            crcupdate(p)
            p.setanno('cif-duration', duration)
            self.log("snd", p, duration=duration)
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


class TestChannel(unittest.TestCase):

    def setUp(self):
        Trace.Global.reset()

    def test_001(self):
        """Test basic `Channel`, `ChannelModel`, and `ChannelInterface`."""
        class Node(Element):
            name="node"
            tracename="NODE"
            def __init__(self, **kwargs):
                Element.__init__(self, **kwargs)
            def configure(self, rate=None, duration=None, **kwargs):
                agt = self.newchild('agent', Agent, rate=rate, duration=duration)
                cif = self.newchild('cif', ChannelInterface)
                # connect ports
                agt.TX.connect(cif.RXU)
                cif.TXU.connect(agt.RX)
                # monitor
                mon = self.newchild('monitor', FSM, tracename=self.tracename+".MON")
                mon.goto(self.MON, cif)
            def MON(self, fsm, cif):
                while fsm.active():
                    c = cif
                    yield waitevent, fsm, (c.txdata,c.txdone,c.rxdata,c.rxdone)
                    if c.txdata in fsm.eventsFired:
                        fsm.log("txdata", c.txdata.signalparam, cif=cif.traceid)
                    if c.txdone in fsm.eventsFired:
                        fsm.log("txdone", c.txdone.signalparam, cif=cif.traceid)
                    if c.rxdata in fsm.eventsFired:
                        fsm.log("rxdata", c.rxdata.signalparam, cif=cif.traceid)
                    if c.rxdone in fsm.eventsFired:
                        fsm.log("rxdone", c.rxdone.signalparam, cif=cif.traceid)
                yield fsm.stop()

        initialize()
        stoptime = 10.0
        verbose = 90
        kwargs = {'verbose':verbose}
        ch = Channel(**kwargs)
        n0 = Node(rate=0.5, **kwargs)
        n1 = Node(rate=0.5, **kwargs)
        n2 = Node(rate=0, **kwargs)
        n3 = Node(rate=0.5, **kwargs)
        ch.add_edge(n0.cif, n2.cif)
        ch.add_edge(n1.cif, n2.cif)
        ch.add_edge(n3.cif, n0.cif)
        simulate(until=stoptime)
        #ch.trace.output()
        # check halfduplex operation
        ndrp, nsent = 0, 0
        for e in ch.trace.events:
            if (e['event']=="DRP") and (e['obj']=="CIF") and (int(e['nid'])==0):
                ndrp += 1
            if (e['event']=="SND") and (e['obj']=="CIF") and (int(e['uid'])==3):
                nsent +=1
        self.assertEqual(ndrp, nsent)   # halfduplex mode error!
        # check for collisions
        ncoll, nsent0, nsent1 = 0, 0, 0
        for e in ch.trace.events:
            if (e['event']=="RCV") and (e['obj']=="CIF") and (int(e['nid'])==2):
                if ('cif-collision' in e) and (len(e['cif-collision'])>0): ncoll += 1
            if (e['event']=="SND") and (e['obj']=="CIF"):
                if   (int(e['uid'])==0): nsent0 +=1
                elif (int(e['uid'])==1): nsent1 +=1
        self.assertEqual(nsent0+nsent1, ncoll)  # collision error!

    def test_radio(self):
        """Test `Propagation` and `Radio`."""
        class MobileNode(Element):
            name = "mobilenode"
            tracename = "MNODE"
            def __init__(self, **kwargs):
                Element.__init__(self, **kwargs)
            def configure(self, rate=None, duration=None, pos=None, **kwargs):
                agt = self.newchild('agent', Agent, rate=rate, duration=duration)
                radio = self.newchild('radio', Radio)
                mobi = self.newchild('motion', Motion, pos=pos)
                # connect ports
                agt.TX.connect(radio.RXU)
                radio.TXU.connect(agt.RX)

        initialize()
        stoptime = 10.0
        verbose = 90
        kwargs = {'verbose':verbose}
        ch = Channel(model=Propagation, **kwargs)
        n0 = MobileNode(rate=0.5, pos=(0,  0), **kwargs)
        n1 = MobileNode(rate=0.5, pos=(0, 50), **kwargs)
        n2 = MobileNode(rate=0.0, pos=(0,100), **kwargs)
        n3 = MobileNode(rate=0.0, pos=(0,200), **kwargs)
        ch.add_edge(n0.radio, n2.radio)
        ch.add_edge(n1.radio, n2.radio)
        #ch.add_edge(n3.radio, n0.radio)
        simulate(until=stoptime)
        #ch.trace.output()
        # check freespace pathloss model
        prop = Propagation()
        d02 = prop.distance(n0.radio, n2.radio)
        d12 = prop.distance(n1.radio, n2.radio)
        c = const.SPEED_OF_LIGHT
        pi = numpy.pi
        epl02 = prop.n * linear2db(4*pi*d02*(prop.fc/c))    # expected pathloss
        epl12 = prop.n * linear2db(4*pi*d12*(prop.fc/c))
        rpl02 = prop.freespace(d02)                         # reported pathloss
        rpl12 = prop.freespace(d12)
        apl02, apl12 = None,None                            # actual pathloss
        RF0, RF1 = n0.radio.traceid, n1.radio.traceid
        rfid0, rfid1 = n0.radio.uid, n1.radio.uid
        nsent0, nsent1, ncoll = 0, 0, 0
        for e in ch.trace.events:
            if (e['event']=="RCV") and (e['obj']=="RF") and (int(e['nid'])==2):
                if ('cif-collision' in e) and (len(e['cif-collision'])>0): ncoll += 1
                if ('cif-src' in e) and (e['cif-src']==RF0):
                    apl02 = float(e['pathloss'].strip("dBm") )
                if ('cif-src' in e) and (e['cif-src']==RF1):
                    apl12 = float(e['pathloss'].strip("dBm") )
            if (e['event']=="SND") and (e['obj']=="RF"):
                if   (int(e['uid'])==rfid0): nsent0 +=1
                elif (int(e['uid'])==rfid1): nsent1 +=1
        self.assertAlmostEquals(epl02, apl02, 2)    # expected != actual
        self.assertAlmostEquals(epl12, apl12, 2)    # expected != actual
        self.assertAlmostEquals(rpl02, apl02, 2)    # reported != actual
        self.assertAlmostEquals(rpl12, apl12, 2)    # reported != actual
        self.assertAlmostEquals(nsent0+nsent1, ncoll)  # collision error!
