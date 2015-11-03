#!  /usr/bin/env python

"""
Runs tests for `Element`.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-05-19 23:11:00 -0500 (Thu, 19 May 2011) $
* $LastChangedRevision: 4986 $

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

class TestElement(unittest.TestCase):

    def setUp(self):
        Trace.Global.reset()

    def test_001(self):
        class MyElement(Element):
            name="element"
            tracename="ELEM"
            def __init__(self, *args, **kwargs):
                Element.__init__(self, *args, **kwargs)
            def configure(self, priority=False, **kwargs):
                self.priority = priority
                self.newchild("IN", Port, tracename=self.tracename+".IN")
                self.newchild("OUT", Port, tracename=self.tracename+".OUT")
                self.newchild("fsm", FSM, tracename=self.tracename+".FSM")
                self.fsm.goto(self.SERVICE, 1.0)
            def SERVICE(self, fsm, delay):
                while fsm.active():
                    yield self.IN.recv(fsm, fn=5)
                    yield hold, fsm, delay
                    if self.priority:
                        for a,b in fsm.got: self.log("fwd", a)
                        S = [(p, p._id) for p in fsm.got]
                    else:
                        for a in fsm.got: self.log("fwd", a)
                        S = [p for p in fsm.got]
                    yield self.OUT.send(fsm, S)
                yield fsm.stop()

        class Producer(Element):
            name="producer"
            tracename="PROD"
            def __init__(self, *args, **kwargs):
                Element.__init__(self, *args, **kwargs)
            def configure(self, **kwargs):
                self.newchild("OUT", Port, tracename=self.tracename+".OUT")
                self.newchild("fsm", FSM, tracename=self.tracename+".FSM")
                self.fsm.goto(self.RUN, 1.0)
            def RUN(self, fsm, delay):
                while fsm.active():
                    yield hold, fsm, delay
                    p = Packet()/"helloworld"
                    p.setanno('parent', self)
                    S = [p]
                    for p in S: self.log("snd", p)
                    yield self.OUT.send(fsm, S)
                yield fsm.stop()

        class Consumer(Element):
            name="consumer"
            tracename="CONS"
            def __init__(self, *args, **kwargs):
                Element.__init__(self, *args, **kwargs)
            def configure(self, priority=False, **kwargs):
                self.newchild("IN", Port, tracename=self.tracename+".IN", \
                              priority=priority)
                self.newchild("fsm", FSM, tracename=self.tracename+".FSM")
                self.fsm.goto(self.SERVICE)
            def SERVICE(self, fsm):
                while not fsm.terminated():
                    yield self.IN.recv(fsm, fn="all")
                    for p in fsm.got:
                        self.log("rcv", p)
                yield fsm.stop()

        initialize()
        stoptime=15.0
        verbose=30
        priority=False
        kwargs = {'verbose':verbose, 'priority':priority}
        # create blocks
        src = Producer(**kwargs)
        el0 = MyElement(**kwargs)
        snk = Consumer(**kwargs)
        # connect blocks
        src.OUT.connect(el0.IN)
        el0.OUT.connect(snk.IN)
        simulate(until=stoptime)
        #src.trace.output()
        # check output
        nsent, nrcvd = 0, 0
        for e in src.trace.events:
            if (e['event']=="RCV") and (e['obj']=="CONS"):
                if (int(e['nid'])==0): nrcvd += 1
            if (e['event']=="SND") and (e['obj']=="PROD"):
                if (int(e['nid'])==0): nsent += 1
        self.assertEqual(nrcvd, 5*int((nsent-1)/5) )
