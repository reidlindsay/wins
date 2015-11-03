#!  /usr/bin/env python

"""
Runs Queue tests.

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

class TestQueue(unittest.TestCase):

    def setUp(self):
        Trace.Global.reset()

    def test_001(self):
        initialize()
        stoptime = 10.0
        verbose = 110
        # test parameters
        priority = False
        initpackets = 5
        S = [Packet() for k in range(initpackets) ]
        if priority: S = [(p, p._id+100) for p in S]
        # create consumer-producer
        q = Queue(initialBuffered=S, priority=priority)
        p = FSM(name="producer", tracename="PROD", verbose=verbose)
        c = FSM(name="consumer", tracename="CONS", verbose=verbose)
        c.addchild('queue', q)
        p.goto(self.PRODUCE, q)
        c.goto(self.CONSUME, q, "all")
        p.start(), c.start()
        simulate(until=stoptime)
        #c.trace.output()
        # error checking
        if priority:
            expected_value = [(p.name, p._id) for p, prio in S[::-1] ]
        else:
            expected_value = [(p.name, p._id) for p in S[:5] ]
        nexpected = len(expected_value)
        result = []
        for e in c.trace.events:
            if (e['event']=="RCV"):
                result.append((e['packet'], e['pid']) )
            if len(result)==nexpected: break
        for k in range(nexpected):
            epkt, eid = expected_value[k]
            rpkt, rid = result[k]
            self.assertEqual(epkt, rpkt)
            self.assertEqual(eid, rid)

    def PRODUCE(self, fsm, queue):
        while fsm.active():
            yield hold, fsm, 1.0
            p = Raw("hellworld")
            fsm.log("snd", p)
            yield queue.insert(fsm, [p], prio=p._id )
            self.assertTrue(fsm.stored(queue), "insert() error!")
        yield fsm.stop()

    def CONSUME(self, fsm, queue, fn=1):
        while fsm.active():
            yield hold, fsm, 1.0
            yield queue.remove(fsm, fn=fn)
            self.assertTrue(fsm.acquired(queue), "remove() error!")
            if isinstance(fn, int):
                self.assertEqual(len(fsm.got), fn, "remove() error! unexpected " + \
                    "number of args dequeued %s != %s"%(len(fsm.got), fn) )
            for p in fsm.got: fsm.log("rcv", p)
        yield fsm.stop()
