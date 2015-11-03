#!  /usr/bin/env python

"""
Runs tests for mobility.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-12-17 11:48:37 -0600 (Sat, 17 Dec 2011) $
* $LastChangedRevision: 5379 $

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

import numpy

class TestMobile(unittest.TestCase):

    def setUp(self):
        Trace.Global.reset()

    def test_motion(self):
        class Tester(Motion):
            def __init__(self, **kwargs):
                Motion.__init__(self, **kwargs)
            def configure(self, **kwargs):
                Motion.configure(self, **kwargs)
                self.newchild('fsm', FSM, **kwargs)

        initialize()
        stoptime = 10.0
        verbose = 100
        kwargs = {'verbose':verbose}
        m = Motion(**kwargs)
        f = m.newchild('fsm', FSM, tracename=m.tracename+".FSM")
        f.goto(self.M0, m)
        f.start()
        simulate(until=stoptime)
        #m.trace.output()

    def M0(self, fsm, mobi):
        p, s = mobi.position, mobi.speed
        mobi.log("POS", position=p, speed=s)
        r = numpy.array((0,0) )
        self.assertTrue(all(p==r), "set_position() error!, %s != %s"%(p,r))
        mobi.move((1,0) )
        yield hold, fsm, 1.0
        yield fsm.goto(self.M1, mobi)

    def M1(self, fsm, mobi):
        p, s = mobi.position, mobi.speed
        mobi.log("POS", position=p, speed=s)
        r = numpy.array((1,0) )
        self.assertTrue(all(p==r), "move() error!, %s != %s"%(p,r))
        self.assertEqual(s, 1.0, "move() error!")
        yield hold, fsm, 1.0
        p, s = mobi.position, mobi.speed
        mobi.log("POS", position=p, speed=s)
        mobi.stop()
        yield hold, fsm, 1.0
        yield fsm.goto(self.M2, mobi)

    def M2(self, fsm, mobi):
        p, s = mobi.position, mobi.speed
        mobi.log("POS", position=p, speed=s)
        r = numpy.array((2,0) )
        self.assertTrue(all(p==r), "stop() error!, %s != %s"%(p,r))
        self.assertEqual(s, 0, "stop() error!")
        yield fsm.stop()

    def test_rwp(self):
        class Tester(Element):
            tracename = "TEST"
            def __init__(self, **kwargs):
                Element.__init__(self, **kwargs)
            def configure(self, mobile=None, **kwargs):
                self.mobile = mobile
                f = self.newchild('run', FSM)
                f.goto(self.RUN)
            def RUN(self, fsm):
                if self.mobile is None: yield fsm.stop()
                yield hold, fsm, 1.0
                self.log("POS", pos=["%.2f"%(x) for x in self.mobile.position])
                yield fsm.goto(self.RUN)

        initialize()
        stoptime = 20.0
        verbose = 50
        border = [0, 10.0, 0, 10.0]
        pause = 1.0
        pos = None  #[5.0, 5.0]
        kwargs = {'verbose':verbose, 'speed':1.0, 'border':border, \
                  'pause':pause, 'pos':pos}
        m = RandomWaypoint(**kwargs)
        t = Tester(mobile=m)
        simulate(until=stoptime)
        #m.trace.output()
