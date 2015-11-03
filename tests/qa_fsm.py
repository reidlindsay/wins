#!  /usr/bin/env python

"""
Test `FSM`.

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
from copy import copy, deepcopy

class TestFSM(unittest.TestCase):

    def setUp(self):
        Trace.Global.reset()

    def tearDown(self):
        Trace.Global.reset()

    def test_goto(self):
        """Test `FSM.goto()`."""
        class Tester(Traceable):
            name = "tester"
            tracename = "TEST"
            def __init__(self, **kwargs):
                Traceable.__init__(self, **kwargs)
                fsm = self.newchild('fsm', FSM)
                self.fsm.goto(self.IDLE)

            def IDLE(self, fsm):
                yield hold, fsm, 1.0
                yield fsm.goto(self.BUSY)

            def BUSY(self, fsm):
                yield hold, fsm, 1.0
                yield fsm.goto(self.IDLE)

        initialize()
        stoptime = 9.0
        verbose = 100
        t = Tester(verbose=verbose)
        t.fsm.start()
        simulate(until=stoptime)
        # check trace
        nidle, nbusy = 0, 0
        for e in t.trace.events:
            if ('obj' in e) and (e['obj']==t.fsm.tracename):
                if ('event' in e) and (e['event']=="IDLE"): nidle += 1
                if ('event' in e) and (e['event']=="BUSY"): nbusy += 1
        idle_expected, busy_expected = int(stoptime/2+0.5), int(stoptime/2+1.0)
        self.assertEqual(idle_expected, nidle, "FSM.goto() error!")
        self.assertEqual(busy_expected, nbusy, "FSM.goto() error!")
        #t.trace.output()

    def test_001(self):
        """Test `FSM` state transitions."""
        initialize()
        stoptime = 10.0
        verbose = 100
        f = FSM(verbose=verbose)
        f.goto(self.S0)
        f.start()
        simulate(until=stoptime)

        #f.trace.output()

    def test_002(self):
        """Test `FSM` state transitions."""
        initialize()
        stoptime = 10.0
        verbose = 100
        f = FSM(verbose=verbose)
        f.goto(self.S0)
        f.start()
        simulate(until=stoptime)

        #f.trace.output()

    def S0(self, fsm):
        yield hold, fsm, 1.0
        t = now()
        self.assertEqual(t, 1.0, "state transition error!")
        yield fsm.goto(self.S1)

    def S1(self, fsm):
        yield hold, fsm, 1.0
        t = now()
        self.assertEqual(t, 2.0, "state transition error!")
        yield fsm.goto(self.S2)

    def S2(self, fsm):
        yield fsm.stop()
        self.assertFalse(True, "sleep() error!")

    def test_timer(self):
        """Test `Timer`."""
        initialize()
        stoptime = 10.0
        verbose = 100
        duration = 5.0
        timer = Timer(duration, start=True, verbose=verbose)
        m = FSM(verbose=verbose, tracename="MON")
        m.goto(self.TMON, timer)
        m.start()
        s = FSM(verbose=verbose, tracename="STOP")
        s.goto(self.TSTOP, timer)
        s.start()
        simulate(until=stoptime)

        #timer.trace.output()

    def TMON(self, fsm, timer):
        yield hold, fsm, 1.0
        d = timer.duration
        self.assertFalse(timer.ispaused)
        self.assertAlmostEquals(timer.timepassed,1.0)   # time passed failed
        self.assertAlmostEquals(timer.timeleft, d-1.0)  # time left failed
        yield timer.pause(fsm)
        yield hold, fsm, 1.0
        tpassed, tleft = timer.timepassed, timer.timeleft
        self.assertAlmostEquals(tpassed,1.0)    # timepassed failed after pause
        self.assertAlmostEquals(tleft, d-1.0)   # timeleft failed after pause
        yield timer.resume(fsm)
        yield hold, fsm, 1.0
        tpassed, tleft = timer.timepassed, timer.timeleft
        self.assertAlmostEquals(tpassed,2.0)    # timepassed failed after resume
        self.assertAlmostEquals(tleft, d-2.0)   # timeleft failed after resume
        yield waitevent, fsm, (timer.done, timer.kill)
        tdone = timer.done in fsm.eventsFired
        tkill = timer.kill in fsm.eventsFired
        self.assertTrue(tdone or tkill)         # events not fired!
        tpassed, tleft = timer.timepassed, timer.timeleft
        if tdone:
            self.assertAlmostEquals(tpassed, d) # timepassed failed after done
            self.assertAlmostEquals(tleft,   0) # timeleft failed after done
        elif tkill:
            fsm.log("KILL", timepassed=timer.kill.signalparam)
            self.assertTrue(tpassed<d)          # stop() failed to signal kill
        yield hold, fsm, 2.0

    def TSTOP(self, fsm, timer):
        d = timer.duration
        yield hold, fsm, 5.0
        timer.halt()

