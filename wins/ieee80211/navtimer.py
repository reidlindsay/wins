#!  /usr/bin/env python

"""
Timer for Network allocation vector (NAV), used with virtual carrier sense;
contains `NAVTimer` class.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-09-28 01:49:15 -0500 (Wed, 28 Sep 2011) $
* $LastChangedRevision: 5168 $

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

:var NAV_VERBOSE: Verbose threshold for `NAVTimer`.
"""
__docformat__ = "restructuredtext en"

from SimPy.Simulation import SimEvent, hold, waitevent
from wins.element import Element
from wins.fsm import FSM, Timer
from wins.helper import time2usec

NAV_VERBOSE = 55

class NAVTimer(Element):
    """Provides API to maintain network allocation vector and associated virtual
    carrier sense timer.

    :CVariables:
     * `timer`: Internal `Timer` object.
     * `done`: SimPy Event signalled when NAV timer successfully expires.
     * `kill`: SimPy Event signalled when NAV timer stops/pauses unexpectedly.
     * `fired`: Check if NAV timer (i.e. `timer`) has fired.
     * `running`: Check if NAV timer is still running.
     * `stopped`: Check if NAV timer has been stopped.
    """
    name = "nav timer"
    tracename = "NAV"
    def __init__(self, **kwargs):
        """Constructor."""
        self.__timer = None
        Element.__init__(self, **kwargs)

    timer = property(fget=lambda self: self.__timer)
    done = property(fget=lambda self: self.timer.done)
    kill = property(fget=lambda self: self.timer.kill)
    fired = property(fget=lambda self: \
                     isinstance(self.timer,Timer) and self.timer.fired)
    running = property(fget=lambda self: \
                       isinstance(self.timer,Timer) and self.timer.running)
    stopped = property(fget=lambda self: \
                       isinstance(self.timer,Timer) and self.timer.stopped)

    def configure(self, **kwargs):
        """Initialize `timer` to None and create monitor."""
        self.__timer = None
        self.__wakeup = SimEvent(name=self.name+".wake")
        mon = self.newchild("mon", FSM, tracename=self.tracename+".MON")
        #mon.goto(self.MON)
        mon.goto(self.MONIDLE)

    def update(self, proc, t=None):
        """Blocking call to cancel active timer, if needed, and starts a new
        timer.

        :param t: New timer value.
        :return: Blocking clause.

        If `t` is None, this method will just cancel any active timer and reset
        `timer` to None.
        """
        # cancel active timer
        if isinstance(self.timer, Timer):
            self.timer.halt()
            #self.delchild("navtimer")
            self.__timer = None
        # start new timer
        if (t>0):
            self.__timer = self.newchild("navtimer", Timer, t, start=True, \
                                         tracename=self.tracename+".TIMER")
            self.__wakeup.signal(self.timer)
        return hold, proc, 0

    def reset(self):
        """Non-blocking call to cancel active timer."""
        return self.update(None, None)

    def MONIDLE(self, fsm):
        """MONIDLE state; monitor `timer` when NAV Timer is in IDLE state."""
        # wait for timer to be set
        self.log_idle()
        yield waitevent, fsm, self.__wakeup
        t = self.__wakeup.signalparam
        assert isinstance(t, Timer)
        yield fsm.goto(self.MONBUSY, t)

    def MONBUSY(self, fsm, t):
        """MONIDLE state; monitor `timer` when NAV Timer is in BUSY state."""
        self.log_busy(duration=time2usec(t.duration) )
        yield waitevent, fsm, (t.done, t.kill)
        # timer fired
        if (t.done in fsm.eventsFired):
            pass
        # timer stopped
        elif (t.kill in fsm.eventsFired) and t.stopped:
            pass
        # otherwise -> raise exception
        else:
            raise RuntimeError, "[NAVTIMER]: Monitor indicates " + \
                                "NAV timer paused unexpectedly!"
        # done monitoring timer
        yield fsm.goto(self.MONIDLE)

    def MON(self, fsm):
        """MON state; monitor `timer` events."""
        while fsm.active():
            # wait for timer to be set
            self.log_idle()
            yield waitevent, fsm, self.__wakeup
            t = self.__wakeup.signalparam
            # monitor timer events
            while isinstance(t, Timer):
                self.log_busy(duration=time2usec(t.duration) )
                yield waitevent, fsm, (t.done, t.kill)
                # timer fired
                if (t.done in fsm.eventsFired):
                    pass
                # timer stopped
                elif (t.kill in fsm.eventsFired) and t.stopped:
                    pass
                # otherwise -> raise exception
                else:
                    raise RuntimeError, "[NAVTIMER]: Monitor indicates " + \
                            "NAV timer paused unexpectedly!"
                # done monitoring timer
                t = None
        return

    def log_busy(self, **kwargs):
        """Log NAV busy."""
        self.log("NAVBUSY", **kwargs)

    def log_idle(self, **kwargs):
        """Log NAV idle."""
        self.log("NAVIDLE", **kwargs)

    def log(self, evt=None, p=None, *args, **kwargs):
        """Overloaded to check verbose level and set common annotations."""
        force = False
        if ('verbose' in kwargs): force = (kwargs['verbose']>NAV_VERBOSE)
        if self.verbose>NAV_VERBOSE or force:
            Element.log(self, evt, p, *args, **kwargs)
