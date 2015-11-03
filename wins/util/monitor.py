#!  /usr/bin/env python

"""
Monitor progress by printing status messages at regular intervals.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-06-26 22:51:38 -0500 (Sun, 26 Jun 2011) $
* $LastChangedRevision: 5032 $

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

from SimPy.Simulation import hold, now
from wins.fsm import FSM
import sys

class Monitor(FSM):
    """Monitor progress of simulation.

    Prints status at regular intervals.

    :ivar period: Duration between each output update.
    :ivar width: Number of updates per line.

    Running Monitor
    ===============
    Start a monitor like you would a normal FSM:

        m = Monitor()
        m.start()
    """
    def __init__(self, period=0, width=25, **kwargs):
        """Constructor.

        :param period: Duration between each output update.
        :param width: Number of updates per line.
        """
        self.period = period
        self.width = width
        self._tic = 0
        kwargs['initstate'] = self.MON
        FSM.__init__(self, **kwargs)

    def MON(self, fsm):
        """MON state; run monitor functionality."""
        if not (self.period>0): fsm.stop()
        self.tic = (self.tic + 1)%(self.width)
        yield hold, fsm, self.period
        sys.stdout.write(".")
        if (self.tic == 0):
            sys.stdout.write(" t = %.5f\n"%(now()) )
        sys.stdout.flush()
        yield fsm.goto(self.MON)
