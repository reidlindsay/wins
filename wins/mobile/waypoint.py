#!  /usr/bin/env python

"""
Random Waypoint mobility model; contains `RandomWaypoint` class.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-12-18 16:54:42 -0600 (Sun, 18 Dec 2011) $
* $LastChangedRevision: 5380 $

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

from SimPy.Simulation import hold
from wins.fsm import FSM
from wins.mobile.motion import Motion
import numpy as np

RWP_VERBOSE = 40

class RandomWaypoint(Motion):
    """Random Waypoint mobility model.

    :cvar DefaultSpeed: Default speed in meters-per-second.
    :cvar DefaultPause: Default pause time in seconds.
    :cvar DefaultBorder: Default boundary for mobility model.
    """
    name = "random waypoint mobility"
    tracename = "RWP"
    DefaultSpeed = 0.0
    DefaultPause = 0.0
    DefaultBorder = (0, 500, 0, 500)
    def __init__(self, *args, **kwargs):
        """Constructor."""
        Motion.__init__(self, *args, **kwargs)

    def configure(self, pos=None, speed=None, pause=None, border=None, **kwargs):
        """Configure motion parameters.

        :param speed: 2-tuple containing (min-speed, max-speed).
        :param pause: 2-tuple containing (min-pause, max-pause).
        :param border: 4-tuple containing (x-min, x-max, y-min, ymax).

        The speed when traveling to the next waypoint is uniformly distributed
        in the specified range. Pause time at a waypoint is uniformly
        distributed in the specified range. The next waypoint is uniformly
        chosen at random from the area specified by the border.
        
        If speed (or pause) are not tuples, they are assumed to be singleton
        values and represent constant speed (or pause) values.
        """
        Motion.configure(self, pos=pos, **kwargs)
        if speed is None:  speed  = self.DefaultSpeed
        if pause is None:  pause  = self.DefaultPause
        if border is None: border = self.DefaultBorder
        # set up parameters
        try:    x = iter(speed)
        except: speed = (speed, speed)
        try:    x = iter(pause)
        except: pause = (pause, pause)
        try:    x = iter(border)
        except: border = (0, border, 0, border)
        # configure RWP
        self._speed = speed
        self._pause = pause
        self.border = border
        # initialize position
        if pos is None:
            B = border
            x = np.random.uniform(B[0],B[1])
            y = np.random.uniform(B[2],B[3])
            self.position = (x,y)
        # create FSM to run RWP
        f = self.newchild("rwp", FSM, tracename=self.tracename+".RWP")
        f.goto(self.INIT)

    def INIT(self, fsm):
        """INIT state; start moving by default."""
        yield fsm.goto(self.MOVE)

    def MOVE(self, fsm):
        """MOVE state; pick a new waypoint and start traveling there."""
        self.stop()
        S = self._speed
        p0 = np.array(self.position)
        # pick new destination
        B = self.border
        x = np.random.uniform(B[0], B[1])
        y = np.random.uniform(B[2], B[3])
        p1 = np.array([x,y])
        # pick new speed
        s = np.random.uniform(S[0], S[1])
        # calculate new velocity
        d = np.linalg.norm(p1-p0)
        if (s>0):
            t = d/s
            self.velocity = (p1-p0)/t
            self.log("MOVE", pos=list(p0), dest=list(p1), v=s, t=t)
            yield hold, fsm, t
        else:
            self.log("STOP", pos=list(p0), v=s)
            yield fsm.stop()    # stop since you will never reach there
        self.position = p1
        # goto PAUSE
        yield fsm.goto(self.PAUSE)

    def PAUSE(self, fsm):
        """PAUSE state; pick pause time and wait at waypoint."""
        self.stop()
        P = self._pause
        # pick pause time
        t = np.random.uniform(P[0], P[1])
        self.log("PAUSE", pos=list(self.position), t=t)
        # pause and then goto MOVE
        yield hold, fsm, t
        yield fsm.goto(self.MOVE)

    def log(self, event=None, p=None, *args, **kwargs):
        """Overloaded to check verbose level and set common annotations."""
        force = False
        if ('verbose' in kwargs): force = (kwargs['verbose']>RWP_VERBOSE)
        if self.verbose>RWP_VERBOSE or force:
            Motion.log(self, event, p, *args, **kwargs)
