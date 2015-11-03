#!  /usr/bin/env python

"""
Base class for simulating motion of `Node` objects; contains `Motion` class.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-05-02 15:01:35 -0500 (Mon, 02 May 2011) $
* $LastChangedRevision: 4961 $

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

:var MOTION_VERBOSE: 
    Constant enumeration to control verbose threshold in this file.

    Setting the verbose level of an `Motion` above this threshold will cause
    the corresponding output in this file to be written (or logged).
"""
__docformat__ = "restructuredtext en"

from SimPy.Simulation import now
from wins.element import Element

import numpy
from copy import copy

class Motion(Element):
    """Base class for simulating motion of `Node` objects.

    :CVariables:
     * `position`: Property to access position information.
     * `velocity`: Property to access velocity information.
    """
    name = "motion"
    tracename = "MOTION"
    def __init__(self, **kwargs):
        """Constructor."""
        self.__tic = None
        self.__pos = None
        self.__vec = None
        Element.__init__(self, **kwargs)

    position = property(fget=lambda self: self.get_position(), \
                        fset=lambda self,p: self.set_position(p) )
    velocity = property(fget=lambda self: self.get_velocity(), \
                        fset=lambda self,v: self.set_velocity(v) )
    speed = property(fget=lambda self: self.get_speed() )

    def configure(self, pos=None, **kwargs):
        """Initialize position and velocity parameters.

        :param pos: Initial position.
        """
        if pos is None: pos = (0, 0)
        self.position = pos
        #self.log("POS", pos=pos)

    def get_velocity(self):
        """Current velocity vector."""
        return copy(self.__vec)

    def set_velocity(self, v):
        """Set velocity vector."""
        p = self.position   # update position
        self.__vec = numpy.array(v)
        return p

    def get_position(self):
        """Get updated position information."""
        t = now()
        tdiff = t - self.__tic
        if tdiff>0:
            # update position
            pdiff = self.velocity * tdiff
            self.__pos += pdiff
        self.__tic = t
        return copy(self.__pos)

    def set_position(self, p):
        """Jump to new position `p`; and resume motion from new location."""
        self.__tic = now()
        self.__pos = numpy.array(p)
        if self.__vec is None: self.__vec = self.__pos * 0

    def get_speed(self):
        """Calculate current speed."""
        return numpy.linalg.norm(self.get_velocity() )

    def move(self, v):
        """Move in a new direction with a new velocity vector."""
        return self.set_velocity(v)

    def stop(self):
        """Stop current motion."""
        v = self.position*0
        return self.set_velocity(v)
