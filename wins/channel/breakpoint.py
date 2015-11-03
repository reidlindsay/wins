#!  /usr/bin/env python

"""
Breakpoint channel model uses freespace propagation until specified break
point distance; contains `BreakPoint` classes.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-06-27 14:46:49 -0500 (Mon, 27 Jun 2011) $
* $LastChangedRevision: 5039 $

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

from wins.channel.propagation import Propagation
from wins.helper import linear2db

class Breakpoint(Propagation):
    """Breakpoint pathloss model uses threshold policy for propagation model.

    This channel model uses freespace pathloss up to a specified breakpoint
    distance. Beyond that distance (`bpdist`), the model uses a higher pathloss
    exponent (`n`):

        PL (db) = PL0 + 10*n*log10(d/d0), (d>d0)

    where PL0 is the pathloss at the breakpoint distance d0 (i.e. `bpdist`).

    :CVariables:
     * `bpdist`: Breakpoint distance (in meters) [default = 5.0].
    """
    name = "breakpoint propagation"
    tracename = "BP"
    bpdist = 5.0
    def __init__(self, bpdist=None, n=3.5, **kwargs):
        """Constructor.

        :param bpdist: Breakpoint distance (in meters) [default = 5.0].
        :param n: Pathloss exponent [default = 3.5].
        :param kwargs: Additional keyword arguments passed to `Propagation`
                       constructor.
        """
        cls = self.__class__
        if bpdist is None: bpdist = cls.bpdist
        self.bpdist = bpdist
        Propagation.__init__(self, n=n, **kwargs)

    def apply(self, p, u, v):
        """Apply breakpoint pathloss model.

        :param p: Packet to modify.
        :param u: Transmitting `ChannelInterface`.
        :param v: Receiving `ChannelInterface`.
        :return: Modified packet.

        This method resets the 'pathloss' annotations after calling
        `Propagation.apply()`.
        """
        pkt = Propagation.apply(self, p, u, v)
        # overwrite pathloss
        dist = self.distance(u,v)
        pl = self.pathloss(dist, bpdist=self.bpdist, fc=self.fc, n=self.n)
        pkt.setanno('pathloss', pl)
        return pkt

    def pathloss(self, dist, **kwargs):
        """Overloaded to implement breakpoint pathloss model; calls
        `calcpathloss()` using local parameters.

        :param dist: Separation distance between transmitter and receiver.
        :param kwargs: Additional keyword arguments passed to `calcpathloss()`.
        :return: Pathloss (in dB).
        """
        for s in ['bpdist', 'fc', 'n']:
            exec("found = (\'%s\' in kwargs)"%(s) )
            if not found: exec("kwargs[\'%s\'] = self.%s"%(s,s) )
        return self.calcpathloss(dist, **kwargs)

    @classmethod
    def calcpathloss(cls, dist, bpdist=None, n=None, **kwargs):
        """Calculate propagation loss using breakpoint pathloss model.

        :param dist: Separation distance between transmitter and receiver.
        :param bpdist: Breakpoint distance (in meters).
        :param n: Pathloss exponent.
        :param kwargs: Additional parameters passed to `freespace()`.
        :return: Pathloss (in dB).

        The breakpoint pathloss model uses freespace propagation up to the
        breakpoint distance and the specified pathloss thereafter.
        """
        if n is None: n = cls.n
        if ((bpdist>0) and (dist>bpdist)):
            bploss = cls.freespace(bpdist, n=2.0, **kwargs)
            PL = bploss + n*linear2db(dist/(1.0*bpdist) )
        else:
            PL = cls.freespace(dist, n=2.0, **kwargs)
        return PL

    def __str__(self):
        return Propagation.__str__(self) + \
               "{bpdist: %.2f, n: %.1f, fc: %.1f GHz}"%(self.bpdist, self.n, self.fc*1e-9)
