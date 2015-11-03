#!  /usr/bin/env python

"""
Reference model for propagation over wireless channel; contains
`ReferencePropagation` class.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2010-10-19 16:34:02 -0500 (Tue, 19 Oct 2010) $
* $LastChangedRevision: 4811 $

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

class ReferencePropagation(Propagation):
    """Reference model for propagation over wireless channel.

    This channel model uses pathloss at reference distance to bias the freespace
    propagation model. The pathloss at separation distance d is:

        PL (db) = PL0 + 10*n*log10(d/d0),

    where PL0 is the pathloss at reference distance d0.

    :CVariables:
     * `refdist`: Reference distance (in meters) [default = 1.0].
     * `refloss`: Pathloss at reference distance (in dB) [default = None].
    """
    name = "reference propagation"
    tracename = "REFPROP"
    refdist = 1.0
    refloss = None
    def __init__(self, refdist=None, refloss=None, **kwargs):
        """Constructor.

        :param refdist: Reference distance (in meters) [default = 1.0].
        :param refloss: Pathloss at reference distance (in dB) [default = None].

        If `refloss` is not specified, the default value of `refloss` is set to
        the freespace pathloss at `refdist` with a pathloss exponent of 2.
        """
        cls = self.__class__
        if refdist is None: refdist = cls.refdist
        if refloss is None: refloss = cls.refloss
        self.refdist = refdist
        self.refloss = refloss
        Propagation.__init__(self, **kwargs)

    def apply(self, p, u, v):
        """Apply reference pathloss model.

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
        pl = self.pathloss(dist, fc=self.fc, n=self.n)
        pkt.setanno('pathloss', pl)
        return p

    def pathloss(self, dist, **kwargs):
        """Overloaded to implement reference pathloss model; calls
        `calcpathloss()` using local parameters.

        :param dist: Separation distance between transmitter and receiver.
        :param kwargs: Additional keyword arguments passed to `calcpathloss()`.
        :return: Pathloss (in dB).
        """
        for s in ['refdist', 'refloss', 'fc', 'n']:
            exec("found = (\'%s\' in kwargs)"%(s) )
            if not found: exec("kwargs[\'%s\'] = self.%s"%(s,s) )
        return self.calcpathloss(dist, **kwargs)

    @classmethod
    def calcpathloss(cls, dist, refdist=None, refloss=None, n=None, **kwargs):
        """Calculate propagation loss using reference pathloss model.

        :param dist: Separation distance between transmitter and receiver.
        :param refdist: Reference distance (in meters).
        :param refloss: Pathloss at reference distance (in dB).
        :param n: Pathloss exponent.
        :return: Pathloss (in dB).

        If reference distance `refdist` is None, and no suitable class variable
        is defined, this method will simply return the `freespace` pathloss. If
        `refdist` is defined and `refloss` is None, this method will compute the
        reference pathloss as the freespace pathloss using loss exponent 2.
        """
        if n is None: n = cls.n
        if refdist>0:
            if refloss is None: refloss = cls.freespace(refdist, n=2, **kwargs)
            PL = refloss + n*linear2db(dist/(1.0*refdist) )
        else:
            PL = cls.freespace(dist, n=n, **kwargs)
        return PL

    def __str__(self):
        return Propagation.__str__(self) + \
               "{n: %.1f, fc: %.1f GHz}"%(self.n, self.fc*1e-9)
