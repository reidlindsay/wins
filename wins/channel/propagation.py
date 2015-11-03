#!  /usr/bin/env python

"""
Basic propagation model for wireless channel; contains `Propagation` class.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-05-19 23:09:01 -0500 (Thu, 19 May 2011) $
* $LastChangedRevision: 4985 $

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

from wins.channel.channelbase import ChannelModel
from wins.channel.interface   import ChannelInterface
from wins.mobile.motion       import Motion
from wins.helper import linear2db, db2linear
from wins.base   import Reference
from wins import const

import numpy as np

class Propagation(ChannelModel):
    """Basic propagation model for wireless channels.

    By default, this class implements the Freespace pathloss model using unit
    antenna gain and no system loss. Overload the `apply()` method to customize
    a `Propagation` subclass. The following parameters are defined in this class
    and may be redefined in dervied classes:

        * `n`:  Pathloss exponent.
        * `fc`: Center frequency for channel.

    Packet Annotations and Propagation
    ==================================
    By default, this class will set the annotation listed below. **Derived
    classes must set 'pathloss', 'cm-delay', and any other needed annotations.**

    ========== ==============================================================
    Name        Description
    ========== ==============================================================
    pathloss    Pathloss between transmitter and receiver.
    ---------- --------------------------------------------------------------
    doppler     Doppler shift from perspective of receiver (in Hz).
    ---------- --------------------------------------------------------------
    cm-delay    Propagation delay annotation marked by `ChannelModel`.
    ========== ==============================================================
    
    Pathloss is defined as a ratio of transmit power to receive power (i.e.
    pathloss is greater than or equal to 1). Thus received power is equal to
    transmit power divided by pathloss, or in decibels:

        Pr (db) = Pt (db) - PL (db).

    See `Link Budget`_ and `Friis transmission equation`_ for more about how
    pathloss is computed.

    Propagation delay is a function of separation `distance()` and the speed of
    light (`const.SPEED_OF_LIGHT`):

        delay = (separation distance)/(speed of light) (meters/sec)

    .. _`Link Budget`: http://en.wikipedia.org/wiki/Link_budget
    .. _`Friis Transmission Equation`: http://en.wikipedia.org/wiki/Friis_transmission_equation#Modifications_to_the_basic_equation

    :CVariables:
     * `n`:  Pathloss exponent.
     * `fc`: Center frequency for channel.
    """
    n = 2
    fc = 2.4e9
    name = "propagation"
    tracename = "PROP"
    def __init__(self, fc=None, n=None, **kwargs):
        """Constructor.

        :param fc: Center frequency of channel [default=2.4e9]
        :param n:  Pathloss exponent.
        """
        # get parameters
        cls = self.__class__
        if fc is None: fc = cls.fc
        if n is None:  n  = cls.n
        # set parameters
        self.fc = fc
        self.n = n
        ChannelModel.__init__(self, **kwargs)

    def apply(self, p, u, v):
        """Apply freespace channel model.

        :param p: Packet to modify.
        :param u: Transmitting `ChannelInterface`.
        :param v: Receiving `ChannelInterface`.
        :return: Modified packet.

        This method set the 'pathloss' and 'cm-delay' annotations. Overload
        this method to change how the channel model applies channel annotations
        (i.e. to implement a new channel model).

        This method also computes the doppler shift between the transmitter and
        receiver and sets the 'doppler' annotation.
        """
        dist = self.distance(u,v)
        pl = self.pathloss(dist, fc=self.fc, n=self.n)
        p.setanno('pathloss', pl)
        delay = dist/const.SPEED_OF_LIGHT
        p.setanno('cm-delay', delay)
        fd = self.doppler(u,v, fc=self.fc)
        p.setanno('doppler', fd)
        return p

    def pathloss(self, dist, fc=None, n=None):
        """Calculate pathloss over specified distance using local (or specified)
        parameters; **overload this method as needed**.

        :param dist: Separation distance (in meters).
        :param fc: Center frequency of channel.
        :param n:  Pathloss exponent.
        :return: Pathloss (in dB).

        By default this method returns the `freespace` pathloss. Overload this
        method to change how pathloss is computed.
        """
        if fc is None: fc = self.fc
        if n is None:  n = self.n
        return self.freespace(dist, fc=fc, n=n)

    @classmethod
    def freespace(cls, dist, fc=None, n=None):
        """Calculate freespace pathloss over distance `dist`.

        :param dist: Separation distance between transmitter and receiver.
        :param fc: Center frequency of channel.
        :param n: Pathloss exponent.
        :return: Freespace pathloss in dB.

        If `fc` or `n` are not provided, this method will simply use the default
        values associated with the class. Freespace pathloss is given by:

            PL (db) = 10n*log10(4*pi*d/wavelength),

        where wavelength = (fc/c) and d = separation distance.

        :note: This method assumes unit antenna gain and no system loss. This
               method also uses `const.SPEED_OF_LIGHT` as the speed of light
               parameter 'c' in meters/second.
        """
        assert (dist>0), "[PROPAGATION]: freespace() pathloss can only " + \
                "be computed for positive separation distance!"
        # set parameters
        if fc is None: fc = cls.fc
        if n is None:  n = cls.n
        c = const.SPEED_OF_LIGHT
        # calculate pathloss in dB
        wavelength = c/fc
        Lfs = n*linear2db(4*np.pi*dist/wavelength)
        return Lfs

    @classmethod
    def distance(cls, u, v):
        """Calculate the distance between `u` and `v`.

        :param u: `ChannelInterface` or location of transmitter.
        :param v: `ChannelInterface` or location of receiver.

        :note: This method calls `getposition()` to determine the position of
               `u` and `v`.
        """
        # extract location information
        if isinstance(u, Reference): u = u._deref
        if isinstance(v, Reference): v = v._deref
        x = cls.position(u)
        y = cls.position(v)
        # calculate distance
        dist = np.linalg.norm(x-y)
        return dist

    @staticmethod
    def position(c):
        """Determine position of a `ChannelInterface`.

        :param c: `ChannelInterface` or position information.
        :return: Location as a `numpy.array`.

        If `c` is a `ChannelInterface` object, this method will attempt to get
        the location of the interface by searching for a `Motion` object
        registered with its `container` under the nickname 'motion'. If no such
        child is found in the container, this method throws an AssertionError.

        If `c` is already position information, this method will convert it to a
        `numpy.array`. Overload this method to change how position is found.

        :note: This method throws an AssertionError if it cannot find a valid
               'motion' in the `ChannelInterface`.
        """
        if isinstance(c, Reference): c = c._deref
        if isinstance(c, ChannelInterface):
            m = c.container.getchild('motion')
            assert isinstance(m, Motion), "[PROPAGATION]: getposition() " + \
                    "could not find `Motion` element in container of %s!"%(c)
            c = m.position
        pos = np.array(c)
        return pos

    @staticmethod
    def velocity(c):
        """Determine velocity for `ChannelInterface` `c`.

        :param c: `ChannelInterface` or velocity vector.
        :return: Velocity as a `numpy.array`.
        """
        if isinstance(c, Reference): c = c._deref
        if isinstance(c, ChannelInterface):
            m = c.container.getchild('motion')
            assert isinstance(m, Motion), "[PROPAGATION]: getvelocity() " + \
                    "could not find `Motion` element in container of %s!"%(c)
            c = m.velocity
        v = np.array(c)
        return v

    @classmethod
    def doppler(cls, u, v, fc=None):
        """Calculate doppler shift between `u` and `v`.

        :param u: `ChannelInterface` of transmitter.
        :param v: `ChannelInterface` of receiver.
        :param fc: Center frequency of channel.
        :return: Doppler shift (w.r.t. receiver) in Hz.

        :note: This method calls `getvelocity()` to determine the position of
               `u` and `v`.
        """
        # get parameters
        if fc is None: fc = cls.fc
        # extract location information
        if isinstance(u, Reference): u = u._deref
        if isinstance(v, Reference): v = v._deref
        # calculate relative velocity and position (w.r.t. receiver)
        w = cls.velocity(u) - cls.velocity(v)
        x = cls.position(u) - cls.position(v)
        # calculate relative speed of transmitter (w.r.t. receiver)
        s = np.linalg.norm(w)
        c = const.SPEED_OF_LIGHT
        assert (s<c), "[PROPAGATION]: Cannot travel faster than speed of light!"
        if (np.dot(w,x)>0):
            sr = s/c    # transmitter moving away from receiver
        elif (np.dot(w,x)<0):
            sr = -s/c   # transmitter moving towards receiver
        else:
            sr = 0      # tangential motion
        # calculate doppler shift
        fd = fc * (1.0/(1-sr) - 1.0)
        return fd

    def log_forward(self, p, *args, **kwargs):
        """Overload to log 'pathloss' annotation."""
        if p.hasanno('pathloss'):
            pathloss = p.getanno('pathloss')
            kwargs['pathloss'] = "%.2f dB"%(pathloss)
        ChannelModel.log_forward(self, p, *args, **kwargs)

