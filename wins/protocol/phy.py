#!  /usr/bin/env python

"""
Base class for implementing physical layer (or Layer 1) protocols; contains
`PHY` class.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-10-04 15:05:26 -0500 (Tue, 04 Oct 2011) $
* $LastChangedRevision: 5183 $

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

from wins.channel.interface import ChannelInterface, strcollision
from wins.channel.radio     import Radio
from wins.element import Element
from wins.packet  import Packet, ANNO
from wins.crc import CRC32

class PHY(Element):
    """Base class for physical layer (i.e. Layer 1) protocols.

    The functions of this layer are to:

        * Encode packets into waveforms.
        * Mark outgoing packets with proper annotations (e.g. 'cif-duration',
          'cif-txpower', etc.)
        * Detect the arrival of packets.
        * Decode received waveforms to recover packet.
        * Mark incoming packets with error annotations and others as needed
          (e.g. 'sinr', 'ber', 'per', etc.)

    This base class only implements helper functions and convenience methods;
    subclasses derived from `PHY` should implement the above functionality.
    Subclasses should also be responsible for simulating appropriate timing.

    :cvar cif:   Property to access pointer to `ChannelInteface`.
    :cvar radio: Property to access `cif`; acts as alias for channel interface.
    :cvar bandwidth: Property to access `Radio.bandwidth` for `radio`.
    :cvar mtu: Default value for maximum transmission unit.
    :cvar rate: Property to access list of valid rate enumerations.

    :ivar mtu: Maximum transmission unit.
    :ivar __cif: Private pointer to `ChannelInterface` for `PHY`; use `cif`
                 property or aliases to access/modify.

    :note: By default, `mtu` is None. This implies that no bound exists.
    """
    name = "phy"
    tracename = "PHY"
    mtu = None
    def __init__(self, mtu=None, **kwargs):
        """Constructor.

        :param mtu: Maximum transmission unit.
        """
        if mtu is None: mtu = self.__class__.mtu
        self.mtu = mtu
        self.__cif = None
        Element.__init__(self, **kwargs)

    cif = property(fget=lambda self: self.__cif, \
                   fset=lambda self,v: self.set_channelif(v) )
    radio = property(fget=lambda self: self.cif, \
                     fset=lambda self,r: self.set_radio(r) )
    rate = property(fget=lambda self: self.get_valid_rates())
    bandwidth = property(fget=lambda self: self.get_bandwidth() )

    def configure(self, cif=None, radio=None, **kwargs):
        """Set up pointer to channel interface."""
        if cif is not None: self.cif = cif
        if radio is not None: self.radio = radio

    def duration(self, p):
        """Calculate duration of packet `p`.

        :return: Duration of packet `p` in seconds.
        
        By default this method returns zero.
        """
        return 0

    def seterror(self, p):
        """Convenience method to set error annotation (or parameters) in `p`.

        :return: Modified packet with updated parameters/annotations.

        By default, if packet `p` has a `CRC32` layer, this method will set the
        'crcerror' field, otherwise this method will set the 'crcerror'
        annotation to 1.

        **Overload this method as needed.**
        """
        hascrc = CRC32.supported(p)
        hasanno = ANNO.supported(p)
        if hascrc:
            p[CRC32].crcerror = 1
        elif hasanno:
            p.setanno('crcerror', 1)
        else:
            raise RuntimeError, "[PHY]: seterror() failed to find packet!"

    def txproctime(self, p=None):
        """Get transmit processing delay for packet `p` from annotations and
        local settings; **overload as needed.**.

        :return: Additional processing time needed to encode the packet, i.e. in
                 excess of transmit duration of packet.

        By default, this method returns zero.
        """
        return 0

    def rxproctime(self, p=None):
        """Get receive processing delay for packet `p` from annotations and
        local settings; **overload as needed.**.

        :return: Additional processing time needed to decode the packet, i.e. in
                 excess of receive duration of packet.

        By default, this method returns zero.
        """
        return 0

    def set_channelif(self, c):
        """Set pointer to `ChannelInterface` for physical layer."""
        assert isinstance(c, ChannelInterface) or (c is None), \
               "[PHY]: set_channelif() cannot set non-ChannelInterface object!"
        self.__cif = c
        return c

    def set_radio(self, r):
        """Set pointer to `Radio` for physical layer; wrapper for
        `set_channelif()` method."""
        assert isinstance(r, Radio) or (r is None), \
               "[PHY]: set_radio() cannot set non-Radio object!"
        return self.set_channelif(r)

    def get_bandwidth(self):
        """Get bandwidth attribute for `radio`.

        :returns: `Radio.bandwidth` in Hz (or None).

        If no `Radio` is associated with the PHY, this method returns None.
        """
        r, bw = self.radio, None
        if isinstance(r, Radio): bw = r.bandwidth
        return bw

    def interval_heap(self, *args, **kwargs):
        """Calls `ChannelIF.interval_heap()` on associated channel interface `cif`."""
        return self.cif.interval_heap(*args, **kwargs)

    def sinr_heap(self, *args, **kwargs):
        """Calls `ChannelIF.sinr_heap()` on associated channel interface `cif`."""
        return self.cif.sinr_heap(*args, **kwargs)

    def get_valid_rates(self):
        """Get list of rate enumerations supported by the PHY.

        :note: This abstract method *must* be defined in subclasses.
        """
        errmsg = "[PHY]: Abstract class PHY does not define get_valid_rates()!"
        raise RuntimeError, errmsg

    def get_datarate(self, r):
        """Get the data rate in bits-per-second (bps).

        :param r: Rate enumeration.

        :note: This abstract method *must* be defined in subclasses.
        """
        errmsg = "[PHY]: Abstract class PHY does not define get_datarate()!"
        raise RuntimeError, errmsg

    def calclength(self, duration, rate):
        """Calculate length of packet in bytes.

        :param duration: Duration of packet in seconds.
        :param rate: Rate enumeration.

        :note: This abstract method *must* be defined in subclasses.
        """
        errmsg = "[PHY]: Abstract class PHY does not define calclength()!"
        raise RuntimeError, errmsg

    def cleananno(self, p):
        """Clean up annotations before passing to upper layers."""
        assert ANNO.supported(p), "[PHY]: ANNO not supported!"
        # replace/remove unwanted annotations
        if p.hasanno('cif-collision'):
            coll = strcollision(p)
            assert (coll is not None)
            p.delanno('cif-collision')
            p.setanno('cif-collision', coll, priv=True)
        return p
