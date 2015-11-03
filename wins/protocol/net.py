#!  /usr/bin/env python

"""
Base class for implementing network (or Layer 3) protocols; contains `NET`
class.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-09-17 12:20:43 -0500 (Sat, 17 Sep 2011) $
* $LastChangedRevision: 5128 $

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

from wins.element import Element
from wins import const

class NET(Element):
    """Base class for implementing network (i.e. Layer 3) protocols.

    This base class just provides addressing support; see `Routing` for base
    class needed to build routing protocols.

    NET Ports
    =========
    Every `NET` object has the following two ports:

        1. "TXD" - sends traffic to a downstream element (usually `ARP`).
        #. "RXD" - receives traffic from a downstream element (usually `ARP`).

    These ports are needed to properly interface `NET` protocols with `ARP`.
    Additional ports should added by subclasses as needed.

    NET and Addressing
    ==================
    The format of network layer addresses depends on the protocol type
    enumeration `ptype`. This parameter is used by `ARP` to determine how to
    properly format messages destined for a network layer or how to interpret
    messages sent by a network layer. By default, `NET` uses IPv4 addressing.

    :cvar address: Property to access/modify address information.
    :cvar ptype: Property to access/modify protocol-type enumeration.
    :cvar addr: Property to access `address`; acts as alias for address info.
    :cvar broadcast: Default value for broadcast address
                     [default=`IP_BROADCAST_ADDR`].
    :cvar nulladdr: Default value of null address [default=`IP_NULL_ADDR`]
    :cvar mtu: Default value for maximum transmission unit.

    :ivar broadcast: Broadcast address of `NET`.
    :ivar nulladdr: Null address (used by `ARP`).
    :ivar mtu: Maximum transmission unit.
    :ivar __addr: Private variable to maintain address information.
    :ivar __ptype: Private variable to maintain protocol-type enumeration.
    """
    name = "network"
    tracename = "NET"
    broadcast = const.IP_BROADCAST_ADDR
    nulladdr  = const.IP_NULL_ADDR
    mtu = None
    def __init__(self, broadcast=None, nulladdr=None, mtu=None, **kwargs):
        """Constructor.

        :param broadcast: Broadcast address.
        :param nulladdr: Null address (used by `ARP`).
        :param mtu: Maximum transmission unit.

        If the `broadcast` address is not specified, then the default broadcast
        address is used.
        """
        if broadcast is None: broadcast = self.__class__.broadcast
        if nulladdr  is None: nulladdr  = self.__class__.nulladdr
        if mtu is None: mtu = self.__class__.mtu
        self.broadcast = broadcast
        self.mtu = mtu
        self.__addr = None
        self.__ptype = None
        Element.__init__(self, **kwargs)
        assert self.hasport("RXD") and self.hasport("TXD"), \
                "[NET]: Invalid Port configuration! " + \
                "Must have 'RXD' and 'TXD' ports!"

    address = property(fget=lambda self: self.__addr, \
                       fset=lambda self,v: self.set_address(v,None) )
    ptype = property(fget=lambda self: self.__ptype, \
                     fset=lambda self: self.set_address(None,v) )
    addr = property(fget=lambda self: self.address, \
                    fset=lambda self,v: setattr(self,'address',v) )

    def configure(self, addr="auto", ptype="auto", **kwargs):
        """Set up ports and addressing.

        :param addr: Address for `NET` (see `NET and Addressing` for more).
        :param ptype: Protocol type enumeration for `ARP`; IPv4 by default.

        By default, this method will use the automatic configuration option for
        `set_address()`. It will also set up the "TXD" and "RXD" ports.
        """
        self.set_address(addr, ptype)
        self.addport("TXD")     # port to send traffic to downstream protocol
        self.addport("RXD")     # port to recv traffic from downstream protocol

    def set_address(self, addr=None, ptype=None):
        """Set address information.

        :param addr: Address for `NET`.
        :param ptype: Protocol-type enumeration for `ARP`.

        If the parameter `addr` (or `ptype`) is set to "auto", the address (or
        protocol type enumeration) will be determined automatically. If `None`
        is passed in, nothing is done.

        Under automatic configuration, this method uses IPv4 addressing and
        automatically configures the address using `uid`.
        """
        # set protocol type
        if isinstance(ptype, str) and (ptype.lower() == "auto"):
            ptype = const.ARP_PTYPE_IP
        if ptype is not None:
            self.__ptype = ptype
        # set address
        if isinstance(addr, str) and (addr.lower() == "auto"):
            # automatically determine address from uid
            if (self.ptype==const.ARP_PTYPE_IP):
                addr, id, alen = "", int(self.uid), 4
                for k in range(alen):
                    addr = "%d"%(id%256) + addr
                    if (k+1)<alen: addr = "." + addr
                    id = id /256
            else:
                raise RuntimeError, \
                        "[NET]: Cannot automatically set_address()! " + \
                        " Unsupported protocol type (%s)!"%(self.ptype)
        if addr is not None:
            self.__addr = addr

    def log(self, event=None, p=None, *args, **kwargs):
        """Overloaded to check verbose level and set common annotations."""
        kwargs.update(self.get_net_anno(p))
        Element.log(self, event, p, *args, **kwargs)

    def get_net_anno(self, *args, **kwargs):
        """Internal method to get relevant annotations."""
        netargs = {}
        netargs['addr'] = self.address
        return netargs
