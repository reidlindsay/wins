#!  /usr/bin/env python

"""
Base class for implementing medium-access control (or Layer 2) protocols;
contains `MAC` class.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-09-19 21:56:57 -0500 (Mon, 19 Sep 2011) $
* $LastChangedRevision: 5134 $

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

from SimPy.Simulation import SimEvent
from wins.element import Element
from wins.protocol.phy import PHY
from wins.packet import ANNO
from wins.crc import CRC32
from wins import const

from scapy.all import Packet

class MAC(Element):
    """Base class for medium-access control (i.e. Layer 2) protocols.

    MAC Ports
    =========
    Every `MAC` object has the following two ports:

        1. "TXU" - sends received packets to upstream protocol (usually `ARP`).
        #. "RXU" - receives traffic from an upstream element (usually `ARP`).

    These ports are needed to properly interface `MAC` protocols with `ARP`.
    Additional ports should added by subclasses as needed.

    MAC and Addressing
    ==================
    The format of MAC layer addresses depends on the hardware type enumeration
    `htype`. This parameter is used by `ARP` to determine how to properly format
    messages destined for a MAC layer or how to interpret messages delivered by
    a MAC layer. By default, `MAC` uses Ethernet addressing.

    :cvar phy: Property to access pointer to `PHY`.
    :cvar address: Property to access/modify address information.
    :cvar htype: Property to access/modify hardware type enumeration.
    :cvar addr: Property to access `address`; acts as alias for address info.
    :cvar broadcast: Default value for broadcast address
                     [default=`ETHERNET_BROADCAST_ADDR`].
    :cvar nulladdr: Default value of null address [default=`ETHERNET_NULL_ADDR`]
    :cvar mtu: Default value for maximum transmission unit.

    :ivar drpdata: SimEvent that should be signalled when a data packet is
                    dropped by the MAC.
    :ivar ackdata: SimEvent that should be signalled when a data packet is
                   successfully acknowledged.
    :ivar promiscuous: Boolean flag; should be used to indicate whether or not
                       MAC is operating in promiscuous mode [default=False].
    :ivar broadcast: Broadcast address of `MAC`.
    :ivar nulladdr: Null address (used by `ARP`).
    :ivar mtu: Maximum transmission unit.
    :ivar __phy: Private pointer to `PHY` for `MAC`; use `phy` property to
                 access/modify.
    :ivar __addr: Private variable to maintain address information.
    :ivar __htype: Private variable to maintain hardware-type enumeration.

    :ivar checkcrc: Boolean flag; if true check CRC for errors before passing to
                    upper layer; otherwise ignore CRC [default=True].
    """
    name = "mac"
    tracename = "MAC"
    broadcast = const.ETHERNET_BROADCAST_ADDR
    nulladdr  = const.ETHERNET_NULL_ADDR
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
        self.__phy = None
        self.__addr = None
        self.__htype = None
        self.checkcrc = None
        self.promiscuous = None
        self.drpdata = SimEvent()
        self.ackdata = SimEvent()
        Element.__init__(self, **kwargs)
        errmsg = "[MAC]: Invalid Port configuration! " + \
                 "Must have 'RXU' and 'TXU' ports!"
        assert self.hasport("RXU") and self.hasport("TXU"), errmsg
        self.drpdata.name = "%s.%s"%(self.name, "drpdata")
        self.ackdata.name = "%s.%s"%(self.name, "ackdata")

    phy = property(fget=lambda self: self.__phy, \
                   fset=lambda self,v: self.set_phy(v) )
    address = property(fget=lambda self: self.__addr, \
                       fset=lambda self,v: self.set_address(v,None) )
    htype = property(fget=lambda self: self.__htype, \
                     fset=lambda self: self.set_address(None,v) )
    addr = property(fget=lambda self: self.address, \
                    fset=lambda self,v: setattr(self,'address',v) )

    def configure(self, phy=None, addr="auto", htype="auto", \
                        checkcrc=True, promiscuous=False, **kwargs):
        """Set up pointers, ports, and addressing.

        :param phy: Pointer to `PHY` for `MAC`.
        :param addr: Address for `MAC` (see `MAC and Addressing` for more).
        :param htype: Hardware type enumeration for `ARP`; Ethernet by default.
        :param checkcrc: Boolean flag; if false `haserror()` will ignore the CRC
                         when determining if an error has occurred.
        :param promiscuous: Boolean flag; should be used to indicate whether or
                            not MAC is operating in promiscuous mode
                            [default=False].

        By default, this method will use the automatic configuration option for
        `set_address()`. It will also set up the "TXU" and "RXU" ports.
        """
        # set parameters
        self.checkcrc = checkcrc
        self.promiscuous = promiscuous
        if phy is not None: self.phy = phy
        self.set_address(addr, htype)
        self.addport("TXU")     # port to send traffic to upstream protocol
        self.addport("RXU")     # port to receive traffic from upstream protocol

    def connect(self, p):
        """Convenience method to connect to a downstream `Element` with 'TXU'
        and 'RXU' ports; **overload this method as needed**."""
        if isinstance(p, Element):
            if p.hasport("TXU") and p.hasport("RXU"):
                self.TXD.connect(p.getport("RXU"))
                p.getport("TXU").connect(self.RXD)
                return 0
        raise RuntimeError, "[MAC]: Unable to connect() to %s! "%(p) + \
                "Could not find necessary 'TXU' and 'RXU' ports!"

    def duration(self, p, *args, **kwargs):
        """Convenience method to get duration of packet `p`.

        :return: Duration of packet `p` in seconds.

        By default, this method calls `PHY.duration()` on the physical layer
        associated with this `MAC` (i.e. `phy`). If `phy` is not a valid `PHY`,
        this method will attempt to find and return the 'cif-duration' annotation.

        Overload this method as needed.

        :note: If duration cannot be determined, this method may raise an
               exception.
        """
        if isinstance(self.phy, PHY):
            return self.phy.duration(p, *args, **kwargs)
        elif ANNO.supports(p, 'cif-duration'):
            return p.getanno('cif-duration')
        else:
            raise RuntimeError, "[MAC]: Could not determine duration! " + \
                    "No valid PHY and 'cif-duration' annotation not found!"

    def haserror(self, p, *args, **kwargs):
        """Convenience method to determine if packet `p` has an error.

        :param args: Additional arguments passed to `CRC32.haserror()`.
        :param kwargs: Keyword arguments passed to `CRC32.haserror()`.
        :return: Boolean; true if error is found; false otherwise.

        By default, this method calls `CRC32.haserror()`. If no `CRC32` layer is
        found, or no 'crcerror' annotation is found the packet is assumed to be
        error-free.
        
        **Overload this method as needed to change this operation.**

        :note: This method ignores the CRC when `checkcrc` is set to false.
        """
        hascrc = isinstance(p, Packet) and p.haslayer(CRC32)
        hasanno = ANNO.supported(p) and p.hasanno('crcerror')
        crcerror = (hascrc or hasanno) and CRC32.haserror(p,*args,**kwargs)
        if not self.checkcrc:
            crcerror = False    # assume no error when checkcrc is false
        return crcerror

    def set_address(self, addr=None, htype=None):
        """Set address information.

        :param addr: Address for `MAC`.
        :param htype: Hardware-type enumeration for `ARP`.

        If the parameter `addr` (or `htype`) is set to "auto", the address (or
        hardware type enumeration) will be determined automatically. If `None`
        is passed in, nothing is done.

        Under automatic configuration, this method uses Ethernet frames and
        automatically configures the address using `uid`.
        """
        # set hardware type
        if isinstance(htype, str) and (htype.lower() == "auto"):
            htype = const.ARP_HTYPE_ETHERNET
        if htype is not None:
            self.__htype = htype
        # set address
        if isinstance(addr, str) and (addr.lower() == "auto"):
            # automatically determine address from uid
            if (self.htype==const.ARP_HTYPE_ETHERNET):
                addr, id, alen = "", int(self.uid), 6
                for k in range(alen):
                    addr = "%02X"%(id%256) + addr
                    if (k+1)<alen: addr = ":" + addr
                    id = id /256
            else:
                raise RuntimeError, \
                        "[MAC]: Cannot automatically set_address()! " + \
                        " Unsupported hardware type (%s)!"%(self.htype)
        if addr is not None:
            self.__addr = addr

    def set_phy(self, p):
        """Set pointer to `PHY` for MAC layer."""
        assert isinstance(p, PHY) or (p is None), \
               "[MAC]: set_phy() cannot set pointer to non-PHY object!"
        self.__phy = p
