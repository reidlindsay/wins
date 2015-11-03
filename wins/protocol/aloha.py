#!  /usr/bin/env python

"""
Implementation of pure Aloha MAC protocol; contains `Aloha` class.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-09-27 22:15:57 -0500 (Tue, 27 Sep 2011) $
* $LastChangedRevision: 5167 $

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

:var ALOHA_VERBOSE:
    Constant enumeration to control verbose thresholds in this file.
  
    Setting the verbose level of a `Aloha` above this threshold will cause the
    corresponding output in this file to be written (or logged).
"""
__docformat__ = "restructuredtext en"

from SimPy.Simulation import hold
from wins.protocol.mac import MAC
from wins.element import Element
from wins.packet  import Packet, ANNO
from wins.helper  import time2usec
from wins.crc import crcupdate, crcremove
from wins.fsm import FSM

from wins import const

from scapy.all import Ether

ALOHA_VERBOSE = 54

class Aloha(MAC):
    """Implementation of pure Aloha MAC protocol.

    The pure Aloha accepts Ethernet packets (or another packet format as
    specified by the `htype` enumeration) from upstream protocols and sends the
    packet out. This upstream protocol is usually (but not necessarily `ARP`).
    This implementation of Aloha does not utilize acknowledgements.

    Aloha and Ports
    ===============
    In addition to the ports created by `MAC`, this class has two ports:

        1. "TXD" - sends traffic to a downstream element/protocol.
        #. "RXD" - receives traffic from a downstream element/protocol.

    Overload `configure()` to change how ports are set up.

    :note: This class does not support carrier sense functionality.
    """
    name = "aloha"
    tracename = "ALOHA"
    def __init__(self, **kwargs):
        """Constructor."""
        MAC.__init__(self, **kwargs)

    def configure(self, **kwargs):
        """Configure `MAC`, add ports, and create `FSM`."""
        MAC.configure(self, **kwargs)
        # add downstream ports
        self.addport("TXD")     # port to send traffic to downstream protocol
        self.addport("RXD")     # port to recv traffic from downstream protocol
        # create FSM to manage send/recv execution of mac
        txfsm = self.newchild("txfsm", FSM, tracename=self.tracename+".TX")
        rxfsm = self.newchild("rxfsm", FSM, tracename=self.tracename+".RX")
        txfsm.goto(self.SEND)
        rxfsm.goto(self.RECV)

    def connect(self, p):
        """Connect to an `Element` with a 'TXU' and 'RXU' port."""
        if isinstance(p, Element):
            if p.hasport("TXU") and p.hasport("RXU"):
                self.TXD.connect(p.getport("RXU"))
                p.getport("TXU").connect(self.RXD)
                return 0
        raise RuntimeError, "[ALOHA]: Unable to connect() to %s! "%(p) + \
                "Could not find necessary 'TXU' and 'RXU' ports!"

    def SEND(self, fsm):
        """Manage downstream traffic received from upper layer.

        Based on the value of `htype`, this method will transition to the
        appropriate message handler state (e.g. `ETHSEND` for
        `const.ARP_HTYPE_ETHERNET`).
        """
        while fsm.active():
            yield self.RXU.recv(fsm, 1)
            assert fsm.acquired(self.RXU) and (len(fsm.got)==1), \
                    "[ALOHA]: Error receiving from 'RXU' port!"
            p = fsm.got[0]
            if (self.htype==const.ARP_HTYPE_ETHERNET):
                yield fsm.goto(self.ETHSEND, p)
            else:
                raise RuntimeError, \
                        "[ALOHA]: Cannot handle incoming packet!" + \
                        " Unsupported hardware type (%s)!"%(self.htype)
        return

    def ETHSEND(self, fsm, p):
        """ETHSEND state; send ethernet frame.

        :param p: Ethernet packet to transmit.

        By default, this state appends a crc to packet `p` and sends it to
        `Port` 'TXD' without any address checking. After simulating the duration
        of the packet as reported by `duration()`, this state execution method
        returns to `SEND`.
        """
        assert (self.htype==const.ARP_HTYPE_ETHERNET), "[ALOHA]: In ETHSEND," + \
                " non-Ethernet hardware type (%s) not allowed!"%(self.htype)
        assert isinstance(p, Packet) and p.haslayer(Ether), \
                "[ALOHA]: ETHSEND cannot handle non-Ether packet!"
        eth = p[Ether]
        addr, src, dst, etype = self.address, eth.src, eth.dst, eth.type
        pkt = crcupdate(eth)
        duration = self.duration(pkt)
        self.log_send(pkt, addr=addr, src=src, dst=dst, type=etype, \
                      duration=time2usec(duration) )
        yield self.TXD.send(fsm, [pkt])
        yield hold, fsm, duration
        yield fsm.goto(self.SEND)

    def RECV(self, fsm):
        """Manage upstream traffic received from lower layer.

        Based on the value of `htype`, this method will transition to the
        appropriate message handler state (e.g. `ETHRECV` for
        `const.ARP_HTYPE_ETHERNET`).
        """
        while fsm.active():
            yield self.RXD.recv(fsm, 1)
            assert fsm.acquired(self.RXD) and (len(fsm.got)==1), \
                    "[ALOHA]: Error receiving from 'RXD' port!"
            p = fsm.got[0]
            if (self.htype==const.ARP_HTYPE_ETHERNET):
                yield fsm.goto(self.ETHRECV, p)
            else:
                raise RuntimeError, \
                        "[ALOHA]: Cannot handle incoming packet!" + \
                        " Unsupported hardware type (%s)!"%(self.htype)
        return

    def ETHRECV(self, fsm, p):
        """ETHRECV state; receive ethernet frames.

        By default, this method will find ethernet frames, check for CRC errors
        (if available), and forward error-free packets upstream. All other
        packets will be dropped.

        Packets forwarded upstream will be Ethernet frames stripped of their
        CRC-32 using `crcremove()`. Upon completion, this state execution method
        returns to the `RECV` state.

        :note: If no `CRC32` layer (or 'crcerror' annotation) is available, the
               packet is assumed to be error-free (see `MAC.haserror()`).
        """
        assert (self.htype==const.ARP_HTYPE_ETHERNET), "[ALOHA]: In ETHRECV," +\
                " non-Ethernet hardware type (%s) not allowed!"%(self.htype)
        # drop non-Ethernet packet
        addr = self.address
        iseth = isinstance(p, Packet) and p.haslayer(Ether)
        if not iseth:
            self.log_drop(p, addr=addr, drop="non-Ethernet packet")
            yield fsm.goto(self.RECV)
        # set condition flags
        eth = p[Ether]
        addr, src, dst, etype = self.address, eth.src, eth.dst, eth.type
        crcerror = self.haserror(eth)   # check for errors
        isforme  = (dst==addr) or (dst==self.broadcast)
        promiscuous = self.promiscuous
        kwargs = {'addr':addr,'src':src,'dst':dst,'type':etype} 
        # handle packet
        if (isforme or promiscuous) and (not crcerror):
            pkt = crcremove(eth)
            self.log_recv(pkt,promiscuous=promiscuous,**kwargs)
            yield self.TXU.send(fsm, [pkt])
        elif isforme and crcerror:
            self.log_drop(eth, drop="CRC error", **kwargs)
        else:
            self.log_drop(eth, drop="not for me", **kwargs)
        # return to RECV
        yield fsm.goto(self.RECV)

    def log_send(self, p, *args, **kwargs):
        """Convenience method for logging send event."""
        if self.verbose>ALOHA_VERBOSE:
            self.log("snd", p, *args, **kwargs)

    def log_recv(self, p, *args, **kwargs):
        """Convenience method for logging receive event."""
        if self.verbose>ALOHA_VERBOSE:
            self.log("rcv", p, *args, **kwargs)

    def log_drop(self, p, *args, **kwargs):
        """Convenience method for logging drop event."""
        if self.verbose>ALOHA_VERBOSE:
            self.log("drp", p, *args, **kwargs)
