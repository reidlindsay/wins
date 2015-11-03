#!  /usr/bin/env python

"""
Base class for implementing routing protocols; contains `Routing` class.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-10-23 04:50:01 -0500 (Sun, 23 Oct 2011) $
* $LastChangedRevision: 5281 $

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

:var USE_NETSRC_SPOOFING:
    Boolean flag; if true, allow 'net-src' annotation to be used by upstream
    protocols to spoof the network address of packets sent by the protocol.

:var ROUTING_VERBOSE:
    Constant enumeration to control verbose thresholds in this file.
  
    Setting the verbose level of a `Routing` above this threshold will cause the
    corresponding output in this file to be written (or logged).
"""
__docformat__ = "restructuredtext en"

from SimPy.Simulation import hold

from wins.protocol.net import NET
from wins.protocol.routetable import RouteTable
from wins.packet import ANNO
from wins.base   import Reference
from wins.fsm    import FSM
from wins import const

from scapy.all import TCP, UDP, IP, ICMP, Packet

from wins.protocol.route_support import *

# FIXME: Should upstream protocols be able to use 'net-src' to spoof the source
#        address of the network protocol?
USE_NETSRC_SPOOFING = 0
ROUTING_VERBOSE = 30

class Routing(NET):
    """Base class for routing protocols.

    The base class provides the methods needed for maintaining a local routing
    table and implements a simple static routing protocol. If a route to the
    destination cannot be found for a packet, this protocol will drop the
    packet. Subclasses can implement custom routing algorithms as needed.

    The routing table contains two entries by default:

        1. A route to itself (i.e. `address` -> `address`).
        #. A broadcast entry (i.e. `broadcast` -> `broadcast`).

    This is implemented in `set_table()`. Overload methods as needed.

    Routing and Ports
    =================
    In addition to the ports created by `NET`, this class has two ports:

        1. "TXU" - sends traffic to an upstream element/protocol.
        #. "RXU" - receives traffic from an upstream element/protocol.

    Overload `configure()` to change how ports are set up.

    Routing and Annotations
    =======================
    The following annotations are used or marked by `Routing`.

    ========== ==========================================================
    Name        Description
    ========== ==========================================================
    net-dst     Network address for destination of packet. If not marked,
                `Routing` assumes packet is for `broadcast`.
                *This should be set by an upstream protocol.*
    ---------- ----------------------------------------------------------
    net-src     `Routing` marks this annotation to indicate the address
                of the sender. An upstream protocol may overload this
                annotation to spoof the send address of a packet.
    ========== ==========================================================

    :cvar table: Property to access/modify `RouteTable`.
    :cvar MaxTTL: Maximum Time-to-live value for packets in the network.

    :ivar __table: Private member to maintina `RouteTable`.

    :note: The routing `table` is also accessible as a child with the nick name
           'routetable'. This is set throught `set_table()`.
    """
    name = "routing"
    tracename = "RTG"
    MaxTTL = 16
    def __init__(self, **kwargs):
        """Constructor."""
        self.__table = None
        NET.__init__(self, **kwargs)
        assert isinstance(self.table, RouteTable), \
                "[ROUTING]: Must have valid RouteTable!"
        assert self.hasroute(self.address), \
                "[ROUTING]: RouteTable has no entry for local address!"
        assert self.hasroute(self.address), \
                "[ROUTING]: RouteTable has no entry for broadcast address!"

    table = property(fget=lambda self: self.__table, \
                     fset=lambda self,t: self.set_table(t) )

    def configure(self, **kwargs):
        """Call `NET.configure()`, initialize routing table, add ports, and
        spawn `FSM` to run static routing algorithm."""
        NET.configure(self, **kwargs)
        self.table = self.newchild('routetable', RouteTable, \
                                   name=self.name+".rtable", \
                                   tracename=self.tracename+".RTABLE")
        self.addport("TXU")     # port to send traffic to upstream protocol
        self.addport("RXU")     # port to recv traffic from upstream protocol
        # create FSM to manage send/recv execution of routing
        txfsm = self.newchild("txfsm", FSM, tracename=self.tracename+".TX")
        rxfsm = self.newchild("rxfsm", FSM, tracename=self.tracename+".RX")
        txfsm.goto(self.SEND)
        rxfsm.goto(self.RECV)

    ##############################
    # TX STATES
    ##############################
    def SEND(self, fsm):
        """SEND state; monitor traffic from upstream protocol and pass to
        routing algorithm.

        This state uses the 'net-src' and 'net-dst' annotations to determine the
        addresses of the source and destination. If not available, these values
        default to the local `address` and `broadcast` address.

        This state then uses `ptype` to determine which packet handler should
        process the new data packet, and spawns a new thread to run the handler.
        
        :note: Currently, IPv4 is the only protocol type supported. Non-IPv4
               packets are passed to `ERRSEND` for handling. Also, this state
               will check for loopback addresses and send them back to the upper
               layer protocol.
        """
        # error messages
        rxuerror = "[ROUTING]: Error occurred receiving from RXU port!"
        pkterror = "[ROUTING]: Invalid packet from upper layer!"
        # get packet to send from upperlayer
        yield self.RXU.recv(fsm, 1)
        assert fsm.acquired(self.RXU) and (len(fsm.got)==1), rxuerror
        p = fsm.got[0]
        if isinstance(p, Reference): p = p._deref
        assert ANNO.supported(p), pkterror
        # process new data packet
        addr = self.address
        src, dst = addr, self.broadcast
        if p.hasanno('net-src'): src = p.getanno('net-src')
        if p.hasanno('net-dst'): dst = p.getanno('net-dst')
        r = self.set_sendanno(p, src, dst)
        # send loopback back to upper layer
        if (dst==addr):
            self.log_loopback(p, src=src, dst=dst)
            yield self.TXU.send(fsm, [p])
            yield fsm.goto(self.SEND)           # continue in SEND
        # otherwise -> use appropriate send function based on ptype
        if (self.ptype == const.ARP_PTYPE_IP):
            f = FSM.launch(self.IPSEND, r, src, dst)        # IP ptype
        else:
            f = FSM.launch(self.ERRSEND, r, src, dst)       # unknown ptype
        # continue in SEND
        yield fsm.goto(self.SEND)

    def ERRSEND(self, fsm, p, src, dst):
        """ERRSEND state; process packets for unknown `ptype`.

        :param p: Packet to send.
        :param src: Network address of source.
        :param dst: Network address of destination.

        By default this method just raises an exception. Overload this state in
        subclasses to add support for new protocol types.
        """
        errmsg = "[ROUTING]: Unsupported protocol type (%s)!"%(self.ptype)
        yield hold, fsm, 0          # yield required to make valid generator
        raise RuntimeError, errmsg  # throw exception

    def IPSEND(self, fsm, p, src, dst, **kwargs):
        """IPSEND state; encapsulate packet using available routing info and
        pass to appropriate handler.

        :param p: Packet to send.
        :param src: Network address of source.
        :param dst: Network address of destination.
        :param kwargs: Additional keywords passed to `encap_data()`.

        This state uses `hasroute()` to determine if a route exists to the
        desired destination. If its known, transition to `IPFORWARD`, otherwise
        go to `IPROUTING`.

        :note: This state assumes an IPv4 protocol type is being used.
        """
        # enacapsulate packet into valid IP packet
        ip = self.encap_data(p, src, dst, **kwargs)
        if not isinstance(ip, IP):
            self.log_drop(p, drop="encap error")
            yield fsm.stop()    # HALT and discard packet
        # is route to next hop known?
        if self.hasroute(dst):
            yield fsm.goto(self.IPFORWARD, ip)
        else:
            yield fsm.goto(self.IPROUTING, ip)

    def IPFORWARD(self, fsm, p):
        """IPFORWARD state; update packet as needed before passing to
        `IPDELIVER`.

        :param p: IP packet to send.

        :note: This state assumes a valid IPv4 packet was went to it.
        """
        # get IP parameters
        ip = p[IP]
        dst = ip.dst
        # deliver to nexthop
        nexthop = self.nexthop(dst)
        if nexthop:
            yield fsm.goto(self.IPDELIVER, ip, nexthop)

    def IPDELIVER(self, fsm, p, nexthop):
        """IPDELIVER state; deliver packet to lower layer.

        :param p: IP packet to send.
        :param nexthop: Network address of next hop destination.

        This state updates the 'net-src' and 'net-dst' annotations using the
        values of `address` and `nexthop` respectively. Then the packet is sent
        to the lower layer (if the TTL has not expired).

        :note: This state assumes a valid IPv4 packet was went to it.
        """
        # error messages
        ttldrop = "TTL expired"
        # update annotations
        addr = self.address
        ip = self.set_sendanno(p, addr, nexthop)
        if isinstance(ip, IP):
            ttl = ip.ttl
            if (ttl<1):
                self.log_drop(ip, drop=ttldrop, nexthop=nexthop)
            else:
                chksum = self.updatechksum(ip, overwrite=True)
                if (ip.src==addr): self.log_send(ip, nexthop=nexthop)
                else:              self.log_forward(ip, nexthop=nexthop)
                yield self.TXD.send(fsm, [ip])

    def IPROUTING(self, fsm, p):
        """IPROUTING state; find route for packet.

        :param p: Packet to send.
        :param src: Network address of source.
        :param dst: Network address of destination.
        :param kwargs: Additional keyword arguments.

        *Overload this method to implement a routing protocol.*

        :note: This state assumes a valid IPv4 packet was went to it.
        """
        # error messages
        droproute = "no route found"
        # get IP parameters
        ip = p[IP]
        dst = ip.dst
        # drop packet
        assert not self.hasroute(dst)
        self.log_drop(ip, drop=droproute)
        yield fsm.stop()    # HALT and discard packet

    ##############################
    # RX STATES
    ##############################
    def RECV(self, fsm):
        """RECV state; monitor traffic from lower layer.

        This state receives packets on `Port` 'RXD'. It then uses `ptype` to
        determine which packet handler should process the new data packet, and
        spawns a new thread to run the handler.
        
        :note: Currently, IPv4 is the only protocol type supported. Non-IPv4
               packets are passed to `ERRRECV` for handling.
        """
        # error messages
        rxderror = "[ROUTING]: Error occurred while receiving " + \
                   "traffic from RXD port!"
        pkterror = "[ROUTING]: Invalid packet from lower layer!"
        # get packet from lower layer
        yield self.RXD.recv(fsm, 1)
        assert fsm.acquired(self.RXD) and (len(fsm.got)==1), rxderror
        p = fsm.got[0]
        if isinstance(p, Reference): p = p._deref
        assert ANNO.supported(p), pkterror
        # use appropriate recv function based on ptype
        if (self.ptype == const.ARP_PTYPE_IP):
            f = FSM.launch(self.IPRECV, p)      # IP ptype
        else:
            f = FSM.launch(self.ERRRECV, p)     # unknown ptype
        # continue in RECV
        yield fsm.goto(self.RECV)

    def ERRRECV(self, fsm, p):
        """ERRRECV state; process packets for unknown `ptype`.

        By default this method just raises an exception. Overload this method in
        subclass to add support for new protocol types.
        """
        errmsg = "[ROUTING]: Unsupported protocol type (%s)!"%(self.ptype)
        yield hold, fsm, 0          # yield required to make valid generator
        raise RuntimeError, errmsg  # throw exception

    def IPRECV(self, fsm, p):
        """IPRECV state; error check and classify IP packet.

        :param p: Received packet.

        This state uses `isctrl()` and `isdata()` to classify packets and send
        them to `IPCTRL` or `IPDATA` for handling.

        :note: This state drops any non-IP packets.
        """
        # error messages
        dropnonip = "non-IP packet in IPRECV"
        dropinvalid = "invalid packet received"
        # check if valid IP packet was received
        self.log_recv(p)
        drop = self.checkiprecv(p)
        if drop:
            self.log_drop(p, drop=drop)
            yield fsm.stop()    # HALT and discard packet
        # classify IP packet
        ip = p[IP]
        yield fsm.goto(self.IPCLASSIFY, ip)

    def IPCLASSIFY(self, fsm, p):
        """IPCLASSIFY state; classify packet received in `IPRECV`.

        :param p: Received IP packet.

        *Overload to change how packets are classified and handled*.

        :note: This state assumes that packet `p` passed all validity checks in
               `IPRECV` (i.e. `checkiprecv()`).
        """
        # error messages
        drop = "not for me|no route found"
        # get IP parameters
        ip = p[IP]
        src, dst = ip.src, ip.dst
        # classify IP packet
        if (dst==self.address) or (dst==self.broadcast):
            # packet for me -> detach payload
            payload = ip.payload
            ip.remove_payload()
            self.log_recv(payload, src=src, dst=dst)
            if (ip.proto==const.IP_PROTO_IP):
                yield fsm.goto(self.IPRECV, payload)    # IPENCAP -> reclassify
            elif (ip.proto!=const.IP_PROTO_NONE):
                pkt = self.set_recvanno(payload, src, dst)
                yield self.TXU.send(fsm, [pkt])         # DATA -> send to upper layer
        elif self.hasroute(dst):
            # has route to dst -> update parameters + forward copy to dst
            pkt = ip.copy()
            pkt.ttl = ip.ttl - 1
            chksum = self.updatechksum(pkt, overwrite=True)
            yield fsm.goto(self.IPFORWARD, pkt)
        else:
            # otherwise -> drop
            self.log_drop(ip, drop=drop)

    ##############################
    # ANNOTATION METHODS
    ##############################
    def set_sendanno(self, p, src, dst):
        """Set relevant annotations for outgoing packet.

        :param p: Packet to modify.
        :param src: Network address of source.
        :param dst: Network address of destination.
        :return: Modified packet.

        This method sets the 'net-src' and 'net-dst' annotations.
        """
        errmsg = "[ROUTING]: Got non-Packet data from upper layer()!"
        assert ANNO.supported(p), errmsg
        # set address annotations
        p.setanno('net-src', src)
        p.setanno('net-dst', dst)
        # set annotations used just for logging
        p.setanno('net-root', str(p.traceid))
        return p

    def set_recvanno(self, p, src, dst):
        """Set relevant annotations for received data to upper layer.

        :param p: Packet to modify.
        :param src: Network address of source.
        :param dst: Network address of destination.
        :return: Modified packet.

        This method sets the 'net-src' and 'net-dst' annotations.
        """
        errmsg = "[ROUTING]: Got non-Packet data from upper layer()!"
        assert ANNO.supported(p), errmsg
        # set address annotations
        p.setanno('net-src', src)
        p.setanno('net-dst', dst)
        return p

    ##############################
    # HELPER METHODS
    ##############################
    def encap_data(self, p, src, dst, force=False, **kwargs):
        """Encapsulate data packet (if needed) for delivery to next hop.

        :param p: Packet to deliver.
        :param src: Source address.
        :param dst: Destination address.
        :param force: If true, force encapsulation into a new IP packet.
        :param kwargs: Additional keywords passed to IP constructor.
        :return: Encapsulated/modified packet for delivery.

        If `p` is already an IP packet and `force` is false, this method will
        update 'src', 'dst', and 'ttl' fields as needed. Otherwise, this method
        encapsulates `p` with an IP header.

        :note: By default this method assumes that `ptype` indicates IPv4.
               *Overload this method to handle other protocol types*.
        """
        # if no encapsulation needed -> update parameters
        isip = isinstance(p, Packet) and p.haslayer(IP)
        if isip and (not force):
            # update IP parameters
            ip = p[IP]
            ip.src, ip.dst = src, dst
            if 'ttl' in kwargs: ip.ttl = kwargs['ttl']
            proto = self.getproto(ip.payload)
            ip.proto = proto
            ip.id = ip._id
            # update chksum
            chksum = self.updatechksum(ip, overwrite=True)
            return ip
        # set up TTL parameter
        ttl = None
        if 'ttl' in kwargs:
            ttl = kwargs['ttl']
            del kwargs['ttl']
        elif (dst==self.broadcast):
            ttl = 1                     # broadcast only use 1 hop
        # correct TTL if needed
        if (ttl is None): ttl = self.MaxTTL
        ttl = min(ttl, self.MaxTTL)
        # create IP packet
        ip = IP(**kwargs)
        ip.add_payload(p)
        # Scapy kludge: build payload to make sure no missing fields are unset
        if isinstance(p, Packet): s = p.do_dissect(p.build())
        # use recursive call to set remaining parameters
        return self.encap_data(ip, src, dst, ttl=ttl, force=False)

    def checkiprecv(self, p, chksum=False, proto=True):
        """Check for valid IP packet and errors.

        :param p: Packet received in `IPRECV`.
        :param chksum: If true, check 'chksum' of `p`.
        :param proto: If true, check consistency of 'proto' field.
        :return: String containing any error condition, or nothing for no error.
        """
        # error messages
        dropnonip  = "non-IP packet in IPRECV"
        dropchksum = "IP chksum failed"
        dropproto  = "invalid IP proto"
        # valid IP packet?
        isip = isinstance(p, Packet) and p.haslayer(IP)
        if not isip: return dropnonip
        # get IP packet/parameters
        ip = p[IP]
        # chksum ok?
        if chksum:
            oldchksum = ip.chksum
            newchksum = self.updatechksum(ip, overwrite=False)
            if (oldchksum!=newchksum): return dropchksum
        # check proto field
        if proto:
            payproto = self.getproto(ip.payload)
            if (payproto!=ip.proto): return dropproto
        return ""   # no error

    def updatechksum(self, p, overwrite=True):
        """Update checksum and len field of IP packet.

        :param p: IP packet to check.
        :param overwrite: If true, overwrite the 'chksum' field of packet `p`
                          [default=True].
        :return: Newly computed checksum.
        """
        # check for valid IP packet
        errmsg = "[ROUTING]: updatechksum() requires IP packet!"""
        assert p.haslayer(IP), errmsg
        # remove ANNO from IP before recomputing checksum
        ip = p[IP]
        anno = None
        if ip.haslayer(ANNO):
            anno = ip[ANNO]
            anno.underlayer.remove_payload()
        # update parameters and recalculate chksum
        ip.len = len(ip)
        pkt = ip.copy()
        del pkt.chksum                  # reset "automatic" field
        s = pkt.do_dissect(pkt.build())
        # reattach ANNO and set chksum
        if anno: ip.add_payload(anno)
        chksum = pkt.chksum
        # overwrite chksum and return new chksum
        if overwrite:
            ip.chksum = chksum
        return chksum

    def getproto(self, p):
        """Convenience method to determine IPv4 protocol opcode for packet `p`.

        :return: IPv4 protocol opcode for `p`; this can be used when
                 encapsulating `p` in an IPv4 packet.

        By default, this method will return `const.IP_PROTO_NONE` if no valid
        payload can be found. Unidentifiable `Packet` types will be classified
        with the proto value `const.IP_PROTO_HOST`.
        """
        # error messages
        errorptype = "[ROUTING]: getproto() only supports " + \
                     "IPv4 ptype (%s)"%(const.ARP_PTYPE_IP)
        # check for IPv4 ptype
        assert (self.ptype==const.ARP_PTYPE_IP), errorptype
        if isinstance(p, Reference): p = p._deref
        # determine protocol type of packet
        proto = const.IP_PROTO_NONE
        if   isinstance(p, TCP):    proto = const.IP_PROTO_TCP
        elif isinstance(p, UDP):    proto = const.IP_PROTO_UDP
        elif isinstance(p, IP):     proto = const.IP_PROTO_IP
        elif isinstance(p, ICMP):   proto = const.IP_PROTO_ICMP
        elif isinstance(p, ANNO):   proto = const.IP_PROTO_NONE
        elif isinstance(p, Packet): proto = const.IP_PROTO_HOST
        return proto

    ##############################
    # ROUTE TABLE METHODS
    ##############################
    def hasroute(self, *args, **kwargs):
        """Wrapper to access `RouteTable.hasroute()` for `table`."""
        return self.table.hasroute(*args, **kwargs)

    def addroute(self, *args, **kwargs):
        """Wrapper to access `RouteTable.addroute()` for `table`."""
        return self.table.addroute(*args, **kwargs)

    def delroute(self, *args, **kwargs):
        """Wrapper to access `RouteTable.delroute()` for `table`."""
        return self.table.delroute(*args, **kwargs)

    def nexthop(self, *args, **kwargs):
        """Wrapper to access `RouteTable.nexthop()` for `table`."""
        return self.table.nexthop(*args, **kwargs)

    def cost(self, *args, **kwargs):
        """Wrapper to access `RouteTable.cost()` for `table`."""
        return self.table.cost(*args, **kwargs)

    def set_table(self, rt):
        """Add new routing table as a child wih nickname 'routetable', and add
        basic routing entries (i.e. `address` and `broadcast`)."""
        assert isinstance(rt, RouteTable) or (rt is None), \
                "[ROUTING]: Cannot set_table() to non-RouteTable!"
        table = self.addchild('routetable', rt)
        table.addroute(self.address)
        table.addroute(self.broadcast)
        self.__table = table

    ##############################
    # LOGGING METHODS
    ##############################
    def log_send(self, p, *args, **kwargs):
        """Convenience method for logging send event."""
        self.log("snd", p, *args, **kwargs)

    def log_recv(self, p, *args, **kwargs):
        """Convenience method for logging receive event."""
        self.log("rcv", p, *args, **kwargs)

    def log_forward(self, p, *args, **kwargs):
        """Convenience method for logging forward event."""
        self.log("fwd", p, *args, **kwargs)

    def log_drop(self, p, *args, **kwargs):
        """Convenience method for logging drop event."""
        self.log("drp", p, *args, **kwargs)

    def log_loopback(self, p, *args, **kwargs):
        """Convenience method for logging loopback."""
        self.log("LO", p, *args, **kwargs)

    def log(self, evt=None, p=None, *args, **kwargs):
        """Overloaded to check verbose level and set common annotations."""
        force = False
        if ('verbose' in kwargs): force = (kwargs['verbose']>ROUTING_VERBOSE)
        if self.verbose>ROUTING_VERBOSE or force:
            kwargs.update(self.get_route_anno(p))
            kwargs.update(self.get_ip_anno(p))
            NET.log(self, evt, p, *args, **kwargs)

    def get_route_anno(self, p):
        """Internal method to get relevant annotations."""
        kwargs = {}
        if not isinstance(p, Packet): return kwargs
        kwargs['addr'] = self.address
        if p.hasanno('net-src'):
            kwargs['net-src'] = p.getanno('net-src')
        if p.hasanno('net-dst'):
            kwargs['net-dst'] = p.getanno('net-dst')
        return kwargs

    def get_ip_anno(self, p):
        """Internal method to get relevant annotations for an IP packet."""
        kwargs = {}
        if not (isinstance(p, Packet) and p.haslayer(IP)): return kwargs
        ip = p[IP]
        kwargs['ip.src'] = ip.src
        kwargs['ip.dst'] = ip.dst
        kwargs['ip.proto'] = ip.proto
        kwargs['ip.ttl'] = ip.ttl
        return kwargs
