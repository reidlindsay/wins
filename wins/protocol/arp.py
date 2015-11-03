#!  /usr/bin/env python

"""
Implementation of address resolution protocol; contains `ARP` class.

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

:var USE_CRC32:
    Boolean; if true, add/remove CRC-32 when sending/receiving packets.

:var USE_ARPCACHE:
    Boolean; if true, cache ARP information from all incoming packets.
"""
__docformat__ = "restructuredtext en"


from scapy.all import IP, Ether

from SimPy.Simulation import hold
from wins.base import Reference
from wins.crc  import crcupdate, crcremove
from wins.fsm  import FSM
from wins.element import Element
from wins.packet  import Packet, ANNO
from wins.protocol.net import NET
from wins.protocol.mac import MAC
from wins.protocol.arp_support import *

from wins import const

USE_CRC32 = 0
USE_ARPCACHE = 1

class ARP(Element):
    """Address resolution protocol.

    This element translates packets from protocol address space to hardware
    address space. That is, based on the value of the `NET.ptype` and
    `MAC.htype` parameters, this module will translate between network and MAC
    addresses.

    Network Traffic
    ===============
    If a packet from a `NET` module has an unknown protocol destination address,
    `ARP` will drop the packet and send an ARP Request.

    MAC Traffic
    ===========
    `ARP` does not perform address filtering on traffic from connected `MAC`
    protocols. This means that traffic from a `MAC` will be forwarded directly
    to its associated `NET` (unless it is an ARP message).

    :CVariables:
     * `listen`: Property to access dictionary containing listener processes to
       monitor elements connected to `ARP`.

     * `table`: Property to access dictionary containing mapping from protocol
       addresses to hardware addresses (use accessor/modifier method
       `hasentry()`, `addentry()`, and `delentry()`).

     * `useshared`: Property to access boolean flag indicating if ARP
       information should be globally shared.

       The flag can only be set during `configure()` and must be passed through
       the constructor.

    :IVariables:
     * `__listen`: Private dictionary accessed through `listen`.
     * `__table`: Private dictionary accessed through `table`.
     * `__useshared`: Private flag accessed through `useshared`.
    """
    name = "arp"
    tracename = "ARP"
    shared = {}
    def __init__(self, **kwargs):
        self.__listen = {}
        self.__table = {}
        self.__useshared = False
        Element.__init__(self, **kwargs)

    listen = property(fget=lambda self: self.__listen)
    table = property(fget=lambda self: self.__table)
    useshared = property(fget=lambda self: self.__useshared)

    def configure(self, useshared=False, **kwargs):
        """Clear `table` and set up parameters.

        :param useshared: Boolean flag indicating if ARP information should be
                          globally shared.

        :note: `useshared` can only be set if passed through the constructor as
               a keyword argument.
        """
        self.__useshared = useshared
        self.table.clear()

    def addport(self, *args, **kwargs):
        """Depricated; do not use directly; use `connect()` instead."""
        raise RuntimeError, "[ARP]: addport() is depricated! " + \
                            "Use connect() instead!"

    def delport(self, *args, **kwargs):
        """Depricated; do not use directly; use `disconnect()` instead."""
        raise RuntimeError, "[ARP]: delport() is depricated! " + \
                            "Use disconnect() instead!"

    def connect(self, net, mac):
        """Connect a `NET` and `MAC` to this `ARP` module.

        :param net: `NET` to be associated with `mac`.
        :param mac: `MAC` associated with `net`.

        This method will disconnect `net` and `mac` prior to connecting the
        elements to the `ARP` module. After creating the connection, this method
        will add the address info to `table` using `addentry()`. It will also
        add a mapping for the broadcast addresses as well.
        """
        if isinstance(net, Reference): net = net._deref
        if isinstance(mac, Reference): mac = mac._deref
        errmsg = "[ARP]: Cannot connect to non-NET (%s)!"%(net)
        assert isinstance(net, NET), errmsg
        errmsg =  "[ARP]: Cannot connect to non-MAC (%s)!"%(mac)
        assert isinstance(mac, MAC), errmsg
        self.disconnect(net)
        self.disconnect(mac)
        # add ports
        addport = lambda p: Element.addport(self, p)
        ptx, prx = addport((net,"TX")), addport((net,"RX"))
        htx, hrx = addport((mac,"TX")), addport((mac,"RX"))
        # connect net <-> mac through ARP
        ptx.connect(net.getport("RXD"))     # connect to net
        net.getport("TXD").connect(prx)
        htx.connect(mac.getport("RXU"))     # connect to mac
        mac.getport("TXU").connect(hrx)
        # create listeners as children
        tracename = "%s.%s"%(self.tracename,net.traceid)
        f = self.newchild((net, "listen"), FSM, tracename=tracename)
        f.goto(self.PSEND, net, mac)
        tracename = "%s.%s"%(self.tracename,mac.traceid)
        g = self.newchild((mac, "listen"), FSM, tracename=tracename)
        g.goto(self.HRECV, net, mac)
        # add listeners and start
        self.listen[net] = (mac, f)
        self.listen[mac] = (net, g)
        f.start(), g.start()
        errmsg =  "[ARP]: Error connecting (%s) <-> (%s)!"%(net, mac)
        assert self.connected(net) and self.connected(mac), errmsg
        # add entry (paddr, ptype) <-> (haddr, htype); and broadcast entry
        self.addentry(net.addr, mac.addr, net.ptype, mac.htype)
        self.addentry(net.broadcast, mac.broadcast, net.ptype, mac.htype)

    def disconnect(self, x):
        """Disconnect a previously connect element.

        :param x: Previously connected `NET` or `MAC` to disconnect.

        If `x` is not connected, this method will do nothing.
        """
        if isinstance(x, Reference): x = x._deref
        if (x not in self.listen): return
        errmsg = "[ARP]: Error! Could not find ports corresponding to (%s)!"%(x)
        assert self.hasport((x,"TX")) and self.hasport((x,"RX")), errmsg
        # disconnect and delete entry for net and associated mac
        y, listener = self.listen[x]
        if isinstance(x, MAC) and isinstance(y, NET):
            x.getport("TXU").disconnect()   # disconnect MAC
        elif isinstance(x, NET) and isinstance(x, MAC):
            self.delentry(x.addr, y.addr, x.ptype, y.htype)
            x.getport("TXD").disconnect()   # disconnect NET
        else:
            errmsg = "[ARP]: Cannot disconnect (%s) non-(MAC or NET)!"%(x)
            raise RuntimeError, errmsg
        # delete ports, child, and listener
        self.delport((x,"TX")), self.delport((x,"RX"))
        self.delchild((x,"listen"))
        del self.listen[x]
        self.disconnect(y)      # disconnect partner
        errmsg =  "[ARP]: disconnect() failed!"
        assert not (self.connected(x) or self.connected(y) ), errmsg

    def connected(self, x):
        """Check if `x` is properly connected to `ARP`."""
        if isinstance(x, Reference): x = x._deref
        # check if MAC
        conn = False
        if isinstance(x, MAC):
            conn = self.hasport((x,"TX")) and self.hasport((x,"RX")) \
                    and (x in self.listen) \
                    and (x.getport("TXU").target is self.getport((x,"RX")) ) \
                    and (self.getport((x,"TX")).target is x.getport("RXU").target)
        elif isinstance(x, NET):
            conn = self.hasport((x,"TX")) and self.hasport((x,"RX")) \
                    and (x in self.listen) \
                    and (x.getport("TXD").target is self.getport((x,"RX")) ) \
                    and (self.getport((x,"TX")).target is x.getport("RXD").target)
        return conn

    def hasentry(self, paddr, haddr, ptype=None, htype=None):
        """Check if mapping from `paddr` to `haddr` exists in `table`.

        :param paddr: Protocol address.
        :param haddr: Hardware address.
        :param ptype: Protocol type [default=`const.ARP_PTYPE_IP`]
        :param htype: Hardware type [default=`const.ARP_HTYPE_ETHERNET`]

        :note: If `useshared` is true, this method will also check `shared` for
               the requested entry.
        """
        if ptype is None: ptype = const.ARP_PTYPE_IP
        if htype is None: htype = const.ARP_HTYPE_ETHERNET
        pinfo = (paddr, ptype)
        hinfo = (haddr, htype)
        found = (pinfo in self.table) and (self.table[pinfo]==hinfo)
        if (not found) and self.useshared:
            # look in shared table
            found = (pinfo in self.shared) and (self.shared[pinfo]==hinfo)
        return found

    def addentry(self, paddr, haddr, ptype=None, htype=None, **kwargs):
        """Add mapping from `paddr` to `haddr` to `table`.

        :param paddr: Protocol address.
        :param haddr: Hardware address.
        :param ptype: Protocol type [default=`const.ARP_PTYPE_IP`]
        :param htype: Hardware type [default=`const.ARP_HTYPE_ETHERNET`]
        :param kwargs: Additional keywords passed to `log()`.

        This method will also log a trace event 'CACHE' for this entry.

        :note: If `useshared` is true, this method will add the entry to
               `shared` as well as the local `table`.
        """
        if ptype is None: ptype = const.ARP_PTYPE_IP
        if htype is None: htype = const.ARP_HTYPE_ETHERNET
        pinfo = (paddr, ptype)
        hinfo = (haddr, htype)
        self.table[pinfo] = hinfo
        if self.useshared:
            # add to shared table
            self.shared[pinfo] = hinfo
        # log cached entry
        self.log_cache(paddr=paddr, haddr=haddr, \
                       ptype=ptype, htype=htype, **kwargs)

    def delentry(self, paddr, haddr, ptype=None, htype=None):
        """Remove mapping from `paddr` to `haddr` from `table`.

        :param paddr: Protocol address.
        :param haddr: Hardware address.
        :param ptype: Protocol type [default=`const.ARP_PTYPE_IP`]
        :param htype: Hardware type [default=`const.ARP_HTYPE_ETHERNET`]

        If (`paddr`, `haddr`) is not found in `table`, this method does nothing.

        :note: If `useshared` is true, this method will also try to delete the
               entry from `shared`.
        """
        if self.hasentry(padd,haddr,ptype,htype):
            pinfo = (paddr, ptype)
            hinfo = (haddr, htype)
            if pinfo in self.table: del self.table[pinfo]
            if self.useshared and (pinfo in self.shared):
                # delete from shared
                del self.shared[pinfo]

    def translate(self, addr, atype=None):
        """Translate a protocol address to corresponding hardware address (or
        vice versa).

        :param addr: Address to translate (e.g. protocol to hardware).
        :param atype: Address type enumeration (e.g. `const.ARP_PTYPE_IP`).
        :return: 2-tuple containing translated (address, address type).

        This method first checks for a matching protocol address. If no match is
        found, the method searches for a matching hardware address. If neither
        is found, this method returns (None, None).
        """
        # check for protocol address
        ptype = atype
        if atype is None: ptype = const.ARP_PTYPE_IP
        pinfo = (addr, ptype)
        table, found = self.table, (pinfo in self.table)
        if self.useshared and (not found):
            table, found = self.shared, (pinfo in self.shared)
        if found: return table[pinfo]
        # check for hardware address
        htype = atype
        if atype is None: htype = const.ARP_HTYPE_ETHERNET
        hinfo = (addr, htype)
        table, found = self.table, (hinfo in self.table.values())
        if self.useshared and (not found):
            table, found = self.shared, (hinfo in self.shared.values())
        if found:
            idx = table.values().index(hinfo)
            pinfo = table.keys()[idx]
            return pinfo
        # otherwise return (None, None)
        return None, None

    def PSEND(self, fsm, net, mac):
        """PSEND state; Manage outgoing traffic from network protocol.

        Based on the value of `ptype`, this method will invoke the appropriate
        message handler to forward traffic from `net` to `mac`.
        """
        if isinstance(net, Reference): net = net._deref
        if isinstance(mac, Reference): mac = mac._deref
        errmsg = "[ARP]: Must be connected to (%s,%s) to SEND!"%(net,mac)
        assert self.connected(net) and self.connected(mac), errmsg
        # get message from rxport
        rxport = self.getport((net,"RX"))
        yield rxport.recv(fsm, 1)
        errmsg = "[ARP]: Error in receiving from net 'RX' port in SEND!"
        assert fsm.acquired(rxport) and (len(fsm.got)==1), errmsg
        # send packet to appropriate message handler
        p = fsm.got[0]
        fname = "ARPTX.%s(%s)"%(self.traceid, p.traceid)
        if (net.ptype==const.ARP_PTYPE_IP):
            f = FSM(name=fname)
            f.goto(self.IPSEND, net, mac, p)
            f.start()
        else:
            f = FSM(name=fname)
            f.goto(self.PERRSEND, net, mac, p)
            f.start()
        # continue in PSEND
        yield fsm.goto(self.PSEND, net, mac)

    def IPSEND(self, fsm, net, mac, p):
        """IPSEND state; determine hardware address of destination and pass
        network protocol packet to hardware handler `HSEND`.

        :param net: Associated `NET`.
        :param mac: Associated `MAC`.
        :param p: IP packet to send.

        This method uses 'net-dst' annotation to determine the hardware address
        of the destination. If this is not available, it will use the 'dst'
        field of the IP packet `p`.
        
        If a 'net-src' annotation is provided, it should match the address of
        the interface sending the packet. If it is not provided, this method
        will set the annotation using the address of the network interface.

        When the hardware address of the destination cannot be determined, the
        FSM will drop the packet `p` and transition to the `ARPREQ` state to
        send an ARP Request.

        :note: This state will intercept loopback messages (i.e. IP packets
               intended for `net`) and forward them back to `net`.
        """
        # check parameters
        errmsg = "[ARP]: In IPSEND with non-IP ptype (%s)!"%(net.ptype)
        assert (net.ptype==const.ARP_PTYPE_IP), errmsg
        errmsg = "[ARP]: IPSEND unexpectedly got non-IP packet!"
        assert isinstance(p,Packet) and p.haslayer(IP), errmsg
        # use annotations (if possible) to get ip information
        ip = p[IP]
        ethertype, ptype = const.ETHERTYPE_IP, const.ARP_PTYPE_IP
        psrc, pdst = net.addr, ip.dst
        if not ip.hasanno('net-src'): ip.setanno('net-src', psrc)
        if not ip.hasanno('net-dst'): ip.setanno('net-dst', pdst)
        psrc, pdst = ip.getanno('net-src'), ip.getanno('net-dst')
        ip.setanno('net-addr', net.addr)
        # check net-src annotation
        errmsg = "[ARP]: IP spoofing not allowed! (%s != %s)"%(psrc, net.addr)
        assert (psrc==net.addr), errmsg
        # translate IP address to hardware address and check if dest is found
        dstfound = True
        hwdst, htype = self.translate(pdst, ptype)
        if (hwdst is None): dstfound = False
        # intercept loopback
        if (pdst==net.addr):
            txport = self.getport((net,"TX"))
            self.log_loopback(ip, netaddr=net.addr, src=psrc, dst=pdst)
            yield txport.send(fsm, [ip])
            yield fsm.stop()    # HALT packet handler
        # send to hardware handler
        if dstfound:
            yield fsm.goto(self.HSEND, net, mac, ip, hwdst)
        else:
            self.log_drop(ip, src=psrc, dst=pdst, ptype=ptype, \
                          drop="dst not found in IPSEND")
            yield fsm.goto(self.ARPREQ, net, mac, ip, psrc, pdst)

    def PERRSEND(self, fsm, net, mac, p):
        """PERRSEND state; handles unknown ptype values."""
        yield hold, fsm, 0  # yield to make valid generator
        errmsg = "[ARP]: Got unsupported ptype (%s) in PSEND!"%(net.ptype)
        raise RuntimeError, errmsg

    def HSEND(self, fsm, net, mac, *args, **kwargs):
        """HSEND state; determine htype of hardware interface and send to
        appropriate handler.

        :param net: Associated `NET`.
        :param mac: Associated `MAC`.
        :param args: Additional arguments passed to handler.
        :param kwargs: Additional keywords passed to handler.
        """
        htype = mac.htype
        if (htype==const.ARP_HTYPE_ETHERNET):
            yield fsm.goto(self.ETHSEND, net, mac, *args, **kwargs)
        else:
            yield fsm.goto(self.HERRSEND, net, mac, *args, **kwargs)

    def ETHSEND(self, fsm, net, mac, p, dst=None):
        """ETHSEND state; Encapsulate `p` in Ethernet frame and send to `mac`.

        :param net: Associated `NET`.
        :param mac: Associated `MAC`.
        :param p: Ethernet packet to send.
        :param dst: Ethernet address of destination [default=`MAC.broadcast`]

        Upon completion this method returns to `returnstate` with the parameters
        `net` and `mac`.
        """
        # get parameters for Ethernet packet
        src, htype = mac.addr, mac.htype
        if dst is None: dst = mac.broadcast
        ethertype = self.get_ethertype(p)
        eth = Ether(src=src, dst=dst, type=ethertype)
        eth.add_payload(p)
        # update CRC
        pkt = eth
        if USE_CRC32: pkt = crcupdate(eth)
        # send to TX port
        self.log_send(pkt, src=src, dst=dst, type=ethertype)
        txport = self.getport((mac, "TX"))
        yield txport.send(fsm, [pkt])

    def HERRSEND(self, fsm, net, mac, *args, **kwargs):
        """HERRSEND state; handles unknown htype values."""
        yield hold, fsm, 0  # yield to make valid generator
        errmsg = "[ARP]: Got unsupported htype (%s) in HSEND!"%(mac.htype)
        raise RuntimeError, errmsg

    def ARPREQ(self, fsm, net, mac, p, src, dst, *args, **kwargs):
        """ARPREQ state; send ARP Request message."""
        psrc, pdst, ptype = src, dst, net.ptype
        hwsrc, hwdst, htype = mac.addr, mac.nulladdr, mac.htype
        arp = ARPRequest(hwtype=htype, ptype=ptype, hwsrc=hwsrc, \
                             psrc=psrc, hwdst=hwdst, pdst=pdst)
        bcast = mac.broadcast
        yield fsm.goto(self.HSEND, net, mac, arp, bcast, *args, **kwargs)

    def HRECV(self, fsm, net, mac):
        """HRECV state; manage incoming traffic from MAC protocol.

        :param net: Associated `NET`.
        :param mac: Associated `MAC`.

        Based on the value of `htype`, this method will invoke the appropriate
        message handler to classify and handle traffic from `mac` to `net`.
        """
        if isinstance(net, Reference): net = net._deref
        if isinstance(mac, Reference): mac = mac._deref
        errmsg = "[ARP]: Must be connected to (%s,%s) to RECV!"%(net,mac)
        assert self.connected(net) and self.connected(mac), errmsg
        # get message from rxport
        rxport = self.getport((mac, "RX"))
        yield rxport.recv(fsm, 1)
        errmsg = "[ARP]: Error in receiving from mac 'RX' port in RECV!"
        assert fsm.acquired(rxport) and (len(fsm.got)==1), errmsg
        # send packet to appropriate message handler
        p = fsm.got[0]
        fname = "ARPTX.%s(%s)"%(self.traceid, p.traceid)
        if (mac.htype==const.ARP_HTYPE_ETHERNET):
            f = FSM(name=fname)
            f.goto(self.ETHRECV, net, mac, p)
            f.start()
        else:
            f = FSM(name=fname)
            f.goto(self.HERRRECV, net, mac, p)
            f.start()
        # continue in HRECV
        yield fsm.goto(self.HRECV, net, mac)

    def ETHRECV(self, fsm, net, mac, p):
        """ETHRECV state; classify incoming Ethernet packets and pass to
        corrseponding message handler based on ethertype.

        :param net: Associated `NET`.
        :param mac: Associated `MAC`.
        :param p: Ethernet packet to handle.

        If `p` is an ARP Request or Reply, this method will forward the packet
        to `ARPRECV` for handling. Otherwise, this method forwards all other
        received packets to their appropriate message handler based on the
        protocol type of `net` and the ethertype of `p`.

        If `USE_ARPCACHE` is turned on this method will cache ARP information
        using incoming packets.

        :note: If `p` is not an Ethernet frame, this state execution method will
               drop the incoming packet and return to `RECV`.
        """
        errmsg = "[ARP]: In ETHRECV with non-Ethernet htype (%s)!"%(mac.htype)
        assert (mac.htype==const.ARP_HTYPE_ETHERNET), errmsg
        iseth =  isinstance(p,Packet) and p.haslayer(Ether)
        if not iseth:
            self.log_drop(p, drop="non-Ethernet packet in ETHRECV")
            yield fsm.stop()    # HALT packet handler
        # get ethernet frame and payload
        eth = p[Ether]
        hsrc, hdst, ethertype = eth.src, eth.dst, eth.type
        # remove CRC and get payload
        pkt = eth
        if USE_CRC32: pkt = crcremove(eth)
        payload = pkt.payload
        pkt.remove_payload()
        # set MAC annotations from ethernet packet
        payload.setanno('mac-src', eth.src)
        payload.setanno('mac-dst', eth.dst)
        payload.setanno('mac-addr', mac.addr)
        # classify and handle
        if (ethertype==const.ETHERTYPE_IP):
            yield fsm.goto(self.IPRECV, net, mac, payload)
        elif (ethertype==const.ETHERTYPE_ARP):
            yield fsm.goto(self.ARPRECV, net, mac, payload)
        else:
            yield fsm.goto(self.HERRRECV, net, mac, payload)

    def IPRECV(self, fsm, net, mac, p):
        """IPRECV state; extract IP Packet and forward to upstream protocol
        `net`; this state does not do address checking.

        :param net: Associated `NET`.
        :param mac: Associated `MAC`.
        :param p: IP packet to receive.

        Upon completion, this state execution method will return to `RECV`.

        :note: **This method does not do any kind of address checking.** If a
               valid IP packet cannot be found in `ip`, this method will drop
               the incoming packet.
        """
        # check PTYPE and IP packet
        if not (net.ptype==const.ARP_PTYPE_IP):
            self.log_drop(p, drop="non-IP ptype in IPRECV")
            yield fsm.stop()    # HALT packet handler
        isip =  isinstance(p,Packet) and p.haslayer(IP)
        if not isip:
            self.log_drop(p, drop="non-IP packet in IPRECV")
            yield fsm.stop()    # HALT packet handler
        # get IP packet and forward to net
        ip = p[IP]
        psrc, pdst = ip.src, ip.dst
        self.log_recv(ip, src=ip.src, dst=ip.dst, netaddr=net.addr)
        txport = self.getport((net,"TX"))
        yield txport.send(fsm, [ip])

    def HERRRECV(self, fsm, net, mac, *args, **kwargs):
        """HERRRECV state; handles unknown htype values."""
        yield hold, fsm, 0  # yield to make valid generator
        errmsg = "[ARP]: Got unsupported htype (%s) in HRECV!"%(mac.htype)
        raise RuntimeError, errmsg

    def ARPRECV(self, fsm, net, mac, arp):
        """ARPRECV state; receive and respond to an ARP message.

        :param net: Associated `NET`.
        :param mac: Associated `MAC`.
        :param arp: ARP packet to receive.

        If `arp` is an ARP request for `net`, this state execution method will
        respond with an ARP reply. If `arp` is an ARP reply, this method will
        cache the ARP information from the reply. Otherwise this state will
        raise an exception.

        ARP replies responding to ARP requests will be sent out using the
        appropriate send method (i.e. depending on the hardware type of `mac`).

        If `USE_ARPCACHE` is turned on, this method will cache ARP information
        from all received messages.

        :note: If `arp` is not a valid ARP request or ARP reply, this state will
               raise an exception.
        """
        # get ARP parameters
        x, isrequest, isreply = None, isarprequest(arp), isarpreply(arp)
        if isrequest: x = get_arprequest(arp)
        elif isreply: x = get_arpreply(arp)
        isforme, addr, ptype, htype = False, net.addr, net.ptype, mac.htype
        kwargs = {'addr': addr, 'ptype': ptype,  'htype': htype}
        if x:
            isforme = (x.pdst==addr) and (x.ptype==ptype) and (x.hwtype==htype)
            kwargs.update({'psrc':x.psrc,'hwsrc':x.hwsrc, \
                           'pdst':x.pdst,'hwdst':x.hwdst, \
                           'ptype':x.ptype, 'hwtype':x.hwtype})
            if USE_ARPCACHE:
                # cache source information
                self.addentry(x.psrc, x.hwsrc, x.ptype, x.hwtype)
        # process ARP message
        self.log_recv(x, **kwargs)
        if isrequest and isforme:
            # ARP Request -> send ARP Reply
            yield fsm.goto(self.ARPREP, net, mac, x)
        elif isreply and isforme:
            # ARP Reply -> cache source information
            self.addentry(x.psrc, x.hwsrc, x.ptype, x.hwtype)
        elif isforme:
            # not a valid ARP message
            errmsg = "[ARP]: ARPRECV cannot find ARP message!"
            raise RuntimeError, errmsg
        elif isrequest or isreply:
            kwargs['drop'] = "ARP message not for me in ARPRECV"
            self.log_drop(x, **kwargs)

    def ARPREP(self, fsm, net, mac, arpreq, *args, **kwargs):
        """ARPREP state; send ARP Reply."""
        x = get_arprequest(arpreq)
        errmsg = "[ARP]: Cannot reply to non-ARP-Request message!"
        assert (x is not None), errmsg
        # form reply
        psrc, pdst, ptype = x.pdst, x.psrc, x.ptype
        hwsrc, hwdst, hwtype = mac.addr, x.hwsrc, x.hwtype
        arp = ARPReply(hwtype=hwtype, ptype=ptype, hwsrc=hwsrc, \
                       psrc=psrc, hwdst=hwdst, pdst=pdst)
        bcast = mac.broadcast
        yield fsm.goto(self.HSEND, net, mac, arp, bcast, *args, **kwargs)

    def get_ethertype(self, p):
        """Get Ethernet type of packet `p`."""
        ethertype = None
        if isinstance(p, Reference): p = p._deref
        # check packet to determine ethertype
        if isinstance(p, IP):
            ethertype = const.ETHERTYPE_IP          # IP packet
        elif isarprequest(p) or isarpreply(p):
            ethertype = const.ETHERTYPE_ARP         # ARP packet
        else:
            errmsg = "[ARP]: Unknown Ethernet type for packet!"
            raise RuntimeError, errmsg
        return ethertype

    def log_cache(self, *args, **kwargs):
        """Convenience method for logging entry cache."""
        self.log("cache", *args, **kwargs)

    def log_loopback(self, p, *args, **kwargs):
        """Convenience method for logging loopback."""
        self.log("LO", p, *args, **kwargs)

    def log_send(self, p, *args, **kwargs):
        """Convenience method for logging send event."""
        self.log("snd", p, *args, **kwargs)

    def log_recv(self, p, *args, **kwargs):
        """Convenience method for logging receive event."""
        self.log("rcv", p, *args, **kwargs)

    def log_drop(self, p, *args, **kwargs):
        """Convenience method for logging drop event."""
        self.log("drp", p, *args, **kwargs)

    def log(self, evt=None, p=None, *args, **kwargs):
        """Overloaded to check verbose level and set common annotations."""
        force = False
        if ('verbose' in kwargs): force = (kwargs['verbose']>ARP_VERBOSE)
        if self.verbose>ARP_VERBOSE or force:
            Element.log(self, evt, p, *args, **kwargs)
