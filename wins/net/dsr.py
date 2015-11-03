#!  /usr/bin/env python

"""
Implementation of Dynamic Source Routing protocol (RFC 4728); contains `DSR` class.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-12-20 23:58:28 -0600 (Tue, 20 Dec 2011) $
* $LastChangedRevision: 5394 $

:author: Ketan Mandke <kmandke@mail.utexas.edu>

:copyright:
    Copyright 2009-2011 The University of Texas at Austin

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

DSRDEBUG = 0

from SimPy.Simulation import SimEvent
from SimPy.Simulation import hold, waitevent
from wins.fsm import FSM
from wins.base import Reference
from wins.helper import time2msec
from wins import const

from wins.protocol.arp import ARP
from wins.protocol.mac import MAC
from wins.protocol.net import NET
from wins.protocol.routing import Routing

from wins.net.dsr_support import *
from wins.net.dsr_support import _dsr_error_types
from wins.net.dsr_rreq import RouteRequestTable, RouteRequestEntry
from wins.net.routecache import RouteCache

from wins.traffic import AGT

from numpy import random

class DSR(Routing):
    """Dynamic Source Routing protocol (RFC 4728).

    This implementation uses link-level acknowledgements to do route
    maintenance (see `MAC.drpdata` and `MAC.ackdata`).

    :ivar rreqtable: `RouteRequestTable` used to manage route requests.
    :ivar rreqrate: Rate annotation applied to all Route Requests.
    :ivar datarate: Rate annotation applied to all unicast messages.
    :ivar maintbuffer: Dictionary containing buffer of packets being maintained
                       by DSR, indexed by next hop addresses.

    :cvar mac: Property to access pointer to `MAC`.
    :cvar rrt: Alias for `rreqtable`.
    :cvar MaxMaintRexmt: Maximum number of retransmission for route maintenance.
    """
    name = "DSR"
    tracename = "DSR"
    MaxMaintRexmt = 0
    DiscoveryHopLimit = 255
    MaxTTL = DiscoveryHopLimit
    BroadcastJitter = 10e-3
    MaintenanceTimeout = 10.0
    RexmtBufferSize = 50
    MaxGratuitousRexmt = 6
    def __init__(self, *args, **kwargs):
        """Constructor."""
        # set up parameters
        self.__mac = None
        self.rreqrate = None
        self.datarate = None
        self.maintbuffer = None
        # additional parameters for signalling
        self.sndrreq = SimEvent()
        self.drprreq = SimEvent()
        self.finrreq = SimEvent()
        self.sndfail = SimEvent()
        # call base constructor
        Routing.__init__(self, *args, **kwargs)

    mac = property(fget=lambda self: self.get_mac(), \
                   fset=lambda self, m: self.set_mac(m) )
    rrt = property(fget=lambda self: self.rreqtable)
    promiscuous = property(fget=lambda self: self.get_promiscuous() )

    def configure(self, mac=None, rreqrate=None, datarate=None, **kwargs):
        """Configure pointers and parameters.

        :param mac: `MAC` module corresponding to this network protocol.
        :param rreqrate: Rate index used for flooding RREQ messages.
        :param datarate: Rate index used for sending other DSR messages.
        :param kwargs: Additional keywords passed to `Routing.configure()`.

        If `mac` is not provided, this module will attempt to automatically
        determine the appropriate `MAC`.

        If `rreqrate` and/or `datarate` are not specified, `DSR` will not
        attempt to set the rate annotation of the outgoing packet. It is
        important that the `rreqrate` be higher than `datarate` to ensure stable
        links for DSR.
        """
        Routing.configure(self, **kwargs)
        self.mac = mac
        self.rreqrate = rreqrate
        self.datarate = datarate
        self.maintbuffer = {}
        self.table = self.newchild('routingtable', RouteCache, \
                                   name=self.name+".rcache", \
                                   tracename=self.tracename+".RCACHE")
        rreqtable  = self.newchild('rreqtable', RouteRequestTable, \
                                   name=self.name+".rreqtable", \
                                   tracename=self.tracename+".RREQTBL")

    ##############################
    # TX STATES
    ##############################
    def IPFORWARD(self, fsm, p):
        """IPFORWARD state; start route maintenance for IP+DSR if needed.

        :param p: IP packet to send.

        :note: This state assumes a valid IP+DSR packet was went to it.
        """
        # get IP parameters
        ip, dsr = p[IP], p[DSRPacket]
        src, dst = ip.src, ip.dst
        # start route maintenance or deliver
        self.checkdsropt(ip, exception=True, incoming=False)
        # send broadcast immediately
        if (dst==self.broadcast):
            yield fsm.goto(self.IPDELIVER, ip, dst)
        # otherwise ...
        srcroute = getdsr_srcroute(ip)
        if srcroute:
            path = srcroute.addresses
            segsleft, naddr = srcroute.segsleft, len(path)
            nexthop = path[-segsleft]
        else:
            nexthop = dst
        # cache route to nexthop and start route maintenance for ip
        self.addroute(nexthop)
        self.maintain(ip, nexthop)      # put in maintbuffer

    def IPROUTING(self, fsm, p):
        """IPROUTING state; start route discovery for `p`."""
        # error messages
        dropbuffer = "sendbuffer overflow"
        # get IP/DSR/Option parameters
        ip, dsr = p[IP], p[DSRPacket]
        src, dst = ip.src, ip.dst
        # has route -> resend
        if self.hasroute(dst):
            yield fsm.goto(self.IPSEND, ip, src, dst)
        # otherwise -> start route discovery
        rre = self.rrt.getentry(dst)
        drop = rre.buffer(p)
        kwargs = {'src':src, 'dst':dst}
        kwargs['nbuffered'] = len(rre.sendbuffer)
        kwargs['timeleft']  = time2msec(rre.timeleft())
        self.debug("SENDBUFF", p, **kwargs)
        for d in drop:
            self.log_drop(d, drop=dropbuffer)
        # send RREQ
        yield fsm.goto(self.TXRREQ, dst)

    ##############################
    # TX DSR OPTIONS
    ##############################
    def TXRREQ(self, fsm, target, options=[], rexmt=0):
        """TXRREQ state; create and send route request.

        :param target: Target for route discovery.
        :param options: Additional DSR options [default=None].

        :note: This state is persistent. It will keep trying to send until
        """
        # error messages
        rreqerror = "[DSR]: Error getting RREQ Table entry!"
        droprexmt = "max rexmt exceeded"
        drophasrt = "route already exists"
        # pause for jitter
        jitter = random.uniform(0,1)*self.BroadcastJitter
        yield hold, fsm, jitter
        # check route and rexmt count
        if self.hasroute(target):
            self.log("RREQSTOP", target=target, rexmt=rexmt, drop=drophasrt)
            rre = self.rrt.getentry(target)
            while rre.sendbuffer:
                p = rre.sendbuffer.pop(0)
                f = FSM.launch(self.IPRECOVER, p, self.address, target)
            self.rrt.delentry(target)
            # signal that RREQ has finished
            self.finrreq.signal(target)
            yield fsm.stop()    # HALT and stop sending RREQ
        if (rexmt>self.rrt.MaxRequestRexmt):
            self.log("RREQDROP", target=target, rexmt=rexmt, drop=droprexmt)
            rre = self.rrt.getentry(target)
            while rre.sendbuffer:
                p = rre.sendbuffer.pop(0)
                self.log_drop(p, drop=droprexmt)
            self.rrt.delentry(target)
            # signal that RREQ has been abandoned
            self.drprreq.signal(target)
            yield fsm.stop()    # HALT and stop sending RREQ
        # get RREQ parameters
        sendrreq = self.rrt.sendrreq(target)
        rre = self.rrt.getentry(target)
        tleft = rre.timeleft()
        # get parameters for logging
        kwargs = {'rexmt':rexmt, 'nbuffered': len(rre.sendbuffer)}
        kwargs['options']  = [o.tracename for o in options]
        kwargs['jitter'] = time2msec(jitter)
        kwargs['timeleft'] = time2msec(tleft)
        # cannot send RREQ? -> RREQ is busy, drop attempt
        if not sendrreq:
            self.debug("RREQBUSY", target=target, **kwargs)
            yield fsm.stop()    # HALT and allow other RREQ to finish
        # otherwise -> send RREQ
        ID, ttl = rre.ID, rre.ttl
        # create DSR+RREQ+options
        nextheader = self.getproto(None)
        dsr = DSRPacket(nextheader=nextheader)
        rreq = DSROPT_RREQ(identification=ID, target=target)
        dsr.options = [rreq] + [o for o in options]
        # create IP+DSR
        proto = self.getproto(dsr)
        src, dst = self.address, self.broadcast
        ip = IP(src=src, dst=dst, proto=proto, ttl=ttl)
        ip.add_payload(dsr)
        # send RREQ -> wait for timeout, then rexmt
        self.debug("TXRREQ", ip, target=target, **kwargs)
        f = FSM.launch(self.IPDELIVER, ip, dst)
        # signal that RREQ has been sent to target
        self.sndrreq.signal(target)
        # wait for send RREQ to timeout before trying again
        yield hold, fsm, tleft
        yield fsm.goto(self.TXRREQ, target, options, rexmt+1)

    def TXRREP(self, fsm, src, dst, rreq, addresses):
        """TXRREP state; create and send route reply to `dst`."""
        self.debug("TXRREP", src=src, dst=dst, addresses=addresses)
        # cache reverse route from rreq?
        opt = getdsr_rreq(rreq)
        if opt and DSR_USE_REVERSE_ROUTES:
            rpath = [a for a in opt.addresses]
            rpath.reverse()
            if rpath: self.addroute(dst, rpath)     # more than one hop away
            else:     self.addroute(dst)            # one hop neighbor
        # create RREP option
        rrep = DSROPT_RREP()
        rrep.addresses = [a for a in addresses]
        # no route to dst -> send RREQ+RREP
        if not self.hasroute(dst):
            yield fsm.goto(self.TXRREQ, dst, options=[rrep])
        # otherwise -> ...
        if opt and DSR_USE_REVERSE_ROUTES:
            path = [a for a in rpath]
        else:
            path = self.srcroute(dst)
        ttl = len(path)+1
        # create DSR
        nextheader = self.getproto(None)
        dsr = DSRPacket(nextheader=nextheader)
        dsr.options = [rrep]
        # create IP
        proto = self.getproto(dsr)
        ip = IP(src=src, dst=dst, proto=proto, ttl=ttl)
        ip.add_payload(dsr)
        # add SRCROUTE?
        pkt = self.updateroute(ip, path)
        #assert (pkt is not None)
        assert self.safe((pkt is not None))
        yield fsm.goto(self.IPFORWARD, pkt)

    def TXRERR(self, fsm, src, dst, errortype, type_specific):
        """TXRERR state; create and send a new route error message."""
        self.debug("TXRERR", src=src, dst=dst, unreachable=type_specific)
        # erro messages
        errornotme = "[DSR]: Cannot send RERR not from me!"
        assert self.safe(src==self.address), errornotme
        # create RERR option
        rerr = DSROPT_RERR()
        rerr.err_src, rerr.err_dst = src, dst
        rerr.errortype, rerr.type_specific = errortype, type_specific
        # create DSR + RERR Option
        nextheader = self.getproto(None)
        dsr = DSRPacket(nextheader=nextheader)
        dsr.options = [rerr]
        # create IP + DSR
        proto = self.getproto(dsr)
        ip = IP(src=src, dst=dst, proto=proto)
        ip.add_payload(dsr)
        # no route -> start route discovery?
        if not self.hasroute(dst):
            yield fsm.goto(self.IPSEND, ip, src, dst)
        # otherwise -> ...
        path = self.srcroute(dst)
        ttl = len(path)+1
        # add SRCROUTE?
        ip.ttl = ttl
        pkt = self.updateroute(ip, path)
        #assert (pkt is not None)
        assert self.safe((pkt is not None))
        yield fsm.goto(self.IPFORWARD, pkt)

    def MAINT(self, fsm, nexthop, rexmt=0, send=True):
        """MAINT state; perform route maintenace for nexthop.

        :note: Assume maintenance buffer contains valid IP+DSR packets.
        """
        yield hold, fsm, 0  # yield to other threads
        #assert (nexthop in self.maintbuffer)
        assert self.safe((nexthop in self.maintbuffer))
        # error messages
        droproute = "broken route"
        # get maintenance parameters
        buff = self.maintbuffer[nexthop]['buffer']
        if (len(buff)<1):
            del self.maintbuffer[nexthop]
            yield fsm.stop()        # HALT and stop maintenance
        # get head of buffer
        p = buff[0]
        ip, dsr = p[IP], p[DSRPacket]
        addr, src, dst = self.address, ip.src, ip.dst
        # check if path is still valid
        opt, path = getdsr_srcroute(ip), None   # one hop away?
        if opt:
            segsleft = opt.segsleft             # more than one hop?
            path = [a for a in opt.addresses[-segsleft:]]
        pathok = self.hasroute(dst, path)
        # if path is broken -> drop or resend
        if not pathok:
            p = buff.pop(0)
            if (addr==src):
                f = FSM.launch(self.IPRECOVER, p, src, dst)
            else:
                self.log_drop(p, drop=droproute)
            yield fsm.goto(self.MAINT, nexthop)     # continue with next packet
        # otherwise -> send head of buffer?
        if send:
            f = FSM.launch(self.IPDELIVER, ip, nexthop)
        # wait for link-level feedback (ACK/DROP)
        mac, feedback = self.mac, None
        yield waitevent, fsm, (mac.ackdata, mac.drpdata)
        if (mac.ackdata in fsm.eventsFired):
            p = mac.ackdata.signalparam
            if self.issame(ip, p): feedback = "ackdata"
        elif (mac.drpdata in fsm.eventsFired):
            p = mac.drpdata.signalparam
            if self.issame(ip, p): feedback = "drpdata"
        # process feedback
        if (feedback=="ackdata"):
            p = buff.pop(0)
            yield fsm.goto(self.MAINT, nexthop)     # continue with next packet
        elif (feedback=="drpdata"):
            rexmt += 1
            norexmt = (rexmt>self.MaxMaintRexmt)
            # rexmt not exceeded -> try again
            if not norexmt:
                self.debug("REXMT%d"%(rexmt), ip)
                yield fsm.goto(self.MAINT, nexthop, rexmt)
            # otherwise -> broken link!!
            etype = DSR_ERROR_NODE_UNREACHABLE
            esrc, unreachable = self.address, nexthop
            self.debug("DROPDATA", src=esrc, unreachable=unreachable)
            self.removelink(esrc, unreachable)
            # signal broken link
            self.sndfail.signal(esrc)
            # clear out maintenance buffer
            errdst = []
            while buff:
                p = buff.pop(0)
                ip, dsr = p[IP], p[DSRPacket]
                src, dst = ip.src, ip.dst
                # send RERR (for non-RREP messages)
                rrep = getdsr_rrep(ip)
                sendrerr = (not rrep) and (src not in errdst)
                # recover packet or send RERR
                if (addr==src):
                    f = FSM.launch(self.IPRECOVER, p, src, dst)
                elif sendrerr:
                    errdst.append(src)
            # send RERR to sources
            for edst in errdst:
                f = FSM.launch(self.TXRERR, esrc, edst, etype, unreachable)
            # continue to allow graceful shutdown
            yield fsm.goto(self.MAINT, nexthop)
        else:
            # feedback not for me -> continue waiting
            yield fsm.goto(self.MAINT, nexthop, rexmt, send=False)

    def IPRECOVER(self, fsm, p, src, dst):
        """IPRECOVER state; attempt to resend packet `p`.

        :note: This method assumes that `p` is a valid IP+DSR packet.
        """
        self.debug("IPRECOVER", p, src=src, dst=dst)
        # error messages
        errornotme = "[DSR]: Cannot recover message not from me!"
        errorbcast = "[DSR]: Cannot recover broadcast messages!"
        droprcount = "recover attempts exceeded"
        assert self.safe(src==self.address), errornotme
        assert self.safe(dst!=self.broadcast), errorbcast
        # drop RREP or RERR packets
        rrep = getdsr_rrep(p)
        rerr = getdsr_rerr(p)
        drop = ""
        if rrep: drop = "do not recover RREP"
        if rerr: drop = "do not recover RERR"
        if drop:
            self.log_drop(p, src=src, dst=dst, drop=drop)
            yield fsm.stop()            # HALT and discard packet
        # keep secret annotation for keeping track of recovery attempts
        rcount = 0
        if p.hasanno('iprecover'): rcount = p.getanno('iprecover') + 1
        p.setanno('iprecover', rcount, priv=True)
        # max retries exceeded?
        if (rcount>1): # FIXME: use (>self.MaxMaintRexmt) instead?
            self.log_drop(p, src=src, dst=dst, rcount=rcount, drop=droprcount)
            yield fsm.stop()
        # otherwise -> remove stale route, update, and resend
        ip, dsr = p[IP], p[DSRPacket]
        dsr.options = [o for o in dsr.options if (not getdsr_srcroute(o))]
        assert (dsr.nextheader!=const.IP_PROTO_NONE)
        yield fsm.goto(self.IPSEND, ip, src, dst)

    ##############################
    # RX STATES
    ##############################
    def IPCLASSIFY(self, fsm, p, cache=True):
        """IPCLASSIFY state; overloaded to implement DSR classification.

        :param p: Received IP packet.
        :param cache: Cache routing information from `p`.

        :note: This state assumes that packet `p` passed all validity checks in
               `IPRECV` (i.e. `checkiprecv()`).
        """
        errmsg = "[DSR]: IPCLASSIFY does not support promiscuous mode."
        #assert (not self.promiscuous), errmsg
        assert self.safe((not self.promiscuous)), errmsg
        # error messages
        errordsropt  = "[DSR]: Unprocessed DSR options still remain!"
        dropnotforme = "not for me"
        # get IP+DSR packets
        ip = p[IP]
        dsr = ip[DSRPacket]
        # get IP/DSR parameters
        addr, bcast = self.address, self.broadcast
        src, dst = ip.src, ip.dst
        # get DSR options
        rreq, rrep     = getdsr_rreq(ip), getdsr_rrep(ip)
        rerr, srcroute = getdsr_rerr(ip), getdsr_srcroute(ip)
        isforme, isbcast = (dst==addr), (dst==bcast)
        unicast, promiscuous = (not isbcast), self.promiscuous
        # cache routing info from packet
        if cache: self.cacheroute(ip)
        # classify DSR options
        if srcroute and unicast:
            yield fsm.goto(self.RXSRCROUTE, ip, dsr, srcroute)
        elif rreq and isbcast:
            yield fsm.goto(self.RXRREQ,     ip, dsr, rreq)
        elif rrep and isforme:
            yield fsm.goto(self.RXRREP,     ip, dsr, rrep)
        elif rerr and isforme:
            yield fsm.goto(self.RXRERR,     ip, dsr, rerr)
        elif dsr.options:
            #assert (self.promiscuous or (not isforme)), errordsropt
            assert self.safe((self.promiscuous or (not isforme))), errordsropt
        # classify payload of packet
        if (dst==addr) or (dst==bcast):
            # packet for me -> detach payload
            payload = dsr.payload
            dsr.remove_payload()
            self.log_recv(payload, src=src, dst=dst)
            if (dsr.nextheader==const.IP_PROTO_IP):
                yield fsm.goto(self.IPRECV, payload)    # IPENCAP -> reclassify
            elif (dsr.nextheader!=const.IP_PROTO_NONE):
                pkt = self.set_recvanno(payload, src, dst)
                yield self.TXU.send(fsm, [pkt])         # DATA -> send to upper layer
        else:
            # otherwise -> drop
            self.log_drop(ip, drop=dropnotforme)

    ##############################
    # RX DSR OPTIONS
    ##############################
    def RXSRCROUTE(self, fsm, ip, dsr, opt):
        """RXSRCROUTE state; process source route option.

        :note: Assumes `checkiprecv()` passed.
        """
        self.debug("RXSRCROUTE", ip)
        # get IP/DSR/Option parameters
        segsleft = opt.segsleft
        path = opt.addresses[-segsleft]
        # forward packet to next hop in source route
        if (segsleft>1):
            # intermediate hop -> update SRCROUTE and send
            ip.ttl = ip.ttl - 1
            opt.segsleft = segsleft - 1
            yield fsm.goto(self.IPFORWARD, ip)
        elif (segsleft==1):
            # last hop -> strip SRCROUTE and send
            ip.ttl = ip.ttl - 1
            dsr.options = [o for o in dsr.options if (not getdsr_srcroute(o))]
            yield fsm.goto(self.IPFORWARD, ip)

    def RXRREQ(self, fsm, ip, dsr, opt):
        """RXRREQ state; process route request option.

        :note: Assumes `checkiprecv()` passed.
        """
        self.debug("RXRREQ", ip)
        # error messages
        droprreq = "in source route or not new RREQ"
        # get IP/DSR/Option parameters
        addr, src, dst = self.address, ip.src, ip.dst
        ID, target = opt.identification, opt.target
        rreqonly = (dsr.nextheader==const.IP_PROTO_NONE)
        # check if target reached
        if (addr==target):
            # send RREP [and process remaining packet]
            path = [a for a in opt.addresses] + [target]
            f = FSM.launch(self.TXRREP, addr, src, opt, path)
            # strip RREQ and set dst=target
            dsr.options = [o for o in dsr.options if (not getdsr_rreq(o))]
            ip.dst = target
            # check for other options
            rerr, rrep = getdsr_rerr(dsr), getdsr_rrep(dsr)
            if rerr:
                yield fsm.goto(self.RXRERR, ip, dsr, rerr)
            elif rrep:
                yield fsm.goto(self.RXRREP, ip, dsr, rrep)
            # RREQ only -> HALT and discard packet
            if rreqonly: yield fsm.stop()
            # otherwise -> remove DSR options and reclassify
            dsr.options = []
            yield fsm.goto(self.IRECV, ip)
        # otherwise -> attempt to forward RREQ
        assert (addr!=target)
        inroute = (src==addr) or (addr in opt.addresses)
        newrreq = self.rrt.addcache(src, ID, target)
        # check if RREQ should be dropped
        if (inroute or (not newrreq)):
            self.log_drop(ip, drop=droprreq)
            yield fsm.stop()    # HALT and discard packet
        # update options and forward RREQ
        ip.ttl = ip.ttl - 1
        opt.addresses = [a for a in opt.addresses] + [addr]
        yield fsm.goto(self.IPFORWARD, ip)

    def RXRREP(self, fsm, ip, dsr, opt):
        """RXRREP state; handle route reply message.

        :note: Assumes `checkiprecv()` passed.
        """
        self.debug("RXRREP", ip)
        # get IP/DSR/Option parameters
        src, dst = ip.src, ip.dst
        addr, target = self.address, opt.addresses[-1]
        rreponly = (dsr.nextheader==const.IP_PROTO_NONE)
        # only keep RREP options and update route cache
        assert (dst==addr)
        dsr.options = [o for o in dsr.options if (getdsr_rrep(o))]
        self.cacheroute(ip)
        #assert (self.hasroute(target))
        assert self.safe((self.hasroute(target)))   # verify route to target
        # get entry from Route Request Table and forward packets
        rre = self.rrt.getentry(target)
        while rre.sendbuffer:
            p = rre.sendbuffer.pop(0)
            f = FSM.launch(self.IPRECOVER, p, addr, target)
        self.rrt.delentry(target)
        # more than RREP?
        if not rreponly:
            # remove DSR options and reclassify
            dsr.options = []
            yield fsm.goto(self.IRECV, ip)

    def RXRERR(self, fsm, ip, dsr, opt):
        """RXRERR state; handle route error message.

        :note: Assumes `checkiprecv()` passed.
        """
        self.debug("RXRRER", ip)
        # error messages
        droptype = "unsupported route error type"
        # get IP/DSR/Option parameters
        addr, src, dst = self.address, ip.src, ip.dst
        errsrc, errdst = opt.err_src, opt.err_dst
        errortype, unreachable = opt.errortype, opt.type_specific
        rerronly = (dsr.nextheader==const.IP_PROTO_NONE)
        #assert (addr==dst)
        assert self.safe((addr==dst))
        # process RERR packet (or drop unsupported errortype)
        if (errortype==DSR_ERROR_NODE_UNREACHABLE):
            self.removelink(errsrc, unreachable)
        else:
            self.log_drop(ip, drop=droptype)
        # more than RERR?
        if not rerronly:
            # remove DSR options and reclassify
            dsr.options = []
            yield fsm.goto(self.IRECV, ip)

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
        isdsr = isip and p.haslayer(DSRPacket)
        if isdsr and (not force):
            ip, dsr = p[IP], p[DSRPacket]
            dsr.nextheader = self.getproto(dsr.payload)
            if self.hasroute(ip.dst):
                p = self.updateroute(ip)
            # update IP parameters
            return Routing.encap_data(self, p, src, dst, force=False, **kwargs)
        # create DSR
        nextheader = self.getproto(p)
        dsr = DSRPacket(nextheader=nextheader)
        dsr.add_payload(p)
        # create IP+DSR
        proto = self.getproto(dsr)
        ip = IP(src=src, dst=dst, proto=proto, **kwargs)
        ip.add_payload(dsr)
        # update and return
        return self.encap_data(ip, src, dst, force=False, **kwargs)

    def maintain(self, p, nexthop):
        """Start or udpate route maintenace on IP+DSR packet.

        :note: This assumes a IP+DSR packet was passed to it.
        """
        # error messages
        errornexthop = "[DSR]: maintain() found invalid nexthop!"
        dropbuffer = "maintenance buffer overflow"
        dropttl = "TTL expired"
        # get IP/DSR/Option parameters
        ip, dsr = p[IP], p[DSRPacket]
        addr, src, dst, ttl = self.address, ip.src, ip.dst, ip.ttl
        # TTL expired?
        if (ttl<1):
            self.log_drop(p, drop=dropttl)
            return
        # start or continue route maintenance
        if nexthop in self.maintbuffer:
            # add packet to end of maintenance buffer
            buff = self.maintbuffer[nexthop]['buffer']
            if (len(buff)<self.RexmtBufferSize):
                buff.append(ip)
                self.debug("MAINTBUFF", ip)
            else:
                self.log_drop(ip, drop=dropbuffer)
        else:
            # add new buffer and launch thread to maintain it
            self.maintbuffer[nexthop] = {'buffer':[ip]}
            f = FSM.launch(self.MAINT, nexthop)

    def cacheroute(self, p):
        """Cache relevant route information from IP+DSR packet.

        :note: This method assumes that the packet has passed through
               `checkiprecv()`.
        """
        errmsg = "[DSR]: cacheroute() does not support promiscuous mode."
        #assert (not self.promiscuous), errmsg
        assert self.safe((not self.promiscuous)), errmsg
        # get ip/dsr packets and parameters
        ip = p[IP]
        dsr = p[DSRPacket]
        addr, bcast = self.address, self.broadcast
        src, dst, ttl = ip.src, ip.dst, ip.ttl
        # classify based on DSR option
        rreq, rrep     = getdsr_rreq(ip), getdsr_rrep(ip)
        rerr, srcroute = getdsr_rerr(ip), getdsr_srcroute(ip)
        # process source route option
        if srcroute:
            # cache remaining route to dst
            segsleft = srcroute.segsleft
            path = srcroute.addresses[-segsleft:]
            if (segsleft>1):     self.addroute(dst, path[1:])
            elif (segsleft==1): self.addroute(dst)
        # process route request option
        if rreq:
            # no routing info to cache, but use target as dst for other options
            dst = rreq.target
        # process route reply option
        if rrep:
            # check if on route or intended dst of RREP
            inroute = (addr in rrep.addresses[:-1])
            isforme = (dst==addr)
            # cache routing info
            idx = None
            if inroute:   idx = rrep.addresses.index(addr)+1
            elif isforme: idx = 0
            if (idx is not None):
                target = rrep.addresses[-1]
                rt = rrep.addresses[idx:-1]
                if rt: self.addroute(target, rt)    # more than one hop away
                else:  self.addroute(target)        # target is neighbor
        # process route error option
        if rerr:
            errsrc, unreachable = rerr.err_src, rerr.type_specific
            if (rerr.errortype==DSR_ERROR_NODE_UNREACHABLE):
                self.removelink(errsrc, unreachable)

    def checkiprecv(self, p, nextheader=True, **kwargs):
        """Overloaded to check for valid IP+DSR packet.

        :param p: Packet received in `IPRECV`.
        :param args: Additional arguments passed to base class method.
        :param kwargs: Additional keywords passed to base class method.
        :return: String containing any error condition, or nothing for no error.
        """
        # error messages
        dropnondsr = "non-IP+DSR packet in IPRECV"
        dropheader = "invalid IP+DSR nextheader"
        dropdsropt = "invalid DSR options!"
        # call base method
        drop = Routing.checkiprecv(self, p, **kwargs)
        if drop: return drop
        # valid IP+DSR?
        isdsr = p.haslayer(IP) and p.haslayer(DSRPacket)
        if not isdsr: return dropnondsr
        dsr = p[DSRPacket]
        # check nextheader
        if nextheader:
            payproto = self.getproto(dsr.payload)
            if (payproto!=dsr.nextheader): return dropheader
        # check DSR options
        if not self.checkdsropt(p, exception=False):
            return dropdsropt
        return ""   # no error

    def checkdsropt(self, p, exception=True, incoming=True):
        """Check DSR options of incoming packet for any violations.

        :param p: Packet to check.
        :param exception: If true, throws exception for violation.
        :param incoming: If true, treat `p` as a received packet.
        :return: Boolean flag; if True, DSR packet is okay.

        :note: Do not use this method to check options for outgoing packets
               without setting `incoming` to False.
        """
        if isinstance(p, Reference): p = p._deref
        isip = isinstance(p, Packet) and p.haslayer(IP)
        isdsr = isip and p.haslayer(DSRPacket)
        if not isdsr:
            if exception: raise RuntimeError, "[DSR]: non-IP+DSR packet!"
            return False
        # get IP+DSR packets
        ip  = p[IP]
        dsr = ip[DSRPacket]
        # get IP/DSR parameters
        addr, bcast = self.address, self.broadcast
        src, dst = ip.src, ip.dst
        # classify DSR options
        rreq, rrep     = getdsr_rreq(ip), getdsr_rrep(ip)
        rerr, srcroute = getdsr_rerr(ip), getdsr_srcroute(ip)
        isforme, isbcast = (dst==addr), (dst==bcast)
        unicast, promiscuous = (not isbcast), self.promiscuous
        # error check options
        dsrok = True
        if srcroute:
            path = srcroute.addresses
            segsleft, naddr = srcroute.segsleft, len(path)
            dsrok &= not rreq
            dsrok &= not isbcast
            dsrok &= (segsleft>0) and (segsleft<naddr+1)
            dsrok &= not (isforme and (not promiscuous))
            if dsrok and (not promiscuous):
                if incoming:
                    # incoming -> I am current hop?
                    currhop = path[-segsleft]
                    dsrok &= (currhop==addr)
                else:
                    # outgoing -> I am right before nexthop?
                    currhop = path[-segsleft]
                    if (currhop==path[0]):
                        dsrok &= segsleft==naddr
                        dsrok &= (addr==src)
                    else:
                        dsrok &= (segsleft<naddr)
                        if dsrok: dsrok &= (addr==path[-(segsleft+1)])
            #if exception: assert dsrok, "[DSR]: Invalid SRCROUTE option!"
            if exception: assert self.safe(dsrok), "[DSR]: Invalid SRCROUTE option!"
        if rreq:
            dsrok &= not srcroute
            dsrok &= not unicast
            #if exception: assert dsrok, "[DSR]: Invalid RREQ option!"
            if exception: assert self.safe(dsrok), "[DSR]: Invalid RREQ option!"
        if rrep:
            dsrok &= not (isbcast and (not rreq))
            dsrok &= not rerr
            dsrok &= (len(rrep.addresses)>0)
            #if exception: assert dsrok, "[DSR]: Invalid RREP option!"
            if exception: assert self.safe(dsrok), "[DSR]: Invalid RREP option!"
        if rerr:
            dsrok &= not (isbcast and (not rreq))
            dsrok &= not rrep
            dsrok &= (rerr.errortype==DSR_ERROR_NODE_UNREACHABLE)
            #if exception: assert dsrok, "[DSR]: Invalid RERR option!"
            if exception: assert self.safe(dsrok), "[DSR]: Invalid RERR option!"
        return dsrok

    def updateroute(self, p, path=None, chksum=True):
        """Update IP+DSR routing info using a particular route.

        :param p: IP+DSR packet.
        :param path: New route; if `None`, get it from route cache.
        :param chksum: If true, update checksum.
        :return: Packet with updated routing info (or `None` if no valid route).
        """
        ip, dsr = p[IP], p[DSRPacket]
        # get IP/DSR/Option parameters
        addr, src, dst = self.address, ip.src, ip.dst
        if path is None: path = self.srcroute(dst)
        # no path found?
        if path is None: return None
        # strip old SRCROUTE
        dsr.options = [o for o in dsr.options if (not getdsr_srcroute(o))]
        # apply new route
        ip.ttl = len(path) + 1
        if path:
            srcroute = DSROPT_SRCROUTE()
            srcroute.segsleft = len(path)
            srcroute.addresses = [a for a in path]
            dsr.options = [srcroute] + [o for o in dsr.options]
        # update chksum?
        if chksum: chksum = self.updatechksum(ip, overwrite=True)
        return ip

    def set_sendanno(self, *args, **kwargs):
        """Overloaded to set rate annotations in outgoing packets."""
        p = Routing.set_sendanno(self, *args, **kwargs)
        isip = isinstance(p, Packet) and p.haslayer(IP)
        isdsr = isip and p.haslayer(DSRPacket)
        if not isdsr: return p
        # remove stale annotations
        if p.hasanno('phy-fixed-rate'): p.delanno('phy-fixed-rate')
        # set rate annotation
        rreq = getdsr_rreq(p)
        rate = None
        if rreq: rate = self.rreqrate   # assume RREQ rate is higher
        else:    rate = self.datarate
        if (rate is not None):
            p.setanno('phy-rate', rate)
            # use fixed-rate anno for broadcast packets
            if rreq: p.setanno('phy-fixed-rate', rate)
        return p

    def getproto(self, p):
        """Overloaded to check for DSR packets."""
        if isinstance(p, Reference): p = p._deref
        # determine protocol type of packet
        if isinstance(p, DSRPacket): proto = const.IP_PROTO_DSR
        else:                        proto = Routing.getproto(self, p)
        return proto

    def issame(self, pa, pb):
        """See if packets have the same IP+DSR header information."""
        isip  = isinstance(pa,Packet) and pa.haslayer(IP)
        isdsr = isip and pa.haslayer(DSRPacket)
        if not isdsr: return False
        isip  = isinstance(pb,Packet) and pb.haslayer(IP)
        isdsr = isip and pb.haslayer(DSRPacket)
        if not isdsr: return False
        ipa, dsra = pa[IP], pa[DSRPacket]
        ipb, dsrb = pb[IP], pb[DSRPacket]
        # check IP/DSR fields
        same = True
        same &= (ipa.src==ipb.src) and (ipa.dst==ipb.dst)
        same &= (ipa.id==ipb.id) and (ipa.proto==ipb.proto)
        if (ipa.chksum is not None) and (ipb.chksum is not None):
            same &= (ipa.chksum==ipb.chksum)
        same &= (dsra.length==dsrb.length)
        same &= (dsra.nextheader==dsrb.nextheader)
        same &= (len(dsra.options)==len(dsrb.options))
        return same

    ##############################
    # ROUTE TABLE METHODS
    ##############################
    def srcroute(self, *args, **kwargs):
        """Get source route from route cache."""
        return self.table.srcroute(*args, **kwargs)

    def removelink(self, a, b):
        """Remove link from route cache."""
        return self.table.removelink(a, b, self.address)

    ##############################
    # PROPERTY METHODS
    ##############################
    def get_mac(self):
        """Get MAC corresponding to this module."""
        if isinstance(self.__mac, MAC): return self.__mac
        # otherwise try to get MAC from ARP
        mac, arp = None, self.TXD.target.parent
        if isinstance(arp, ARP):
            if self in arp.listen: mac, f = arp.listen[self]
        # final check on MAC
        if not isinstance(mac, MAC): mac = None
        return mac

    def set_mac(self, mac):
        """Set MAC corresponding to this module."""
        m = None
        if isinstance(mac, MAC): m = mac
        self.__mac = m

    def get_promiscuous(self):
        """Get `MAC.promiscuous` flag."""
        p = None
        if self.mac: p = self.mac.promiscuous
        return p

    ##############################
    # LOGGING METHODS
    ##############################
    def get_ip_anno(self, p):
        """Internal method to get relevant annotations for an IP+DSR packet."""
        kwargs = Routing.get_ip_anno(self, p)
        kwargs.update(self.get_dsr_anno(p) )
        kwargs.update(self.get_agt_anno(p) )
        return kwargs

    def get_agt_anno(self, p):
        """Internal method to get relevant annotations for a AGT packet."""
        kwargs = {}
        isagt = isinstance(p, Packet) and p.haslayer(AGT)
        if not isagt: return kwargs
        agt = p[AGT]
        kwargs['agt-root'] = agt.traceid
        return kwargs

    def get_dsr_anno(self, p):
        """Internal method to get relevant annotations for a DSR packet."""
        kwargs = {}
        isdsr = isinstance(p, Packet) and p.haslayer(DSRPacket)
        if not isdsr: return kwargs
        dsr = p[DSRPacket]
        kwargs['dsr.nextheader'] = dsr.nextheader
        # classify based on DSR option
        rreq, rrep     = getdsr_rreq(dsr), getdsr_rrep(dsr)
        rerr, srcroute = getdsr_rerr(dsr), getdsr_srcroute(dsr)
        # get parameters from options
        if rreq:
            kwargs['rreq.ID'] = rreq.identification
            kwargs['rreq.target'] = rreq.target
            kwargs['rreq.addresses'] = rreq.addresses
        if rrep:
            kwargs['rrep.addresses'] = rrep.addresses
        if srcroute:
            kwargs['srcroute.segsleft'] = srcroute.segsleft
            kwargs['srcroute.addresses'] = srcroute.addresses
        if rerr:
            kwargs['rerr.errortype'] = _dsr_error_types[rerr.errortype]
            kwargs['rerr.err_src'] = rerr.err_src
            kwargs['rerr.err_dst'] = rerr.err_dst
            kwargs['rerr.type_specific'] = rerr.type_specific
        # get other parameters
        if p.hasanno('phy-rate'):
            kwargs['phy-rate'] = p.getanno('phy-rate')
        return kwargs

    def debug(self, *args, **kwargs):
        """Log debug statements to trace."""
        if DSRDEBUG:
            self.log(*args, **kwargs)

    def safe(self, x):
        """Print trace to standard out before throwing exception."""
        if not x:
            self.trace.output()
        return x
