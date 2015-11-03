#!  /usr/bin/env python

"""
Distributed Control Function of IEEE 802.11; contains `DCF` class.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-10-22 18:58:45 -0500 (Sat, 22 Oct 2011) $
* $LastChangedRevision: 5271 $

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

from SimPy.Simulation import SimEvent, waitevent, hold, now
from wins.protocol.csmac import CSMAC
from wins.packet import ANNO, Packet
from wins.helper import time2usec
from wins.fsm import FSM, Timer
from wins.crc import crcupdate, crcremove
from wins import const

from wins.ieee80211.dot11_support import *
from wins.ieee80211.dot11a import Dot11APHY
from wins.ieee80211.dot11n import Dot11NPHY
from wins.trace import Traceable

from scapy.all import Ether
from numpy import inf, random

from wins.ieee80211.navtimer import NAVTimer

class DCF(CSMAC):
    """Distributed Control Function of IEEE 802.11 MAC protocol.

    DCF and Ports
    =============
    The DCF module has two downstream port ('TXD' and 'RXD'). These ports are
    used to interact with the associated physical layer.

    CSMA/CA
    =======
    Carrier sense multiple access with collision avoidance (CSMA/CA) may be used
    for the IEEE 802.11 MAC protocol. This mode of operation will not use
    RTS/CTS reservation messages. Instead, after carrier sense backoff, the
    protocol will immediately send the data message. Use the boolean flag
    `usecsma` to enable this mode of operation.

    :CVariables:
     * `pifs`: Property to access point-coordination interframe spacing.
     * `difs`: Property to access distributed interframe spacing.
     * `acktimeout`: Property to access ACK timeout.
     * `ctstimeout`: Property to access CTS timeout.
     * `ackduration`: Property to access duration of ACK.
     * `ctsduration`: Property to access duration of CTS.

    :IVariables:
     * `sifs`: Short interframe spacing.
     * `slottime`: Duration of a contention slot.
     * `cwmin`: Minimum contention window size.
     * `cwmax`: Maximum contention window size.
     * `cslot`: Contention slot value.
     * `datatosend`: Packet currently being handled; `None` indicates not busy.
     * `retrycount`: Retry counter.
     * `retrylimit`: Maximum number of retries allowed.
     * `usecsma`: Boolean flag; if true, use CSMA/CA.
     * `rxdata`: SimEvent signalled when a packet arrives on `Port` 'RXD'.
     * `_ctsduration`: Internal member used to cache duration of CTS.
     * `_ackduration`: Internal member used to cache duration of ACK.

    :note: `DCF` only handles Ethernet packets. (i.e. `htype` must be
           `const.ARP_HTYPE_ETHERNET`).
    """
    name = "distributed coordination function"
    tracename = "DCF"
    cwmin = DOT11_CWMIN
    cwmax = DOT11_CWMAX
    retrylimit = DOT11_RETRYLIMIT
    def __init__(self, cwmin=None, cwmax=None, retrylimit=None, \
                       usecsma=False, **kwargs):
        """Constructor.

        :param cwmin: Minimum contention window size.
        :param cwmax: Maximum contention window size.
        :param retrylimit: Maximum number of retries allowed.
        :param usecsma: Boolean flag; if true, use CSMA/CA without RTS-CTS
                        reservation messages.
        :param kwargs: Additional keywords passed to `configure()`.

        The default parameters are specified by the class.
        """
        if cwmin is None: cwmin = self.__class__.cwmin
        if cwmax is None: cwmax = self.__class__.cwmax
        if retrylimit is None: retrylimit = self.__class__.retrylimit
        # timing parameters
        self.sifs, self.slottime = None, None
        self.cwmin, self.cwmax   = cwmin, cwmax
        self.cslot = None
        # events and other members
        self.datatosend = None
        self.retrycount = None
        self.retrylimit = retrylimit
        self.usecsma = usecsma
        self.rxdata = SimEvent(name=".rxdata")
        self._ctsduration = None
        self._ackduration = None
        # call CSMAC constructor
        CSMAC.__init__(self, **kwargs)
        self.rxdata.name = "%s.rxdata"%(self.name)

    pifs = property(fget=lambda self: self.sifs + self.slottime)
    difs = property(fget=lambda self: self.pifs + self.slottime)
    ctstimeout = property(fget=lambda self: self.get_ctstimeout() )
    acktimeout = property(fget=lambda self: self.get_acktimeout() )
    navbusy = property(fget=lambda self: self.nav.running)
    navidle = property(fget=lambda self: not self.navbusy)
    ctsduration = property(fget=lambda self: self.get_ctsduration() )
    ackduration = property(fget=lambda self: self.get_ackduration() )

    def configure(self, **kwargs):
        """Configure downstream ports and `FSM`."""
        CSMAC.configure(self, **kwargs)
        # add downstream and control ports
        self.addport("TXD", tracename=self.tracename+".TXD")
        self.addport("RXD", tracename=self.tracename+".RXD")
        self.addport("TXC", tracename=self.tracename+".TXC") # control port
        # create FSM to manage send/recv execution of DCF
        txfsm = self.newchild("txfsm", FSM, tracename=self.tracename+".TX")
        rxfsm = self.newchild("rxfsm", FSM, tracename=self.tracename+".RX")
        txfsm.goto(self.IDLE)
        rxfsm.goto(self.RECV)
        # set up timing parameters
        self.set_timing()
        nav = self.newchild("nav", NAVTimer, tracename=self.tracename+".NAV")

    def encapsulate(self, p, src=None, dst=None, **kwargs):
        """Convenience method to encapsulate an packet in an IEEE 802.11 data
        header (i.e. `Dot11Data`).

        :param p: Packet to encapsulate.
        :param src: Source address [default=`address`]
        :param dst: Destination address [default=`broadcast`]
        :param kwargs: Additional keywords passed to `Dot11Data` constructor.
        :return: Newly created `Dot11Data` packet.

        :note: This method adds/updates a CRC using `crcupdate()`.
        """
        if src is None: src = self.address
        if dst is None: dst = self.broadcast
        addr1, addr2 = dst, src
        pargs = {'addr1': addr1, 'addr2': addr2}
        kwargs.update(pargs)
        data = self.dot11data(**kwargs)
        data.add_payload(p)
        pkt = crcupdate(data)
        return pkt

    def retry(self, count=None):
        """Update retry count and set backoff parameters.

        :param count: If specified, `retrycount` will be set to this value;
                      otherwise, increment `retrycount`.

        Slot value `cslot` is set to a random integer between [0, CW), where
        the contention window CW is defined as:

            CW = CWmin + 2^`retrycount`
        """
        if count is None: count = self.retrycount + 1
        cwsize = min(self.cwmax, self.cwmin * (2**count) )
        self.cslot = random.randint(cwsize)
        self.retrycount = count
        # update RETRY flag if first retry
        if ((count==1) and self.isdot11data(self.datatosend)):
            pkt = self.get_dot11data(self.datatosend)
            dst = pkt.addr1
            errmsg = "[DCF]: Cannot retry() with broadcast packet!"
            assert (dst != self.broadcast), errmsg
            # update retry field
            pkt.FCfield |= DOT11_FC_RETRY
            self.datatosend = crcupdate(pkt)
        if count>0:
            self.log("retry%d"%(count), self.datatosend, retrycount=count, retrylimit=self.retrylimit)

    def IDLE(self, fsm):
        """IDLE state; reset parameters and check for new data from 'RXU'."""
        assert (self.htype==const.ARP_HTYPE_ETHERNET), \
                "[DCF]: Unsupported hardware type (%s)!"%(self.htype)
        # reset parameters and check RXU
        self.cslot = None
        self.datatosend = None
        self.retrycount = self.retrylimit + 1
        # csbusy -> go to RXBUSY
        if self.isbusy: yield fsm.goto(self.RXBUSY)
        yield self.RXU.recv(fsm, 1, \
                renege=(self.csbusy, self.csidle, self.rxdata) )
        # RXU -> new data to transmit -> go to TXDATA
        if fsm.acquired(self.RXU):
            assert (len(fsm.got)==1), \
                    "[DCF]: Error receiving from RXU in IDLE!"
            p = fsm.got[0]
            yield fsm.goto(self.TXDATA, p)
        # csbusy -> go to RXBUSY
        if (self.csbusy in fsm.eventsFired):
            p = self.csbusy.signalparam
            yield fsm.goto(self.RXBUSY)
        # rxdata -> go to RXPKT
        if (self.rxdata in fsm.eventsFired):
            p = self.rxdata.signalparam
            yield fsm.goto(self.RXPKT, p)
        # otherwise -> ignore
        ignore = self.csidle in fsm.eventsFired
        # continue in IDLE
        yield fsm.goto(self.IDLE)

    def TXDATA(self, fsm, pkt):
        """TXDATA state; initialize `datatosend` and associated parameters
        before transitioning to `BACKOFF`."""
        assert (self.datatosend is None), \
                "[DCF]: 'datatosend' already set in TXDATA!"
        assert (self.htype==const.ARP_HTYPE_ETHERNET), \
                "[DCF]: Unsupported hardware type (%s)!"%(self.htype)
        assert isinstance(pkt, Ether) and pkt.haslayer(Ether), \
                "[DCF]: Got non-Ether packet in TXDATA!"
        # process Ethernet frame
        eth = pkt[Ether]
        addr, src, dst = self.addr, eth.src, eth.dst
        pkt = crcupdate(eth)
        # initialize datatosend and other parameters
        self.datatosend = self.encapsulate(pkt, src=src, dst=dst)
        isbroadcast = (dst==self.broadcast)
        self.retry(count=0)
        if isbroadcast: self.retrycount = self.retrylimit
        self.datatosend.setanno('mac-txts', now())
        # go to BACKOFF
        yield fsm.goto(self.BACKOFF)

    def BACKOFF(self, fsm):
        """BACKOFF state; perform backoff operation."""
        assert self.isdot11data(self.datatosend), \
                "[DCF]: Cannot determine 'datatosend' in BACKOFF!"
        # retry limit exceeded -> DROP PACKET! -> go to IDLE
        if self.retrycount>self.retrylimit:
            self.log_drop(self.datatosend, drop="retry limit exceeded")
            pkt = self.datatosend.payload
            self.datatosend.remove_payload()
            p = crcremove(pkt)
            self.drpdata.signal(p)
            yield fsm.goto(self.IDLE)
        # csbusy -> go to RXBUSY
        if self.isbusy:
            yield fsm.goto(self.RXBUSY)
        # check for nav timer -> start nav backoff
        if self.navbusy: yield fsm.goto(self.NAVBACKOFF)
        # start backoff timer
        backoff = self.difs + self.cslot*self.slottime
        timer = self.newchild("backofftimer", Timer, backoff, start=True, \
                              tracename=self.tracename+".BACKOFF")
        self.log("BACKOFF", self.datatosend, backoff=time2usec(backoff), cslot=self.cslot)
        yield waitevent, fsm, (timer.done, timer.kill, self.csbusy, self.rxdata)
        csbusy = (self.csbusy in fsm.eventsFired)
        rxdata = (self.rxdata in fsm.eventsFired)
        # timer done -> go to TXRTS or TXBCAST
        if (timer.done in fsm.eventsFired):
            isbroadcast = (self.datatosend.addr1==self.broadcast)
            if isbroadcast:    ns = self.TXBCAST
            elif self.usecsma: ns = self.TXUCAST
            else:              ns = self.TXRTS
            yield fsm.goto(ns)
        # timer kill -> raise exception
        elif (timer.kill in fsm.eventsFired):
            raise RuntimeError, "[DCF]: Unexpected kill signal " + \
                                "from timer in BACKOFF!"
        # csbusy/rxdata -> halt timer -> update cslot -> go to RXBUSY/RXPKT
        elif csbusy or rxdata:
            yield timer.pause(fsm)
            if timer.timepassed>self.difs:
                rslot = int(timer.timeleft/self.slottime)
                self.cslot = rslot      # update cslot
            timer.halt()
            if rxdata:
                p = self.rxdata.signalparam
                yield fsm.goto(self.RXPKT, p)
            else:
                yield fsm.goto(self.RXBUSY)
        # otherwise -> raise error!
        else:
            raise RuntimeError, "[DCF]: Unexpected interruption in BACKOFF!"

    def NAVBACKOFF(self, fsm):
        """NAVBACKOFF state; defer access for NAV and virtual carrier sense."""
        nav = self.nav
        assert nav.running, "[DCF]: NAV timer not running in NAVBACKOFF!"
        yield waitevent, fsm, (nav.done, nav.kill, self.csbusy, self.rxdata)
        # nav done -> go to RESUME
        if (nav.done in fsm.eventsFired):
            yield fsm.goto(self.RESUME)
        # nav kill -> raise exception
        elif (nav.kill in fsm.eventsFired):
            raise RuntimeError, "[DCF]: Unexpected kill signal " + \
                                "from NAV timer in NAVBACKOFF!"
        # csbusy -> go to RXBUSY
        elif (self.csbusy in fsm.eventsFired):
            p = self.csbusy.signalparam
            yield fsm.goto(self.RXBUSY)
        # rxdata -> go to RXPKT
        elif (self.rxdata in fsm.eventsFired):
            p = self.rxdata.signalparam
            yield fsm.goto(self.RXPKT, p)
        # otherwise -> raise error!
        else:
            raise RuntimeError, "[DCF]: Interrupted in NAVBACKOFF!"
        # go to RESUME
        yield fsm.goto(self.RESUME)

    def TXBCAST(self, fsm):
        """TXBCAST state; broadcast `datatosend`."""
        assert self.isdot11data(self.datatosend), \
                "[DCF]: Cannot determine 'datatosend' in TXBCAST!"
        assert (self.datatosend.addr1==self.broadcast), \
                "[DCF]: Non-broadcast 'datatosend' in TXBCAST!"
        data = self.datatosend
        self.send_data(data)
        # send and hold for duration
        duration = self.duration(data)
        src, dst = data.addr2, data.addr1
        self.log_send(data, src=src, dst=dst, duration=time2usec(duration) )
        yield self.TXD.send(fsm, [data])
        yield hold, fsm, duration
        # go back to IDLE
        yield fsm.goto(self.IDLE)

    def TXRTS(self, fsm):
        """TXRTS state; send RTS for `datatosend`."""
        assert self.isdot11data(self.datatosend), \
                "[DCF]: Cannot determine 'datatosend' in TXRTS!"
        assert not (self.datatosend.addr1==self.broadcast), \
                "[DCF]: Cannot send broadcast 'datatosend' in TXRTS!"
        # create RTS
        src, dst = self.datatosend.addr2, self.datatosend.addr1
        rts = self.dot11rts(addr1=dst, addr2=src)
        if (self.retrycount>0):
            rts.FCfield |= DOT11_FC_RETRY
        # calculate NAV
        self.send_rts(rts, self.datatosend)
        nav = self.rtsnav(self.datatosend)
        rts.ID = nav
        pkt = crcupdate(rts)
        # send and hold for duration
        duration = self.duration(pkt)
        self.log_send(pkt, src=src, dst=dst, nav=nav, \
                      duration=time2usec(duration), retry=self.retrycount)
        yield self.TXD.send(fsm, [pkt])
        yield hold, fsm, duration
        # go to RXCTS
        yield fsm.goto(self.RXCTS)

    def RXCTS(self, fsm):
        """RXCTS state; wait for CTS response."""
        assert self.isdot11data(self.datatosend), \
                "[DCF]: Cannot determine 'datatosend' in RXCTS!"
        # start timeout timer
        timeout = self.ctstimeout
        timer = self.newchild("ctstimeout", Timer, timeout, start=True, \
                              tracename=self.tracename+".CTSTIMEOUT")
        yield waitevent, fsm, (self.rxdata, timer.done, timer.kill)
        # rxdata -> check for CTS
        if (self.rxdata in fsm.eventsFired):
            timer.halt()
            p, cts = self.rxdata.signalparam, None
            if self.isdot11cts(p):
                cts = self.get_dot11cts(p)
                addr, addr1 = self.address, cts.addr1
                crcerror = self.haserror(cts)
                # CTS for me -> update NAV -> pause for SIFS -> go to TXUCAST
                if (addr1==addr) and (not crcerror):
                    self.recv_cts(cts)     # process CTS
                    self.log_recv(cts, addr=addr, addr1=addr1, nav=cts.ID)
                    yield self.navupdate(fsm, cts.ID*1e-6)
                    yield hold, fsm, self.sifs
                    yield fsm.goto(self.TXUCAST)
            # otherwise -> retry -> go to RXPKT
            self.retry()
            yield fsm.goto(self.RXPKT, p)
        # timer done -> retry -> go to RESUME
        elif (timer.done in fsm.eventsFired):
            self.retry()
            yield fsm.goto(self.RESUME)
        # timer kill -> raise exception
        elif (timer.kill in fsm.eventsFired):
            raise RuntimeError, "[DCF]: Unexpected kill signal " + \
                                "from timer in RXCTS!"
        # otherwise -> raise error!
        else:
            raise RuntimeError, "[DCF]: Unexpected interruption in RXCTS!"

    def TXUCAST(self, fsm):
        """TXUCAST state; transmit unicast `datatosend` packet."""
        assert self.isdot11data(self.datatosend), \
                "[DCF]: Cannot determine 'datatosend' in TXUCAST!"
        assert not (self.datatosend.addr1==self.broadcast), \
                "[DCF]: Cannot send broadcast 'datatosend' in TXUCAST!"
        # FIXME: Assume all packet formatting has been done. No need to worry
        # about setting parameters or updating CRC.
        data = self.get_dot11data(self.datatosend)
        self.send_data(data)
        # send and hold for duration
        duration = self.duration(data)
        self.log_send(data, addr1=data.addr1, addr2=data.addr2)
        yield self.TXD.send(fsm, [data])
        yield hold, fsm, duration
        # go to RXACK
        yield fsm.goto(self.RXACK)

    def RXACK(self, fsm):
        """RXACK state; wait on ACK for `datatosend`."""
        assert self.isdot11data(self.datatosend), \
                "[DCF]: Cannot determine 'datatosend' in RXACK!"
        data = self.get_dot11data(self.datatosend)
        isbroadcast = (data.addr1==self.broadcast)
        assert not isbroadcast, \
                "[DCF]: Broadcast 'datatosend' cannot get ACK in RXACK!"
        # start timeout timer
        timeout = self.acktimeout
        timer = self.newchild("acktimeout", Timer, timeout, start=True, \
                              tracename=self.tracename+".ACKTIMEOUT")
        yield waitevent, fsm, (self.rxdata, timer.done, timer.kill)
        # rxdata -> check for ACK
        if (self.rxdata in fsm.eventsFired):
            timer.halt()
            p, ack = self.rxdata.signalparam, None
            if self.isdot11ack(p):
                ack = self.get_dot11ack(p)
                src, dst = data.addr2, data.addr1
                addr, addr1 = self.address, ack.addr1
                crcerror = self.haserror(ack)
                # ACK for me -> success -> go to IDLE
                if (addr1==src) and (not crcerror):
                    self.recv_ack(ack)     # process ACK
                    self.log_recv(ack, addr1=addr1, src=src, dst=dst)
                    yield fsm.goto(self.IDLE)
            # otherwise -> failure -> retry -> go to RXPKT
            self.retry()
            yield fsm.goto(self.RXPKT, p)
        # timer done -> timeout -> retry -> go to RESUME
        elif (timer.done in fsm.eventsFired):
            self.retry()
            yield fsm.goto(self.RESUME)
        # timer kill -> raise exception
        elif (timer.kill in fsm.eventsFired):
            raise RuntimeError, "[DCF]: Unexpected kill signal " + \
                                "from timer in RXACK!"
        # otherwise -> raise error!
        else:
            raise RuntimeError, "[DCF]: Unexpected interruption in RXACK!"
        return

    def RXBUSY(self, fsm):
        """RXBUSY state; check for `rxdata` and `csidle`."""
        yield waitevent, fsm, (self.rxdata, self.csidle)
        # rxdata -> go to RXPKT to classify
        if (self.rxdata in fsm.eventsFired):
            p = self.rxdata.signalparam
            yield fsm.goto(self.RXPKT, p)
        # csidle -> go back to RESUME
        elif (self.csidle in fsm.eventsFired):
            p = self.csidle.signalparam
            yield fsm.goto(self.RESUME)
        else:
            raise RuntimeError, \
                    "[DCF]: Unexpected interruption in RXBUSY!"
        return

    def RXPKT(self, fsm, pkt):
        """RXPKT state; classify incoming packets."""
        if (self.htype==const.ARP_HTYPE_ETHERNET):
            yield fsm.goto(self.RXETH, pkt)
        else:
            raise RuntimeError, "[DCF]: Unsupported hardware " + \
                    "type (%s) in RXPKT!"%(self.htype)

    def RXETH(self, fsm, pkt):
        """RXETH state; classify incoming Ethernet `pkt`.

        All packets that have an error will be dropped.
        """
        # drop packets with errors
        addr = self.address
        crcerror = self.haserror(pkt)
        if crcerror:
            self.log_drop(pkt, drop="CRC error")
            yield fsm.goto(self.RESUME)
        # classify Dot11 packet
        isdata = self.isdot11data(pkt)
        isrts  = self.isdot11rts(pkt)
        iscts  = self.isdot11cts(pkt)
        isack  = self.isdot11ack(pkt)
        if isdata:
            # receive DATA? -> go to RXDATA
            data = self.get_dot11data(pkt)
            yield fsm.goto(self.RXDATA, data)
        elif isrts:
            # receive RTS? -> go to RXRTS
            rts = self.get_dot11rts(pkt)
            yield fsm.goto(self.RXRTS, rts)
        elif iscts:
            # drop unsolicited CTS
            cts = self.get_dot11cts(pkt)
            addr1 = cts.addr1
            # update NAV -> log drop
            nav = cts.ID
            drop = "unsolicited CTS"
            self.log_drop(cts, addr=addr, addr1=addr1, drop=drop, nav=nav)
            yield self.navupdate(fsm, nav*1e-6)
            # process unsolicited CTS
            self.recv_cts(cts)
        elif isack:
            # unsolicited ack?
            ack = self.get_dot11ack(pkt)
            addr1, drop = ack.addr1, "unsolicited ack"
            self.log_drop(pkt, addr=addr, addr1=addr1, drop=drop)
            # process unsolicited ACK
            self.recv_ack(ack)
        else:
            raise RuntimeError, "[DCF]: Got unexpected message in RXETH!"
        # go to RESUME
        yield fsm.goto(self.RESUME)

    def RXDATA(self, fsm, pkt):
        """RXDATA state; process received DATA message."""
        crcerror = self.haserror(pkt)
        if crcerror:
            self.log_drop(pkt, drop="CRC error")
            yield fsm.goto(self.RESUME)
        assert self.isdot11data(pkt), "[DCF]: Cannot find DATA in RXDATA!"
        data = self.get_dot11data(pkt)
        isbroadcast = (data.addr1==self.broadcast)
        isforme     = (data.addr1==self.address)
        promiscuous = self.promiscuous
        src, dst    = data.addr2, data.addr1
        iseth =  isinstance(data, Packet) and data.haslayer(Ether)
        # isbroadcast or isforme -> send to RXU
        if (isforme or isbroadcast or promiscuous) and iseth:
            self.recv_data(data)
            eth = data[Ether]
            p = crcremove(eth)
            self.log_recv(data, src=src, dst=dst, promiscuous=promiscuous)
            data.remove_payload()
            yield self.TXU.send(fsm, [p])
            assert fsm.stored(self.TXU), \
                    "[DCF]: Error sending packet to 'TXU'!"
            if isforme: yield fsm.goto(self.TXACK, data)
        # non-Ethernet packets or not for me -> drop
        else:
            drop = "non-Ethernet packet"
            if not isforme: drop = "not for me"
            self.log_drop(pkt, src=src, dst=dst, drop=drop)
        # go to RESUME
        yield fsm.goto(self.RESUME)

    def TXACK(self, fsm, pkt):
        """TXACK state; transmit ACK message in response to `data`."""
        assert self.isdot11data(pkt), "[DCF]: Cannot find Dot11Data in TXACK!"
        data = self.get_dot11data(pkt)
        assert not (data.addr1==self.broadcast), \
                "[DCF]: Cannot send ACK for broadcast data!"
        addr1 = data.addr2
        ack = self.dot11ack(addr1=addr1)
        self.send_ack(ack)
        pkt = crcupdate(ack)
        # pause for SIFS
        yield hold, fsm, self.sifs
        # send and hold duration
        duration = self.duration(pkt)
        self.log_send(pkt, addr1=addr1, duration=duration)
        yield self.TXD.send(fsm, [pkt])
        yield hold, fsm, duration
        yield fsm.goto(self.RESUME)

    def RXRTS(self, fsm, pkt):
        """RXRTS state; update NAV and process RTS message."""
        crcerror = self.haserror(pkt)
        if crcerror:
            self.log_drop(pkt, drop="CRC error")
            yield fsm.goto(self.RESUME)
        assert self.isdot11rts(pkt), "[DCF]: Cannot find RTS in RXRTS!"
        addr, rts = self.address, self.get_dot11rts(pkt)
        addr1, addr2 = rts.addr1, rts.addr2
        isforme = (addr1==addr)
        # check NAV? -> update NAV or drop RTS
        nav = rts.ID
        navbusy = self.navbusy
        navidle = self.navidle
        # process incoming RTS
        self.recv_rts(rts)
        # NAV IDLE -> send CTS
        if (isforme and navidle):
            self.log_recv(rts, addr1=addr1, addr2=addr2, nav=nav)
            yield fsm.goto(self.TXCTS, rts)
        # NAV BUSY -> drop RTS
        elif (isforme and navbusy):
            drop = "NAV busy"
            self.log_drop(rts, addr1=addr1, addr2=addr2, drop=drop, nav=nav)
        # not for me -> drop RTS
        else:
            drop = "not for me"
            self.log_drop(rts,addr=addr,addr1=addr1,addr2=addr2,drop=drop)
            yield self.navupdate(fsm, nav*1e-6)
        # go to RESUME
        yield fsm.goto(self.RESUME)

    def TXCTS(self, fsm, rts):
        """TXCTS state; send CTS response message."""
        assert self.isdot11rts(rts), "[DCF]: Cannot find RTS in TXCTS!"
        addr, addr1 = self.address, rts.addr2
        # create CTS
        pkt = self.dot11cts(addr1=addr1)
        cts = crcupdate(pkt)
        # update nav
        self.send_cts(cts, rts)
        nav = self.ctsnav(rts)
        cts.ID = nav
        pkt = crcupdate(cts)
        # pause for SIFS
        yield hold, fsm, self.sifs
        # send and hold duration
        duration = self.duration(pkt)
        self.log_send(pkt, addr=addr, addr1=addr1, nav=nav)
        yield self.TXD.send(fsm, [pkt])
        yield hold, fsm, duration
        # set NAV and resume
        yield self.navupdate(fsm, nav*1e-6)
        yield fsm.goto(self.RESUME)

    def RESUME(self, fsm):
        """RESUME state; resume operation in `IDLE` or `BACKOFF`."""
        if self.datatosend:
            yield fsm.goto(self.BACKOFF)
        else:
            yield fsm.goto(self.IDLE)

    def RECV(self, fsm):
        """RECV state; check for receive data from downstream element."""
        yield self.RXD.recv(fsm, 1)
        assert fsm.acquired(self.RXD) and (len(fsm.got)==1), \
                "[DCF]: Error receiving from RXD in RECV state!"
        p = fsm.got[0]
        self.rxdata.signal(p)
        # continue in RECV
        yield fsm.goto(self.RECV)

    def send_rts(self, rts, data):
        """Additional processing for outgoing RTS."""
        errmsg = "[DCF]: Cannot process non-RTS in send_rts()!"
        assert self.isdot11rts(rts), errmsg
        errmsg = "[DCF]: Cannot process non-DATA in send_rts()!"
        assert self.isdot11data(data), errmsg
        rts.setanno('mac-root', str(data.traceid))
        if data.hasanno('net-root'):
            rts.setanno('net-root', data.getanno('net-root'))

    def recv_rts(self, rts):
        """Additional processing for incoming RTS."""
        errmsg = "[DCF]: Cannot process non-RTS in recv_rts()!"
        assert self.isdot11rts(rts), errmsg

    def send_cts(self, cts, rts):
        """Additional processing for outgoing CTS."""
        errmsg = "[DCF]: Cannot process non-CTS in send_cts()!"
        assert self.isdot11cts(cts), errmsg
        errmsg = "[DCF]: Cannot process non-RTS in send_cts()!"
        assert self.isdot11rts(rts), errmsg

    def recv_cts(self, cts):
        """Additional processing for incoming CTS."""
        errmsg = "[DCF]: Cannot process non-CTS in recv_cts()!"
        assert self.isdot11cts(cts), errmsg

    def send_data(self, data):
        """Additional processing for outgoing DATA."""
        errmsg = "[DCF]: Cannot process non-DATA in send_data()!"
        assert self.isdot11data(data), errmsg
        data.setanno('mac-root', str(data.traceid))

    def recv_data(self, data):
        """Additional processing for incoming DATA."""
        errmsg = "[DCF]: Cannot process non-DATA in recv_data()!"
        assert self.isdot11data(data), errmsg
        data.setanno('mac-rxts', now())

    def send_ack(self, ack):
        """Additional processing for outgoing ACK."""
        errmsg = "[DCF]: Cannot process non-ACK in send_ack()!"
        assert self.isdot11ack(ack), errmsg

    def recv_ack(self, ack):
        """Additional processing for incoming ACK."""
        errmsg = "[DCF]: Cannot process non-ACK in recv_ack()!"
        assert self.isdot11ack(ack), errmsg
        # Expecting ACK?
        if self.isdot11data(self.datatosend):
            data = self.get_dot11data(self.datatosend)
            dst, src = data.addr1, data.addr2
            # ACK for me? -> received ACK -> signal ackdata
            if (src==ack.addr1):
                pkt = self.datatosend.payload
                self.datatosend.remove_payload()
                p = crcremove(pkt)
                self.ackdata.signal(p)

    def rtsnav(self, p, *args, **kwargs):
        """Calculate NAV value for RTS given DATA `p`.

        :param p: DATA packet.
        :param args: Additional arguments passed to compute `duration` of DATA.
        :param kwargs: Keywords passed to compute `duration` of DATA.
        :return: Integer; representing NAV value.

        This method uses `duration()` to compute the NAV as follows:

            NAV = SIFS + CTS + SIFS + DATA + SIFS + ACK

        :note: It is assumed that `p` is a valid packet.
        """
        d = self.sifs + self.ctsduration
        d += self.sifs + self.duration(p, *args, **kwargs)
        d += self.sifs + self.ackduration
        nav = int(d*1e6)
        return nav

    def ctsnav(self, rts):
        """Calculate NAV value for CTS.

        :param rts: RTS packet.
        :return: Integer; representing NAV value.

        This method computes CTS NAV as follows:

            NAV = RTSNAV - SIFS - CTS

        :note: Assumes NAV for `rts` is greater than a SIFS + CTS duration.
        """
        errmsg = "[DCF]: Invalid non-RTS message in ctsnav()!"
        assert self.isdot11rts(rts), errmsg
        # compute CTS NAV
        nav = rts.ID - int((self.sifs+self.ctsduration)*1e6)
        return nav

    def navupdate(self, proc, t=None):
        """Update NAVTimer.

        :param proc: Process that blocks on NAV update.
        :param t: New timer value for virtual carrier sense busy.

        If `t` is None, this method will reset the NAVTimer `nav`.
        """
        return self.nav.update(proc, t)

    def get_ctstimeout(self):
        """Calculate timeout for CTS messages.

        CTSTIMEOUT = SIFS + CTS duration + SLOTTIME
        """
        timeout = self.sifs + self.ctsduration + self.slottime
        return timeout

    def get_acktimeout(self):
        """Calculate timeout for ACK messages.

        ACKTIMEOUT = SIFS + ACK duration + SLOTTIME
        """
        timeout = self.sifs + self.ackduration + self.slottime
        return timeout

    def get_ctsduration(self, force=False):
        """Calculate duration of CTS message.

        :param force: If true, ignore any cached value.
        """
        if force or (self._ctsduration is None):
            cts = self.dot11cts()
            pkt = crcupdate(cts)
            self._ctsduration = self.duration(pkt)
        return self._ctsduration

    def get_ackduration(self, force=False):
        """Calculate duration of ack message.

        :param force: If true, ignore any cached value.
        """
        if force or (self._ackduration is None):
            ack = self.dot11ack()
            pkt = crcupdate(ack)
            self._ackduration = self.duration(pkt)
        return self._ackduration

    def set_timing(self, phymode=None):
        """Set up timing parameters (e.g. `sifs`, `slottime`, etc.).

        :param phymode: Enumeration for physical layer mode of operation
                        [default=None].

        This method sets up timing parameters based on the `phymode` enumeration
        which enumerates the physical layer mode of operation.

        :note: This method is called from `configure()`.
        """
        m = DOT11A_PHY_MODE 
        if phymode is None:
            p = self.phy
            if isinstance(p, Dot11APHY):    m = DOT11A_PHY_MODE
            elif isinstance(p, Dot11NPHY):  m = DOT11N_PHY_MODE
        else:
            m = phymode
        # check PHY mode
        assert (m in DOT11_TIMING), \
                "[DCF]: Cannot set timing from invalid PHY mode (%s)!"(m)
        # set up timing
        self.sifs = DOT11_TIMING[m]['sifs']
        self.slottime = DOT11_TIMING[m]['slottime']

    def get_datarate(self, r):
        """Get the data rate in bits-per-second (bps).

        :param r: Rate enumeration.

        This method calls `get_datarate()` for `phy`.
        """
        return self.phy.get_datarate(r)


    def calclength(self, duration, rate):
        """Calculate length of packet in bytes.

        :param duration: Duration of packet in seconds.
        :param rate: Rate enumeration.

        This method calls `calclength()` for `phy`.
        """
        return self.phy.calclength(duration, rate)

    def isdot11data(self, p):
        """Check if packet is DATA; *overload as needed*."""
        return isdot11data(p)

    def isdot11rts(self, p):
        """Check if packet is RTS; *overload as needed*."""
        return isdot11rts(p)

    def isdot11cts(self, p):
        """Check if packet is CTS; *overload as needed*."""
        return isdot11cts(p)

    def isdot11ack(self, p):
        """Check if packet is ACK; *overload as needed*."""
        return isdot11ack(p)

    def get_dot11data(self, p):
        """Extract DATA from `p`; *overload as needed*."""
        return get_dot11data(p)

    def get_dot11rts(self, p):
        """Extract RTS from `p`; *overload as needed*."""
        return get_dot11rts(p)

    def get_dot11cts(self, p):
        """Extract CTS from `p`; *overload as needed*."""
        return get_dot11cts(p)

    def get_dot11ack(self, p):
        """Extract ACK from `p`; *overload as needed*."""
        return get_dot11ack(p)

    def dot11data(self, *args, **kwargs):
        """Create new `Dot11Data` packet."""
        return Dot11Data(*args, **kwargs)

    def dot11rts(self, *args, **kwargs):
        """Create new `Dot11RTS` packet."""
        return Dot11RTS(*args, **kwargs)

    def dot11cts(self, *args, **kwargs):
        """Create new `Dot11CTS` packet."""
        return Dot11CTS(*args, **kwargs)

    def dot11ack(self, *args, **kwargs):
        """Create new `Dot11ACK` packet."""
        return Dot11Ack(*args, **kwargs)

    def connect(self, p):
        """Overloaded to connect and call `set_phy()`."""
        self.set_phy(p)
        return CSMAC.connect(self, p)

    def log_send(self, p, *args, **kwargs):
        """Convenience method for logging send event."""
        if self.verbose>DOT11_VERBOSE:
            kwargs['addr'] = self.address
            kwargs['retrycount'] = self.retrycount
            kwargs['retrylimit'] = self.retrylimit
            if p.hasanno('cif-duration'):
                kwargs['cif-duration'] = time2usec(p.getanno('cif-duration') )
            if p.hasanno('phy-rate'):
                kwargs['phy-rate'] = p.getanno('phy-rate')
            self.log("snd", p, *args, **kwargs)

    def log_recv(self, p, *args, **kwargs):
        """Convenience method for logging receive event."""
        if self.verbose>DOT11_VERBOSE:
            kwargs['addr'] = self.address
            if p.hasanno('phy-sinr'):
                kwargs['phy-sinr'] = "%.4f dB"%(p.getanno('phy-sinr') )
            self.log("rcv", p, *args, **kwargs)

    def log_drop(self, p, *args, **kwargs):
        """Convenience method for logging drop event."""
        if self.verbose>DOT11_VERBOSE:
            kwargs['addr'] = self.address
            self.log("drp", p, *args, **kwargs)

    def log(self, event=None, p=None, *args, **kwargs):
        """Overloaded to check verbose level and set common annotations."""
        force = False
        if ('verbose' in kwargs): force = (kwargs['verbose']>DOT11_VERBOSE)
        if self.verbose>DOT11_VERBOSE or force:
            kwargs.update(self.get_dcf_anno(p))
            CSMAC.log(self, event, p, *args, **kwargs)

    def get_dcf_anno(self, p):
        """Convenience method to extract annotations and convert to strings."""
        kwargs = {}
        if not isinstance(p, Packet): return kwargs
        kwargs['addr'] = self.address
        if p.hasanno('cif-duration'):
            kwargs['cif-duration'] = time2usec(p.getanno('cif-duration') )
        if p.hasanno('phy-rate'):
            kwargs['phy-rate'] = p.getanno('phy-rate')
        if p.hasanno('phy-sinr'):
            kwargs['phy-sinr'] = "%.4f dB"%(p.getanno('phy-sinr') )
        if p.hasanno('net-root'):
            kwargs['net-root'] = p.getanno('net-root')
        if p.hasanno('mac-root'):
            kwargs['mac-root'] = p.getanno('mac-root')
        if p.hasanno('mac-txts'):
            kwargs['mac-txts'] = p.getanno('mac-txts')
        if p.hasanno('mac-rxts'):
            kwargs['mac-rxts'] = p.getanno('mac-rxts')
        return kwargs
