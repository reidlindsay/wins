#!  /usr/bin/env python

"""
Implements IEEE 802.11a PHY protocol; contains `Dot11APHY` and `Dot11ARadio`
classes.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-10-05 02:17:53 -0500 (Wed, 05 Oct 2011) $
* $LastChangedRevision: 5186 $

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

from SimPy.Simulation import SimEvent, hold, now
from wins.protocol.csphy import CSPHY
from wins.channel.radio  import Radio
from wins.digital.mqam   import MQAM
from wins.digital.rcpc   import RCPC
from wins.ieee80211.dot11a_support import *
from wins.helper import time2usec, linear2db, db2linear, monitor_events
from wins.packet import ANNO
from wins.fsm    import FSM
from wins import const

from numpy import random, inf, ceil

class Dot11APHY(CSPHY):
    """Implementation of IEEE 802.11a physical layer protocol.

    Modulation and Coding Parameters
    ================================

    The table below describes modulation and convolutional coding parameters for
    IEEE 802.11a.

    =========== ============= ============ ====== ======= =========
    Rate Index    Data Rate    Modulation   Nbps     M    Code rate
    ----------- ------------- ------------ ------ ------- ---------
         0          6 Mbps       BPSK         1      2        1/2
         1          9 Mbps       BPSK         1      2        3/4
         2         12 Mbps       QPSK         2      4        1/2
         3         18 Mbps       QPSK         2      4        3/4
         4         24 Mbps      16-QAM        4     16        1/2
         5         36 Mbps      16-QAM        4     16        3/4
         6         48 Mbps      64-QAM        6     64        2/3
         7         54 Mbps      64-QAM        6     64        3/4
    =========== ============= ============ ====== ======= =========

    :note: Rate index is used in Packets when referring to rate. Use
           `DOT11A_DATARATE` to convert this to a bitrate value.

    Dot11APHY and Ports
    ===================
    Every `Dot11APHY` has the following ports:

        1. "TXU" - sends decoded packets to an upstream element.
        #. "RXU" - receives traffic to send from an upstream element.
        #. "TXD" - sends encoded packets to a downstream element.
        #. "RXD" - receives traffic from a downstream element.

    The upstream element of a `Dot11APHY` is usually a `CSMAC` (or subclass),
    and the downstream element is usually a `ChannelInterface` (or subclass).

    Dot11APHY and Annotations
    =========================
    This class uses/sets the following annotations:

    ============== ==========================================================
    Name            Description
    ============== ==========================================================
    cif-txpower     Power of transmitted packets (in dBm).
    -------------- ----------------------------------------------------------
    rate            Rate index used to specify coding and modulation scheme.
    -------------- ----------------------------------------------------------
    cif-duration    Duration of packet is marked on outgoing packets using
                    the `duration()` method during `SEND`.
    -------------- ----------------------------------------------------------
    cif-collision   List of packets that arrived at the same time as the
                    current packet (i.e. collided with the marked packet).
    -------------- ----------------------------------------------------------
    dot11a-detect   Indicates whether or not the physical layer packet
                    detection was successful (see `framedetect()`).
    -------------- ----------------------------------------------------------
    dot11a-header   Indicates whether or not the physical layer header
                    decoding was successful (see `decode_header()`).
    -------------- ----------------------------------------------------------
    dot11a-per      Packet-error rate (PER) modeled by the physical layer
                    (the `decode_header()` will initially set this, but
                    `decode_data()` will overwrite it as necessary).
    ============== ==========================================================

    :IVariables:
     * `detect`: Internal SimEvent signalled when a packet is detected.
     * `mod`: `MQAM`, child member for modulation.
     * `coder`: `RCPC`, child member for coding.
     * `detectdelay`: Duration required for packet detection.

    :CVariables:
     * `detectdelay`: Duration required for packet detection.
    """
    name = "IEEE 802.11a"
    tracename = "80211A"
    detectdelay = DOT11A_TDETECT
    def __init__(self, **kwargs):
        """Constructor."""
        self.detect = SimEvent(name=".detect")
        self.detectdelay = None
        CSPHY.__init__(self, **kwargs)
        self.detect.name = "%s%s"%(self.name,".detect")
        #monitor_events(self.detect)

    def configure(self, detectdelay=None, **kwargs):
        """Call `CSPHY.configure()`; add ports, `FSM`, and other members."""
        if detectdelay is None: self.detectdelay = self.__class__.detectdelay
        CSPHY.configure(self, **kwargs)
        # add ports and FSM
        self.addport("RXU"), self.addport("TXU")
        self.addport("TXD"), self.addport("RXD")
        # create FSM to manage send/recv execution of phy
        txfsm = self.newchild("txfsm", FSM, tracename=self.tracename+".TX")
        rxfsm = self.newchild("rxfsm", FSM, tracename=self.tracename+".RX")
        txfsm.goto(self.SEND)
        rxfsm.goto(self.LISTEN)
        # set up other members
        mod = self.newchild("mod", MQAM, tracename=self.tracename+".MQAM")
        coder = self.newchild("coder", RCPC, tracename=self.tracename+".RCPC")

    def connect(self, radio):
        """Convenience method to connect PHY to an RF front end (i.e. `Radio`).

        **Overload this method to change how a connection is made.**

        :note: This method also calls `set_radio()` to set the appropriate
               pointer for the physical layer.
        """
        assert isinstance(radio, Radio), \
                "[DOT11A]: Cannot connect to non-Radio!"
        self.set_radio(radio)
        self.TXD.connect(radio.getport("RXU"))
        radio.getport("TXU").connect(self.RXD)

    def duration(self, p, rate=None):
        """Calculate duration of packet `p` using `calcduration()` method.

        :param p: Packet to compute duration for (or packet length in bytes).
        :param rate: Optional rate index to denote modulation/coding scheme.
        :return: Duration of packet in seconds.

        This method checks for the 'phy-rate' annotation, which is a rate index
        to denote the coding and modulation scheme to be used. If a 'phy-rate'
        annotation is not found (or specified as a parameter), this method uses
        the base rate (i.e. rate 0) to calculate the duration of the waveform.
        """
        plen = p
        if isinstance(p, Packet):
            if p.hasanno("rate") and (rate is None):
                rate = p.getanno("rate")
            plen = len(p)
        if (rate is None): rate = 0
        return self.calcduration(plen, rate)

    def calcper_header(self, p, **kwargs):
        """Calculate probability of error for header decoding.

        :param p: Packet being decoded.
        :param kwargs: Additional keywords arguments passed to `sinr_heap()`
                       (or `sinr()`).
        :return: PER for header decoding.

        This method sets the 'dot11a-sinr' and 'dot11a-per' annotations. The
        operation of this method depends on `DOT11A_USE_PIECEWISE_PER`.
        """
        for a in ['cif-rxts']:
            errmsg = "[DOT11APHY]: calcper_header() cannot find " + \
                     "'%s' annotation!"%(a)
            assert ANNO.supports(p, a), errmsg
        # calculate PER using appropriate method
        plen = len(p.payload)
        if DOT11A_USE_PIECEWISE_PER:
            sinrheap = self.sinr_heap(p, **kwargs)
            t0 = p.getanno('cif-rxts') + DOT11A_TSHORT + DOT11A_TLONG
            t1 = t0 + DOT11A_TSIGNAL
            xheap = [(max(ta,t0), min(tb,t1), sinr) for (ta,tb,sinr) \
                        in sinrheap if (ta<t1) and (tb>t0)]
            errmsg = "[DOT11APHY]: Unable to find valid data from SINR heap!"
            assert (len(xheap)>0), errmsg
            # calculate piecewise PER and average SINR
            psuccess, stot = 1.0, 0.0
            for ta, tb, sinr in xheap:
                alpha = (tb-ta)/(t1-t0)
                hlen = len(Dot11A())*alpha
                stot += db2linear(sinr)*alpha
                psuccess *= 1.0 - self.calcper(hlen, 0, sinr)
            per = 1.0 - psuccess
            sinr = linear2db(stot)
        else:
            sinr, hlen = self.sinr(p, **kwargs), len(p)-plen
            # configure modulation and coding to calculate PER
            per = self.calcper(hlen, 0, sinr)
        # set annotations and return PER
        p.setanno('dot11a-sinr', sinr)
        p.setanno('dot11a-per', per)
        return per

    def calcper_data(self, p, **kwargs):
        """Calculate probability of error for data decoding.

        :param p: Packet being decoded.
        :param kwargs: Additional keywords arguments passed to `sinr_heap()`
                       (or `sinr()`).
        :return: PER for decoding packet payload.

        This method sets the 'dot11a-sinr' and 'dot11a-per' annotations. The
        operation of this method depends on `DOT11A_USE_PIECEWISE_PER`.
        """
        for a in ['cif-rxts','cif-duration']:
            errmsg = "[DOT11APHY]: calcper_data() cannot find " + \
                     "'%s' annotation!"%(a)
            assert ANNO.supports(p, a), errmsg
        # verify header parameters
        plen = len(p.payload)
        rate, length = p.rate, p.length
        assert (0<=rate<len(DOT11A_DATARATE) ), \
                "[DOT11A]: Invalid rate option (%s)!"%(rate)
        assert (p.length==plen), "[DOT11A]: Header length reported " + \
                "does not equal payload length; %s!=%s"%(p.length, plen)
        # calculate PER using appropriate method
        if DOT11A_USE_PIECEWISE_PER:
            sinrheap = self.sinr_heap(p, **kwargs)
            t1 = p.getanno('cif-rxts') + p.getanno('cif-duration')
            t0 = t1 - self.calcnofdm(plen, rate)*DOT11A_TSYM
            xheap = [(max(ta,t0), min(tb,t1), sinr) for (ta,tb,sinr) \
                        in sinrheap if (ta<t1) and (tb>t0)]
            errmsg = "[DOT11APHY]: Unable to find valid data from SINR heap!"
            assert (len(xheap)>0), errmsg
            # calculate piecewise PER and average SINR
            psuccess, stot = 1.0, 0.0
            for ta, tb, sinr in xheap:
                alpha = (tb-ta)/(t1-t0)
                dlen = plen*alpha
                stot += db2linear(sinr)*alpha
                psuccess *= 1.0 - self.calcper(dlen, rate, sinr)
            per = 1.0 - psuccess
            sinr = linear2db(stot)
        else:
            # configure modulation and coding to calculate PER
            sinr, plen = self.sinr(p, **kwargs), length
            per = self.calcper(plen, rate, sinr)
        # set annotations and return PER
        p.setanno('dot11a-sinr', sinr)
        p.setanno('dot11a-per', per)
        return per

    def calcper(self, p, rate, sinr):
        """Calculate packet-error rate for given parameters.

        :param p: Packet or packet length in bytes.
        :param rate: Rate enumeration to indicate modulation and coding scheme.
        :param sinr: Signal-to-interference-and-noise ratio.
        """
        # configure modulation and coding to calculate PER
        assert (0<=rate<len(DOT11A_DATARATE) ), \
                "[DOT11A]: Invalid rate enumeration!"
        plen = p
        if isinstance(p, Packet): plen = len(p)
        mod, coder = self.mod, self.coder
        mod.mtype = DOT11A_MTYPE[rate]
        coder.rate = DOT11A_CODERATE[rate]
        uber = mod.ber(sinr)
        per = coder.per(plen, uber)
        return per

    def encode(self, p, rate=None, txpower=None):
        """Encapsulate MPDU in a `Dot11A` packet; set annotations.

        :param p: MPDU to encode.
        :param rate: Optional rate index to denote modulation/coding scheme.
        :param txpower: Optional transmit power (in dBm).
        :return: `Dot11A` packet with MPDU as payload.

        This method will make sure that the following annotations are set:
            * rate [default=0]
            * txpower [default=`DOT11A_MAXPOWER`]
            * duration [from `duration()`]

        If the `rate` parameter is not specified and no 'phy-rate' annotation is
        found in packet `p`, this method will use the base rate (i.e. zero) to
        encode packet `p`.

        If the `txpower` parameter is not specified and no 'cif-txpower'
        annotation is found in packet `p`, this method will use
        `DOT11A_MAXPOWER` as the default transmit power.
        """
        # check parameters
        if p.hasanno('phy-rate') and (rate is None):
            rate = p.getanno('phy-rate')
        if p.hasanno('cif-txpower') and (txpower is None):
            txpower = p.getanno('cif-txpower')
        if (rate is None): rate = 0
        if (txpower is None): txpower = DOT11A_MAXPOWER
        duration = self.duration(p, rate=rate)
        # set annotations
        p.setanno('phy-rate', rate)
        p.setanno('cif-txpower', txpower)
        p.setanno('cif-duration', duration)
        # encap in Dot11A
        length = len(p)
        w = Dot11A(rate=rate, length=length)/p
        return w

    def framedetect(self, p, thresh=None):
        """Apply packet detection model for detecting training sequence; based
        on signal-to-interference-and-noise ratio (SINR).

        :param p: `Dot11A` packet being received.
        :param thresh: Packet detection SINR threshold (in dB)
                       [default=`DOT11A_FDTHRESHOLD`].
        :return: Boolean flag; if true, packet detection was successful.

        This method checks to make sure `p` is a `Dot11A` packet, and that the
        received SINR is greater than the receiver detection threshold `thresh`.
        This method will also mark the 'dot11a-detect' annotation to indicate
        the success (or failure) of the frame detection.

        **Overload this method to change how frame detection works.**

        :note: If `p` is not a `Dot11A` packet, this method will set the
               'dot11a-detect' annotation to false and return false.
        """
        if not isinstance(p, Dot11A):
            if ANNO.supported(p): p.setanno('dot11a-detect', False)
            return False
        # check SINR
        if thresh is None: thresh = DOT11A_FDTHRESHOLD
        sinr = self.sinr(p)
        detect = True
        if sinr<thresh: detect = False
        # mark annotations
        p.setanno('dot11a-detect', detect)
        return detect

    def decode_header(self, p):
        """Apply physical layer model for header decoding.

        :param p: `Dot11A` packet being received.
        :return: Boolean flag; if true, header decoding was successful.

        This method uses `mod` and `coder` to determine the error
        characteristics of the header decoding process. The following conditions
        must be met in order to have successful header decoding:

            * packet `p` must be a `Dot11A` packet,
            * the 'dot11a-detect' annotation must be true,
            * the header decoding must succeed,
            * and all header parameters must be valid.
        
        This method marks the 'dot11a-header' annotation to mark success (or
        reason for failure) of header decoding, and sets the 'dot11a-per'
        annotation to indicate the probability of error in decoding the header.
        """
        header, danno = "", 'dot11a-detect'
        isdot11a = isinstance(p, Dot11A)
        detected = isdot11a and p.hasanno(danno) and p.getanno(danno)
        if not (isdot11a and detected):
            if not isdot11a: header = "not Dot11A packet"
            if not detected: header = "not detected"
            p.setanno('dot11a-header', header)
            return False
        # check header parameters
        plen = len(p.payload)
        rate, length = p.rate, p.length
        okrate = (0<=rate<len(DOT11A_DATARATE) )
        oklen  = (p.length==plen)
        okpar  = self.parity(p)
        if not (okrate and oklen and okpar):
            header = "header parameters failed"
            p.setanno('dot11a-header', header)
            return False
        # decode header of Dot11A
        per = self.calcper_header(p, force=False)
        sinr = p.getanno('dot11a-sinr')
        # simulate packet header decoding errors
        header = "success"
        if (random.uniform(0,1)<per): header = "decoding failed"
        # mark annotations
        p.setanno('dot11a-header', header+" (SINR = %.2f dB)"%(sinr) )
        p.setanno('dot11a-per', per)
        return (header=="success")

    def decode(self, p):
        """Apply physical layer model for decoding a `Dot11A` packet.

        :param p: `Dot11A` packet to decode.
        :return: Decoded payload (or `None` if header decoding fails).

        This method uses `mod` and `coder` to simulate packet errors. This
        method decodes the `Dot11A` header and then decodes the payload. If
        `decode_header()` fails, this method will return `None` and skip payload
        decoding. Otherwise, this method will simulate errors and call
        `seterror()` if any packet errors occur.

        This method also overwrites the 'dot11a-per' annotation with the packet
        error probability for the payload.

        :note: This method will log a drop event if `decode_header()` fails.
        """
        header, hanno = "header failed", "dot11a-header"
        # update SINR (and SINR heap)
        sinr = self.sinr(p, force=True)
        if not self.decode_header(p):
            if p.hasanno(hanno): header = p.getanno(hanno)
            self.log_drop(p, header=header)
            return None
        # verify header decoding
        detectanno = 'dot11a-detect'
        errmsg =  "[DOT11APHY]: Cannot decode payload of non-Dot11A packet!"
        assert isinstance(p, Dot11A), errmsg
        errmsg = "[DOT11APHY]: Cannot decode payload that has not been detected!"
        assert p.hasanno(detectanno) and p.getanno(detectanno), errmsg
        # calculate PER
        per = self.calcper_data(p, force=False)
        sinr = p.getanno('dot11a-sinr')
        # simulate packet header decoding errors
        error = False
        if (random.uniform(0,1)<per): error = True
        # mark annotations
        pkt = p.payload
        if error: self.seterror(pkt)
        pkt.setanno('dot11a-per', per)
        return pkt

    def sinr(self, p, **kwargs):
        """Calculate signal-to-interference-and-noise ratio (SINR).

        :param p: Packet to compute SINR for.
        :param kwargs: Additional keyword arguments passed to `sinr_heap()`.

        This method uses the 'rxpower', 'noisepower', and 'cif-collision'
        annotations to calculate the SINR of the received packet p.

        :note: This method sets the 'phy-sinr' annotation indicating the SINR (in dB).
        """
        for a in ['rxpower', 'noisepower', 'cif-collision']:
            errmsg = "[DOT11APHY]: sinr() cannot find '%s' annotation!"%(a)
            assert ANNO.supports(p, a), errmsg
        # get SINR heap
        sinrheap = self.sinr_heap(p,**kwargs)
        minsinr = +inf
        # find minimum SINR
        for ta,tb,sinr in sinrheap:
            if sinr<minsinr: minsinr = sinr
        # set annotations
        p.setanno('phy-sinr', minsinr)
        return minsinr

    def parity(self, p):
        """Check parity of `Dot11A` packet header.

        :param p: `Dot11A` packet.
        :return: Boolean; true if parity check passes.

        :note: If `DOT11A_USEPARITY` is false, this method returns true.
        """
        assert isinstance(p, Dot11A), \
               "[DOT11APHY]: Cannot check parity of non-Dot11A packet!"
        return p.checkpar()

    def SEND(self, fsm):
        """SEND state; simulate encoding and send process.

        This state performs the following tasks:

            1. Get packet from 'RXU' port.
            2. Call `encode()` to generate waveform for outgoing packet.
            3. Mark 'cif-duration' annotation with value returned by `duration()`.
            4. Simulate `txproctime` for outgoing packet.
            5. Send outgoing waveform to 'TXD'.
            6. Simulate `duration` of waveform.
        """
        while fsm.active():
            yield self.RXU.recv(fsm, 1)
            assert fsm.acquired(self.RXU) and (len(fsm.got)==1), \
                    "[DOT11APHY]: Error receiving from RXU port in SEND()!"
            p = fsm.got[0]
            # simulate encoding and Tx processing time
            w = self.encode(p)
            duration = self.duration(p)
            if ANNO.supported(w): w.setanno('cif-duration', duration)
            yield hold, fsm, self.txproctime(p)
            # send waveform and simulate duration
            self.log_send(w, duration=time2usec(duration) )
            yield self.TXD.send(fsm, [w])
            assert fsm.stored(self.TXD), \
                    "[DOT11APHY]: Error sending to TXD in SEND!"
            yield hold, fsm, duration
        return

    def LISTEN(self, fsm):
        """LISTEN state; monitor `radio` and manage packet detection."""
        r = self.radio
        assert isinstance(r, Radio), "[DOT11A]: Cannot find radio in LISTEN!"
        while fsm.active():
            # check rxenergy -> set csbusy?
            rxenergy = r.rxenergy()
            rxhigh = r.inreceive and (r.rxenergy()>DOT11A_CSTHRESHOLD)
            if rxhigh: self.set_csbusy(rxenergy="high, %.2f dBm"%(rxenergy), \
                        rxbuffer=[x.traceid for x in r.rxbuffer]  )
            else:      self.set_csidle(rxenergy="%.2f dBm"%(rxenergy) )
            # monitor events and RXD port
            yield self.RXD.recv(fsm, 1, \
                  renege=(r.rxdata, r.rxdone, r.txdata, r.txdone, self.detect) )
            # RXD -> ignore incoming packets in LISTEN
            if fsm.acquired(self.RXD):
                assert (len(fsm.got)==1), "[DOT11A]: Received unexpected " + \
                        "number of packets from 'RXD' port in LISTEN state!"
                p = fsm.got[0]
                self.log_drop(p, drop="not detected in LISTEN")
            # rxdata -> start DETECT thread
            if r.rxdata in fsm.eventsFired:
                p = r.rxdata.signalparam
                fname = "detect(%s)"%(p._id)
                ### XXX ####
                f = FSM()
                #f = self.newchild(fname, FSM, tracename=fname.upper() )
                f.goto(self.DETECT, p)
                f.start()
            # detect -> set csbusy -> goto DECODE
            if self.detect in fsm.eventsFired:
                p = self.detect.signalparam
                rxenergy = "%.2f dBm"%(r.rxenergy() )
                sinr  = "%.2f dB"%(self.sinr(p) )
                self.set_csbusy(p, detect=True, rxenergy=rxenergy)
                danno = 'dot11a-detect'
                errmsg = "[DOT11A]: Cannot find 'dot11a-detect' " + \
                         "annotation in detected packet!"
                assert ANNO.supports(p, danno) and p.getanno(danno), errmsg
                #yield hold, fsm, 0
                self.log("detect", p, rxenergy=rxenergy, sinr=sinr, \
                        rxbuffer=[x.traceid for x in r.rxbuffer] )
                yield fsm.goto(self.DECODE, p)
            # ignore otherwise
            ignore = r.txdata in fsm.eventsFired
            ignore = ignore or (r.txdone in fsm.eventsFired)
            ignore = ignore or (r.rxdone in fsm.eventsFired)
            if ignore: pass
        return

    def DECODE(self, fsm, pkt):
        """DECODE state; monitor `radio` and manage packet decoding."""
        r = self.radio
        assert isinstance(r, Radio), "[DOT11A]: Cannot find radio in DECODE!"
        assert self.isbusy, "[DOT11A]: Carrier sense *not* busy in DECODE!"
        while fsm.active():
            # monitor events and RXD port
            yield self.RXD.recv(fsm, 1, \
                  renege=(r.rxdata, r.rxdone, r.txdata, r.txdone, self.detect) )
            # receive pkt -> apply error model and forward to upper layers
            if fsm.acquired(self.RXD):
                assert (len(fsm.got)==1), "[DOT11A]: Received unexpected " + \
                        "number of packets from 'RXD' port in DECODE state!"
                p = fsm.got[0]
                if (p is pkt):
                    payload = self.decode(p)
                    if payload:
                        self.log_recv(p)
                        self.cleananno(p)   # replace/remove unwanted annotations
                        p.remove_payload()
                        yield self.TXU.send(fsm, [payload])
                        yield hold, fsm, const.EPSILON  # pause before resuming
                    pkt = None
            # rxdone received before RXD -> interface dropped packet?
            if (r.rxdone in fsm.eventsFired):
                p = r.rxdone.signalparam
                if (p is pkt):
                    qlen = self.RXD.length
                    drop = ANNO.supports(p,'cif-drp') and p.getanno('cif-drp')
                    errmsg = "[DOT11A]: Unexpected rxdone received in DECODE!"
                    assert (qlen==0) and drop, errmsg
                    self.log_drop(p, drop="interface dropped packet")
                    pkt = None
            # rxdata -> drop new packets
            if (r.rxdata in fsm.eventsFired):
                p = r.rxdata.signalparam
                assert (p is not pkt), \
                        "[DOT11A]: Unexpected rxdata for pkt in DECODE!"
                self.log_drop(p, drop="ignore rxdata in DECODE")
            # detect -> drop detected packets
            if (self.detect in fsm.eventsFired):
                p = r.rxdata.signalparam
                assert (p is not pkt), \
                        "[DOT11A]: Unexpected detect for pkt in DECODE!"
                self.log_drop(p, drop="ignore detect in DECODE", \
                        decode=pkt.traceid)
            # ignore otherwise
            ignore = r.txdata in fsm.eventsFired
            ignore = ignore or (r.txdone in fsm.eventsFired)
            if ignore: pass
            # check if DECODE is done
            if pkt is None:
                yield fsm.goto(self.LISTEN)
        return

    def DETECT(self, fsm, pkt, thresh=None):
        """DETECT state; simulate physical layer frame detection.

        :param pkt: `Dot11A` packet being detected.
        :param thresh: SINR threshold used by `framedetect()`.

        This method signals the `detect` SimEvent when `pkt` is detected. Before
        calling `framedetect()`, this state will pause for the packet detection
        duration `detectdelay`.

        :note: Upon completion this state method sleeps the calling `FSM`.
        """
        r = self.radio
        assert isinstance(r, Radio), "[DOT11A]: Cannot find radio in DETECT!"
        assert (pkt in r.rxbuffer), \
                "[DOT11A]: Current pkt is not in radio rxbuffer!"
        yield hold, fsm, self.detectdelay
        if r.inreceive and self.framedetect(pkt, thresh=thresh):
            self.detect.signal(pkt)
        return

    def get_valid_rates(self):
        """Get list of rate enumerations supported by PHY.

        This method uses the PHY configuration to determine the list of valid
        rate enumerations supported by the PHY (i.e. `ntx`).
        """
        ntx = self.radio.ntx
        rates = range(len(DOT11A_DATARATE))
        return rates

    def get_datarate(self, r):
        """Get the data rate in bits-per-second (bps).

        :param r: Rate enumeration.

        This method gets data rate from `DOT11A_DATARATE`.
        """
        errmsg = "[DOT11APHY]: Invalid rate index (%s)!"%(r)
        assert (0<=r<len(DOT11A_DATARATE) ), errmsg
        return DOT11A_DATARATE[r]

    def log_send(self, p, *args, **kwargs):
        """Convenience method for logging send event."""
        if self.verbose>DOT11A_VERBOSE:
            if isinstance(p, Dot11A):
                kwargs['phy-rate'] = p.rate
                kwargs['length'] = p.length
            if p.hasanno('cif-txpower'):
                kwargs['cif-txpower'] = "%.2f dBm"%(p.getanno('cif-txpower') )
            if p.hasanno('cif-duration'):
                kwargs['cif-duration'] = time2usec(p.getanno('cif-duration') )
            self.log("snd", p, *args, **kwargs)

    def log_recv(self, p, *args, **kwargs):
        """Convenience method for logging receive event."""
        if self.verbose>DOT11A_VERBOSE:
            if isinstance(p, Dot11A):
                kwargs['phy-rate'] = p.rate
                kwargs['length'] = p.length
            if p.hasanno('phy-sinr'):
                kwargs['phy-sinr'] = "%.2f dB"%(p.getanno('phy-sinr') )
            if p.hasanno('rxpower'):
                kwargs['rxpower'] = "%.2f dBm"%(p.getanno('rxpower') )
            if p.hasanno('noisepower'):
                kwargs['noisepower'] = "%.2f dBm"%(p.getanno('noisepower') )
            if p.hasanno('cif-duration'):
                kwargs['cif-duration'] = time2usec(p.getanno('cif-duration') )
            if p.hasanno('dot11a-per'):
                kwargs['dot11a-per'] = "%.5g"%(p.getanno('dot11a-per') )
            if p.hasanno('crcerror'):
                crcerror = p.getanno('crcerror')
                if crcerror: kwargs['crc'] = "FAIL"
                else:        kwargs['crc'] = "OK"
            self.log("rcv", p, *args, **kwargs)

    def log_drop(self, p, *args, **kwargs):
        """Convenience method for logging drop event."""
        if self.verbose>DOT11A_VERBOSE:
            if isinstance(p, Dot11A):
                kwargs['phy-rate'] = p.rate
                kwargs['length'] = p.length
            if p.hasanno('phy-sinr'):
                kwargs['phy-sinr'] = "%.2f dBm"%(p.getanno('phy-sinr') )
            if p.hasanno('crcerror'):
                crcerror = p.getanno('crcerror')
                if crcerror: kwargs['crc'] = "FAIL"
                else:        kwargs['crc'] = "OK"
            self.log("drp", p, *args, **kwargs)
        self.cleananno(p)

    def calcduration(self, plen, rate=None):
        """Calculate duration of a packet of length `plen`.

        :param plen: Packet length (in bytes).
        :param rate: Optional rate index to denote modulation/coding scheme.
        :return: Duration of packet in seconds.

        If `rate` is not specified this method will use the base rate (i.e. rate
        index 0) to calculate the duration of the waveform.
        """
        if (rate is None): rate = 0
        # calculate duration
        nofdm = self.calcnofdm(plen, rate)  # number of data OFDM symbols
        d = DOT11A_TSHORT + DOT11A_TLONG + \
            DOT11A_TSIGNAL + DOT11A_TSYM*nofdm
        return d

    def calcnofdm(self, plen, rate):
        """Calculate the number of OFDM symbols in the data payload of a packet.

        :param plen: Length of packet in bytes).
        :param rate: Rate index to denote modulation/coding scheme.
        :return: Number of OFDM symbols in payload of 802.11n waveform.
        """
        # check packet length and rate
        assert not (plen<0), "[DOT11APHY]: Cannot compute " + \
                             "duration of negative length packet!"
        assert 0<=rate<len(DOT11A_NDBPS), \
               "[DOT11APHY]: Invalid rate index (%s)!"%(rate)
        blen = 8*plen
        nbits = 16 + blen + 6           # service bits + PSDU + tail-bits
        ndbps = DOT11A_NDBPS[rate]
        nofdm = ceil(1.0*nbits/ndbps)
        return int(nofdm)

class Dot11ARadio(Radio):
    """Radio interface for IEEE 802.11a physical layer.
    :CVariables:
     * `fc`: Carrier frequency.
     * `Lrx`: Analog system loss in receiver [default=`DOT11A_ANALOGLOSS`].
     * `bandwidth`: System bandwidth.
    """
    name = "IEEE 802.11a radio"
    fc = DOT11A_CARRIER
    Lrx = DOT11A_ANALOGLOSS
    bandwidth = DOT11A_BANDWIDTH
    def __init__(self, **kwargs):
        """Constructor."""
        Radio.__init__(self, **kwargs)
    def set_recvanno(self, p):
        """Overloaded to modify noisepower with noise figure.

        :return: Modified packet `p`.

        This method adds the `DOT11A_NOISEFIGURE` to the noisepower specified by
        `Radio`.
        """
        r = Radio.set_recvanno(self, p)
        if ANNO.supports(r, 'noisepower'):
            noisepower = r.getanno('noisepower')
            noisepower += DOT11A_NOISEFIGURE
            r.setanno('noisepower', noisepower)
        return r
