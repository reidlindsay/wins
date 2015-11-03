#!  /usr/bin/env python

"""
Implements IEEE 802.11n PHY protocol; contains `Dot11NPHY` and `Dot11NRadio`
classes.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-11-15 18:33:21 -0600 (Tue, 15 Nov 2011) $
* $LastChangedRevision: 5340 $

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

from SimPy.Simulation import SimEvent, hold, now, waitevent
from wins.protocol.csphy import CSPHY
from wins.channel.radio  import Radio
from wins.digital.mqam   import MQAM
from wins.digital.rcpc   import RCPC
from wins.ieee80211.dot11n_support import *
from wins.ieee80211.dot11n_dsp import Dot11N_DSP
from wins.channel.interface import strcollision
from wins.helper import time2usec, linear2db, db2linear, monitor_events
from wins.packet import ANNO, Packet
from wins.fsm    import FSM, Timer
from wins.crc    import CRC32
from wins import const

from numpy import random, inf, floor, ceil, sqrt
from numpy import exp, pi, sin, sinc
from numpy import arange

iorj = complex(0,1)

import sys

class Dot11NPHY(CSPHY):
    """Implementation of IEEE 802.11n physical layer protocol.

    Packet Format
    =============
     * STF: short training sequence (8.0 usec)
     * LTF: long training sequence (8.0 usec)
     * SIGNAL: PHY header (8.0 usec)
     * E-LTF: extension LTF (only used for added spatial streams)
     * PAYLOAD: data symbols

    ------------------------------------------
    | STF | LTF | SIGNAL | [E-LTF] | PAYLOAD |
    ------------------------------------------

    **See `Dot11N` for more on composition of header, i.e. SIGNAL field.**
    

    Modulation and Coding Parameters
    ================================

    The table below describes modulation and convolutional coding parameters for
    IEEE 802.11n.

    =========== ============= ============ ====== ======= ======= =========
    Rate Index    Data Rate    Modulation   Nbps     M     Nsts   Code rate
    ----------- ------------- ------------ ------ ------- ------- ---------
         0        6.5 Mbps       BPSK         1      2       1        1/2
         1         13 Mbps       QPSK         2      2       1        1/2
         2       19.5 Mbps       QPSK         2      4       1        3/4
         3         26 Mbps      16-QAM        4      4       1        1/2
         4         39 Mbps      16-QAM        4     16       1        3/4
         5         52 Mbps      64-QAM        6     16       1        2/3
         6       58.5 Mbps      64-QAM        6     64       1        3/4
         7         65 Mbps      64-QAM        6     64       1        5/6
         8         13 Mbps       BPSK         1      2       2        1/2
         9         26 Mbps       QPSK         2      2       2        1/2
        10         39 Mbps       QPSK         2      4       2        3/4
        11         52 Mbps      16-QAM        4      4       2        1/2
        12         78 Mbps      16-QAM        4     16       2        3/4
        13        104 Mbps      64-QAM        6     16       2        2/3
        14        117 Mbps      64-QAM        6     64       2        3/4
        15        130 Mbps      64-QAM        6     64       2        5/6
        ...         ...          ...         ...    ...     ...       ...
        23        195 Mbps      64-QAM        6     64       3        3/4
        ...         ...          ...         ...    ...     ...       ...
        31        260 Mbps      64-QAM        6     64       4        3/4
    =========== ============= ============ ====== ======= ======= =========

    :note: Rate index is used in Packets when referring to rate. Use
           `DOT11N_DATARATE` to convert this to a bitrate value.

    Dot11NPHY and Ports
    ===================
    Every `Dot11NPHY` has the following ports:

        1. "TXU" - sends decoded packets to an upstream element.
        #. "RXU" - receives traffic to send from an upstream element.
        #. "TXD" - sends encoded packets to a downstream element.
        #. "RXD" - receives traffic from a downstream element.

    The upstream element of a `Dot11NPHY` is usually a `CSMAC` (or subclass),
    and the downstream element is usually a `ChannelInterface` (or subclass).

    Dot11NPHY and Annotations
    =========================
    This class uses/sets the following annotations:

    ============== ===========================================================
    Name            Description
    ============== ===========================================================
    cif-txpower     Power of transmitted packets (in dBm).
    -------------- -----------------------------------------------------------
    cif-duration    Duration of packet is marked on outgoing packets using
                    the `duration()` method during `SEND`.
    -------------- -----------------------------------------------------------
    cif-collision   List of packets that arrived at the same time as the
                    current packet (i.e. collided with the marked packet).
    -------------- -----------------------------------------------------------
    phy-rate        Rate index used to specify coding and modulation scheme.
    -------------- -----------------------------------------------------------
    phy-sinr        SINR metric used to indicate channel quality
                    (for upper layers). 
    -------------- -----------------------------------------------------------
    dot11n-detect   Indicates whether or not the physical layer packet
                    detection was successful (see `framedetect()`).
    -------------- -----------------------------------------------------------
    dot11n-start    Timestamp indicating what time the PHY *thinks* the
                    packet starts (see `framedetect()`).
    -------------- -----------------------------------------------------------
    dot11n-header   Indicates whether or not the physical layer header
                    decoding was successful (see `decode_header()`).
    -------------- -----------------------------------------------------------
    dot11n-per      Packet-error rate (PER) modeled by the physical layer
                    (the `decode_header()` will initially set this, but
                    `decode_data()` will overwrite it as necessary).
    -------------- -----------------------------------------------------------
    dot11n-fdr      Frame detection rate (FDR) modeled by the physical layer
                    (this annotation is set during `framedetect()`).
    -------------- -----------------------------------------------------------
    dot11n-sinr     Signal-to-interference-plus-noise ratio (SINR) computed   
                    in `calcper_header()` and overwritten by `calcper_data()`.
    -------------- -----------------------------------------------------------
    dot11n-model-*  Model-prediction annotations. Duplicate annotations kept
                    for tracking model-predicted outcomes. These annotations
                    include: 'dot11n-model-detect', 'dot11n-model-fdr',
                    'dot11n-model-header', 'dot11n-model-her' (header error
                    rate), 'dot11n-model-error' (payload decoding error),
                    'dot11n-model-per' (payload error rate).
    ============== ===========================================================

    :IVariables:
     * `detect`: Internal SimEvent signalled when a packet is detected.
     * `header`: Internal SimEvent signalled when a packet header is decoded.
     * `mod`: `MQAM`, child member for modulation.
     * `coder`: `RCPC`, child member for coding.

    :CVariables:
     * `usewaveform`: If flag is true, use waveform level simulation.
     * `detectdelay`: Duration required for packet detection.
    """
    name = "IEEE 802.11n"
    tracename = "80211N"
    detectdelay = DOT11N_TDETECT
    usewaveform = True
    def __init__(self, **kwargs):
        """Constructor."""
        self.detect = SimEvent(name=".detect")
        self.header = SimEvent(name=".header")
        self.usecapture = False
        # call superclass constructor
        CSPHY.__init__(self, **kwargs)
        self.detect.name = "%s%s"%(self.name,".detect")
        self.header.name = "%s%s"%(self.name,".header")

    def configure(self, detectdelay=None, usecapture=False, **kwargs):
        """Call `CSPHY.configure()`; add ports, `FSM`, and other members.

        :param detectdelay: Delay for detecting new frame.
        :param usecapture: Enable capture effect [default=False].
        :param kwargs: Additional keywords passed to `Dot11N_DSP` constructor.
        """
        CSPHY.configure(self, **kwargs)
        # set parameters
        cls = self.__class__
        if detectdelay is None: self.detectdelay = cls.detectdelay
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
        self.usecapture = usecapture
        # set up waveform members
        dsp = self.newchild("dsp", Dot11N_DSP, tracename=self.tracename+".DSP", **kwargs)

    def connect(self, radio):
        """Convenience method to connect PHY to an RF front end (i.e. `Radio`).

        **Overload this method to change how a connection is made.**

        :note: This method also calls `set_radio()` to set the appropriate
               pointer for the physical layer.
        """
        assert isinstance(radio, Radio), \
                "[DOT11N]: Cannot connect to non-Radio!"
        self.set_radio(radio)
        self.TXD.connect(radio.getport("RXU"))
        radio.getport("TXU").connect(self.RXD)

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
        yield self.RXU.recv(fsm, 1)
        assert fsm.acquired(self.RXU) and (len(fsm.got)==1), \
                "[DOT11N]: Error receiving from RXU port in SEND()!"
        p = fsm.got[0]
        # simulate encoding and Tx processing time
        w = self.encode(p)
        duration = p.getanno('cif-duration')
        yield hold, fsm, self.txproctime(p)
        # send waveform and simulate duration
        self.log_send(w, duration=time2usec(duration) )
        yield self.TXD.send(fsm, [w])
        errmsg = "[DOT11N]: Error sending to TXD in SEND!"
        assert fsm.stored(self.TXD), errmsg
        yield hold, fsm, duration
        # continue in SEND
        yield fsm.goto(self.SEND)

    def encode(self, p, **kwargs):
        """Encapsulate MPDU in a `Dot11N` packet.

        :param p: MPDU to encode.
        :param kwargs: Keywords arguments passed to `set_sendanno()`.
        :return: `Dot11N` packet with MPDU as payload.
        """
        # set any missing annotations
        self.set_sendanno(p, **kwargs)
        # get parameters
        plen = len(p)
        ntx = self.radio.Ntx
        rate = p.getanno('phy-rate')
        txpower = p.getanno('cif-txpower')
        # for waveform level simulation ...
        if self.usewaveform:
            self.dsp.encode(p, rate, ntx)
        # encap in Dot11N
        w = Dot11N(rate=rate, length=plen)
        w.add_payload(p)
        w.updatecrc()
        return w

    def set_sendanno(self, p, rate=None, txpower=None):
        """Set default annotations for missing PHY annotations.

        :param rate: Optional rate index to denote modulation/coding scheme.
        :param txpower: Optional transmit power (in dBm).

        This method will make sure that the following annotations are set:
            * rate [default=0]
            * txpower [default=`DOT11N_MAXPOWER`]
            * duration [from `duration()`]

        If the `rate` parameter is not specified and no 'phy-rate' annotation is
        found in packet `p`, this method will use the base rate (i.e. zero) to
        encode packet `p`.

        If the `txpower` parameter is not specified and no 'cif-txpower'
        annotation is found in packet `p`, this method will use
        `DOT11N_MAXPOWER` as the default transmit power.
        """
        assert ANNO.supported(p), "[DOT11NPHY]: ANNO not supported!"
        # get parameters from p (or use defaults)
        if p.hasanno('phy-rate'):
            if (rate is None): rate = p.getanno('phy-rate')
        if p.hasanno('cif-txpower'):
            if (txpower is None): txpower = p.getanno('cif-txpower')
        if (rate is None): rate = 0
        if (txpower is None): txpower = DOT11N_MAXPOWER
        # check parameters
        txpower = min(txpower, DOT11N_MAXPOWER)
        ntx = self.radio.Ntx
        errmsg = "[DOT11N]: Invalid Ntx %d not in [%d, %d]!"%(ntx,1,4)
        assert (0<ntx<=4), errmsg
        errmsg = "[DOT11N]: Invalid rate %d, ntx = %d "%(rate, ntx)
        assert (0<=rate<8*ntx), errmsg
        # set annotations
        plen = len(p)
        duration = self.duration(plen, rate=rate)
        p.setanno('phy-rate', rate)
        p.setanno('cif-txpower', txpower)
        p.setanno('cif-duration', duration)
        # remove some annotations from outgoing packets
        rannolist = ['dot11n-detect', 'dot11n-fdr', 'dot11n-start', \
                     'dot11n-header', 'dot11n-sinr', 'dot11n-per', \
                     'dot11n-per', 'dot11n-error', \
                     'dot11n-model-detect', 'dot11n-model-fdr', \
                     'dot11n-model-header', 'dot11n-model-her', \
                     'dot11n-model-per', 'dot11n-model-error']
        for a in rannolist:
            if p.hasanno(a): p.delanno(a)
        return p

    def duration(self, p, rate=None):
        """Calculate duration of packet `p` using `calcduration()` method.

        :param p: Packet to compute duration for (or packet length in bytes).
        :param rate: Optional rate index to denote modulation/coding scheme.
        :return: Duration of packet in seconds.

        This method checks for the 'phy-rate' annotation, which is a rate
        index to denote the coding and modulation scheme to be used. If a
        'phy-rate' annotation is not found (or specified as a parameter),
        this method uses the base rate (i.e. rate 0) to calculate the duration
        of the waveform.
        """
        plen = p
        if isinstance(p, Packet): plen = len(p)
        if ANNO.supports(p, 'phy-rate'):
            if (rate is None): rate = p.getanno('phy-rate')
        if (rate is None): rate = 0
        return self.calcduration(plen, rate)

    def calcduration(self, plen, rate=None):
        """Calculate duration of a packet of length `plen`.

        :param plen: Packet length (in bytes).
        :param rate: Optional rate index to denote modulation/coding scheme.
        :return: Duration of packet in seconds.

        If `rate` is not specified this method will use the base rate (i.e. rate
        index 0) to calculate the duration of the waveform.
        """
        if (rate is None): rate = 0
        # check packet length and rate
        errmsg = "[DOT11N]: Packet has negative length!"
        assert not (plen<0), errmsg
        errmsg = "[DOT11N]: Invalid rate index (%s)!"%(rate)
        assert 0<=rate<len(DOT11N_NDBPS), errmsg
        # calculate duration
        nofdm = self.calcnofdm(plen, rate)  # number of data OFDM symbols
        ndltf = 4                           # number of extension LTF symbols
        if (rate<16): ndltf = 2
        if (rate<8):  ndltf = 1
        neltf = self.radio.Ntx          # add non-standard extension LTF symbols
        if (neltf<2): neltf = 0         # Ntx=0 -> Neltf=0
        d = DOT11N_TSHORT + DOT11N_TLONG + DOT11N_TSIGNAL + \
            DOT11N_TSYM*(ndltf + neltf - 1) + DOT11N_TSYM*nofdm
        return d

    def calcnofdm(self, plen, rate):
        """Calculate the number of OFDM symbols in the data payload of a packet.

        :param plen: Length of packet in bytes).
        :param rate: Rate index to denote modulation/coding scheme.
        :return: Number of OFDM symbols in payload of 802.11n waveform.
        """
        # check packet length and rate
        assert not (plen<0), "[DOT11N]: Cannot compute " + \
                             "duration of negative length packet!"
        assert 0<=rate<len(DOT11N_NDBPS), \
               "[DOT11N]: Invalid rate index (%s)!"%(rate)
        blen = 8*plen
        nbits = 16 + blen + 6           # service bits + PSDU + tail-bits
        ndbps = DOT11N_NDBPS[rate]
        nofdm = ceil(1.0*nbits/ndbps)   # number of data OFDM symbols
        return int(nofdm)

    def calclength(self, duration, rate):
        """Calculate length of packet in bytes.

        :param duration: Duration of packet in seconds.
        :param rate: Rate enumeration.

        This method uses the local configuration to compute the packet length.

        :note: This method does not compute the actual number of bytes, rather
               the maximum number of bytes that might be in the packet.
        """
        errmsg = "[DOT11NPHY]: Invalid rate index (%s)!"%(r)
        assert (0<=r<len(DOT11N_DATARATE) ), errmsg
        # calculate duration of payload
        ndltf = 4                           # number of extension LTF symbols
        if (rate<16): ndltf = 2
        if (rate<8):  ndltf = 1
        neltf = self.radio.Ntx          # add non-standard extension LTF symbols
        if (neltf<2): neltf = 0         # Ntx=0 -> Neltf=0
        d = duration
        d -= DOT11N_TSHORT + DOT11N_TLONG + DOT11N_TSIGNAL + \
             DOT11N_TSYM*(ndltf + neltf - 1)
        # calculate number of bytes from number of OFDM symbols in payload
        nofdm = int(d/DOT11N_TSYM)
        ndbps = DOT11N_NDBPS[rate]      # number of data bits per OFDM symbol
        nbits = nofdm*ndbps
        blen = nbits - (16 + 6)         # strip service and tail bits
        plen = int(blen/8)
        return plen

    def LISTEN(self, fsm, header=False):
        """LISTEN state; monitor `radio` and manage packet detection.

        :param header: Boolean flag; if true, indicates that LISTEN is currently
                       decoding the header of a packet.

        The `header` argument is used to allow LISTEN to be repurposed for
        decoding of the header or signal field of an 802.11n packet.
        """
        r = self.radio
        assert isinstance(r, Radio), "[DOT11N]: Cannot find radio in LISTEN!"
        # check rxenergy -> set csbusy?
        Erx = r.rxenergy()
        rxhigh = r.inreceive and (r.rxenergy()>DOT11N_CSTHRESHOLD)
        rxenergy = "%.2f dBm"%(Erx)
        if rxhigh:
            self.set_csbusy(rxenergy="high, "+rxenergy, \
                            rxbuffer=[x.traceid for x in r.rxbuffer]  )
        elif (not header):
            # only set idle when not decoding header
            self.set_csidle(rxenergy=rxenergy )
        # monitor events and RXD port
        evlist =  (r.rxdata, r.rxdone, r.txdata, r.txdone)
        evlist += (self.detect, self.header)
        yield self.RXD.recv(fsm, 1, renege=evlist)
        # RXD -> ignore incoming packets in LISTEN
        if fsm.acquired(self.RXD):
            errmsg = "[DOT11N]: 'RXD' port received wrong number of packets!"
            assert (len(fsm.got)==1), errmsg
            p = fsm.got[0]
            drop = "not detected in LISTEN"
            if ANNO.supports(p, 'dot11n-drop'):
                drop = p.getanno('dot11n-drop')
            self.log_drop(p, drop=drop)
        # r.rxdata -> start DETECT thread
        if r.rxdata in fsm.eventsFired:
            p = r.rxdata.signalparam
            if self.usewaveform: self.dsp.addrxinput(p)
            self.start_detect(p)
        # detect -> set csbusy -> start HEADER thread
        if (self.detect in fsm.eventsFired) and (not header):
            p = self.detect.signalparam
            # check for annotations
            for a in ['dot11n-detect', 'dot11n-fdr', 'dot11n-start']:
                errmsg = "[DOT11N]: Cannot find '%s' annotation!"%(a)
                assert ANNO.supports(p, a), errmsg
            errmsg = "[DOT11N]: Error! 'dot11n-detect' " + \
                     "annotation not set in detected packet!"
            assert p.getanno('dot11n-detect'), errmsg
            # log short-training SINR
            sinr  = "%.2f dB"%(self.sinr_part(p, DOT11N_STF) )
            rxenergy = "%.2f dBm"%(r.rxenergy() )   # current value in radio
            self.set_csbusy(p, detect=True, rxenergy=rxenergy)
            # log detect event
            kwargs = {'dot11n-fdr': "%.5g"%(p.getanno('dot11n-fdr')), \
                      'dot11n-start': "%s"%(p.getanno('dot11n-start')) }
            kwargs.update(self.get_dot11n_anno(p))
            self.log("detect", p, rxenergy=rxenergy, sinr=sinr, \
                    rxbuffer=[x.traceid for x in r.rxbuffer], **kwargs)
            self.start_header(p)
            header = True   # continue to LISTEN(header=True)
        # header -> goto DECODE?
        if self.header in fsm.eventsFired:
            assert (header), "[DOT11N]: Unexpected header signalled!"
            assert self.isbusy, "[DOT11N]: Carrier sense *not* busy!"
            p = self.header.signalparam
            # check for annotations
            for a in ['dot11n-header']:
                errmsg = "[DOT11N]: Cannot find '%s' annotation!"%(a)
                assert ANNO.supports(p, a), errmsg
            # did header decoding succeed? -> goto DECODE
            hstatus = p.getanno('dot11n-header')
            if (hstatus=="success"):
                yield fsm.goto(self.DECODE, p)
            else:
                self.log_drop(p, drop=hstatus)
                p.setanno('dot11n-drop', hstatus)
                header = False  # continue to LISTEN(header=False)
        # ignore other events
        ignore = r.txdata in fsm.eventsFired
        ignore = ignore or (r.txdone in fsm.eventsFired)
        ignore = ignore or (r.rxdone in fsm.eventsFired)
        if ignore: pass
        # remain in LISTEN
        yield fsm.goto(self.LISTEN, header)

    def start_detect(self, p, **kwargs):
        """Spawn thread to run packet detection."""
        fname = "detect(%s)"%(p._id)
        f = FSM()
        f.goto(self.DETECT, p)
        f.start(**kwargs)
        return f

    def start_header(self, p, **kwargs):
        """Spawn thread to run header decoding."""
        fname = "header(%s)"%(p._id)
        f = FSM()
        f.goto(self.HEADER, p)
        f.start(**kwargs)
        return f

    def DETECT(self, fsm, p, *args, **kwargs):
        """DETECT state; simulate physical layer frame detection.

        :param p: `Dot11N` packet being detected.
        :param args: Additional arguments passed to `framedetect()`.
        :param kwargs: Additional keywords passed to `framedetect()`.

        This method signals the `detect` SimEvent when `p` is detected. Before
        calling `framedetect()`, this state will pause for the packet detection
        duration `detectdelay`.

        :note: Upon completion this state method sleeps the calling `FSM`.
        """
        r = self.radio
        assert isinstance(r, Radio), "[DOT11N]: Cannot find radio in DETECT!"
        errmsg = "[DOT11N]: %s is not in radio rxbuffer! @ %.8f"%(p.traceid, now())
        assert (p in r.rxbuffer), errmsg
        yield hold, fsm, self.detectdelay
        # update SINR (and SINR heap)
        sinr = self.sinr_part(p, DOT11N_STF, force=True)
        # check if packet is detected
        if r.inreceive and self.framedetect(p, *args, **kwargs):
            self.detect.signal(p)
        else:
            ignore = "detection failed"
            if not r.inreceive: ignore += "|in transmit mode"
            self.log("ignore", p, ignore=ignore)
        yield fsm.stop()

    def framedetect(self, p, threshmin=None, threshmax=None):
        """Apply packet detection model for detecting training sequence; based
        on signal-to-interference-and-noise ratio (SINR).

        :param p: `Dot11N` packet being received.
        :param threshmin: Below this packet detection SINR threshold (in dB),
                          the incoming frame cannot be detected
                          [default=`DOT11N_FDTHRESHOLD0`].
        :param threshmin: Above this packet detection SINR threshold (in dB),
                          the incoming frame will alway be detected
                          [default=`DOT11N_FDTHRESHOLD1`].
        :return: Boolean flag; if true, packet detection was successful.

        This method checks to see if the packet can be detected by the packet
        detector. The method will determine the frame detection rate (FDR) and
        perform coin flipping to determine if the packet is detected. This
        methods sets the 'dot11n-fdr' annotation to indicate the FDR (which is
        determined as a function of SINR, `threshmin`, and `threshmax`).
        Finally, this method will mark the 'dot11n-detect' annotation to
        indicate the success (or failure) of the frame detection.

        If the packet is detected, this method also sets the 'dot11n-start'
        annotation to indicate what time the PHY thinks the packet starts.

        **Overload this method to change how frame detection works.**

        :note: If `p` is not a `Dot11N` packet, this method will set the
               'dot11n-detect' annotation to false and return false.
        """
        # check if packet is valid
        if not isinstance(p, Dot11N):
            p.setanno('dot11n-fdr', 0.0)
            p.setanno('dot11n-detect', False)
            return False
        # Calculate SINR over SHORT training sequence
        detect = False
        sinr = self.sinr_part(p, DOT11N_STF)
        # determine detection thresholds
        try:    threshmin = float(threshmin)
        except: threshmin = DOT11N_FDTHRESHOLD0
        try:    threshmax = float(threshmax)
        except: threshmax = DOT11N_FDTHRESHOLD1
        errmsg ="[DOT11N]: Invalid frame detection " + \
                "thresholds (%.2f dB, %.2f dB)"%(threshmin, threshmax)
        assert not (threshmin>threshmax), errmsg
        # calculate frame detection rate (FDR)
        if (threshmin<>threshmax):
            fdr = (sinr - threshmin)/(threshmax - threshmin)
        else:
            fdr = float(sinr>threshmin)
        fdr = max(0.0, min(fdr, 1.0))
        tstart = p.getanno('cif-rxts')
        # make model-prediction of detection outcome
        #  -> use FDR to determine if packet was detected
        if not ((fdr<>0) and (fdr<>1)): detect = (fdr>0)
        else: detect = (random.uniform(0,1)<fdr)
        # save model-prediction
        model_detect, model_fdr = detect, fdr
        # waveform-level frame detection?
        if self.usewaveform:
            # check number of collider before using detection gate
            ncollision = 0
            if p.hasanno('cif-collision'): ncollision = len(p.getanno('cif-collision'))
            # use SINR as a final gate for detection when collisions are present
            if (sinr<0) and (ncollision>0):
                detect, fdr = False, 0.0
            else:
                # use waveform-level decoding to decode header
                param = self.dsp.decode_header(p, False)  # don't check header crc
                # -> get header parameters from receiver after decoding
                detect = param['detect']
                tstart = param['start']
        # mark annotations
        p.setanno('dot11n-sinr', sinr)
        p.setanno('dot11n-detect', detect)
        p.setanno('dot11n-fdr', fdr)
        p.setanno('dot11n-start', tstart)
        # mark model-prediction annotations
        p.setanno('dot11n-model-detect', model_detect)
        p.setanno('dot11n-model-fdr', model_fdr)
        return detect

    def HEADER(self, fsm, p, *args, **kwargs):
        """HEADER state; simulate header decoding.

        :param p: `Dot11N` packet being decoding.

        This method signals the `header` SimEvent after completing the header
        decoding for `pkt`. This state will also set the 'dot11n-header'
        annotation to indicate the success or failure of header decoding.

        :note: This state will throw an assertion error if packet `p` was not
               detected prior to be passed to the SEM.
        """
        # check annotations
        for a in ['cif-rxts', 'dot11n-detect']:
            errmsg = "[DOT11N]: Cannot find '%s' annotation!"%(a)
            assert ANNO.supports(p, a), errmsg
        assert p.getanno('dot11n-detect')
        assert isinstance(p, Dot11N)
        # check current time
        tnow = now()
        t0 = p.getanno('cif-rxts')
        t1 = t0 + self.detectdelay  # HEADER only called after detectdelay
        errmsg = "[DOT11N]: HEADER started at invalid time! " + \
                 "t1 = %.7f, now = %.7f, diff = %.6g"%(t1, tnow, abs(t1-tnow))
        assert (abs(t1-tnow)<2*const.EPSILON), errmsg
        # pause to wait until the end of the header (SIGNAL field)
        hdelay = DOT11N_TSHORT + DOT11N_TLONG + DOT11N_TSIGNAL - (t1-t0)
        # simulate capture effect?
        if self.usecapture and (hdelay>0):
            timer = Timer(hdelay, start=True)
            yield waitevent, fsm, (timer.done, self.detect)
            if self.detect in fsm.eventsFired:
                # capture new packet
                pkt = self.detect.signalparam
                timer.halt()
                # drop current packet
                p.setanno('dot11n-header', "capture new packet")
                self.header.signal(p)
                yield hold, fsm, const.EPSILON  # pause to let other threads run
                # re-fire detect with new packet
                self.detect.signal(pkt)
                return
            else:
                # continue normal processing
                assert (timer.done in fsm.eventsFired)
        else:
            # pause
            yield hold, fsm, hdelay
        # update SINR (and SINR heap)
        sinr = self.sinr_part(p, DOT11N_SIGNAL, force=True)
        # decode header (SIGNAL field)
        header = self.decode_header(p)
        p.setanno('dot11n-header', header)
        self.header.signal(p)
        return

    def decode_header(self, p):
        """Apply physical layer model for header decoding.

        :param p: `Dot11N` packet being received.
        :return: Status string; "success" if header was successfully decode,
                 otherwise indicates error.

        This method uses `mod` and `coder` to determine the error
        characteristics of the header decoding process. The following conditions
        must be met in order to have successful header decoding:

            * packet `p` must be a `Dot11N` packet,
            * the 'dot11n-detect' annotation must be true,
            * the header decoding must succeed,
            * and all header parameters must be valid.
        
        This method marks the 'dot11n-header' annotation to mark success (or
        reason for failure) of header decoding, and sets the 'dot11n-per'
        annotation to indicate the probability of error in decoding the header.

        :note: This method assumes that the packet `p` is a valid `Dot11N`
               packet and that it has been detected, i.e. it has the
               'dot11n-detect' annotation set properly.
        """
        assert isinstance(p, Dot11N)
        assert p.getanno('dot11n-detect')
        # get header parameters
        plen = len(p.payload)
        rate, length = p.rate, p.length
        okpar = self.parity(p)
        # check header parameters
        okrate = (0<=rate<len(DOT11N_DATARATE) )
        oklen  = (length == plen)
        if not (okrate and oklen and okpar):
            header = "header parameters failed"
            p.setanno('dot11n-header', header)
            return header
        # calculate PER for header decode of Dot11N
        per = self.calcper_header(p, force=False)
        sinr = p.getanno('dot11n-sinr')
        # make model-prediction of header decoding outcome
        #  -> simulate header decoding: use packet-level abstraction (coin flipping)
        hsuccess = "success"
        hfail    = "header decoding failed"
        header   = hsuccess
        if (random.uniform(0,1)<per): header = hfail
        # save model-prediction
        model_header, model_her = header, per
        # waveform-level header decoding?
        if self.usewaveform:
            # use waveform-level decoding to decode header
            param = self.dsp.decode_header(p, True)
            # -> get header parameters from receiver after decoding
            header = param['status']
            rate   = param['rate']
            length = param['length']
            okpar  = param['crcok']
            if okpar: header = hsuccess
            else:     header = hfail
        # mark annotations
        p.setanno('dot11n-header', header)
        p.setanno('dot11n-sinr', sinr)
        p.setanno('dot11n-per', per)
        # mark model-prediction annotations
        p.setanno('dot11n-model-header', model_header)
        p.setanno('dot11n-model-her', model_her)
        return header

    def DECODE(self, fsm, pkt):
        """DECODE state; monitor `radio` and manage packet decoding."""
        r = self.radio
        assert isinstance(r, Radio), "[DOT11N]: Cannot find radio in DECODE!"
        assert self.isbusy, "[DOT11N]: Carrier sense *not* busy in DECODE!"
        # monitor events and RXD port
        evlist =  (r.rxdata, r.rxdone, r.txdata, r.txdone)
        evlist += (self.detect, self.header)
        yield self.RXD.recv(fsm, 1, renege=evlist)
        # receive pkt -> apply error model and forward to upper layers
        if fsm.acquired(self.RXD):
            errmsg = "[DOT11N]: 'RXD' port received wrong number of packets!"
            assert (len(fsm.got)==1), errmsg
            p = fsm.got[0]
            if (p is pkt):
                # update SINR (and SINR heap)
                sinr = self.sinr_part(p, DOT11N_PAYLOAD, force=True)
                # decode data (payload)
                payload = self.decode_data(p)
                if payload:
                    self.log_recv(p)
                    yield hold, fsm, self.rxproctime(p)
                    self.set_recvanno(payload)
                    p.remove_payload()
                    yield self.TXU.send(fsm, [payload])
                pkt = None
        # r.rxdata -> drop new packets
        if r.rxdata in fsm.eventsFired:
            p = r.rxdata.signalparam
            if self.usewaveform: self.dsp.addrxinput(p)
            errmsg = "[DOT11N]: Unexpected rxdata for %s "%(p.traceid) + \
                     "while in DECODE for %s!"%(pkt.traceid)
            assert (p is not pkt), errmsg
            if self.usecapture:
                self.start_detect(p)    # start detection
            else:
                self.log_drop(p, drop="ignore rxdata in DECODE")    #ignore
        # r.rxdone received before RXD -> interface dropped packet?
        if (r.rxdone in fsm.eventsFired):
            p = r.rxdone.signalparam
            if (p is pkt):
                qlen = self.RXD.length
                drop = p.hasanno('cif-drp') and p.getanno('cif-drp')
                errmsg = "[DOT11N]: Unexpected rxdone received in DECODE!"
                assert (qlen==0) and drop, errmsg
                self.log_drop(p, drop="interface dropped packet")
                pkt = None
        # detect -> ignore? -> drop detected packets
        if (self.detect in fsm.eventsFired):
            p = self.detect.signalparam
            try:    pktid = pkt.traceid
            except: pktid = None
            errmsg = "[DOT11N]: Unexpected detect from %s "%(p.traceid) + \
                     " while in decode for %s in DECODE!"%(pktid)
            assert (p is not pkt), errmsg
            if self.usecapture:
                # capture effect -> drop current packet
                drop = "capture new packet in DECODE"
                self.log_drop(pkt, drop=drop, capture=p.traceid)
                pkt = None
                # re-fire detect signal
                self.detect.signal(p)   # FIXME: Might need to pause first?
            else:
                drop = "ignore detect in DECODE"
                self.log_drop(p, drop=drop, decode=pkt.traceid)
        # header -> ignore in DECODE
        if self.header in fsm.eventsFired:
            p = self.header.signalparam
            drop = "ignore header in DECODE"
            self.log_drop(p, drop=drop, decode=pkt.traceid)
        # ignore other events
        ignore = r.txdata in fsm.eventsFired
        ignore = ignore or (r.txdone in fsm.eventsFired)
        if ignore: pass
        # check if DECODE is done
        if pkt is None: yield fsm.goto(self.LISTEN)
        yield fsm.goto(self.DECODE, pkt)

    def set_recvanno(self, p, sinr=None):
        """Set relevant PHY annotations for upper layers to use.

        :param sinr: Signal-to-interference-and-noise ratio (in dB).

        By default this method will copy the SINR from the 'dot11n-sinr'
        annotation to set the 'phy-sinr' annotation.
        """
        assert ANNO.supported(p), "[DOT11NPHY]: ANNO not supported!"
        # check for SINR annotation
        if sinr is None:
            assert p.hasanno('dot11n-sinr')
            sinr = p.getanno('dot11n-sinr')
        # set annotations
        p.setanno('phy-sinr', sinr)
        # replace/remove unwanted annotations
        self.cleananno(p)
        return p

    def cleananno(self, p):
        """Overload to clean up waveform annotations as well."""
        pkt = CSPHY.cleananno(self, p)
        return self.dsp.delrxwaveform(pkt)

    def decode_data(self, p):
        """Apply physical layer model for decoding data.

        :param p: `Dot11N` packet being received.
        :return: Decoded packet payload marked with appropriate annotations.

        This method uses `mod` and `coder` to simulate packet errors. After
        verifying the header decoding and header parameters, this method applies
        the appropriate error model to simulate the decoding of the packet
        payload. If `usewaveform` is set, this method will call the
        decode_data() method on `dsp`.

        This method also overwrites the 'dot11n-per' annotation with the packet
        error probability for the payload.
        """
        # verify detection and header decoding succeeded
        errmsg = "[DOT11N]: Cannot decode payload of non-Dot11N packet!"
        assert isinstance(p, Dot11N), errmsg
        for a in ['dot11n-header', 'dot11n-detect']:
            errmsg = "[DOT11N]: Cannot find '%s' annotation!"%(a)
            assert ANNO.supports(p, a), errmsg
        errmsg = "[DOT11N]: Cannot decode payload that has not been detected!"
        assert p.getanno('dot11n-detect'), errmsg
        header = p.getanno('dot11n-header')
        errmsg = "[DOT11N]: Header decoding failed before decode_data()!"
        assert (header=="success"), errmsg
        # get packet payload and continue decoding
        pkt = p.payload
        error = True    # assume error by default
        # determine if data decoding error occurred and return pkt
        per = self.calcper_data(p, force=False)
        sinr = p.getanno('dot11n-sinr')
        # make model-prediction of data decoding outcome
        #  -> simulate data decoding: use packet-level abstraction (coin flipping)
        if (random.uniform(0,1)>per): error = False
        # save model-prediction
        model_error, model_per = error, per
        # waveform-level data decoding?
        if self.usewaveform:
            # use waveform-level decoding
            param = self.dsp.decode_data(p)
            data  = param['data']
            error = param['error']
        # mark errors and annotations
        if error: self.seterror(pkt)
        pkt.setanno('dot11n-per', per)
        pkt.setanno('dot11n-sinr', sinr)
        pkt.setanno('dot11n-error', error)
        # mark model-prediction annotations
        pkt.setanno('dot11n-model-error', model_error)
        pkt.setanno('dot11n-model-per', model_per)
        return pkt

    def txproctime(self, p):
        """Calculate processing time for encoding waveform."""
        return const.EPSILON

    def rxproctime(self, p):
        """Calculate processing time for decoding waveform."""
        return const.EPSILON

    def calcper_header(self, p, **kwargs):
        """Calculate probability of error for header decoding.

        :param p: Packet being decoded.
        :param kwargs: Additional keywords arguments passed to `sinr_heap()`
                       (or `sinr()`).
        :return: PER for header decoding.

        This method sets the 'dot11n-sinr' and 'dot11n-per' annotations. The
        operation of this method depends on `DOT11N_USE_PIECEWISE_PER`.
        """
        for a in ['cif-rxts']:
            errmsg = "[DOT11N]: Cannot find '%s' annotation!"%(a)
            assert ANNO.supports(p, a), errmsg
        # get calcper arguments
        rmsdelay, bandwidth = 0, 1
        if p.hasanno('dot11n-rmsdelay'):
            bandwidth = self.radio.bandwidth
            rmsdelay  = p.getanno('dot11n-rmsdelay')
        calcperargs = {'rmsdelay':rmsdelay, 'bandwidth':bandwidth}
        # calculate PER using appropriate method
        plen = len(p.payload)
        hlen = len(p)-plen
        assert (hlen == len(Dot11N())), "[DOT11N]: Invalid header found!"
        if DOT11N_USE_PIECEWISE_PER:
            sinrheap = self.sinr_heap(p, **kwargs)
            t0 = p.getanno('cif-rxts') + DOT11N_TSHORT + DOT11N_TLONG
            t1 = t0 + DOT11N_TSIGNAL
            xheap = [(max(ta,t0), min(tb,t1), sinr) for (ta,tb,sinr) \
                        in sinrheap if (ta<t1) and (tb>t0)]
            errmsg = "[DOT11N]: Unable to find valid data from SINR heap!"
            assert (len(xheap)>0), errmsg
            # calculate piecewise PER and average SINR
            psuccess, stot = 1.0, 0.0
            for ta, tb, sinr in xheap:
                alpha = (tb-ta)/(t1-t0)
                dlen = hlen*alpha
                stot += db2linear(sinr)*alpha
                psuccess *= 1.0 - self.calcper(dlen, 0, sinr, **calcperargs)
            per = 1.0 - psuccess
            sinr = linear2db(stot)
        else:
            # use SINR over the duration of the SIGNAL field
            sinr, hlen = self.sinr_part(p, DOT11N_SIGNAL, **kwargs), hlen
            # configure modulation and coding to calculate PER
            per = self.calcper(hlen, 0, sinr, **calcperargs)
        # set annotations and return PER
        p.setanno('dot11n-sinr', sinr)
        p.setanno('dot11n-per', per)
        return per

    def calcper_data(self, p, **kwargs):
        """Calculate probability of error for data decoding.

        :param p: Packet being decoded.
        :param kwargs: Additional keywords arguments passed to `sinr_heap()`
                       (or `sinr()`).
        :return: PER for decoding packet payload.

        This method sets the 'dot11n-sinr' and 'dot11n-per' annotations. The
        operation of this method depends on `DOT11N_USE_PIECEWISE_PER`.
        """
        # verify annotations
        for a in ['cif-rxts','cif-duration']:
            errmsg = "[DOT11N]: calcper_data() cannot find " + \
                     "'%s' annotation!"%(a)
            assert ANNO.supports(p, a), errmsg
        # get calcper arguments
        rmsdelay, bandwidth = 0, 1
        if p.hasanno('dot11n-rmsdelay'):
            bandwidth = self.radio.bandwidth
            rmsdelay  = p.getanno('dot11n-rmsdelay')
        calcperargs = {'rmsdelay':rmsdelay, 'bandwidth':bandwidth}
        # verify header parameters
        plen = len(p.payload)
        rate, length = p.rate, p.length
        assert (0<=rate<len(DOT11N_DATARATE) ), \
                "[DOT11N]: Invalid rate option (%s)!"%(rate)
        assert (p.length==plen), "[DOT11N]: Header length reported " + \
                "does not equal payload length; %s!=%s"%(p.length, plen)
        # calculate PER using appropriate method
        if DOT11N_USE_PIECEWISE_PER:
            sinrheap = self.sinr_heap(p, **kwargs)
            t1 = p.getanno('cif-rxts') + p.getanno('cif-duration')
            t0 = t1 - self.calcnofdm(plen, rate)*DOT11N_TSYM
            xheap = [(max(ta,t0), min(tb,t1), sinr) for (ta,tb,sinr) \
                        in sinrheap if (ta<t1) and (tb>t0)]
            errmsg = "[DOT11N]: Unable to find valid data from SINR heap!"
            assert (len(xheap)>0), errmsg
            # calculate piecewise PER and average SINR
            psuccess, stot = 1.0, 0.0
            for ta, tb, sinr in xheap:
                alpha = (tb-ta)/(t1-t0)
                dlen = plen*alpha
                stot += db2linear(sinr)*alpha
                psuccess *= 1.0 - self.calcper(dlen, rate, sinr, **calcperargs)
            per = 1.0 - psuccess
            sinr = linear2db(stot)
        else:
            # configure modulation and coding to calculate PER
            plen = length
            sinr = self.sinr_part(p, DOT11N_PAYLOAD, **kwargs)
            per = self.calcper(plen, rate, sinr, **calcperargs)
        # set annotations and return PER
        p.setanno('dot11n-sinr', sinr)
        p.setanno('dot11n-per', per)
        return per

    def calcper(self, p, rate, sinr, rmsdelay=0, \
                      rfo=DOT11N_USE_RFO, bandwidth=DOT11N_BANDWIDTH):
        """Calculate packet-error rate for given parameters.

        :param p: Packet or packet length in bytes.
        :param rate: Rate enumeration to indicate modulation and coding scheme.
        :param sinr: Signal-to-interference-and-noise ratio (in dB).
        :param rmsdelay: Root-mean square delay of channel (in seconds).
        :param rfo: Boolean flag; if true, use effective SINR by accounting for
                    the impact of residual frequency offset [default=`DOT11N_USE_RFO`].
        :param bandwidth: System bandwidth (in Hz).
        """
        # configure modulation and coding to calculate PER
        assert (0<=rate<len(DOT11N_DATARATE) ), \
                "[DOT11N]: Invalid rate enumeration!"
        plen = p
        if isinstance(p, Packet): plen = len(p)
        mod, coder = self.mod, self.coder
        mod.mtype = DOT11N_MTYPE[rate]
        coder.rate = DOT11N_CODERATE[rate]
        # XXX FIXME XXX: update SINR to include effect of Residual Frequency Offset
        if rfo:
            sinr_eff = self.calc_sinr_rfo(sinr)
        else:
            sinr_eff = sinr
        if DOT11N_USE_SOVA:
            sinr_eff += 2.0     # soft decision decoding has 2 dB gain
        if (rmsdelay>0):
            sinr_eff = self.calc_sinr_rmsdelay(sinr_eff, rmsdelay, bandwidth)
        # calculate PER
        uber = mod.ber(sinr_eff)
        per = coder.per(plen, uber)
        return per

    def calc_sinr_rfo(self, sinr, rfo=None):
        """Calculate effective SINR in presence of RFO.

        :param sinr: Base-SINR (in dB).
        :param rfo: Residual frequency offset normalized by sampling frequency
                    [default = None].

        If `rfo` is not specified this method will model RFO as a Gaussian R.V.
        with a variance of:

                1.0 / ((2*pi)^2 * N_FFT * SINR_B)

        where N_FFT is the FFT size and SINR_B is the base-SINR.
        """
        sinr_lin = db2linear(sinr)
        sinr_inv = db2linear(-sinr)
        if rfo is None:
            mu_rfo, var_rfo = 0, sinr_inv/(4 * (pi**2) * DOT11N_NFFT)
            phi = random.normal(mu_rfo, sqrt(var_rfo) )
        else:
            phi = rfo
        numer = sinr_lin*(sinc(phi)**2)
        denom = 1 + 0.5947*sinr_lin*(sin(pi*phi)**2)
        sinr_eff = linear2db(numer/denom)
        return sinr_eff

    def calc_sinr_rmsdelay(self, sinr, rmsdelay, bandwidth=1):
        """Calculate effective SINR for a multipath channel.

        :param sinr: Base-SINR (in dB).
        :param rmsdelay: Root-mean square delay spread (in seconds).
        :param bandwidth: System bandwidth (in Hz).
        """
        tau = rmsdelay*bandwidth
        if tau>0:
            sinr_eff = sinr + linear2db(1-exp(-1.0/tau))
        else:
            sinr_eff = sinr
        return sinr_eff

    def sinr(self, p, tstart=None, tstop=None, method=DOT11N_SINR_METHOD, **kwargs):
        """Calculate signal-to-interference-and-noise ratio (SINR).

        :param p: Packet to compute SINR for.
        :param tstart: Specify start time for looking at data from
                       `sinr_heap()` [default=-inf].
        :param tstop: Do not include information past `tstop` in SINR
                      computation [default=+inf].
        :param method: Method for computing SINR [default='min'].
        :param kwargs: Additional keyword arguments passed to `sinr_heap()`.

        This method uses the 'rxpower', 'noisepower', and 'cif-collision'
        annotations to calculate the SINR of the received packet p.

        :note: This method sets the 'dot11n-sinr' annotation indicating the SINR
               (in dB). The only SINR computation methods currently supported
               are 'min' and 'piecewise'.
        """
        # check for annotations
        for a in ['rxpower', 'noisepower']:
            errmsg = "[DOT11N]: sinr() cannot find '%s' annotation!"%(a)
            assert ANNO.supports(p, a), errmsg
        # set default parameters
        if tstart is None: tstart = -inf
        if tstop is None:  tstop  = +inf
        # calculate SINR from sinrheap
        sinrheap = self.sinr_heap(p,**kwargs)
        if (method=='min'):
            # find minimum SINR
            minsinr = +inf
            for ta,tb,sinr in sinrheap:
                # check tstart
                if (tb<tstart): continue
                if (ta<=tstart) and (tb>tstart): ta = tstart
                # check tstop
                if tb>tstop:  tb = tstop
                if sinr<minsinr: minsinr = sinr
                # quit if tstop is reached
                if (tb>=tstop): break
            rval = minsinr
        elif (method=='piecewise'):
            # compute weighted SINR
            wsinr, totaltime = 0, 0
            for ta,tb,sinr in sinrheap:
                # check tstart
                if (tb<tstart): continue
                if (ta<=tstart) and (tb>tstart): ta = tstart
                # check tstop
                if tb>tstop: tb = tstop
                totaltime += tb - ta
                wsinr += db2linear(sinr)*(tb - ta)
                # quit if tstop is reached
                if (tb>=tstop): break
            wsinr = linear2db(wsinr/totaltime)
            rval = wsinr
        else:
            errmsg = "[DOT11N]: Invalid method (%s) in sinr()!"%(method)
            raise RuntimeError, errmsg
        # set annotations
        p.setanno('dot11n-sinr', rval)
        return rval

    def sinr_part(self, p, mode=0, **kwargs):
        """Compute SINR over a particular partition of a `Dot11N` packet.

        :param p: `Dot11N` packet for which SINR will be computed.
        :param mode: Logical bit mask used to represent partition of packet.

        Use `mode` to select one or more contiguous segments of the `Dot11N`
        packet (e.g. DOT11N_STF, DOT11N_STF+DOT11N_LTF, DOT11N_SIGNAL, etc.)
        """
        # check annotations
        assert isinstance(p, Dot11N)
        for a in ['cif-rxts', 'cif-duration']:
            errmsg = "[DOT11N]: Cannot find '%s' annotation!"%(a)
            assert ANNO.supports(p, a), errmsg
        rxts = p.getanno('cif-rxts')
        duration = p.getanno('cif-duration')
        rate, plen = p.rate, p.length
        # determine start and end of partition
        t0, t1 = +inf, -inf
        if (mode & DOT11N_STF):
            t0 = min(t0, rxts)
            t1 = max(t1, t0 + DOT11N_TSHORT)
        if (mode & DOT11N_LTF):
            t0 = min(t0, rxts + DOT11N_TSHORT)
            t1 = max(t1, t0 + DOT11N_TLONG)
        if (mode & DOT11N_SIGNAL):
            t0 = min(t0, rxts + DOT11N_TSHORT + DOT11N_TLONG)
            t1 = max(t1, t0 + DOT11N_TSIGNAL)
        if (mode & DOT11N_PAYLOAD):
            t1 = max(t1, rxts + duration)
            t0 = min(t0, rxts + duration - self.calcnofdm(plen,rate)*DOT11N_TSYM)
        if (mode==0):
            t0, t1 = -inf, +inf
        return self.sinr(p, t0, t1, **kwargs)

    def parity(self, p):
        """Check parity of `Dot11N` packet header.

        :param p: `Dot11N` packet.
        :return: Boolean; true if parity check passes.

        :note: If `USE_DOT11NPARITY` is false, this method returns true.
        """
        assert isinstance(p, Dot11N), \
               "[DOT11N]: Cannot check parity of non-Dot11N packet!"
        return p.checkcrc()

    def get_valid_rates(self):
        """Get list of rate enumerations supported by PHY.

        This method uses the PHY configuration to determine the list of valid
        rate enumerations supported by the PHY (i.e. `ntx`).
        """
        ntx = self.radio.Ntx
        rates = range(len(DOT11N_DATARATE))[:8*ntx]
        return rates

    def get_datarate(self, r):
        """Get the data rate in bits-per-second (bps).

        :param r: Rate enumeration.

        This method gets data rate from `DOT11N_DATARATE`.
        """
        errmsg = "[DOT11NPHY]: Invalid rate index (%s)!"%(r)
        assert (0<=r<len(DOT11N_DATARATE) ), errmsg
        return DOT11N_DATARATE[r]

    def log_send(self, p, *args, **kwargs):
        """Convenience method for logging send event."""
        if self.verbose>DOT11N_VERBOSE:
            if isinstance(p, Dot11N):
                kwargs['length'] = p.length
                kwargs['phy-rate'] = p.rate
            if p.hasanno('cif-txpower'):
                kwargs['cif-txpower'] = "%.2f dBm"%(p.getanno('cif-txpower') )
            if p.hasanno('cif-duration'):
                kwargs['cif-duration'] = time2usec(p.getanno('cif-duration') )
            self.log("snd", p, *args, **kwargs)

    def log_recv(self, p, *args, **kwargs):
        """Convenience method for logging receive event."""
        if self.verbose>DOT11N_VERBOSE:
            if isinstance(p, Dot11N):
                kwargs['length'] = p.length
                kwargs['phy-rate'] = p.rate
            if p.hasanno('tx-cfo') and p.hasanno('rx-cfo'):
                txcfo, rxcfo = p.getanno('tx-cfo'), p.getanno('rx-cfo')
                cfo = txcfo - rxcfo
                kwargs['cfo'] = "%.2f kHz"%((txcfo-rxcfo)/1e3)
            if p.hasanno('rxpower'):
                kwargs['rxpower'] = "%.2f dBm"%(p.getanno('rxpower') )
            if p.hasanno('noisepower'):
                kwargs['noisepower'] = "%.2f dBm"%(p.getanno('noisepower') )
            if p.hasanno('cif-duration'):
                kwargs['cif-duration'] = time2usec(p.getanno('cif-duration') )
            kwargs.update(self.get_dot11n_anno(p))
            kwargs.update(self.get_model_anno(p))
            if p.hasanno('crcerror'):
                crcerror = p.getanno('crcerror')
                if crcerror: kwargs['crc'] = "FAIL"
                else:        kwargs['crc'] = "OK"
            self.log("rcv", p, *args, **kwargs)

    def log_drop(self, p, *args, **kwargs):
        """Convenience method for logging drop event."""
        if self.verbose>DOT11N_VERBOSE:
            if isinstance(p, Dot11N):
                kwargs['length'] = p.length
                kwargs['phy-rate'] = p.rate
            kwargs.update(self.get_dot11n_anno(p))
            kwargs.update(self.get_model_anno(p))
            errmsg = "[DOT11N]: Cannot log drop event without reason!"
            assert ('drop' in kwargs), errmsg
            self.log("drp", p, *args, **kwargs)
        self.cleananno(p)

    def log(self, event=None, p=None, *args, **kwargs):
        """Overloaded to check verbose level and set common annotations."""
        force = False
        if ('verbose' in kwargs): force = (kwargs['verbose']>DOT11N_VERBOSE)
        if self.verbose>DOT11N_VERBOSE or force:
            kwargs.update(self.get_dot11n_anno(p))
            kwargs.update(self.get_model_anno(p))
            CSPHY.log(self, event, p, *args, **kwargs)

    def get_dot11n_anno(self, p):
        """Internal method to extract annotations and convert to strings."""
        kwargs = {}
        if not isinstance(p, Packet): return kwargs
        if p.hasanno('dot11n-sinr'):
            kwargs['dot11n-sinr'] = "%.4f dB"%(p.getanno('dot11n-sinr') )
        if p.hasanno('dot11n-start'):
            kwargs['dot11n-start'] = "%s"%(p.getanno('dot11n-start') )
        if p.hasanno('dot11n-fdr'):
            kwargs['dot11n-fdr'] = "%.5g"%(p.getanno('dot11n-fdr') )
        if p.hasanno('dot11n-start-index'):
            kwargs['dot11n-start-index'] = "%d"%(p.getanno('dot11n-start-index') )
        if p.hasanno('dot11n-cfo'):
            kwargs['dot11n-cfo'] = "%.3g kHz"%(p.getanno('dot11n-cfo')/1e3 )
        if p.hasanno('dot11n-per'):
            kwargs['dot11n-per'] = "%.5g"%(p.getanno('dot11n-per') )
        if p.hasanno('dot11n-detect'):
            kwargs['dot11n-detect'] = "%s"%(p.getanno('dot11n-detect') )
        if p.hasanno('dot11n-header'):
            kwargs['dot11n-header'] = "%s"%(p.getanno('dot11n-header') )
        if p.hasanno('dot11n-error'):
            kwargs['dot11n-error'] = "%s"%(p.getanno('dot11n-error') )
        if p.hasanno('dot11n-channel-fading'):
            kwargs['dot11n-channel-fading'] = "%.4f dB"%(p.getanno('dot11n-channel-fading') )
        if p.hasanno('dot11n-rmsdelay'):
            kwargs['dot11n-rmsdelay'] = "%.5g"%(p.getanno('dot11n-rmsdelay') )
        if p.hasanno('dot11n-tgn-model'):
            kwargs['dot11n-tgn-model'] = "%s"%(p.getanno('dot11n-tgn-model') )
        if p.hasanno('cif-collision'):
            kwargs['cif-collision'] = strcollision(p)
            assert (kwargs['cif-collision'] is not None)
        if p.hasanno('net-root'):
            kwargs['net-root'] = p.getanno('net-root')
        return kwargs

    def get_model_anno(self, p):
        """Internal method to extract model-precition annotations."""
        kwargs = {}
        if not isinstance(p, Packet): return kwargs
        if p.hasanno('dot11n-model-detect'):
            kwargs['dot11n-model-detect'] = "%s"%(p.getanno('dot11n-model-detect') )
        if p.hasanno('dot11n-model-fdr'):
            kwargs['dot11n-model-fdr'] = "%.5g"%(p.getanno('dot11n-model-fdr') )
        if p.hasanno('dot11n-model-her'):
            kwargs['dot11n-model-her'] = "%.5g"%(p.getanno('dot11n-model-her') )
        if p.hasanno('dot11n-model-header'):
            kwargs['dot11n-model-header'] = "%s"%(p.getanno('dot11n-model-header') )
        if p.hasanno('dot11n-model-error'):
            kwargs['dot11n-model-error'] = "%s"%(p.getanno('dot11n-model-error') )
        if p.hasanno('dot11n-model-per'):
            kwargs['dot11n-model-per'] = "%.5g"%(p.getanno('dot11n-model-per') )
        return kwargs

class Dot11NRadio(Radio):
    """Radio interface for IEEE 802.11n physical layer.
    :CVariables:
     * `fc`: Carrier frequency.
     * `Lrx`: Analog system loss in receiver [default=`DOT11N_ANALOGLOSS`].
     * `Ntx`: Number of transmit antennas [default=2].
     * `Nrx`: Number of receive antennas  [default=2].
     * `fomax`: Maximum frequency offset of local oscillator (in ppm)
     * `bandwidth`: System bandwidth.
    """
    name = "IEEE 802.11n radio"
    fc = DOT11N_CARRIER
    Lrx = DOT11N_ANALOGLOSS
    bandwidth = DOT11N_BANDWIDTH
    txpower = DOT11N_MAXPOWER
    Nrx = 1
    Ntx = 1
    fomax = 10.0
    def __init__(self, **kwargs):
        """Constructor."""
        Radio.__init__(self, **kwargs)

    def noisepower(self):
        """Get effective noise power of radio interface.

        :return: Noise power (in dBm).

        The effective noise power is the `thermalnoise` of the system plus the
        noise figure contribution from digital noise.
        """
        noise = Radio.noisepower(self)
        return noise + DOT11N_NOISEFIGURE
