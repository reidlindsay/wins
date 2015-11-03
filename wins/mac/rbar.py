#!  /usr/bin/env python

"""
Receiver-Based AutoRate (RBAR) protocol; contains `RBAR` class.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-10-28 03:08:28 -0500 (Fri, 28 Oct 2011) $
* $LastChangedRevision: 5316 $

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

from SimPy.Simulation import hold, now
from wins.packet import Packet, ANNO
from wins.helper import time2usec
from wins.crc import crcupdate, crcremove

from wins.ieee80211.dcf import DCF
from wins.ieee80211.dot11_support import *
from wins.mac.rateadaptation import RateAdapt, RAI
from wins.mac.rbar_support import *

class RBAR(DCF):
    """Receiver-Based AutoRate (RBAR) protocol.

    Implementation based on description from Holland, Vaidya, and Bahl (2001).

    :CVariables:
     * `base`: Property to access base rate enumeration from `ra`.

    :IVariables:
     * `ra`: `RateAdapt` module.
     * `thresh`: Dictionary containing minimum SNR threshold corresponding to
                 rate specific index [default=`RBAR_THRESHOLD`].
     * `_rshduration`: Internal member used to cache duration of RSH.
    """
    name = "Receiver-Based AutoRate"
    tracename = "RBAR"
    def __init__(self, usecsma=None, **kwargs):
        """Constructor."""
        self.thresh = None
        self._rshduration = None
        DCF.__init__(self, usecsma=False, **kwargs)

    base = property(fget=lambda self: self.ra.base)
    rshduration = property(fget=lambda self: self.get_rshduration() )

    def configure(self, base=None, thresh=None, **kwargs):
        """Configure rate adaptation parameters.

        :param base: Base rate used to initialize `RateAdapt` component [default=0].
        :param thresh: Dictionary containing minimum SNR threshold corresponding
                       to rate specific index [default=`RBAR_THRESHOLD`].
        """
        DCF.configure(self, **kwargs)
        if base is None: base = 0
        if thresh is None: thresh = RBAR_THRESHOLD
        ra = self.newchild("ra", RateAdapt, base=base, tracename=self.tracename+".RA")
        ra.set_rate(self.broadcast) # use default base rate for broadcast
        self.thresh = thresh

    def RXDATA(self, fsm, pkt):
        """RXDATA state; overloaded to support RSH."""
        crcerror = self.haserror(pkt)
        if crcerror:
            self.log_drop(pkt, drop="CRC error")
            yield fsm.goto(self.RESUME)
        assert self.isdot11data(pkt), "[RBAR]: Cannot find DATA in RXDATA!"
        rsh = self.get_dot11rsh(pkt)
        if rsh:
            # process RSH
            yield fsm.goto(self.RXRSH, rsh)
        else:
            # otherwise use original RXDATA() state
            OLDSTATE = lambda *args, **kwargs: DCF.RXDATA(self, *args, **kwargs)
            yield fsm.goto(OLDSTATE, pkt)

    def RXRSH(self, fsm, pkt):
        """RXRSH state; update NAV using RSH."""
        crcerror = self.haserror(pkt)
        if crcerror:
            self.log_drop(pkt, drop="CRC error")
            yield fsm.goto(self.RESUME)
        assert self.isdot11rsh(pkt), "[RBAR]: Cannot find RSH in RXDATA!"
        rsh = self.get_dot11rsh(pkt)
        if rsh:
            # update NAV
            nav = rsh.ID
            yield self.navupdate(fsm, nav*1e-6)
        # go to RESUME
        yield fsm.goto(self.RESUME)

    def TXUCAST(self, fsm):
        """TXUCAST state; overloaded to go to `TXRSH`."""
        yield fsm.goto(self.TXRSH)

    def TXRSH(self, fsm):
        """TXRSH state; overloaded to send RSH prior to DATA."""
        assert DCF.isdot11data(self, self.datatosend), \
                "[RBAR]: Cannot determine 'datatosend' in TXUCAST!"
        assert not (self.datatosend.addr1==self.broadcast), \
                "[RBAR]: Cannot send broadcast 'datatosend' in TXUCAST!"
        # update rate annotation
        data = self.get_dot11data(self.datatosend)
        dst, src = data.addr1, data.addr2
        rate = self.ra.get_rate(dst)
        data.setanno('phy-rate', rate)
        if not self.usecsma:
            # send RSH before sending DATA
            rsh = self.dot11rsh(addr1=dst, addr2=src)
            self.send_rsh(rsh, self.datatosend)
            # recalculate NAV from new rate information
            nav = self.rshnav(self.datatosend)
            rsh.ID = nav
            pkt = crcupdate(rsh)
            # send and hold for duration
            rate, length = pkt[RAI].rate, pkt[RAI].length
            self.log("RSH", pkt, nav=nav, rate=rate, length=length)
            duration = self.duration(pkt)
            self.log_send(pkt, nav=nav)
            yield self.TXD.send(fsm, [pkt])
            yield hold, fsm, duration
            # pause for SIFS before continuing
            yield hold, fsm, self.sifs
        # go to DCF.TXUCAST state
        OLDSTATE = lambda *args, **kwargs: DCF.TXUCAST(self, *args, **kwargs)
        yield fsm.goto(OLDSTATE)

    def send_rts(self, rts, data):
        """Additional processing for outgoing RTS."""
        DCF.send_rts(self, rts, data)
        # get rate and length information
        plen = len(data)
        dst, src = data.addr1, data.addr2
        rate = self.ra.get_rate(dst)
        data.setanno('phy-rate', rate)
        # copy info to RTS/RAI
        rai = rts[RAI]
        rai.rate = rate
        rai.length = plen

    def recv_cts(self, cts):
        """Additional processing for incoming CTS.

        :param cts: Received CTS packet.

        If the CTS is "intended for me", this method will log the reported rate
        adaptation information with the `ra` module.
        """
        errmsg = "[RBAR]: Non-CTS in recv_cts()!"
        assert self.isdot11cts(cts), errmsg
        errmsg = "[RBAR]: Invalid CTS! No RAI found in recv_cts()!"
        assert cts.haslayer(RAI), errmsg
        # check if CTS is for me
        addr, ctsdst = self.address, cts.addr1
        if (addr==ctsdst):
            errmsg = "[RBAR]: Received unsolicited CTS " + \
                     "addressed to me (%s)!"%(addr)
            assert self.isdot11data(self.datatosend), errmsg
            ctssrc = self.datatosend.addr1     # source of CTS
            # log RAI information
            rate = cts[RAI].rate
            self.ra.set_rate(ctssrc, rate)
            self.log("SETRATE", cts, addr=addr, src=ctssrc, dst=ctsdst, rate=rate)

    def send_cts(self, cts, rts):
        """Additional processing for outgoing CTS.

        :param cts: Outgoing CTS packet.
        :param rts: Received RTS packet.

        This method uses the 'phy-sinr' annotation to apply the rate adaptation
        algorithm used by `RBAR`.
        """
        DCF.send_cts(self, cts, rts)
        errmsg = "[RBAR]: Invalid RTS! No RAI found in send_cts()!"
        assert rts.haslayer(RAI), errmsg
        errmsg = "[RBAR]: Invalid CTS! No RAI found in send_cts()!"
        assert cts.haslayer(RAI), errmsg
        errmsg = "[RBAR]: Cannot find 'phy-sinr' annotation in send_cts()!"
        assert rts.hasanno('phy-sinr'), errmsg
        # extract info from RTS+RAI
        sinr = rts.getanno('phy-sinr')
        newrate = self.calcrate(sinr)
        cts[RAI].rate = newrate
        cts[RAI].length = rts[RAI].length

    def send_data(self, data):
        """Additional processing for outgoing DATA."""
        DCF.send_data(self, data)
        dst  = data.addr1
        rate = self.ra.get_rate(dst)
        # allow fixed rate for broadcast packets
        fixedrate = data.hasanno('phy-fixed-rate')
        if (data.addr1==self.broadcast) and fixedrate:
            rate = data.getanno('phy-fixed-rate')
        data.setanno('phy-rate', rate)

    def send_rsh(self, rsh, data):
        """Additional processing for outgoing RSH.

        :param rsh: Outgoing RSH packet.
        :param data: Corresponding DATA packet.
        """
        errmsg = "[RBAR]: Non-RSH found in send_rsh()!"
        assert self.isdot11rsh(rsh), errmsg
        errmsg = "[RBAR]: Non-DATA found in send_rsh()!"
        assert self.isdot11data(data), errmsg
        errmsg = "[RBAR]: Invalid RSH! No RAI found in send_rsh()!"
        assert rsh.haslayer(RAI), errmsg
        # get rate and length parameters
        dst, src = rsh.addr1, rsh.addr2
        assert (dst==data.addr1)
        assert (src==data.addr2)
        rate = self.ra.get_rate(dst)
        length = len(data)
        # set parameters in RSH+RAI
        rsh[RAI].rate = rate
        rsh[RAI].length = length

    def rtsnav(self, p, *args, **kwargs):
        """Overloaded to add RSH to NAV."""
        nav = DCF.rtsnav(self, p, *args, **kwargs)
        # add RSH + SIFS
        d = self.sifs + self.rshduration
        nav += int(d*1e6)
        return nav

    def ctsnav(self, rts):
        """Overloaded to add RSH to NAV and apply rate adaptation..

        :param rts: RTS packet.
        :return: Integer; representing NAV value.

        This method uses `DCF.ctsnav()`, `calcrate()`, and `duration()` to
        calculate the new value of the CTS NAV.
        """
        ctsnav = DCF.ctsnav(self, rts)
        errmsg = "[RBAR]: Invalid RTS! No RAI found in ctsnav()!"
        assert rts.haslayer(RAI), errmsg
        errmsg = "[RBAR]: Cannot find 'phy-sinr' annotation in ctsnav()!"
        assert rts.hasanno('phy-sinr'), errmsg
        # calculate new rate from SINR
        sinr = rts.getanno('phy-sinr')
        newrate = self.calcrate(sinr)
        # extract info from RTS+RAI
        rai = rts[RAI]
        oldrate, length = rai.rate, rai.length
        # recompute NAV using length and new rate
        dold = self.duration(length, oldrate)
        dnew = self.duration(length, newrate)
        nav = ctsnav - int(dold*1e6) + int(dnew*1e6)
        return nav

    def rshnav(self, data, *args, **kwargs):
        """Calculate NAV value for RSH.

        :param data: DATA packet.
        :param args: Arguments passed on to `DCF.rtsnav()`.
        :param kwargs: Keywords passed on to `DCF.rtsnav()`.
        :return: Integer; representing NAV value.

        This method computes the RSH NAV as follows:

            NAV = SIFS + DATA + SIFS + ACK

        :note: This method uses `DCF.rtsnav()` to calculate the new NAV value.
        """
        # Use regular RTS NAV (i.e. from DCF) to compute new NAV value
        rtsnav = DCF.rtsnav(self, data, *args, **kwargs)
        d = self.sifs + self.ctsduration
        nav = rtsnav - int(d*1e6)
        return nav

    def calcrate(self, sinr):
        """Calculate rate enumeration for given SINR using `thresh`."""
        errmsg = "[RBAR]: Cannot find valid dictionary for `thresh`!"
        assert isinstance(self.thresh, dict), errmsg
        # iterate over thresh to find maxrate
        maxrate = None
        for (r,t) in self.thresh.items():
            rbps = self.get_datarate(r)
            if sinr>t:
                if (maxrate is not None):
                    mbps = self.get_datarate(maxrate)
                    if (rbps>mbps): maxrate = r
                else:
                    maxrate = r
        # return maxrate
        if maxrate is None: maxrate = self.base
        return maxrate

    def get_rshduration(self, force=False):
        """Calculate duration of RSH message.

        :param force: If true, ignore any cached value.
        """
        if force or (self._rshduration is None):
            rsh = self.dot11rsh()
            pkt = crcupdate(rsh)
            self._rshduration = self.duration(pkt)
        return self._rshduration

    def isdot11data(self, p):
        """Overloaded to check for DATA or RSH."""
        isdata = DCF.isdot11data(self, p)
        isrsh = self.isdot11rsh(p)
        return (isdata or isrsh)

    def isdot11rsh(self, p):
        """Is packet a Reservation Subheader (RSH)?"""
        return isdot11rsh(p)

    def get_dot11data(self, p):
        """Overloaded to extract DATA or RSH."""
        isrsh = self.isdot11rsh(p)
        if isrsh:
            return self.get_dot11rsh(p)
        else:
            return DCF.get_dot11data(self, p)

    def get_dot11rsh(self, p):
        """Extract RSH from `p`."""
        return get_dot11rsh(p)

    def dot11data(self, *args, **kwargs):
        """Overloaded to add 'phy-rate' annotation to new packet."""
        p = DCF.dot11data(self, *args, **kwargs)
        dst = p.addr1
        rate = self.ra.get_rate(dst)
        p.setanno('phy-rate', rate)
        return p

    def dot11rts(self, *args, **kwargs):
        """Overloaded to add 'phy-rate' annotation to new packet."""
        rargs = {}
        for s in ['rate', 'length']:
            if s in kwargs:
                rargs[s] = kwargs[s]
                del kwargs[s]
        # create RTS + RAI payload
        p = DCF.dot11rts(self, *args, **kwargs)/RAI(**rargs)
        p.setanno('phy-rate', self.base)
        return p

    def dot11cts(self, *args, **kwargs):
        """Overloaded to add 'phy-rate' annotation to new packet."""
        rargs = {}
        for s in ['rate', 'length']:
            if s in kwargs:
                rargs[s] = kwargs[s]
                del kwargs[s]
        # create CTS + RAI payload
        p = DCF.dot11cts(self, *args, **kwargs)/RAI(**rargs)
        p.setanno('phy-rate', self.base)
        return p

    def dot11ack(self, *args, **kwargs):
        """Overloaded to add 'phy-rate' annotation to new packet."""
        p = DCF.dot11ack(self, *args, **kwargs)
        p.setanno('phy-rate', self.base)
        return p

    def dot11rsh(self, *args, **kwargs):
        """Create new `Dot11RSH` packet."""
        rargs = {}
        for s in ['rate', 'length']:
            if s in kwargs:
                rargs[s] = kwargs[s]
                del kwargs[s]
        # create RSH + RAI payload
        p = Dot11RSH(*args, **kwargs)/RAI(**rargs)
        p.setanno('phy-rate', self.base)
        return p

    def log_send(self, p, *args, **kwargs):
        """Updated to print `RAI` parameters."""
        if ANNO.supported(p):
            if p.haslayer(RAI):
                kwargs['rai-rate'] = p[RAI].rate
                kwargs['rai-length'] = p[RAI].length
            if p.hasanno('phy-rate'):
                kwargs['phy-rate'] = p.getanno('phy-rate')
        DCF.log_send(self, p, *args, **kwargs)

    def log_recv(self, p, *args, **kwargs):
        """Updated to print `RAI` parameters."""
        if ANNO.supported(p):
            if p.haslayer(RAI):
                kwargs['rai-rate'] = p[RAI].rate
                kwargs['rai-length'] = p[RAI].length
            if p.hasanno('phy-rate'):
                kwargs['phy-rate'] = p.getanno('phy-rate')
        DCF.log_recv(self, p, *args, **kwargs)
