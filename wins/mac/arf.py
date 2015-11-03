#!  /usr/bin/env python

"""
Auto Rate Fallback (ARF) protocol; contains `ARF` class.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-09-25 00:33:20 -0500 (Sun, 25 Sep 2011) $
* $LastChangedRevision: 5147 $

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

from SimPy.Simulation import SimEvent, waitevent
from wins.fsm import FSM
from wins.packet import ANNO

from wins.ieee80211.dcf import DCF
from wins.ieee80211.dot11_support import *
from wins.mac.rateadaptation import RateAdapt
from wins.mac.arf_support import *

class ARF(DCF):
    """Auto Rate Fallback (ARF) protocol.

    Implementation based on description from Kamerman and Monteban (1997). By
    default ARF operates using CSMA/CA (i.e. not MACAW). To enable RTS/CTS
    messages, set `usecsma` to False.

    ARF Timeout
    ===========
    The ARF timer is used to indicate that the rate of transmissions should be
    increased prior to having `nsuccess` consecutive ACKs. This mechanism is
    meant to allow ARF to explore higher data rates when traffic is sparse.

    The timer can be implemented as a time-based or packet-based timer. The
    time-based timer is a straight forward timeout timer that waits for a
    duration before firing and indicating that the rate of the transmission
    should be increased. The packet-based timer uses the "virtual timer"
    approach described by Qiao, Choi, and Shin (2002). This "timer" counts the
    number of transmission attempts made rather than using a time-based timer.

    This module implements the ARF timer as a time-based timer.

    :CVariables:
     * `base`: Property to access base rate enumeration from `ra`.
     * `rates`: Property to access ordered list of valid rate enumerations. This
                is ordered based on the data rates reported by `get_datarate()`.

    :IVariables:
     * `ra`: `RateAdapt` module.
     * `ackcount`: Dictionary containing number of consecutive successful ACKs
                   received (negative values indicate consecutive dropped ACKs).
                   Values are keyed by destination address.
     * `arftimer`: Dictionary containing ARF timers used for updating transmit
                   data rate. Values are keyed by destination address. This
                   module implements ARF timers as packet-based timers.
     * `probation`: Dictionary containing probation flags indicating if
                    probation period has been entered; values are keyed by
                    destination address.
     * `nsuccess`: Threshold for number of consecutive ACKs that must be
                   received prior to increasing the data rate.
     * `nfailure`: Threshold for number of consecutive ACK failures that
                   will trigger a decrease in the data rate.
     * `timeout`: Timeout value for ARF timer; if exceeded, then increase rate.

    """
    name = "Auto Rate Fallback"
    tracename = "ARF"
    def __init__(self, usecsma=True, **kwargs):
        """Constructor."""
        # create ARF parameters
        self.nsuccess = None
        self.nfailure = None
        self.timeout  = None
        self.ackcount = {}
        self.arftimer = {}
        self.probation = {}
        self._rates = None
        # create ARF events
        self.dropack = SimEvent()
        self.recvack = SimEvent()
        DCF.__init__(self, usecsma=usecsma, **kwargs)
        # update event names
        self.dropack.name = "%s%s"%(self.name, ".dropack")
        self.recvack.name = "%s%s"%(self.name, ".recvack")

    base = property(fget=lambda self: self.ra.base)
    rates = property(fget=lambda self: self.get_rates() )

    def configure(self, base=None, nsuccess=None, nfailure=None, timeout=None, **kwargs):
        """Configure rate adaptation parameters.

        :param base: Base rate used to initialize `RateAdapt` component [default=0].
        :param nsuccess: Threshold for number of consecutive ACKs that must be
                         received prior to increasing the data rate.
        :param nfailure: Threshold for number of consecutive ACK failures that
                         will trigger a decrease in the data rate.
        :param timeout: Timeout that can automatically trigger
        """
        DCF.configure(self, **kwargs)
        if base is None: base = 0
        if nsuccess is None: nsuccess = ARF_NSUCCESS
        if nfailure is None: nfailure = ARF_NFAILURE
        if timeout  is None: timeout  = ARF_TIMEOUT
        ra = self.newchild("ra", RateAdapt, base=base, tracename=self.tracename+".RA")
        ra.set_rate(self.broadcast) # use default base rate for broadcast
        # create FSM to manage ARF
        fsm = self.newchild("manager", FSM, tracename=self.tracename+".MGR")
        fsm.goto(self.MGR)
        # set other parameters
        self.nsuccess = nsuccess
        self.nfailure = nfailure
        self.timeout  = timeout
        self.ackcount = {}
        self.arftimer = {}
        self.probation = {}

    def MGR(self, fsm):
        """MGR state; manage ARF protocol."""
        yield waitevent, fsm, (self.dropack, self.recvack)
        dropack = (self.dropack in fsm.eventsFired)
        recvack = (self.recvack in fsm.eventsFired)
        # get destination parameter from event
        if dropack: dst = self.dropack.signalparam
        if recvack: dst = self.recvack.signalparam
        errmsg = "[ARF]: Invalid state transition!"
        assert (dropack ^ recvack), errmsg
        # initialize ackcount if needed
        if (dst not in self.ackcount):
            self.ackcount[dst] = 0
            self.probation[dst] = False
            self.arftimer[dst] = 0
        nack = self.ackcount[dst]
        probation = self.probation[dst]
        timer = 0
        # check probation condition
        errmsg = "[ARF]: ACK count must be zero in probation!"
        assert not ((nack<>0) and probation), errmsg
        # increment arftimer for every transmission attempt (packet-based timer)
        timer += 1
        # process drop ACK or recv ACK
        if (nack>0):
            if dropack:     nack = -1
            elif recvack:   nack += 1
        elif (nack<0):
            if dropack:     nack -= 1
            elif recvack:   nack = +1
        else:
            if dropack:     nack = -1
            elif recvack:   nack = +1
            # probation? dropack -> decrease rate and continue using lower rate
            if probation and dropack:
                nack, timer, rate = 0, 0, self.decrate(dst)
                self.log("DECRATE", rate=rate, dst=dst, ackcount=nack, timer=timer, probation=probation)
            probation =  False      # reset flag
        # log events
        if dropack: self.log("DROPACK", dst=dst, ackcount=nack, arftimer=timer)
        if recvack: self.log("RECVACK", dst=dst, ackcount=nack, arftimer=timer)
        # check if Nsuccess or Nfailure condition is met
        timerexpired = (timer>self.timeout)
        if (self.nsuccess-nack<1) or timerexpired:
            # (Nack >= +Nsuccess) OR (Timer expired)
            #   -> increase rate, enter probation period
            rate = self.incrate(dst)
            nack, timer = 0, 0      # reinitialize ackcount and timer
            probation = True        # set probation flag
            self.log("INCRATE", rate=rate, dst=dst, ackcount=nack, probation=True)
        elif (nack+self.nfailure<1):
            # Nack <= -Nfailure -> decrease rate
            rate = self.decrate(dst)
            nack, timer = 0, 0     # set ackcount to 0, reset timer
            self.log("DECRATE", rate=rate, dst=dst, ackcount=nack, probation=probation)
        # set ackcount/probation value
        self.ackcount[dst]  = nack
        self.arftimer[dst]  = timer
        self.probation[dst] = probation
        # continue in MGR
        yield fsm.goto(self.MGR)

    def retry(self, *args, **kwargs):
        """Overloaded to manage dropped ACKs."""
        rc = self.retrycount
        DCF.retry(self, *args, **kwargs)
        # Actual retry?
        if (self.retrycount>rc) and self.isdot11data(self.datatosend):
            # retry count was increased -> dropped ACK
            data = self.get_dot11data(self.datatosend)
            dst, src  = data.addr1, data.addr2
            #self.log("DROPACK", dst=dst, src=src)
            self.dropack.signal(dst)

    def send_data(self, data):
        """Additional processing for outgoing DATA."""
        DCF.send_data(self, data)
        dst  = data.addr1
        rate = self.ra.get_rate(dst)
        data.setanno('phy-rate', rate)

    def recv_ack(self, ack):
        """Additional processing for incoming ACK."""
        DCF.recv_ack(self, ack)
        # Expecting ACK?
        if self.isdot11data(self.datatosend):
            data = self.get_dot11data(self.datatosend)
            dst, src = data.addr1, data.addr2
            # ACK for me? -> received ACK
            if (src==ack.addr1):
                #self.log("RECVACK", ack, dst=dst, src=src)
                self.recvack.signal(dst)

    def incrate(self, dst):
        """Increase data rate corresponding to destination `dst`."""
        rate = self.ra.get_rate(dst)
        errmsg = "[ARF]: Invalid rate (%s) in incrate()!"%(rate)
        assert (rate in self.rates), errmsg
        idx = self.rates.index(rate)
        if (idx+1<len(self.rates)):
            rate = self.rates[idx+1]
        self.ra.set_rate(dst, rate)
        return rate

    def decrate(self, dst):
        """Decrease data rate corresponding to destination `dst`."""
        rate = self.ra.get_rate(dst)
        errmsg = "[ARF]: Invalid rate (%s) in incrate()!"%(rate)
        assert (rate in self.rates), errmsg
        idx = self.rates.index(rate)
        if (idx>0):
            rate = self.rates[idx-1]
        self.ra.set_rate(dst, rate)
        return rate

    def get_rates(self, force=False):
        """Get list of valid rate enumerations.

        :param force: If true, recalculate rates; otherwise use cached value
                      [default=False].
        :return: List of valid rate enumerations.
        """
        if force: self._rates = None
        # initialize ordered rates
        if self._rates is None:
            orates = []     # ordered rates
            rates = [r for r in self.phy.rate]
            while rates:
                maxrate, maxidx = None, None
                for k in range(len(rates)):
                    r = rates[k]
                    rbps = self.get_datarate(r)
                    if maxrate is None: maxrate, maxidx = r, k
                    mbps = self.get_datarate(maxrate)
                    if (rbps>mbps):  maxrate, maxidx = r, k
                orates.insert(0, rates.pop(maxidx))     # add max rate to orates
            self._rates = orates
        return self._rates

    def log_send(self, p, *args, **kwargs):
        """Updated to print `ARF` related parameters."""
        if p.hasanno('phy-rate'):
            kwargs['phy-rate'] = p.getanno('phy-rate')
        DCF.log_send(self, p, *args, **kwargs)

    def log_recv(self, p, *args, **kwargs):
        """Updated to print `ARF` related parameters."""
        if p.hasanno('phy-rate'):
            kwargs['phy-rate'] = p.getanno('phy-rate')
        DCF.log_recv(self, p, *args, **kwargs)
