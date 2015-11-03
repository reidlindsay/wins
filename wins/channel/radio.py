#!  /usr/bin/env python

"""
Interface for connecting nodes to a wireless channel; contains `Radio` class.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-10-03 17:37:57 -0500 (Mon, 03 Oct 2011) $
* $LastChangedRevision: 5180 $

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

from wins.channel.interface import ChannelInterface
from wins.helper import linear2db, db2linear
from wins.packet import ANNO

from numpy.random import uniform

class Radio(ChannelInterface):
    """Interface to connect to a wireless channel.

    Packet Annotations and Radio
    ============================
    This class sets the following annotations:

    ============== ==============================================================
    Name            Description
    ============== ==============================================================
    cif-txpower     Power of transmitted packets (in dBm).
    -------------- --------------------------------------------------------------
    rxpower         Power of received packet (in dBm).
                  
                    This value is computed from the 'cif-txpower' and 'pathloss'
                    annotations as well as `Lrx` and `Grx` parameters. 
    -------------- --------------------------------------------------------------
    noisepower      Power of additive thermal noise at receiver (in dBm).
    -------------- --------------------------------------------------------------
    tx-cfo          Local frequency offset at transmitter (in Hz).
    -------------- --------------------------------------------------------------
    rx-cfo          Local frequency offset at receiver (in Hz).
    ============== ==============================================================

    Overload `set_sendanno()` (or `set_recvanno()`) to change how annotations are
    applied to outgoing (or incoming) packets.

    :CVariables:
     * `fc`: Default for center frequency (in Hz).
     * `Gtx`: Default for transmit antenna gain (in dB).
     * `Grx`: Default for receive antenna gain (in dB).
     * `Ltx`: Default for system loss associated with transmitter (in dB).
     * `Lrx`: Default for system loss associated with receiver (in dB).
     * `bandwidth`: Default for system bandwidth (in Hz).
     * `txpower`: Default for transmit power (in dBm).
     * `Ntx`: Default for number of transmit antennas.
     * `Nrx`: Default for number of receive antennas.
     * `fomax`: Default for max frequency offset of local oscillator (in ppm).
     * `cfo`: Property to access frequency offset of local oscillator (in Hz).

    :IVariables:
     * `fc`: Center frequency (in Hz).
     * `Gtx`: Transmit antenna gain (in dB).
     * `Grx`: Receive antenna gain (in dB).
     * `Ltx`: System loss associated with transmitter (in dB).
     * `Lrx`: System loss associated with receiver (in dB).
     * `bandwidth`: System bandwidth (in Hz).
     * `txpower`: Transmit power (in dBm).
     * `Ntx`: Number of transmit antennas.
     * `Nrx`: Number of receive antennas.
     * `fomax`: Maximum frequency offset of local oscillator (in ppm)

    Local frequency offset is set every time `set_sendanno()` or
    `set_recvanno()` is called.

    :note: The base `Radio` class does not consider directionallity. All
           parameters are considered omnidirectional.
    """
    fc = 2.4e9
    Gtx = 0
    Grx = 0
    Ltx = 0
    Lrx = 0
    bandwidth = 20e6
    txpower = 20
    Ntx = 1
    Nrx = 1
    fomax = 0.0      # in ppm
    name = "radio"
    tracename = "RF"
    def __init__(self, **kwargs):
        """Constructor."""
        self.fc = None
        self.Gtx, self.Grx = None, None
        self.Ltx, self.Lrx = None, None
        self.Ntx, self.Nrx = None, None
        self.fomax = None
        self.__cfo = None
        self.bandwidth = None
        ChannelInterface.__init__(self, **kwargs)

    cfo = property(fget=lambda self: self.__cfo)

    def set_cfo(self, fo=None):
        """Set new value for frequency offset of local oscillator.

        :param fo: New value of frequency offset (in Hz) [default=None].

        If `fo` is None, this method will determine a new `cfo` based on a
        uniform distribution based on `fomax`.
        """
        if fo is None:
            fo = uniform(-1.0,1.0)*self.fomax*1e-6*self.fc
        self.__cfo = fo

    def configure(self, fc=None, Gtx=None, Grx=None, \
                  Ltx=None, Lrx=None, Ntx = None, Nrx = None, fomax=None, \
                  txpower = None, bandwidth=None, **kwargs):
        """Set up radio parameters.

        :param fc: Center frequency of `Radio` (in Hz).
        :param Gtx: Transmit antenna gain (in dB).
        :param Grx: Receive antenna gain (in dB).
        :param Ltx: System loss associated with transmitter (in dB).
        :param Lrx: System loss associated with receiver (in dB).
        :param Ntx: Number of transmit antennas.
        :param Nrx: Number of receive antennas.
        :param fomax: Maximum frequency offset of local oscillator (in ppm).
        :param txpower: Transmit power (in dBm).
        :param bandwidth: System bandwidth (in Hz).
        :param kwargs: Keywords passed to `ChannelInterface.configure()`.

        All parameters not specified will be set to their default values as
        specified by the class.
        """
        # get parameters
        cls = self.__class__
        if fc is None: fc = cls.fc
        if Gtx is None: Gtx = cls.Gtx
        if Grx is None: Grx = cls.Grx
        if Ltx is None: Ltx = cls.Ltx
        if Lrx is None: Lrx = cls.Lrx
        if bandwidth is None: bw = cls.bandwidth
        if txpower is None: txpower = cls.txpower
        if Ntx is None: Ntx = cls.Ntx
        if Nrx is None: Nrx = cls.Nrx
        if fomax is None: fomax = cls.fomax
        # set parameters
        self.fc = fc
        self.Gtx, self.Grx = Gtx, Grx
        self.Ltx, self.Lrx = Ltx, Lrx
        self.bandwidth = bw
        self.txpower = txpower
        self.Ntx, self.Nrx = Ntx, Nrx
        self.fomax = fomax
        self.set_cfo()
        ChannelInterface.configure(self, **kwargs)

    def set_sendanno(self, p):
        """Set radio annotations for transmitted packet.

        :return: Modified packet.

        If the 'cif-txpower' annotation has not already been set by another
        protocol, this method will use the local `txpower` value. The transmit
        power is then adjusted by the transmit antenna gain `Grx` and the
        transmitter system loss `Lrx`, such that the new transmit power is:

            Ptx (dBm) = Pold (dBm) + Gtx (dB) - Ltx (dB)

        The 'cif-txpower' annotation is then set to this new value. Finally, the
        method calls `ChannelInterface.set_sendanno()`. Overload this method as
        necessary.
        """
        Gtx = self.Gtx
        Ltx = self.Ltx
        Ptx = self.txpower
        if p.hasanno('cif-txpower'): Ptx = p.getanno('cif-txpower')
        # adjust transmit power
        txpower = Ptx + Gtx - Ltx
        p.setanno('cif-txpower', txpower)
        # update and set local cfo
        self.set_cfo()
        p.setanno('tx-cfo', self.cfo)
        return ChannelInterface.set_sendanno(self, p)

    def set_recvanno(self, p):
        """Set radio annotations for received packet.

        :return: Modified packet.

        This method sets the 'rxpower' annotation to the value returned by
        `rxpower()`, sets the 'noisepower' annotation using `thermalnoise()` and
        `bandwidth`, and then calls ChannelInterface.set_recvanno()`. Overload
        this method as needed.

        :note: This method is called immediately after the packet `p` has been
               inserted in the `rxbuffer`.
        """
        r = ChannelInterface.set_recvanno(self, p)
        rxpower = self.rxpower(r, self.Grx, self.Lrx)
        noisepower = self.noisepower()
        r.setanno('rxpower', rxpower)
        r.setanno('noisepower', noisepower)
        # update and set local cfo
        self.set_cfo()
        r.setanno('rx-cfo', self.cfo)
        return r

    def rxenergy(self):
        """Calculate total receive energy of packets in `rxbuffer`.

        :return: Receive power (in dBm), or -Inf if nothing is in `rxbuffer`.

        This method uses the 'rxpower' annotation of packets in `rxbuffer` to
        determine the total receive power. Any packets that do not support this
        annotation will not contribute to the total power calculated.
        """
        tot = 0     # sum total in milliwatts
        for p in self.rxbuffer:
            if ANNO.supports(p, 'rxpower'):
                prx = p.getanno('rxpower')  # rxpower in dBm
                tot += db2linear(prx)
        tot_dbm = linear2db(tot)
        return tot_dbm

    def noisepower(self):
        """Get noise power for radio interface.

        :return: Noise power (in dBm).

        By default, this returns the `thremalnoise` of the system.
        """
        return self.thermalnoise(self.bandwidth)

    def snr(self, p, **kwargs):
        """Get signal-to-noise-ratio (SNR) for packet.

        :param kwargs: Additional keywords passed to `rxpower`.
        :return: Preprocessing SNR (in dB).

        This method only uses the `rxpower()` and `noisepower()` functions to
        compute the SNR of the packet in question. This means the calculated SNR
        is a pre-processing value.
        """
        rxpower = self.rxpower(p, self.Grx, self.Lrx, **kwargs)
        noisepower = self.noisepower()
        snr = rxpower - noisepower
        return snr

    @classmethod
    def rxpower(cls, p, Grx=None, Lrx=None, txpower=None, pathloss=None):
        """Calculate receive power of a received packet.

        :param Grx: Receive antenna gain (in dB).
        :param Lrx: System loss associated with receiver (in dB).
        :param txpower: Transmit power (in dBm) [default=None].
        :param pathloss: Pathloss (in dB) [defulat=None].
        :return: Receive power (in dBm).

        If `txpower` and `pathloss` are not provided, this method utilizes the
        'cif-txpower' and 'pathloss' annotations along with the `Lrx` and `Grx`
        parameters to compute the receive power. If `Lrx` and `Grx` are not
        provided, these parameters will use the default values specified by the
        class. Overload this method as needed.

        Received power is a function of transmit power, pathloss, system loss at
        the receiver, and receiver antenna gain:

            Prx (dBm) = Ptx (dBm) - PL (dB) + Grx (dB) - Lrx (dB)

        :note: This method throws an AssertionError if either 'cif-txpower' or
               'pathloss' annotations are not available.
        """
        if (txpower is None):
            errmsg = "[RADIO]: Cannot calculate Prx without 'cif-txpower' annotation!"
            assert p.hasanno('cif-txpower'), errmsg
            txpower = p.getanno('cif-txpower')
        if (pathloss is None):
            errmsg = "[RADIO]: Cannot calculate Prx without 'pathloss' annotation!"
            assert p.hasanno('pathloss'), errmsg
            pathloss = p.getanno('pathloss')
        Ptx = txpower
        PL  = pathloss
        if Grx is None: Grx = cls.Grx
        if Lrx is None: Lrx = cls.Lrx
        # calculate Prx
        Prx = Ptx - PL + Grx - Lrx
        return Prx

    @staticmethod
    def thermalnoise(bandwidth):
        """Compute thermal noise floor for communication systems with given
        bandwidth.

        :param bandwidth: Bandwidth of system (in Hz).
        :return: Thermal noise noise power (in dBm).

        See `Thermal Noise`_ for more on how the thermal noise floor of a
        communications system is computed in this method.

        .. _`Thermal Noise`: http://en.wikipedia.org/wiki/Thermal_noise#Noise_in_decibels
        """
        noise_dbm = -174.01 + linear2db(bandwidth)
        return noise_dbm

    def log_send(self, p, *args, **kwargs):
        """Overloaded to log 'cif-txpower' annotation."""
        if p.hasanno('cif-txpower'):
            txpower = p.getanno('cif-txpower')
            kwargs['cif-txpower'] = "%.2f dBm"%(txpower)
        ChannelInterface.log_send(self, p, *args, **kwargs)

    def log_recv(self, p, *args, **kwargs):
        """Overload to log 'cif-txpower', 'rxpower', and 'pathloss' annotation."""
        if p.hasanno('cif-txpower'):
            txpower = p.getanno('cif-txpower')
            kwargs['cif-txpower'] = "%.2f dBm"%(txpower)
        if p.hasanno('rxpower'):
            rxpower = p.getanno('rxpower')
            kwargs['rxpower'] = "%.2f dBm"%(rxpower)
        if p.hasanno('noisepower'):
            noisepower = p.getanno('noisepower')
            kwargs['noisepower'] = "%.2f dBm"%(noisepower)
        if p.hasanno('pathloss'):
            pathloss = p.getanno('pathloss')
            kwargs['pathloss'] = "%.2f dB"%(pathloss)
        ChannelInterface.log_recv(self, p, *args, **kwargs)

    def log_drop(self, p, *args, **kwargs):
        """Overload to log 'cif-txpower', 'rxpower', and 'pathloss' annotation."""
        if p.hasanno('cif-txpower'):
            txpower = p.getanno('cif-txpower')
            kwargs['cif-txpower'] = "%.2f dBm"%(txpower)
        if p.hasanno('rxpower'):
            rxpower = p.getanno('rxpower')
            kwargs['rxpower'] = "%.2f dBm"%(rxpower)
        if p.hasanno('noisepower'):
            noisepower = p.getanno('noisepower')
            kwargs['noisepower'] = "%.2f dBm"%(noisepower)
        if p.hasanno('pathloss'):
            pathloss = p.getanno('pathloss')
            kwargs['pathloss'] = "%.2f dB"%(pathloss)
        ChannelInterface.log_drop(self, p, *args, **kwargs)
