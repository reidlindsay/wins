#!  /usr/bin/env python

"""
Model for quadrature amplitude modulation (QAM); contains `MQAM` class.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2010-10-19 16:34:02 -0500 (Tue, 19 Oct 2010) $
* $LastChangedRevision: 4811 $

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

from wins.digital.modulation import Modulation
from wins.helper import nchoosek, db2linear

from scapy.all import Packet
from numpy import sqrt, ceil, log2
from scipy.stats import erfc

class MQAM(Modulation):
    """Models M-order quadrature amplitude modulation.

    `Quadrature amplitude modulation`_ is a common digital modulation technique
    to convert a bitstream into a complex waveform.

    :CVariables:
     * `Nbps`: Dictionary containing number of bits per symbol indexed by
       supported modulation type (e.g. BPSK, QPSK, 4QAM, etc.). The table below
       lists the supported modulation type enumerations:

            =============== ====== ===========
            Modulation type  Nbps  M (`order`)
            =============== ====== ===========
                 'BPSK'       1        2
            --------------- ------ -----------
             'QPSK','4QAM'    2        4
            --------------- ------ -----------
                '16QAM'       4       16
            --------------- ------ -----------
                '64QAM'       6       64
            =============== ====== ===========

     * `mtype`: Property to access/modify modulation type enumeration (see
       `Nbps` for supported modulation types).

     * `nbps`: Property to access/modify number of bits per symbol.

     * `order`: Property to access/modify modulation order, which is the number
       of symbols in M-QAM constellation.

    :IVariables:
     * `__mtype`: Private variable maintains modulation type; use property
       `mtype` to access/modify.

    .. _`Quadrature amplitude modulation`: http://en.wikipedia.org/wiki/Quadrature_amplitude_modulation
    """
    name = "modulation"
    tracename = "MOD"
    Nbps = {"BPSK":1, "QPSK":2, "4QAM":2, "16QAM":4, "64QAM":6}
    def __init__(self, mtype=None, nbps=None, order=None, **kwargs):
        """Constructor.

        :param mtype: Modulation type enumeration (see `mtype`).
        :param nbps: Number of bits per symbol (see `nbps`).
        :param order: Modulation order (see `order`).

        By default, this class uses BPSK modulation. Either `mtype`, `nbps`, or
        `order` can be used to configure this object.
        """
        self.__mtype = None
        if (order is not None): nbps = log2(order)
        if (mtype is not None):
            assert (mtype.upper() in self.Nbps.keys() ), \
                    "[MQAM]: Modulation type (%s) is not supported!"%(mtype)
            mtype = mtype.upper()
            nbps  = self.Nbps[mtype]
        elif (nbps is not None):
            mtype = self.nbps2mtype(nbps)
        else:
            mtype, nbps = "BPSK", 1
        assert (mtype in self.Nbps) and (nbps==self.Nbps[mtype]), \
               "[MQAM]: Unsupported modultation (%s, nbps=%s)!"(mtype,nbps)
        Modulation.__init__(self, **kwargs)
        # set parameters
        self.mtype = mtype
        self.nbps  = nbps

    mtype = property(fget=lambda self: self.__mtype, \
                     fset=lambda self,m: self.set_mtype(m) )
    nbps  = property(fget=lambda self: self.Nbps[self.mtype], \
                     fset=lambda self,n: self.set_nbps(n) )
    order = property(fget=lambda self: 2**self.nbps, \
                     fset=lambda self,v: self.set_nbps(log2(v)) )

    def ber(self, snr):
        """Calculate bit-error rate (BER) corresponding to the signal-to-noise
        ratio (SNR) provided.

        :param snr: Signal-to-noise ratio (in dB).
        :return: Bit-error rate (BER) in [0,1].

        This method calls `calcber()` with local parameters to compute the
        bit-error rate.

        :note: If `snr` is a `numpy` array, this method will return an array of
               values corresponding to the specified parameters.
        """
        return self.calcber(snr, nbps=self.nbps)

    def per(self, p, snr):
        """Calculate packet-error rate (PER) corresponding to the
        signal-to-noise ratio (SNR) provided.

        :param snr: Signal-to-noise ratio (in dB).
        :return: Packet-error rate (PER) in [0,1].

        This method calls `calcper()` with local parameters to compute the
        packet-error rate.
        """
        return self.calcper(p, snr, nbps=self.nbps)

    def set_mtype(self, m):
        """Set modulation type.

        :param m: Modulation type enumeration; must be supported by `Nbps`.
        """
        m = m.upper()
        assert (m in self.Nbps), \
                "[MQAM]: Modulation type (%s) is not supported!"%(m)
        self.__mtype = m

    def set_nbps(self, n):
        """Set number of bits per symbol.

        :param n: Number of bits per symbol; must be supported by `Nbps`.
        """
        assert (n in self.Nbps.values() ), \
                "[MQAM]: Nbps (%s) is not supported!"%(n)
        self.mtype = self.nbps2mtype(n)

    @classmethod
    def nbps2mtype(cls, nbps):
        """Convenience method for converting number of bits per symbol to a
        proper modulation type enumeration, i.e. a reverse mapping of `Nbps`.

        :param nbps: Number of bits per symbol; must be supported by `Nbps`.
        :return: Modulation type enumeration.
        """
        NBPS = cls.Nbps
        assert (nbps in NBPS.values() ), \
                "[MQAM]: Nbps (%s) is not supported!"%(n)
        mtype = NBPS.keys()[NBPS.values().index(nbps)]
        return mtype

    @classmethod
    def calcber(cls, snr, mtype=None, nbps=None, order=None):
        """Calculate bit-error rate (BER) corresponding to the signal-to-noise
        ratio (SNR) provided in an AWGN channel.

        :param snr: Signal-to-noise ratio (in dB).
        :param mtype: Modulation type enumeration (see `mtype`).
        :param nbps: Number of bits per symbol (see `nbps`).
        :param order: Modulation order (see `order`).
        :return: Bit-error rate (BER) in [0,1].

        By default, this class uses BPSK modulation. Use any of the optional
        parameters `mtype`, `nbps`, or `order` to specify a different modulation
        type. The bit-error rate (BER) calculated in this method assumes an
        additive white gaussian noise (AWGN) channel. See `M-QAM Performance`_
        for more on this calculation.

        :note: If `snr` is a `numpy` array, this method will return an array of
               values corresponding to the specified parameters.

        .. _`M-QAM Performance`: http://en.wikipedia.org/wiki/Quadrature_amplitude_modulation#Quantized_QAM_performance
        """
        snrdb = snr
        # determine modulation type
        if (order is not None): nbps = log2(order)
        if (mtype is not None):
            assert (mtype.upper() in cls.Nbps.keys() ), \
                    "[MQAM]: Modulation type (%s) is not supported!"%(mtype)
            mtype = mtype.upper()
            nbps  = cls.Nbps[mtype]
        elif (nbps is not None):
            mtype = cls.nbps2mtype(nbps)
        else:
            mtype, nbps = "BPSK", 1
        assert (mtype in cls.Nbps) and (nbps==cls.Nbps[mtype]), \
               "[MQAM]: Unsupported modultation (%s, nbps=%s)!"(mtype,nbps)
        # calculate BER
        snr = db2linear(snrdb)
        M = 2**nbps
        if (M == 2):
            x = sqrt(2.0*snr)
            ser = 0.5*erfc(x/sqrt(2) )
            ber = ser
        else:
            x = sqrt(3.0*snr/(M-1))
            ter = 2.0*(1 - sqrt(1.0/M))*0.5*erfc(x/sqrt(2) )
            ser = 1 - (1-ter)**2.0
            ber = ser/log2(M)
        return ber

    @classmethod
    def calcper(cls, p, *args, **kwargs):
        """Calculate packet-error rate; **overload as needed**.

        :param p: Packet (packet length in bytes) to compute PER for.
        :param args: Arguments passed to `calcber()`.
        :param kwargs: Additional keywords also passed to `calcber()`.
        :return: Packet-error rate (PER) in [0,1].

        By default this method uses the following relationship to determine the
        packet-error rate:

            PER = 1 - (1 - BER)^Nbits

        where BER is the bit-error rate (as determined by calling `calcber()`
        with args and kwargs) and Nbits is the number of bits in packet `p`.

        **Overload this method to change how packet-error rate is computed.**
        """
        plen = p
        ber = cls.calcber(*args, **kwargs)
        if isinstance(p, Packet): plen = len(p)
        nbits = 8*plen
        per = 1 - (1-ber)**nbits
        return per
