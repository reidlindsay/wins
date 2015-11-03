#!  /usr/bin/env python

"""
Model for Reed-Solomon error-correcting codes; contains `ReedSolomon` class.

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

from wins.trace import Traceable
from wins.helper import nchoosek

from numpy import log2, floor, ceil, iterable

from scapy.all import Packet

class ReedSolomon(Traceable):
    """Model for Reed-Solomon (R-S) error-correcting codes.

    This error-correcting code is used a variety of applications including CDs,
    wireless protocols, and satellite systems.

    :IVariables:
     * `blocklength`: Block length in symbols; *n*.
     * `messagelength`: Uncoded message length in symbols; *k*.
     * `bitspersymbol`: Number of bits per symbol; *m*.
     * `ncorrectable`: The number of correctable symbol errors; *t*.

    The FEC parameters *(n,k,m,t)* must satisfy the following:

        1. n <= 2^m-1
        #. n-k >= 2t

    The default value of 'm' and 't' are ceil(log2(n+1)) and (n-k) div 2,
    respectively. RS(255,239) and RS(32,38) codes have been used in various
    applications such as deep space communication, digital video broadcast by
    satellites (DVB-S), and data storage.
    """
    name = "reed-solomon"
    tracename = "RS"
    def __init__(self, n, k, m=None, t=None, **kwargs):
        """Constructor.

        :param n: Block length in symbols.
        :param k: Uncoded message length in symbols.
        :param m: Number of bits per symbol [default=ceil(log2(n+1))]
        :param t: Number of correctable errors [default=(n-k) div 2].
        """
        Traceable.__init__(self, **kwargs)
        # set parameters
        self.blocklength   = None
        self.messagelength = None
        self.bitspersymbol = None
        self.ncorrectable  = None
        self.set_rate(n, k, m, t)

    rate = property(fget=lambda self: self.get_rate(), \
                    fset=lambda self, r: self.set_rate(r) )

    def get_rate(self):
        """Get code rate (k/n)."""
        n = self.blocklength
        k = self.messagelength
        return (1.0*k)/n

    def set_rate(self, n, k=None, m=None, t=None):
        """Set rate parameters for code.

        :param n: Block length in symbols.
        :param k: Uncoded message length in symbols.
        :param m: Number of bits per symbol [default=ceil(log2(n+1))]
        :param t: Number of correctable errors [default=(n-k) div 2].
        :return: Code rate (k/n).

        The FEC parameters *(n,k,m,t)* must satisfy the following:

            1. n <= 2^m-1
            #. n-k >= 2t

        This method will assert an error if parameters are invalid.

        :note: This method may also be called with a 2-tuple (n,k) as the first
               argument. This is primarily to enable setting parameters using
               the '=' operator (see `rate` property).
        """
        if k is None:
            assert iterable(n) and (len(n)>1), "[R-S]: Invalid argument! " + \
                                "n must be a tuple if k is not specified!"
            x0, x1 = n[0], n[1]
            n, k = x0, x1
        else:
            assert not iterable(n), "[R-S]: Invalid argument! " + \
                                    "n cannot be a tuple if k is specified!"
        if m is None: m = int(ceil(log2(n+1) ) )
        if t is None: t = int(floor((n-k)/2) )
        self.blocklength   = n
        self.messagelength = k
        self.bitspersymbol = m
        self.ncorrectable  = t
        assert (n<=(2**m-1)) and (n-k>=2*t), \
                "[R-S]: Invalid parameters (n,k,m,t) = (%s,%s,%s,%s)!"%(n,k,m,t)
        return (1.0*k)/n

    def per(self, p, ber):
        """Calculate packet-error rate (PER).

        :param p: Packet (packet length in bytes) to compute PER for.
        :param ber: BER of uncoded system.
        :return: Packet-error rate (PER) in [0,1].

        :note: If `ber` is a `numpy` array, this method will return an array of
               values corresponding to the specified parameters.
        """
        n, k = self.blocklength, self.messagelength
        m, t = self.bitspersymbol, self.ncorrectable
        plen = p
        if isinstance(p, Packet): plen = len(p)
        # calculate uncorrectable error rate for symbols
        ser = 1 - (1 - ber)**m
        # calculate uncorrectable error rate
        s = 0
        for b in range(t+1):
            s +=nchoosek(n,b)*(ser**b)*((1-ser)**(n-b))
        uer = 1 - s
        # calculate PER
        nblocks = int(ceil(plen/k) )
        per = 1 - (1-uer)**nblocks
        return per

