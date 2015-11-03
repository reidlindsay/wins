#!  /usr/bin/env python

"""
Model for (133,171) rate-compatible punctured convolutional (RCPC) code;
contains `RCPC` class.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-04-28 14:37:55 -0500 (Thu, 28 Apr 2011) $
* $LastChangedRevision: 4954 $

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

:var USE_A_D:
    Constant specifying if a_d parameter should be used in computation of
    convolutional coded error-rate (instead of c_d).
"""
__docformat__ = "restructuredtext en"

from wins.trace import Traceable
from wins.helper import nchoosek

from numpy import sqrt, ceil, array
from numpy import log2
from scipy.stats import erfc

from scapy.all import Packet

##  FIXME: 
#   Using c_d instead of a_d seems to produce a more pessimistic result (i.e. a
#   more conservative guess of the error-probability.
USE_A_D = 0

class RCPC(Traceable):
    """Model for (133,171) rate-compatible punctured convolutional (RCPC) code.

    This RCPC code is commonly used in many wireless standards (i.e. IEEE
    802.11a/g/n). See `Punctured convolutional codes`_ for more.

    :CVariables:
     * `params`: Dictionary containing parameters used for characterizing
       bit-error (BER) performance of the (133,171) RCPC code.
       
       Each dictionary key is a string representation of the code rate (e.g.
       "1/2", "2/3", etc.), and each entry is another dictionary containing
       three fields:

        * 'a_d': The distance coefficients used in calculating BER.
        * 'c_d': Complementary coefficients also used in calculating BER.
        * 'dfree': `Free distance`_ of code rate.

     * `rates`: Property to access code rates supported by this class.

     * `rate`: Property to access/modify coding rate of `RCPC`.
     * `coderate`: Property to alias `rate`.

    :IVariables:
     * `__rate`: Private variable to maintain coding rate.

    .. _`Punctured convolutional codes`: http://en.wikipedia.org/wiki/Convolutional_code#Popular_convolutional_codes
    .. _`Free distance`: http://en.wikipedia.org/wiki/Convolutional_code#Free_distance_and_error_distribution
    """
    name = "rcpc code"
    tracename = "RCPC"
    params = {'1/2':{'dfree':10,
                     'a_d':[11,0,38,0,193,0,1331,0,7275,0,40406,0,234969,0,1337714,0,7594819,0,43375588,0],
                     'c_d':[36,0,211,0,1404,0,11633,0,77433,0,502690,0,3322763,0,21292910,0,134365911,0,843425871,0]},
              '2/3':{'dfree':6,
                     'a_d':[1,16,48,158,642,2435,9174,34705,131585,499608],
                     'c_d':[3,70,285,1276,6160,27128,117019,498860,2103891,8784123]},
              '3/4':{'dfree':5,
                     'a_d':[8,31,160,892,4512,23307,121077,625059,3234886,16753077],
                     'c_d':[42,201,1492,10469,62935,379644,2253373,13073811,75152755,428005675]},
              '4/5':{'dfree': 4,
                     'a_d':[3,24,172,1158,7409,48729,319861,2097971,13765538,90315667],
                     'c_d':[12,188,1732,15256,121372,945645,7171532,53399130,392137968,2846810288]},
              '5/6':{'dfree':4,
                     'a_d':[14,69,654,4996,39699,315371,2507890,19921920,158275483,1257455600],
                     'c_d':[92,528,8694,79453,792114,7375573,67884974,610875423,5427275376,47664215639]},
              '6/7':{'dfree':3,
                     'a_d':[1,20,223,1961,18093,169175,1576108,14656816,136394365,1269388805,11812597193],
                     'c_d':[5,169,2725,32233,370861,4169788,45417406,483171499,5063039490,52394710081,536645404278]},
              '7/8':{'dfree':3,
                     'a_d':[2,46,499,5291,56179,599557,6387194,68117821,726098696,7741086706],
                     'c_d':[9,500,7437,105707,1402743,17909268,222292299,2706822556,32434972565,383973015410]},
              }
    def __init__(self, rate=None, **kwargs):
        """Constructor.

        :param rate: Coding rate [default='1/2']
        
        If specified, `rate` **must** be supported by `params`. By default,
        `rate` is set to the mother code rate (i.e. '1/2').
        """
        self.__rate = None
        if rate is None: rate = "1/2"
        Traceable.__init__(self, **kwargs)
        # set parameters
        self.rate = rate

    rate = property(fget=lambda self: self.__rate, \
                    fset=lambda self,r: self.set_rate(r) )
    coderate = property(fget=lambda self: self.rate, \
                        fset=lambda self,r: setattr(self,'rate',r) )
    rates = property(fget=lambda self: self.params.keys() )

    def ber(self, uber):
        """Calculate bit-error rate of coded system (or coded BER) from BER of
        uncoded system (or uncoded BER).

        :param uber: BER of uncoded system.
        :return: BER of coded system; in [0,1].

        This method calls `calcber()` with local parameters to compute the
        bit-error rate of the coded system.

        :note: If `uber` is a `numpy` array, this method will return an array of
               values corresponding to the specified parameters.
        """
        return self.calcber(uber, coderate=self.rate)

    def per(self, p, uber):
        """Calculate packet-error rate (PER); **overload as needed**.

        :param p: Packet (packet length in bytes) to compute PER for.
        :param uber: BER of uncoded system.
        :return: Packet-error rate (PER) in [0,1].

        This method calls `calcper()` with local parameters to compute the
        packet-error rate.

        **Overload this method to change how packet-error rate is computed.**
        """
        return self.calcper(p, uber, coderate=self.rate)

    def set_rate(self, r):
        """Set coding rate `r`."""
        assert (r in self.params), "[RCPC]: Unsupported rate (%s)!"%(r)
        self.__rate = r

    @classmethod
    def hard_decode(cls, uber, dfree):
        """Calculate bit-error probability for hard decision decoded RCPC.

        :param uber:  Probability of bit-error for uncoded system.
        :param dfree: Free distance of convolutional code.
        :return: BER for hard decision decoded RCPC coded system.
        """
        pd, p = 0, uber
        dmin = int(ceil((dfree+1.0)/2) )
        q = 1 - p
        for k in range(dmin, dfree+1):
            pd += nchoosek(dfree,k)*(p**k)*(q**(dfree-k) )
        if (int(dfree)%2<1):
            pd += 0.5*nchoosek(dfree,int(dfree/2))*(p**(dfree/2))*(q**(dfree/2))
        return pd

    @classmethod
    def calcber(cls, uber, coderate=None):
        """Calculate bit-error rate of coded system (or coded BER) from BER of
        uncoded system (or uncoded BER).

        :param uber: BER of uncoded system; in [0, 1].
        :param coderate: Optional code rate; **must be supported by `params`**.
        :return: BER of coded system; in [0,1].

        If `coderate` is not specified this method will use the rate of the
        mother code (i.e. "1/2") when calculating BER.

        :note: If `uber` is a `numpy` array, this method will return an array of
               values corresponding to the specified parameters.

        :note: This method computes an upper bound on the probability of
               bit-error for a coded system utilizing hard-decision decoding.
               This is modeled by calling `hard_decode()`.
        """
        global USE_A_D
        if coderate is None: coderate = "1/2"
        assert (coderate in cls.params), \
                "[RCPC]: Unsupported code rate (%s) in calcber()!"%(coderate)
        # get parameters
        dfree = cls.params[coderate]['dfree']
        a_d   = cls.params[coderate]['a_d']
        c_d   = cls.params[coderate]['c_d']
        numd = len(c_d)
        if USE_A_D: numd = len(a_d)
        # compute coded BER
        cber = 0*uber
        for k in range(numd):
            if USE_A_D: cber += a_d[k]*cls.hard_decode(uber, dfree+k)
            else:       cber += c_d[k]*cls.hard_decode(uber, dfree+k)
        # apply ber ceiling
        blim = 0.5
        try:
            for k in range(len(uber)):
                #blim = ber[k]
                if cber[k]>blim: cber[k] = blim
        except TypeError:
            if cber>blim: cber = blim
        return cber

    @classmethod
    def calcper(cls, p, *args, **kwargs):
        """Calculate packet-error rate of coded system; **overload as needed**.

        :param p: Packet (packet length in bytes) to compute PER for.
        :param args: Arguments passed to `calcber()`.
        :param kwargs: Additional keywords also passed to `calcber()`.
        :return: PER of coded system; in [0,1].

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
