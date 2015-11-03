#!  /usr/bin/env python

"""
Base class for implementing digital modulation; contains `Modulation` class.

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
from scapy.all import Packet

class Modulation(Traceable):
    """Abstract base class for implementing digital modulation.

    Concrete subclasses must implement the `ber()` method.
    """
    name = "modulation"
    tracename = "MOD"
    def __init__(self, **kwargs):
        """Constructor."""
        Traceable.__init__(self, **kwargs)

    def ber(self):
        """Abstract method to determine bit-error rate; **must be implemented to
        create concrete subclass**."""
        raise RuntimeError, "[MODULATION]: Abstract method ber not defined!"

    def per(self, p, *args, **kwargs):
        """Calculate packet-error rate (PER); **overload as needed**.

        :param p: Packet (packet length in bytes) to compute PER for.
        :param args: Arguments passed to `ber()`.
        :param kwargs: Additional keywords also passed to `ber()`.
        :return: Packet-error rate (PER) in [0,1].

        By default this method uses the following relationship to determine the
        packet-error rate:

            PER = 1 - (1 - BER)^Nbits

        where BER is the bit-error rate (as determined by calling `ber()` with
        args and kwargs) and Nbits is the number of bits in packet `p`.

        **Overload this method to change how packet-error rate is computed.**
        """
        plen = p
        ber = self.ber(*args, **kwargs)
        if isinstance(p, Packet): plen = len(p)
        nbits = 8*plen
        per = 1 - (1-ber)**nbits
        return per
