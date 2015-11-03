#!  /usr/bin/env python

"""
Packet definitions, enumerations, and helper functions for RBAR protocol.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-08-29 14:09:23 -0500 (Mon, 29 Aug 2011) $
* $LastChangedRevision: 5117 $

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

:var RBAR_TARGET_BER:
    Target bit-error rate (BER) threshold for adaptation strategy.

:var RBAR_THRESHOLD:
    Dictionary containing rate index and threshold value pairs. These values
    have been set based on a target BER of 10^-6 for a single antenna system
    over an AWGN channel (`RBAR_TARGET_BER`).

:var DOT11_FC_SUBTYPE_RSH: Subtype enumeration for Reservation Subheader (RSH).
"""
__docformat__ = "restructuredtext en"

from scapy.all import Dot11, Packet
from wins.ieee80211.dot11_support import Dot11Data, \
        DOT11_FC_TYPE_DATA, DOT11_FC_SUBTYPE_DATA

RBAR_TARGET_BER = 1e-6
RBAR_THRESHOLD = {0:4.2, 1:7.2, 2:10.1, 3:13.8, 4:16.9, 5:21.7, 6:22.9, 7:24.7}

DOT11_FC_SUBTYPE_RSH = 0x08

class Dot11RSH(Dot11):
    """IEEE 802.11 Data.

    This class inherits from the `Dot11` class from Scapy.

    -------------------------------------------------------------------
    | FC | Duration/ID | Addr1 | Addr2 | Addr3 | SC | Addr4 | Payload |
    -------------------------------------------------------------------
    """
    name = "802.11 RSH"
    def __init__(self, *args, **kwargs):
        """Constructor; initialize header parameters for Data."""
        Dot11.__init__(self, *args, **kwargs)
        self.type = "Data"
        self.subtype = DOT11_FC_SUBTYPE_RSH

def isdot11rsh(p):
    """Is packet a Reservation Subheader (RSH)?"""
    data, isrsh = None, False
    if isinstance(p, Packet):
        if p.haslayer(Dot11): data = p[Dot11]
        elif p.haslayer(Dot11RSH): data = p[Dot11RSH]
    if data:
        isrsh = (data.type==DOT11_FC_TYPE_DATA) \
                and (data.subtype==DOT11_FC_SUBTYPE_RSH)
    return isrsh

def get_dot11rsh(p):
    """Extract Reservation Subheader (RSH) from packet `p`.

    :return: IEEE 802.11 Data packet or `None` if no matching layer found.
    """
    data = None
    if isdot11rsh(p):
        if p.haslayer(Dot11): data = p[Dot11]
        elif p.haslayer(Dot11RSH): data = p[Dot11RSH]
        else: raise RuntimeError, \
                "[RBAR]: get_dot11rsh() failed to find matching layer!"
    return data
