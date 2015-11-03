#!  /usr/bin/env python

"""
Packet definitions, enumerations, and helper functions for IEEE 802.11 protocol.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-09-27 22:15:57 -0500 (Tue, 27 Sep 2011) $
* $LastChangedRevision: 5167 $

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

:var DOT11_VERBOSE:
    Constant enumeration to control verbose threshold of `Dot11MAC`.

    Setting the verbose level of a `Dot11MAC` above this threshold will cause
    the corresponding output in this class to be written (or logged).

:var DOT11_CWMIN: Minimum contention window size for IEEE 802.11 MAC.

:var DOT11_CWMAX: Maximum contention window size for IEEE 802.11 MAC.

:var DOT11_RETRYLIMIT: Maximum number of retries allowed.

:var DOT11_FC_TYPE_DATA: Type enumeration for data.

:var DOT11_FC_TYPE_CONTROL: Type enumeration for control.

:var DOT11_FC_SUBTYPE_DATA: Subtype enumeration for data.

:var DOT11_FC_SUBTYPE_ACK: Subtype enumeration for acknowledgement.

:var DOT11_FC_SUBTYPE_RTS: Subtype enumeration for request-to-send message.

:var DOT11_FC_SUBTYPE_CTS: Subtype enumeration for clear-to-send message.

:var DOT11_FC_TODS: Frame control bit mask for to-DS field.

:var DOT11_FC_FROMDS: Frame control bit mask for from-DS field.

:var DOT11_FC_MF: Frame control bit mask for more fragment (MF) field.

:var DOT11_FC_RETRY: Frame control bit mask for retry field.

:var DOT11_FC_PWMGT: Frame control bit mask for power management (pw-mgt) field.

:var DOT11_FC_MD: Frame control bit mask for more data (MD) field.

:var DOT11_FC_WEP: Frame control bit mask for wep (WEP) field.

:var DOT11_FC_ORDER: Frame control bit mask for order field.

:var DOT11A_PHY_MODE: Enumeration for IEEE 802.11a physical layer.

:var DOT11B_PHY_MODE: Enumeration for IEEE 802.11b physical layer.

:var DOT11G_PHY_MODE: Enumeration for IEEE 802.11g physical layer.

:var DOT11N_PHY_MODE: Enumeration for IEEE 802.11n physical layer.

:var DOT11_TIMING: Dictionary containing IEEE 802.11 timing parameters for
                   various physical layer modes.
"""
__docformat__ = "restructuredtext en"

from scapy.all import Dot11, Packet

DOT11_VERBOSE = 53

DOT11_CWMIN = 32
DOT11_CWMAX = 1024
DOT11_RETRYLIMIT = 7

DOT11_FC_TYPE_DATA    = 0x02
DOT11_FC_TYPE_CONTROL = 0x01
DOT11_FC_SUBTYPE_DATA = 0x00
DOT11_FC_SUBTYPE_RTS  = 0x0B
DOT11_FC_SUBTYPE_CTS  = 0x0C
DOT11_FC_SUBTYPE_ACK  = 0x0D
DOT11_FC_TODS   = 0x01
DOT11_FC_FROMDS = 0x02
DOT11_FC_MF     = 0x04
DOT11_FC_RETRY  = 0x08
DOT11_FC_PWMGT  = 0x10
DOT11_FC_MD     = 0x20
DOT11_FC_WEP    = 0x40
DOT11_FC_ORDER  = 0x80

DOT11A_PHY_MODE = "IEEE 802.11a"
DOT11B_PHY_MODE = "IEEE 802.11b"
DOT11G_PHY_MODE = "IEEE 802.11g"
DOT11N_PHY_MODE = "IEEE 802.11n"
DOT11_TIMING = { \
        DOT11A_PHY_MODE: {'sifs': 16e-6, 'slottime':  9e-6}, \
        DOT11B_PHY_MODE: {'sifs': 10e-6, 'slottime': 20e-6}, \
        DOT11G_PHY_MODE: {'sifs': 10e-6, 'slottime':  9e-6}, \
        DOT11N_PHY_MODE: {'sifs': 16e-6, 'slottime':  9e-6} \
        }

class Dot11Data(Dot11):
    """IEEE 802.11 Data.

    This class inherits from the `Dot11` class from Scapy.

    -------------------------------------------------------------------
    | FC | Duration/ID | Addr1 | Addr2 | Addr3 | SC | Addr4 | Payload |
    -------------------------------------------------------------------
    """
    name = "802.11 DATA"
    def __init__(self, *args, **kwargs):
        """Constructor; initialize header parameters for Data."""
        Dot11.__init__(self, *args, **kwargs)
        self.type = "Data"
        self.subtype = DOT11_FC_SUBTYPE_DATA

class Dot11Ack(Dot11):
    """IEEE 802.11 Acknowledgement.

    This class inherits from the `Dot11` class from Scapy.

    -------------------------
    | FC | Duration/ID | RA |
    -------------------------
    """
    name = "802.11 ACK"
    fields_desc = []
    def __init__(self, *args, **kwargs):
        """Constructor; initialize header parameters for ACK."""
        Dot11.__init__(self, *args, **kwargs)
        self.type = "Control"
        self.subtype = DOT11_FC_SUBTYPE_ACK

    @staticmethod
    def _init_fields_desc():
        """Internal method to redefine `Dot11Ack.fields_desc`; this removes
        unnecessary fields from `Dot11`."""
        Dot11Ack.fields_desc = Dot11.fields_desc[0:6]

Dot11Ack._init_fields_desc()

class Dot11RTS(Dot11):
    """IEEE 802.11 Request-To-Send message.

    This class inherits from the `Dot11` class from Scapy.

    ------------------------------
    | FC | Duration/ID | RA | TA |
    ------------------------------
    """
    name = "802.11 RTS"
    fields_desc = []
    def __init__(self, *args, **kwargs):
        """Constructor; initialize header parameters for RTS."""
        Dot11.__init__(self, *args, **kwargs)
        self.type = "Control"
        self.subtype = DOT11_FC_SUBTYPE_RTS

    @staticmethod
    def _init_fields_desc():
        """Internal method to redefine `Dot11RTS.fields_desc`; this removes
        unnecessary fields from `Dot11`."""
        Dot11RTS.fields_desc = Dot11.fields_desc[0:7]

Dot11RTS._init_fields_desc()

class Dot11CTS(Dot11):
    """IEEE 802.11 Clear-To-Send message.

    This class inherits from the `Dot11` class from Scapy.

    -------------------------
    | FC | Duration/ID | RA |
    -------------------------
    """
    name = "802.11 CTS"
    fields_desc = []
    def __init__(self, *args, **kwargs):
        """Constructor; initialize header parameters for CTS."""
        Dot11.__init__(self, *args, **kwargs)
        self.type = "Control"
        self.subtype = DOT11_FC_SUBTYPE_CTS

    @staticmethod
    def _init_fields_desc():
        """Internal method to redefine `Dot11CTS.fields_desc`; this removes
        unnecessary fields from `Dot11`."""
        Dot11CTS.fields_desc = Dot11.fields_desc[0:6]

Dot11CTS._init_fields_desc()

def isdot11data(p):
    """Is packet an IEEE 802.11 Data packet?"""
    data, isdata = None, False
    if isinstance(p, Packet):
        if p.haslayer(Dot11): data = p[Dot11]
        elif p.haslayer(Dot11Data): data = p[Dot11Data]
    if data:
        isdata = (data.type==DOT11_FC_TYPE_DATA) \
                 and (data.subtype==DOT11_FC_SUBTYPE_DATA)
    return isdata

def isdot11ack(p):
    """Is packet an IEEE 802.11 Acknowledgement packet?"""
    ack, isack = None, False
    if isinstance(p, Packet):
        if p.haslayer(Dot11): ack = p[Dot11]
        elif p.haslayer(Dot11Ack): ack = p[Dot11Ack]
    if ack:
        isack = (ack.type==DOT11_FC_TYPE_CONTROL) \
                and (ack.subtype==DOT11_FC_SUBTYPE_ACK)
    return isack

def isdot11rts(p):
    """Is packet an IEEE 802.11 Request-to-send packet?"""
    rts, isrts = None, False
    if isinstance(p, Packet):
        if p.haslayer(Dot11): rts = p[Dot11]
        elif p.haslayer(Dot11RTS): rts = p[Dot11RTS]
    if rts:
        isrts = (rts.type==DOT11_FC_TYPE_CONTROL) \
                and (rts.subtype==DOT11_FC_SUBTYPE_RTS)
    return isrts

def isdot11cts(p):
    """Is packet an IEEE 802.11 Clear-to-send packet?"""
    cts, iscts = None, False
    if isinstance(p, Packet):
        if p.haslayer(Dot11): cts = p[Dot11]
        elif p.haslayer(Dot11CTS): cts = p[Dot11CTS]
    if cts:
        iscts = (cts.type==DOT11_FC_TYPE_CONTROL) \
                and (cts.subtype==DOT11_FC_SUBTYPE_CTS)
    return iscts

def get_dot11data(p):
    """Extract IEEE 802.11 Data from packet `p`.

    :return: IEEE 802.11 Data packet or `None` if no matching layer found.
    """
    data = None
    if isdot11data(p):
        if p.haslayer(Dot11): data = p[Dot11]
        elif p.haslayer(Dot11Data): data = p[Dot11Data]
        else: raise RuntimeError, \
                "[DOT11]: get_dot11data() failed to find matching layer!"
    return data

def get_dot11ack(p):
    """Extract IEEE 802.11 Acknowledgement (ACK) from packet `p`.

    :return: IEEE 802.11 ACK packet or `None` if no matching layer found.
    """
    ack = None
    if isdot11ack(p):
        if p.haslayer(Dot11): ack = p[Dot11]
        elif p.haslayer(Dot11Ack): ack = p[Dot11Ack]
        else: raise RuntimeError, \
                "[DOT11]: get_dot11ack() failed to find matching layer!"
    return ack

def get_dot11rts(p):
    """Extract IEEE 802.11 Request-to-send (RTS) from packet `p`.

    :return: IEEE 802.11 RTS packet or `None` if no matching layer found.
    """
    rts = None
    if isdot11rts(p):
        if p.haslayer(Dot11): rts = p[Dot11]
        elif p.haslayer(Dot11RTS): rts = p[Dot11RTS]
        else: raise RuntimeError, \
                "[DOT11]: get_dot11rts() failed to find matching layer!"
    return rts

def get_dot11cts(p):
    """Extract IEEE 802.11 Clear-to-send (CTS) from packet `p`.

    :return: IEEE 802.11 CTS packet or `None` if no matching layer found.
    """
    cts = None
    if isdot11cts(p):
        if p.haslayer(Dot11): cts = p[Dot11]
        elif p.haslayer(Dot11CTS): cts = p[Dot11CTS]
        else: raise RuntimeError, \
                "[DOT11]: get_dot11cts() failed to find matching layer!"
    return cts
