#!  /usr/bin/env python

"""
Helper methods and support classes for `ARP` implementation.

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

:var ARP_VERBOSE:
    Constant enumeration to control verbose thresholds in this file.
  
    Setting the verbose level of a `ARP` above this threshold will cause the
    corresponding output in this file to be written (or logged).

:var ARP_OPCODE_REQUEST:
    Enumeration for operation code corresponding to ARP request.

:var ARP_OPCODE_REPLY:
    Enumeration for operation code corresponding to ARP reply.
"""
__docformat__ = "restructuredtext en"

from SimPy.Simulation import hold
from wins.protocol.mac import MAC
from wins.packet import Packet, ANNO

from scapy.all import ARP as ARPpacket

ARP_VERBOSE = 45
ARP_OPCODE_REQUEST = 1
ARP_OPCODE_REPLY   = 2

class ARPRequest(ARPpacket):
    name = "ARP Request"
    def __init__(self, *args, **kwargs):
        """Contructor; initialize 'op' parameter to `ARP_OPCODE_REQUEST`."""
        ARPpacket.__init__(self, *args, **kwargs)
        self.op = ARP_OPCODE_REQUEST

class ARPReply(ARPpacket):
    name = "ARP Reply"
    def __init__(self, *args, **kwargs):
        """Contructor; initialize 'op' parameter to `ARP_OPCODE_REPLY`."""
        ARPpacket.__init__(self, *args, **kwargs)
        self.op = ARP_OPCODE_REPLY

def isarprequest(p):
    """Is packet an ARP request?"""
    arp, isreq = None, False
    if isinstance(p, Packet):
        if p.haslayer(ARPpacket): arp = p[ARPpacket]
        elif p.haslayer(ARPRequest): arp = p[ARPRequest]
    if arp:
        isreq = (arp.op == ARP_OPCODE_REQUEST)
    return isreq

def isarpreply(p):
    """Is packet an ARP reply?"""
    arp, isrep = None, False
    if isinstance(p, Packet):
        if p.haslayer(ARPpacket): arp = p[ARPpacket]
        elif p.haslayer(ARPReply): arp = p[ARPReply]
    if arp:
        isrep = (arp.op == ARP_OPCODE_REPLY)
    return isrep

def get_arprequest(p):
    """Extract ARP request from packet `p`.

    :return: ARP Request message or `None` if no matching layer found.
    """
    arp = None
    if isarprequest(p):
        if p.haslayer(ARPpacket): arp = p[ARPpacket]
        elif p.haslayer(ARPRequest): arp = p[ARPRequest]
        else: raise RuntimeError, \
                "[ARP]: get_arprequest() failed to find matching layer!"
    return arp

def get_arpreply(p):
    """Extract ARP reply from packet `p`.

    :return: ARP Reply message or `None` if no matching layer found.
    """
    arp = None
    if isarpreply(p):
        if p.haslayer(ARPpacket): arp = p[ARPpacket]
        elif p.haslayer(ARPReply): arp = p[ARPReply]
        else: raise RuntimeError, \
                "[ARP]: get_arpreply() failed to find matching layer!"
    return arp
