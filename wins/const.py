#!  /usr/bin/env python

"""
Defines all constants used in `wins`.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-09-23 12:36:44 -0500 (Fri, 23 Sep 2011) $
* $LastChangedRevision: 5140 $

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

To access/modify the constants in this module, ``import const`` and then directly
access/modify ``const.VARNAME``:

    from wins import const
    ...
    print const.VARNAME
    const.VARNAME = ...

:var EPSILON:
    Epsilon time duration.

:var SPEED_OF_LIGHT:
    Constant for speed of light in meters/second (3.0 x 10^8).

:var ETHERTYPE_IP:
    Enumeration for EtherType that indicates the payload of a Ethernet frame is
    an IP Packet.

:var ETHERTYPE_ARP:
    Enumeration for EtherType that indicates the payload of a Ethernet frame is
    an ARP Packet.

:var IP_NULL_ADDR:
    String representation of IP (or layer 3) null address.

:var IP_BROADCAST_ADDR:
    String representation of IP (or layer 3) broadcast address.

:var ETHERNET_NULL_ADDR:
    String representation of Ethernet (or layer 2) null address.

:var ETHERNET_BROADCAST_ADDR:
    String representation of Ethernet (or layer 2) broadcast address.

:var ARP_PTYPE_IP:
    Enumeration for protocol type in ARP Packets to indicate that IP formatted
    addresses are used in protcol.

:var ARP_HTYPE_ETHERNET:
    Enumeration for hardware type in ARP Packets to indicate Ethernet formatted
    addressing is used in hardware.

:var IP_PROTO_UDP:
    Enumeration for IP protocol code corresponding to User Datagram Protocol.

:var IP_PROTO_IP:
    Enumeration for IP protocol code corresponding to IP in IP (encapsulation).

:var IP_PROTO_TCP:
    Enumeration for IP protocol code corresponding to Transmission Control
    Protocol.

:var IP_PROTO_ICMP:
    Enumeration for IP protocol code corresponding to Internet Control Message
    Protocol.

:var IP_PROTO_NONE:
    Enumeration for IP protocol code indicating no additional headers or payload.

:var IP_PROTO_DSR:
    Enumeration for IP protocol code indicating Dynamic Source Routing packet.

:var IP_PROTO_HOST:
    Enumeration for IP protocol code indicating an internal host protocol. This
    is used as a default when some unidentifiable packet is being used.

"""
__docformat__ = "restructuredtext en"

EPSILON = 1e-12
SPEED_OF_LIGHT = 3.0e8

ETHERTYPE_IP  = 0x0800
ETHERTYPE_ARP = 0x0806

IP_NULL_ADDR = "0.0.0.0"
IP_BROADCAST_ADDR = "255.255.255.255"
ETHERNET_NULL_ADDR = "00:00:00:00:00:00"
ETHERNET_BROADCAST_ADDR = "FF:FF:FF:FF:FF:FF"

ARP_PTYPE_IP = 0x0800
ARP_HTYPE_ETHERNET = 1

IP_PROTO_IP   = 4
IP_PROTO_TCP  = 6
IP_PROTO_UDP  = 17
IP_PROTO_ICMP = 1

IP_PROTO_NONE = 59
IP_PROTO_DSR  = 48
IP_PROTO_HOST = 61

ANNO_BLACKLIST = ['dot11n-rxwaveform', \
                  'dot11n-rxnoise', 'dot11n-rxadded', \
                  'dot11n-channel', 'cif-iheap', \
                  'cif-collision']
