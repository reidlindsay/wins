#!  /usr/bin/env python

"""
Packet definitions, enumerations, and helper functions for `DSR` protocol.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-10-23 04:50:01 -0500 (Sun, 23 Oct 2011) $
* $LastChangedRevision: 5281 $

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

:var DSR_MAX_SALVAGE_COUNT: Maximum number of route errors that can be salvaged.

:var DSR_OPTION_TYPE_RREQ: Enumeration for DSR Route Request Option Type.

:var DSR_OPTION_TYPE_RREP: Enumeration for DSR Route Reply Option Type.

:var DSR_OPTION_TYPE_RERR: Enumeration for DSR Route Error Option Type.

:var DSR_OPTION_TYPE_ACKREQ: Enumeration for DSR Acknowledgement Request Option Type.

:var DSR_OPTION_TYPE_ACK: Enumeration for DSR Acknowledgement Option Type.

:var DSR_OPTION_TYPE_SRCROUTE: Enumeration for DSR Source Route Option Type.

:var _dsr_types_names: Internal dictionary for DSR option type enumerations.

:var DSR_ERROR_NODE_UNREACHABLE: DSR error type.

:var DSR_ERROR_FLOW_STATE_NOT_SUPPORTED: DSR error type.

:var DSR_ERROR_OPTION_NOT_SUPPORTED: DSR error type.

:var _dsr_error_types: Internal dictionary for DSR error types.

:var DSR_CACHE_ALL_HOPS: If true, cache all intermediate hops in a source route.

:var DSR_USE_CACHE_REPLY: Allow intermediate node to reply to route request
                          using its rout cache (if possible).

:var DSR_USE_GRATUITOUS_RREQ: Send gratuitous RREQ+RERR for routes broken by an
                              incoming route error message.
"""
__docformat__ = "restructuredtext en"

import struct
from scapy.all import Packet
from scapy.all import BitField, BitEnumField, \
                      ByteField, ByteEnumField, \
                      ShortField
from scapy.all import FieldLenField, StrLenField, \
                      FieldListField, PacketListField
from scapy.all import IPField, IP_PROTOS
from scapy.all import IP, UDP, TCP, ICMP, bind_layers

from wins import const

import sys

# Protocol Configuration
DSR_CACHE_ALL_HOPS = 1
DSR_USE_REVERSE_ROUTES = 1
DSR_USE_CACHE_REPLY = 0
DSR_USE_GRATUITOUS_RREQ = 0

# Protocol Definitions
DSR_MAX_SALVAGE_COUNT = 15
DSR_OPTION_TYPE_RREQ = 1
DSR_OPTION_TYPE_RREP = 2
DSR_OPTION_TYPE_RERR = 3
DSR_OPTION_TYPE_ACKREQ = 160
DSR_OPTION_TYPE_ACK = 32
DSR_OPTION_TYPE_SRCROUTE = 96
DSR_OPTION_TYPE_PAD1 = 224
DSR_OPTION_TYPE_PADN = 0

_dsr_types_names  = { DSR_OPTION_TYPE_RREQ:   "route_request",
                      DSR_OPTION_TYPE_RREP:   "route_reply",
                      DSR_OPTION_TYPE_RERR:   "route_error",
                      DSR_OPTION_TYPE_ACKREQ: "acknowledgement_request",
                      DSR_OPTION_TYPE_ACK:    "acknowledgement",
                      DSR_OPTION_TYPE_SRCROUTE: "source_route",
                      DSR_OPTION_TYPE_PAD1:   "pad1",
                      DSR_OPTION_TYPE_PAD1:   "padN",
                     }

_dsr_option_names = { DSR_OPTION_TYPE_RREQ:   "RREQ",
                      DSR_OPTION_TYPE_RREP:   "RREP",
                      DSR_OPTION_TYPE_RERR:   "RERR",
                      DSR_OPTION_TYPE_ACKREQ: "ACKREQ",
                      DSR_OPTION_TYPE_ACK:    "ACK",
                      DSR_OPTION_TYPE_SRCROUTE: "SRCROUTE",
                      DSR_OPTION_TYPE_PAD1:   "PAD1",
                      DSR_OPTION_TYPE_PAD1:   "PADN",
                     }

DSR_ERROR_NODE_UNREACHABLE = 1
DSR_ERROR_FLOW_STATE_NOT_SUPPORTED = 2
DSR_ERROR_OPTION_NOT_SUPPORTED = 3

_dsr_error_types = {DSR_ERROR_NODE_UNREACHABLE: "node_unreachable",
                    DSR_ERROR_FLOW_STATE_NOT_SUPPORTED: "flow_state_not_supported",
                    DSR_ERROR_OPTION_NOT_SUPPORTED: "option_not_supported"}

class _DSROPT_HDR(Packet):
    """Internal packet structure for fields common among DSR options.

    Frame Format
    ============
    Contains two fields:
        * 'type': DSR Option Type (8 bits)
        * 'optdatalen: DSR Option Length (8 bits)
    """
    fields_desc = [ ByteEnumField("type", 0, _dsr_types_names),
                    ByteField("optdatalen", None) ]

class DSROPT(Packet):
    """Base class for DSR Option packet."""
    fields_desc = [ _DSROPT_HDR,
                    FieldLenField("length", None, fmt="B",
                                  length_of="value", adjust=lambda pkt,l:l+2),
                    StrLenField("value", "",length_from=lambda pkt:pkt.length-2) ]

    def extract_padding(self, p):
        return "",p

    registered_dsr_types = {}
    @classmethod
    def register_variant(cls):
        cls.registered_dsr_types[cls.type.default] = cls
    @classmethod
    def dispatch_hook(cls, pkt=None, *args, **kargs):
        if pkt:
            opt = ord(pkt[0])
            if opt in cls.registered_dsr_types:
                return cls.registered_dsr_types[opt]
        return cls

    def post_build(self, p, pay):
        if self.optdatalen is None:
            l = len(p)-2
            p = p[:1]+struct.pack("B", l)+p[2:]
        return p+pay
    
    def get_tracename(self):
        """Overloaded to use options to determine tracename."""
        name = "DSROPT"
        if self.type in _dsr_option_names:
            name = "%s"%(_dsr_option_names[self.type])
        return name

class DSROPT_RREQ(DSROPT):
    """DSR Route Request Option.

    Frame Format
    ============
    Contains following fields:
        * 'identification': Unique identifier generated by original sender (16 bits).
        * 'target': Address of target node (4 bytes).
        * 'addresses': Array of addresses in the Route Request (4*n bytes).
    """
    type = DSR_OPTION_TYPE_RREQ
    fields_desc = [ _DSROPT_HDR,
                    ShortField("identification", 0),
                    IPField("target", "127.0.0.1"),
                    FieldListField("addresses",[],IPField("","127.0.0.1"), 
                                   length_from=lambda pkt:pkt.optdatalen-6) ]
    def add_addresses(self, *args):
        if self.addresses: self.addresses += list(args)
        else:              self.addresses = list(args)

    def set_addresses(self, addr):
        self.addresses = addr

class DSROPT_RREP(DSROPT):
    """DSR Route Reply Option.

    Frame Format
    ============
    Contains following fields:
        * 'lasthop': 1-bit flag; indicates if last hop in route is external to
                     DSR network (1 bit).
        * 'rsvd': Reserved field (15 bits).
        * 'addresses': Array of addresses in the selected route (4*n bytes).
    """
    type = DSR_OPTION_TYPE_RREP
    fields_desc = [ _DSROPT_HDR,
                    BitEnumField("lasthop", 0, 1, {0:'internal', 1:'external'}),
                    BitField("rsvd", 0, 15),
                    FieldListField("addresses",[],IPField("","127.0.0.1"), 
                                   length_from=lambda pkt:pkt.optdatalen-2) ]
    def add_addresses(self, *args):
        if self.addresses: self.addresses += list(args)
        else:              self.addresses = list(args)

    def set_addresses(self, addr):
        self.addresses = addr

class DSROPT_RERR(DSROPT):
    """DSR Route Error Option.

    Frame Format
    ============
    Contains following fields:
        * 'errortype': DSR error type; see error enumerations for more (1 byte).
        * 'rsvd': Reserved field (4 bits).
        * 'salvage': Copied from salvage field of DSR Source Route (4 bits).
        * 'err_src': Address of node that discovered the error (4 bytes).
        * 'err_dst': Address of node Route Error must be delivered to (4 bytes).
        * 'type_specific': Type specific information for the error (4 bytes).
    """
    type = DSR_OPTION_TYPE_RERR
    fields_desc = [ _DSROPT_HDR,
                    ByteEnumField("errortype", 1, _dsr_error_types),
                    BitField("rsvd", 0, 4),
                    BitField("salvage", 0, 4),
                    IPField("err_src", "127.0.0.1"),
                    IPField("err_dst", "127.0.0.1"),
                    IPField("type_specific", "127.0.0.1") ]

class DSROPT_ACKREQ(DSROPT):
    """DSR Acknowledgement Request Option.

    Frame Format
    ============
    Contains following fields:
        * 'identification': Unique value generated by sender (16 bits).
    """
    type = DSR_OPTION_TYPE_ACKREQ
    fields_desc = [ _DSROPT_HDR,
                    ShortField("identification", 0) ]

class DSROPT_ACK(DSROPT):
    """DSR Acknowledgement Option.

    Frame Format
    ============
    Contains following fields:
        * 'identification': Copied from 'identification' field of ACK request (16 bits).
        * 'src': Address of node originating acknowledgement (4 bytes).
        * 'dst': Address of destination for acknowledgement (4 bytes).
    """
    type = DSR_OPTION_TYPE_ACKREQ
    fields_desc = [ _DSROPT_HDR,
                    ShortField("identification", 0),
                    IPField("src", "0.0.0.0"),
                    IPField("dst", "0.0.0.0") ]

class DSROPT_SRCROUTE(DSROPT):
    """DSR Source Route Option.

    Frame Format
    ============
    Contains following fields:
        * 'firsthop': 1-bit flag; indicates if first hop in route is external to
                     DSR network (1 bit).
        * 'lasthop': 1-bit flag; indicates if last hop in route is external to
                     DSR network (1 bit).
        * 'rsvd': Reserved field (4 bits).
        * 'salvage': Count of number of times packet has been salvaged (4 bits).
        * 'segsleft': Number of route segments remaining (6 bits).
        * 'addresses': Array of addresses in the source route (4*n bytes).
    """
    type = DSR_OPTION_TYPE_SRCROUTE
    fields_desc = [ _DSROPT_HDR,
                    BitEnumField("firsthop", 0, 1, {0:'internal', 1:'external'}),
                    BitEnumField("lasthop", 0, 1, {0:'internal', 1:'external'}),
                    BitField("rsvd", 0, 4),
                    BitField("salvage", 0, 4),
                    BitField("segsleft", 0, 6),
                    FieldListField("addresses",[],IPField("","127.0.0.1"), 
                                   length_from=lambda pkt:pkt.optdatalen-2) ]
    def add_addresses(self, *args):
        if self.addresses: self.addresses += list(args)
        else:              self.addresses = list(args)

    def set_addresses(self, addr):
        self.addresses = addr

class DSRPacket(Packet):
    """DSR packet; use specific subclasses for different options.

    Frame Format
    ============
    DSR header contains following fields:
        * 'nextheader': Indicates protocol following this header; same as
                        'proto' field in IP header (1 byte).
        * 'fsh': 1-bit flag; indicates flow state header is being used (1 bit).
        * 'rsvd': Reserved field (7 bits).
        * 'length': Length of options field after DSR header (2 bytes).
        * 'options': Variable length options.

    :cvar opt: Property to access option field.
    """
    name = "DSR"
    fields_desc = [ ByteEnumField("nextheader", const.IP_PROTO_NONE, IP_PROTOS), 
                    BitEnumField("fsh", 0, 1, {0:"dsr_option", 1:"dsr_flow_state"}),
                    BitField("rsvd", 0, 7),
                    ShortField("length", None),
                    PacketListField("options", [], DSROPT, length_from=lambda p:p.length) ]

    opt = property(fget=lambda self: self.options[0])
    
    def get_tracename(self):
        """Overloaded to use options to determine tracename."""
        name = "DSR"
        if self.options:
            for opt in self.options:
                if opt.type in _dsr_option_names:
                    name += "|%s"%(_dsr_option_names[opt.type])
        return name

    def post_build(self, p, pay):
        """Set variable 'length' field."""
        l = self.length
        if l is None:
            l = len(p)-4
            p = p[:2] + struct.pack("!H", l) + p[4:]
        return p+pay

    def add_options(self, *args):
        if self.options: self.options += list(args)
        else:            self.options = list(args)

    def set_options(self, opt):
        self.options = opt

def _update_IP_with_DSR():
    """Internal method to update `IP` for DSR."""
    proto = None
    for f in IP.fields_desc:
        if (f.name=='proto'): proto = f
    assert (proto is not None)
    proto.i2s[const.IP_PROTO_DSR] = 'dsr'
    # bind protocol layers so that proto field gets set automatically
    bind_layers(IP,        DSRPacket, frag=0, proto=const.IP_PROTO_DSR)
    bind_layers(DSRPacket, UDP,       nextheader=const.IP_PROTO_UDP)
    bind_layers(DSRPacket, TCP,       nextheader=const.IP_PROTO_TCP)
    bind_layers(DSRPacket, ICMP,      nextheader=const.IP_PROTO_ICMP)
    bind_layers(DSRPacket, IP,        nextheader=const.IP_PROTO_IP)

_update_IP_with_DSR()

def DSR_RREQ(*args, **kwargs):
    """Create DSR Route Request Packet."""
    kwargs['options'] = DSROPT_RREQ()
    p = DSRPacket(*args, **kwargs)
    return p

def DSR_RREP(*args, **kwargs):
    """Create DSR Route Reply Packet."""
    kwargs['options'] = DSROPT_RREP()
    p = DSRPacket(*args, **kwargs)
    return p

def DSR_RERR(*args, **kwargs):
    """Create DSR Route Error Packet."""
    kwargs['options'] = DSROPT_RERR()
    p = DSRPacket(*args, **kwargs)
    return p

def DSR_ACKREQ(*args, **kwargs):
    """Create DSR Acknowledgement Request Packet."""
    kwargs['options'] = DSROPT_ACKREQ()
    p = DSRPacket(*args, **kwargs)
    return p

def DSR_ACK(*args, **kwargs):
    """Create DSR Acknowledgement Packet."""
    kwargs['options'] = DSROPT_ACK()
    p = DSRPacket(*args, **kwargs)
    return p

def DSR_SRCROUTE(*args, **kwargs):
    """Create DSR Source Route Packet."""
    kwargs['options'] = DSROPT_SRCROUTE()
    p = DSRPacket(*args, **kwargs)
    return p

def getdsr_rreq(p):
    """Get `DSROPT_RREQ` option from packet `p`."""
    opt = p
    if isinstance(p, Packet) and p.haslayer(DSRPacket):
        dsr, opt = p[DSRPacket], None
        for x in dsr.options:
            # return first option that matches
            if (opt is None) and isinstance(x, DSROPT_RREQ): opt = x
    if not isinstance(opt, DSROPT_RREQ): opt = None
    if opt and (opt.type!=DSR_OPTION_TYPE_RREQ): opt = None
    return opt

def getdsr_rrep(p):
    """Get `DSROPT_RREP` option from packet `p`."""
    opt = p
    if isinstance(p, Packet) and p.haslayer(DSRPacket):
        dsr, opt = p[DSRPacket], None
        for x in dsr.options:
            # return first option that matches
            if (opt is None) and isinstance(x, DSROPT_RREP): opt = x
    if not isinstance(opt, DSROPT_RREP): opt = None
    if opt and (opt.type!=DSR_OPTION_TYPE_RREP): opt = None
    return opt

def getdsr_rerr(p):
    """Get `DSROPT_RERR` option from packet `p`."""
    opt = p
    if isinstance(p, Packet) and p.haslayer(DSRPacket):
        dsr, opt = p[DSRPacket], None
        for x in dsr.options:
            # return first option that matches
            if (opt is None) and isinstance(x, DSROPT_RERR): opt = x
    if not isinstance(opt, DSROPT_RERR): opt = None
    if opt and (opt.type!=DSR_OPTION_TYPE_RERR): opt = None
    return opt

def getdsr_ackreq(p):
    """Get `DSROPT_ACKREQ` option from packet `p`."""
    opt = p
    if isinstance(p, Packet) and p.haslayer(DSRPacket):
        dsr, opt = p[DSRPacket], None
        for x in dsr.options:
            # return first option that matches
            if (opt is None) and isinstance(x, DSROPT_ACKREQ): opt = x
    if not isinstance(opt, DSROPT_ACKREQ): opt = None
    if opt and (opt.type!=DSR_OPTION_TYPE_ACKREQ): opt = None
    return opt

def getdsr_ack(p):
    """Get `DSROPT_ACK` option from packet `p`."""
    opt = p
    if isinstance(p, Packet) and p.haslayer(DSRPacket):
        dsr, opt = p[DSRPacket], None
        for x in dsr.options:
            # return first option that matches
            if (opt is None) and isinstance(x, DSROPT_ACK): opt = x
    if not isinstance(opt, DSROPT_ACK): opt = None
    if opt and (opt.type!=DSR_OPTION_TYPE_ACK): opt = None
    return opt

def getdsr_srcroute(p):
    """Get `DSROPT_SRCROUTE` option from packet `p`."""
    opt = p
    if isinstance(p, Packet) and p.haslayer(DSRPacket):
        dsr, opt = p[DSRPacket], None
        for x in dsr.options:
            # return first option that matches
            if (opt is None) and isinstance(x, DSROPT_SRCROUTE): opt = x
    if not isinstance(opt, DSROPT_SRCROUTE): opt = None
    if opt and (opt.type!=DSR_OPTION_TYPE_SRCROUTE): opt = None
    return opt

def hasdsr_rreq(p):
    """Is `p` a DSR Route Request?"""
    opt = p
    if isinstance(p, Packet) and p.haslayer(DSRPacket):
        dsr = p[DSRPacket]
        return any([isinstance(x, DSROPT_RREQ) for x in dsr.options])
    if not isinstance(opt, DSROPT_RREQ): opt = None
    return opt and (opt.type!=DSR_OPTION_TYPE_RREQ)

def hasdsr_rrep(p):
    """Is `p` a DSR Route Reply?"""
    opt = p
    if isinstance(p, Packet) and p.haslayer(DSRPacket):
        dsr = p[DSRPacket]
        return any([isinstance(x, DSROPT_RREP) for x in dsr.options])
    if not isinstance(opt, DSROPT_RREP): opt = None
    return opt and (opt.type!=DSR_OPTION_TYPE_RREP)

def hasdsr_rerr(p):
    """Is `p` a DSR Route Error?"""
    opt = p
    if isinstance(p, Packet) and p.haslayer(DSRPacket):
        dsr = p[DSRPacket]
        return any([isinstance(x, DSROPT_RERR) for x in dsr.options])
    if not isinstance(opt, DSROPT_RERR): opt = None
    return opt and (opt.type!=DSR_OPTION_TYPE_RERR)

def hasdsr_ackreq(p):
    """Is `p` a DSR Acknowledgement Request?"""
    opt = p
    if isinstance(p, Packet) and p.haslayer(DSRPacket):
        dsr = p[DSRPacket]
        return any([isinstance(x, DSROPT_ACKREQ) for x in dsr.options])
    if not isinstance(opt, DSROPT_ACKREQ): opt = None
    return opt and (opt.type!=DSR_OPTION_TYPE_ACKREQ)

def hasdsr_ack(p):
    """Is `p` a DSR Acknowledgement?"""
    opt = p
    if isinstance(p, Packet) and p.haslayer(DSRPacket):
        dsr = p[DSRPacket]
        return any([isinstance(x, DSROPT_ACK) for x in dsr.options])
    if not isinstance(opt, DSROPT_ACK): opt = None
    return opt and (opt.type!=DSR_OPTION_TYPE_ACK)

def hasdsr_srcroute(p):
    """Is `p` a DSR Source Route?"""
    opt = p
    if isinstance(p, Packet) and p.haslayer(DSRPacket):
        dsr = p[DSRPacket]
        return any([isinstance(x, DSROPT_SRCROUTE) for x in dsr.options])
    if not isinstance(opt, DSROPT_SRCROUTE): opt = None
    return opt and (opt.type!=DSR_OPTION_TYPE_SRCROUTE)
