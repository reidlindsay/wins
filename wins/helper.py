#!  /usr/bin/env python

"""
Helper classes and methods for `wins`.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-09-19 04:45:05 -0500 (Mon, 19 Sep 2011) $
* $LastChangedRevision: 5132 $

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

from SimPy.Simulation import Process, SimEvent, waitevent, activate
from wins import const

def monitor_events(*args, **kwargs):
    """Monitor event(s) to keep them up to date.

    :param args: All arguments are SimEvents to wait on.
    :param kwargs: All keywords are passed to `activate()` Process that monitors
                   events in `args`.
    :return: Monitoring Process.
    """
    proc = Process(name="monitor events")
    activate(proc, _monitor_events(proc, *args), **kwargs)
    return proc

def _monitor_events(proc, *args):
    """Internal process execution method (PEM) to monitor event(s)."""
    assert isinstance(proc, Process), \
           "[monitor_events]: PEM must have Process as first argument!"
    assert all([isinstance(a,SimEvent) for a in args]), \
           "[monitor_events]: PEM arguments must be SimEvents!"
    while proc.active():
        if len(args)<2: yield waitevent, proc, args[0]
        else:           yield waitevent, proc, args

def isnum(x):
    """Check if an object is a float or int."""
    num = False
    num = num or isinstance(x, int)
    num = num or isinstance(x, float)
    return num

from scipy.misc import comb

def nchoosek(*args):
    """Wrapper for combinatorial function (i.e. n choose k)."""
    return int(comb(*args))

import numpy

def linear2db(xlinear):
    """Convert a linear value to decibels.

    :param xlinear: Value(s) in linear scale.

    :returns: 10*log10(`xlinear`)
    """
    yval = []
    xval = xlinear
    isiter = hasattr(xlinear, '__iter__')
    if not isiter: xval = [xlinear]
    for x in xval:
        if   (x>0): y = 10*numpy.log10(x)
        elif (x<0): y = numpy.nan
        else: y = -numpy.inf
        yval.append(y)
    y = numpy.array(yval)
    if not isiter: y = yval[0]
    return y

def db2linear(xdb):
    """Convert a decibel value to a linear scale.

    :param xdb: Value(s) in decibels.

    :returns: 10^(`xdb`/10)
    """
    return 10.0**(xdb/10.0)

def time2usec(t, fmt=None):
    """Convert a float time into a readable time in microseconds.

    :param t: Time value.
    :param fmt: Optional format string [default="%.4f"].
    :return: Formatted string representation of time in microseconds.
    """
    if fmt is None: fmt= "%.2f"
    exec("s = \"%s usec\"%%(t*1e6)"%(fmt) )
    return s

def time2msec(t, fmt=None):
    """Convert a float time into a readable time in milliseconds.

    :param t: Time value.
    :param fmt: Optional format string [default="%.4f"].
    :return: Formatted string representation of time in milliseconds.
    """
    if fmt is None: fmt= "%.3f"
    exec("s = \"%s msec\"%%(t*1e3)"%(fmt) )
    return s

from scapy.all import IP, ICMP, TCP, UDP, ByteField, bind_layers
import struct

class PseudoIPv4(IP):
    """Pseudo header that mimics IPv4 header for computing checksum in transport
    layer protocols (i.e. UDP and TCP)."""
    name = "PseudoIPv4"
    def post_build(self, p, pay):
        p += pay
        l = self.len
        if l is None:
            l = len(p)
            p = p[:10] + struct.pack("!H", l) + p[12:]
        return p+pay
    def extract_padding(self, s):
        """Seperate header and payload."""
        l = self.len - 12
        return s[:l], s[l:]
    def mysummary(self):
        """Create readable summary string."""
        s = self.sprintf("%PseudoIPv4.src% > %PseudoIPv4.dst% %PseudoIPv4.proto%")
        return s
    @staticmethod
    def _init_fields_desc():
        """Internal method to initialize `PseudoIPv4.fields_desc`; this removes
        unnecessary fields from the header."""
        PseudoIPv4.fields_desc = [IP.fields_desc[10],       # src   \
                                  IP.fields_desc[11],       # dst   \
                                  ByteField("zero", 0),     # zero  \
                                  IP.fields_desc[8],        # proto \
                                  IP.fields_desc[3] ]       # len
        bind_layers(PseudoIPv4, IP,   proto=const.IP_PROTO_IP)
        bind_layers(PseudoIPv4, ICMP, proto=const.IP_PROTO_ICMP)
        bind_layers(PseudoIPv4, TCP,  proto=const.IP_PROTO_TCP)
        bind_layers(PseudoIPv4, UDP,  proto=const.IP_PROTO_UDP)

PseudoIPv4._init_fields_desc()
