#!  /usr/bin/env python

"""
Packet definitions, enumerations, and helper functions for `Routing` protocol.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-09-17 12:21:22 -0500 (Sat, 17 Sep 2011) $
* $LastChangedRevision: 5129 $

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
"""
__docformat__ = "restructuredtext en"

from scapy.all import IP, IP_PROTOS
from wins import const

def _update_IP_with_HOST():
    """Internal method to update `IP` for unknown internal host protocols."""
    proto = None
    for f in IP.fields_desc:
        if (f.name=='proto'): proto = f
    assert (proto is not None)
    proto.i2s[const.IP_PROTO_HOST] = 'host'

_update_IP_with_HOST()
