#!  /usr/bin/env python

"""
Package to implement base classes for protocols.

Most protocols should include four ports, namely:

    1. "TXU" - sends received packets to upstream element.
    #. "RXU" - receives traffic from an upstream element.
    #. "TXD" - sends traffic to a downstream element.
    #. "RXD" - receives traffic from a downstream element.

Maintaining this convention will make accessing `Port` configurations for
various protocols more straight forward. Some protocols may strictly enforce
all or part of this convention.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-10-15 01:32:10 -0500 (Sat, 15 Oct 2011) $
* $LastChangedRevision: 5196 $

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

from phy import PHY
from mac import MAC
from csphy import CSPHY, CSPHY_BUSY, CSPHY_IDLE
from csmac import CSMAC
from net import NET
from route_support import *
from routetable import RouteTable, ROUTE_TABLE_VERBOSE
from routing    import Routing, ROUTING_VERBOSE
from aloha import Aloha
from arp import ARP
from arp_support import *
