#!  /usr/bin/env python

"""
**WiNS**: Wi-reless N-etwork S-imulator.

The `wins` package is a network simulator built on top of the event-driven
simulation environment of SimPy_. The emphasis of this design is on flexibly
modeling the wireless channel and the network, MAC, and PHY layers. The
simulator will be used to investigate cooperative communication and cross-layer
design.

Dependencies:

    * Python_ >= 2.4.4 (not version 3.0)
    * SimPy_ >= 1.8
    * Numpy_ >= 1.0.4
    * NetworkX_ >= 0.99

.. _Python: http://www.python.org
.. _SimPy : http://simpy.sourceforge.net/
.. _Numpy : http://numpy.scipy.org
.. _NetworkX : http://networkx.lanl.gov/

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-10-08 10:39:45 -0500 (Sat, 08 Oct 2011) $
* $LastChangedRevision: 5189 $

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
__version__   = 0.1

# import from main module
import const
from helper import *
from base   import Base, Reference
from trace  import Trace, Traceable
from fsm    import FSM, Timer
from packet import Packet, ANNO, HiddenField, FuncField
from queue  import Queue
from element import Port, Element
from crc    import CRC32, crcupdate, crcremove

# import other modules
from util import *
from mobile import *
from channel import *
from digital import *
from protocol import *
import ieee80211
import mac
import net
from traffic import *
