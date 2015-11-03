#!  /usr/bin/env python

"""
Package to implement channel.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-10-04 14:07:37 -0500 (Tue, 04 Oct 2011) $
* $LastChangedRevision: 5181 $

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

from interface import ChannelInterface, strcollision
from interface import CHANNELIF_RX, CHANNELIF_TX
from channelbase import Channel, ChannelModel
from propagation import Propagation
from refprop import ReferencePropagation
from breakpoint import Breakpoint
from radio import Radio

from fixdelay import FixedDelay

