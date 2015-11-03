#!  /usr/bin/env python

"""
Package implements extension modules for WiNS. Includes all SWIG-wrapped objects
and packages for interacting with Hydra PHY.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-06-17 14:29:04 -0500 (Fri, 17 Jun 2011) $
* $LastChangedRevision: 5004 $

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

from crc import CRC
from itpp import *
from waveform import Waveform
#from channel import MimoChannel
from dot11n_backend import Dot11N_Transmitter, Dot11N_Receiver, \
                           Dot11N_Channel, MimoChannel, RNG_init

# import constants for TGn channel models
from dot11n_backend import DOT11N_LOS_MODEL, \
                           DOT11N_TGN_MODEL_A, DOT11N_TGN_MODEL_B, \
                           DOT11N_TGN_MODEL_C, DOT11N_TGN_MODEL_D, \
                           DOT11N_TGN_MODEL_E, DOT11N_TGN_MODEL_F

# import flags for configuring TGn channel models
from dot11n_backend import DOT11N_USE_LOS, \
                           DOT11N_USE_NLOS, \
                           DOT11N_LOS_RICE, \
                           DOT11N_NLOS_FADING, \
                           DOT11N_NLOS_EQUALGAIN, \
                           DOT11N_NLOS_DOPPLER, \
                           DOT11N_TGN_DEFAULT
