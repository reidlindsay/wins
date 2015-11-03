#!  /usr/bin/env python

"""
Packet definitions, enumerations, and helper functions for `ARF` protocol.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-09-13 15:45:28 -0500 (Tue, 13 Sep 2011) $
* $LastChangedRevision: 5125 $

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

:var ARF_TIMEOUT:
    Timeout value used by `ARF` timer to determine when to prematurely increase
    the data rate of the transmission.

:var ARF_NSUCCESS:
    Number of consecutive successful ACKs used that must be received to increase
    the data rate of the transmission.

:var ARF_NFAILURE:
    Number of failed (or unacknowledged) ACKs that will trigger `ARF` to use a
    lower data rate for the next transmission.
"""
__docformat__ = "restructuredtext en"

ARF_TIMEOUT = 15
ARF_NSUCCESS = 10
ARF_NFAILURE = 2
