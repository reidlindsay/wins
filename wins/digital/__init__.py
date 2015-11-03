#!  /usr/bin/env python

"""
Package for digital communications algorithm and models.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2010-03-03 11:05:36 -0600 (Wed, 03 Mar 2010) $
* $LastChangedRevision: 4425 $

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

from modulation import Modulation
from mqam   import MQAM
from rcpc   import RCPC
from reedsolomon import ReedSolomon
