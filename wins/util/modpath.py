#!  /usr/bin/env python

"""
Modify current PATH.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-03-24 11:29:29 -0500 (Thu, 24 Mar 2011) $
* $LastChangedRevision: 4942 $

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

from os.path import splitext, basename, join as pjoin, walk
import os
import sys

def addpath(relpath):
    """Add relative path to current PATH variable.

    :param relpath: Relative path to add to PATH.
    :returns: Added path.
    """
    assert isinstance(relpath, str), "[addpath]: new path must be string!"
    d = os.getcwd()
    rpath = pjoin(d, relpath)
    sys.path.insert(0, rpath)
    return rpath
