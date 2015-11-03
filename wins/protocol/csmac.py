#!  /usr/bin/env python

"""
Base class for implementing carrier sense MAC control protocols such as CSMA/CA;
contains `CSMAC` class.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2010-10-19 16:34:02 -0500 (Tue, 19 Oct 2010) $
* $LastChangedRevision: 4811 $

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

from wins.element import Element
from wins.protocol.mac import MAC
from wins.protocol.csphy import CSPHY

class CSMAC(MAC):
    """Base class for implementing carrier sense MAC protocols such as CSMA/CA.

    This base class contains members to interact with carrier sense properties
    in a `CSPHY`. The actual implementation of this functionality is left to
    subclasses.

    :CVariables:
     * `csmode`: Property to access current carrier sense state of `phy`.
     * `isbusy`: Property to determine if current `csmode` is `CSPHY_BUSY`.
     * `isidle`: Property to determine if current `csmode` is `CSPHY_IDLE`.
     * `csbusy`: Property to access `CSPHY.csbusy` of `phy`.
     * `csidle`: Property to access `CSPHY.csidle` of `phy`.

    :note: In order to implement a `CSMAC`, a valid `CSPHY` must be set as its
           `phy` prior to starting. The `set_phy()` method has been overloaded
           to check that this condition has been satisfied.
    """
    name = "carrier sense MAC"
    tracename = "CSMAC"
    def __init__(self, **kwargs):
        """Constructor."""
        MAC.__init__(self, **kwargs)

    csmode = property(fget=lambda self: self.phy.csmode)
    isbusy = property(fget=lambda self: self.phy.isbusy)
    isidle = property(fget=lambda self: self.phy.isidle)
    csbusy = property(fget=lambda self: self.phy.csbusy)
    csidle = property(fget=lambda self: self.phy.csidle)

    def set_phy(self, p):
        """Check for valid `CSPHY` for MAC layer."""
        assert isinstance(p, CSPHY) or (p is None), \
               "[CSMAC]: set_phy() cannot set pointer to non-CSPHY object!"
        MAC.set_phy(self, p)
