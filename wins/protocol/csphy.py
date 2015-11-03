#!  /usr/bin/env python

"""
Base class for implementing physical layer with carrier sense support which is
needed for implementing `CSMAC` MAC protocols; contains `CSPHY` class.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-09-27 22:15:57 -0500 (Tue, 27 Sep 2011) $
* $LastChangedRevision: 5167 $

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

:var CSPHY_BUSY: Enumeration for carrier sense busy state in `CSPHY`.

:var CSPHY_IDLE: Enumeration for carrier sense idle state in `CSPHY`.

:var CSPHY_VERBOSE:
    Constant enumeration to control verbose thresholds in this file.
  
    Setting the verbose level of a `CSPHY` above this threshold will cause the
    corresponding output in this file to be written (or logged).
"""
__docformat__ = "restructuredtext en"

from SimPy.Simulation import SimEvent, now
from wins.protocol.phy import PHY
from wins.helper import monitor_events

CSPHY_BUSY = "csbusy"
CSPHY_IDLE = "csidle"
CSPHY_VERBOSE = 60

class CSPHY(PHY):
    """Base class for implementing physical layer with carrier sense support,
    which is needed by `CSMAC` MAC protocols.

    The functionality of this class, which should be implemented by subclasses,
    includes the following:

        * `set_csbusy()` when the medium is "busy" (e.g. the received power
          level exceeds some threshold, a packet is detected, etc.)
        * `set_csidle()` when the medium becomes "idle" (e.g. there are no more
          packets being received, the received power falls below some threshold
          for some time, etc.)

    :CVariables:
     * `csmode`: Property to access current carrier sense state. Use
       `set_csbusy()`/`set_csidle()` to modify.
     * `isbusy`: Property to determine if current `csmode` is `CSPHY_BUSY`.
     * `isidle`: Property to determine if current `csmode` is `CSPHY_IDLE`.

    :IVariables:
     * `csbusy`: SimEvent signalled if carrier sense state transitions to busy.
     * `csidle`: SimEvent signalled if carrier sense state transitions to idle.
     * `__csmode`: Private variable to maintain carrier sense state; initialized
       to `CSPHY_IDLE`. Use `set_csbusy()`/`set_csidle()` to modify.

    :note: This element has `csmode` initialized to `CSPHY_IDLE`.
    """
    name = "carrier sense PHY"
    tracename = "CSPHY"
    def __init__(self, **kwargs):
        """Constructor."""
        self.__csmode = CSPHY_IDLE
        self.csbusy = SimEvent(name="csbusy")
        self.csidle = SimEvent(name="csidle")
        PHY.__init__(self, **kwargs)
        # rename events
        self.csbusy.name = "%s(%s).csbusy"%(self.name, self.uid)
        self.csidle.name = "%s(%s).csidle"%(self.name, self.uid)
        # monitor events -> keep up to date
        monitor_events(self.csbusy, self.csidle)

    csmode = property(fget=lambda self: self.__csmode)
    isbusy = property(fget=lambda self: (self.csmode == CSPHY_BUSY) )
    isidle = property(fget=lambda self: (self.csmode == CSPHY_IDLE) )

    def set_csbusy(self, *args, **kwargs):
        """Transition to busy state.

        :param args: Additional arguments passed when signalling `csbusy`.
        :param kwargs: Additional keywords used with `log()`.

        The SimEvent `csbusy` will only be signalled if there is a `state`
        transition from idle to busy.
        """
        ostate = self.csmode
        self.__csmode = CSPHY_BUSY
        if (ostate != self.csmode):
            # IDLE -> BUSY!
            if self.verbose>CSPHY_VERBOSE: self.log(self.csmode,*args,**kwargs)
            self.csbusy.signal(*args)

    def set_csidle(self, *args, **kwargs):
        """Transition to idle state.

        :param args: Additional arguments passed when signalling `csidle`.
        :param kwargs: Additional keywords used with `log()`.

        The SimEvent `csidle` will only be signalled if there is a `csmode`
        transition from busy to idle.
        """
        ostate = self.csmode
        self.__csmode = CSPHY_IDLE
        if (ostate != self.csmode):
            # BUSY -> IDLE!
            if self.verbose>CSPHY_VERBOSE: self.log(self.csmode,*args,**kwargs)
            self.csidle.signal(*args)
