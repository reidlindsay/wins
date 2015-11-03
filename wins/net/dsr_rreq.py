#!  /usr/bin/env python

"""
Implementation of Route Request Table for DSR; contains `RouteRequestTable` and
`RouteRequestEntry` classes.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-12-18 16:54:42 -0600 (Sun, 18 Dec 2011) $
* $LastChangedRevision: 5380 $

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

from SimPy.Simulation import now
from wins.trace import Traceable

from numpy import log2

class RouteRequestTable(Traceable):
    """DSR Route Request Table.

    Manage route requests sent and forwarded by the network protocol. All FIFOs
    are implemented as lists where new objects are appended to the back and old
    objects are popped from the front.

    :ivar cache: Maintains information for RREQ overheard by DSR. Maps source
                 addresses to (ID, target) pairs.
    :ivar sendbuffer: Maintain information for route discovery initiated by the
                      parent `DSR`. Maps target address to `RouteRequestEntry`.

    :note: `RouteRequestTable` *must* be a child of a `DSR` module.
    """
    name = "route request table"
    tracename = "RREQTBL"
    MaxRequestRexmt = 10
    RequestTableIds = 16
    RequestTableSize = 64
    def __init__(self, **kwargs):
        """Constructor."""
        Traceable.__init__(self, **kwargs)
        self.cache = {}
        self.sendbuffer = {}

    def addcache(self, src, ID, target):
        """Add new entry to `cache`.
        
        :param src: Source address that initiated route discovery.
        :param ID: Identification field of RREQ.
        :param target: Target of RREQ.
        :return: Boolean; if true, indicates new entry into cache, otherwise
                 entry already existed (no new information).
        """
        # create entry in cache (if needed)
        if (src not in self.cache): self.cache[src] = []
        # check if (ID, target) is in cache
        key = (ID, target)
        isnew = (key not in self.cache[src])
        if isnew: self.cache[src].append(key)    # FIFO (newest in back)
        # trim old data - pop oldest from front
        while (len(self.cache[src])>self.RequestTableIds): self.cache[src].pop(0)
        return isnew

    def sendrreq(self, target):
        """Attempt to send route request based on table data.

        :param target: Target of RREQ.
        :return: `RouteRequestEntry` for target (or None if attempt fails).

        This method limits the rate at which route discovery occurs. It also
        assumes the parent `DSR` is initiating the route discovery.
        """
        rre = self.getentry(target)
        # check if route discovery can be initiated
        if rre.isalive(): rre = None
        else:             rre.sendrreq()
        return rre

    def getentry(self, target):
        """Return entry into `sendbuffer` corresponding to `target`."""
        if target in self.sendbuffer:
            rre = self.sendbuffer[target]
        else:
            rre = self.sendbuffer[target] = RouteRequestEntry()
        return rre

    def delentry(self, target):
        """Remove entry from `sendbuffer` corresponding to `target`."""
        rre = None
        if target in self.sendbuffer:
            rre = self.sendbuffer[target]
            del self.sendbuffer[target]
        return rre

class RouteRequestEntry:
    """DSR Route Request Entry used by `RouteRequestTable`.

    Entry in routing table corresponding to a particular target.
    """
    MaxID = 4096
    RequestPeriod = 500e-3
    MaxRequestPeriod = 10
    NonpropRequestTimeout = 30e-3       # initial timeout for Route Discovery
    DiscoveryHopLimit = 255
    def __init__(self, ttl=None, ts=None, ID=None):
        """Constructor."""
        self.sendbuffer = []
        self.ttl = ttl
        self.ts  = ts
        self.ID  = ID

    def isalive(self):
        """See if previously initiated route discovery is still alive."""
        return (self.timeleft()>0)

    def timeleft(self):
        """Compute time left before next route discovery can be initiated."""
        tleft = 0
        if self.ttl is None: return tleft
        if self.ts is None: return tleft
        t0 = self.NonpropRequestTimeout
        #if (self.ttl==1): t0 = self.NonpropRequestTimeout
        #else:             t0 = self.RequestPeriod
        maxtimeout = self.MaxRequestPeriod
        timeout = min(self.ttl*t0, maxtimeout)
        tpassed = now() - self.ts
        tleft = max(timeout - tpassed, 0)
        return tleft

    def numxmit(self):
        """Number of transmission attempts on to this target."""
        if self.ttl is None: return 0
        nxmit = int(log2(self.ttl))+1
        return nxmit

    def sendrreq(self):
        """Update parameters for next RREQ."""
        ttl, ts, ID = self.ttl, self.ts, self.ID
        maxttl = self.DiscoveryHopLimit
        if ttl is None:
            # first RREQ
            ttl, ts = 1, now()
            ID = (id(self)+int(ts*1e3))%self.MaxID  # random ID seed
        else:
            # not first RREQ
            ttl, ts = ttl*2, now()
            ID = (ID+1)%(self.MaxID)                # iterate ID
        # set parameters
        self.ttl = min(ttl, maxttl)
        self.ts, self.ID = ts, ID

    def buffer(self, p):
        """Insert packet into `sendbuffer`.

        :param p: Packet to insert into buffer.
        :return: List of packets that were dropped by the buffer.
        """
        self.sendbuffer.append(p)       # treat as FIFO (newest in back)
        # trim buffer as needed
        drop = []
        while (len(self.sendbuffer)>RouteRequestTable.RequestTableSize):
            d = self.sendbuffer.pop(0)      # pop oldest from front
            drop.append(d)
        return drop
