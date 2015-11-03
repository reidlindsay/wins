#!  /usr/bin/env python

"""
Implementation of Route Cache for `DSR`; contains `RouteCache` class.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-12-19 12:26:09 -0600 (Mon, 19 Dec 2011) $
* $LastChangedRevision: 5386 $

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

RCACHEDEBUG = 0

from SimPy.Simulation import now

from wins.base import Reference
from wins.protocol.routetable import RouteTable

from wins.net.dsr_support import *

import sys

class RouteCache(RouteTable):
    """Routing table/cache for source routes in DSR.

    Entries in table are lists of addresses (i.e. source route), instead of a
    single value (or destination). Each route consists is a list of the
    intermediate nodes that must be visited to get to the destination. An empty
    list implies the destination is a neighbor.
    """
    name = "source route cache"
    tracename = "RCACHE"
    DefaultTimeout = None
    def __init__(self, *args, **kwargs):
        """Constructor.

        :param addr: Address for parent network interface [default=None].
        :param args: Additional arguments passed to `RouteTable` constructor.
        :param kwargs: Additional keywords passed to `RouteTable` constructor.
        """
        RouteTable.__init__(self, *args, **kwargs)
        # set default parameters
        if self.timeout is None:
            self.timeout = self.DefaultTimeout

    def hasroute(self, dest, path="ignore"):
        """Overloaded to allow checking for a specific path."""
        r = RouteTable.hasroute(self, dest)
        if r and (path!="ignore"):
            nhlist = [h for (c,t,h) in self.data[dest]['list']]
            r = (path in nhlist)
        return r

    def addroute(self, dest, nexthop=None, cost=None):
        """Cache intermediate destinations."""
        rval = RouteTable.addroute(self, dest, nexthop, cost)
        if DSR_CACHE_ALL_HOPS and nexthop:
            # check if nexthop is valid source route
            isiter = True
            try:    x = iter(nexthop)
            except: isiter = False
            isiter &= not isinstance(nexthop, str)
            errmsg = "[ROUTECACHE]: Invalid source route (%s)"%(nexthop)
            assert (isiter), errmsg
            # cache intermediate hops
            for k in range(len(nexthop)):
                dst = nexthop[k]    # intermediate node
                path = nexthop[:k]  # truncated path to intermediate node
                if not path: path = None
                RouteTable.addroute(self, dst, path, None)
        return rval

    def nexthop(self, dest):
        """Find address of first hop in source route to `dest`."""
        nexthop = None
        if self.hasroute(dest):
            path = self.srcroute(dest)
            nexthop = (path+[dest])[0]
            self.debug("addr = %s\n"%(self.parent.address))
            self.debug("dest = %s\n"%(dest))
            self.debug("path = %s\n"%(path))
            self.debug("nexthop =%s\n"%(nexthop))
            self.debug("\n")
        return nexthop

    def srcroute(self, dest):
        """Find source route to `dest`."""
        if isinstance(dest, Reference): dest = dest._deref
        path = None
        if self.hasroute(dest):
            idx = self._currindex(dest)
            c, ts, path = self.data[dest]['list'][idx]
            if path is None: path = []  # empty source route
        return path

    def getcost(self, dest, nexthop, ts):
        """Overloaded to compute cost according to route length."""
        cost = 1
        if nexthop: cost += len(nexthop)
        return cost

    def removelink(self, a, b, addr=None):
        """Remove all routes that traverse the link between `a` and `b`."""
        rlist = []
        for dst in self.data:
            for path in [nh for c,ts,nh in self.data[dst]['list']]:
                remove = False
                if (path is None):
                    remove = (a==addr) and (b==dst)      # one hop neighbor
                elif (len(path)==0):
                    remove = (a==addr) and (b==dst)      # invalid path
                elif (a in path) and (b==dst):
                    remove = (a==path[-1])               # last hop of route
                elif (a in path) and (b in path):
                    remove = ((path.index(b)-path.index(a))==1) # middle of route
                elif (b==path[0]):
                    remove = (a==addr)                   # start of route
                # add to remove list
                if remove:
                    rlist.append((dst, path))
        # clean up routes
        for d,p in rlist: self.delroute(d, p)
        # return list of affected destinations
        droplist = []
        for d,p in rlist:
            if (d not in droplist): droplist.append(d)
        return droplist

    def debug(self, s, stderr=True):
        """Internal method to print debug statements."""
        fout = sys.stderr
        if not stderr: fout = sys.stdout
        if (RCACHEDEBUG>0): fout.write(s)

