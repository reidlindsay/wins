#!  /usr/bin/env python

"""
Base class for maintaining route table; contains `RouteTable` class.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-12-19 12:26:09 -0600 (Mon, 19 Dec 2011) $
* $LastChangedRevision: 5386 $

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

:var ROUTE_TABLE_VERBOSE:
    Constant enumeration to control verbose thresholds in this file.
  
    Setting the verbose level of a `Routing` above this threshold will cause the
    corresponding output in this file to be written (or logged).
"""
__docformat__ = "restructuredtext en"

from SimPy.Simulation import now
from wins.trace  import Traceable
from wins.base   import Reference

ROUTE_TABLE_VERBOSE = 31

class RouteTable(Traceable):
    """Table to maintain routing information.

    Each routing entry consists of a destination, route, and cost. This class
    supports multiple paths to the same destination. It uses cost to order the
    paths so that `nexthop` and `route` will select the lowest cost path.

    :ivar data: Data structure used to maintain routing information; use
                accessor/modifier methods `addroute()`, `hasroute()`, etc.
    :ivar maxweight: If true, this module will use the maximum weight path.
    :ivar quiet: If true, operate in quiet mode (do not log events).
    :ivar timeout: Timeout used to invalidate stale entries.
    """
    name = "routing table"
    tracename = "RT"
    def __init__(self, maxweight=False, quiet=False, timeout=None, **kwargs):
        """Constructor.

        :param maxweight: If true, this module will use the maximum weight path.
        :param quiet: If true, operate in quiet mode (do not log events).
        """
        Traceable.__init__(self, **kwargs)
        self.data = {}
        self.maxweight = maxweight
        self.quiet = quiet
        self.timeout = timeout

    def addroute(self, dest, nexthop=None, cost=None):
        """Add new route to `dest`.

        :param dest: Address of destination.
        :param nexthop: Address of nexthop [default=`None`].
        :param cost: Cost of `route` [default=None].
        :return: Number of paths to `dest` in routing table.

        If `nexthop` is `None`, it is assumed that `dest` is a one hop neighbor.
        If cost is None, `RouteTable()` uses `cost()` to determine the cost of
        the route through `nexthop`.
        """
        if isinstance(dest, Reference): dest = dest._deref
        if isinstance(nexthop, Reference): nexthop = nexthop._deref
        idx, ts = None, now()  # index, timestamp
        if dest in self.data:
            # add path to list
            nhlist = [h for (c,t,h) in self.data[dest]['list']]
            if nexthop in nhlist:
                idx = nhlist.index(nexthop)
            # new route?
            if (idx is None):
                self.data[dest]['list'].append((cost, ts, nexthop) )
            # overwrite existing route?
            else:
                self.data[dest]['list'][idx] = (cost, ts, nexthop)
            # reset index to indicate entry was modified
            self.data[dest]['index'] = None
        else:
            # add first path to dest
            self.data[dest] = {'list':[(cost, ts, nexthop)], 'index': 0}
        # log action
        if not self.quiet:
            if (idx is None):  self.log("add", dest=dest, nexthop=nexthop, cost=cost)
            else:              self.log("replace", dest=dest, nexthop=nexthop, cost=cost)
        return len(self.data[dest])

    def hasroute(self, dest):
        """Check if route to `dest` exists.

        :return: Boolean; true if `dest` is found, false otherwise.
        """
        self._removestale(dest)
        if isinstance(dest, Reference): dest = dest._deref
        return (dest in self.data) and (len(self.data[dest]['list'])>0)

    def _removestale(self, dest):
        """Internal method to remove stale routing data from table.

        This method is called from `hasroute`.
        """
        if not (self.timeout>0): return   # do nothing
        # Remove all stale routes to dest
        if isinstance(dest, Reference): dest = dest._deref
        if dest in self.data:
            tnow, stale = now(), []
            ndata = len(self.data[dest]['list'])
            # find stale routes
            for k in range(ndata):
                c, ts, nh = self.data[dest]['list'][k]
                if ((tnow-ts)>self.timeout): stale.append(k)
            # log stale routes to be removed
            if not self.quiet:
                for s in stale:
                    nexthop = self.data[dest]['list'][s]
                    self.log("stale", dest=dest, nexthop=nexthop)
            # remove stale routes
            if stale:
                keep = [k for k in range(ndata) if (k not in stale)]
                newdata = [self.data[dest]['list'][k] for k in keep]
                self.data[dest]['list'] = newdata
                self.data[dest]['index'] = None     # reset index
            # remove any empty routing entries
            if (len(self.data[dest]['list'])<1):
                del self.data[dest]                 # remove empty entry

    def delroute(self, dest, nexthop="auto", reset=False):
        """Remove least desirable route to `dest`.

        :param dest: Destination to modify.
        :param nexthop: Specify which route should be removed based on nexthop
                        value [optional].
        :param reset: If true, remove all routes associated with `dest`.
        """
        if isinstance(dest, Reference): dest = dest._deref
        hasroute = self.hasroute(dest)
        if hasroute and reset:
            # remove all
            del self.data[dest]
        elif hasroute and (nexthop=="auto"):
            # remove least desirable
            idx = self._index(dest, not self.maxweight) # least desirable
            self.data[dest]['list'].pop(idx)
            self.data[dest]['index'] = None
        elif hasroute:
            # remove nexthop
            idx = None
            for k in range(len(self.data[dest]['list'])):
                c, ts, nh = self.data[dest]['list'][k]
                if (nh==nexthop): idx = k
            if idx is not None:
                self.data[dest]['list'].pop(idx)
                self.data[dest]['index'] = None     # reset index
        if not self.quiet:
            self.log("del", dest=dest, nexthop=nexthop, reset=reset)

    def nexthop(self, dest):
        """Find address of next hop to destination `dest`.

        :return: Address of next hop, or `None` if no route to `dest` is found.
        """
        if isinstance(dest, Reference): dest = dest._deref
        nh = None
        if self.hasroute(dest):
            idx = self._currindex(dest)
            c, ts, nh = self.data[dest]['list'][idx]
            if nh is None: nh = dest
        return nh

    def cost(self, dest):
        """Find cost of route to destination `dest`.

        :return: Cost of route, or `None` if no route to `dest` is found.
        """
        if isinstance(dest, Reference): dest = dest._deref
        c = None
        if self.hasroute(dest):
            idx = self._currindex(dest)
            c, ts, nh = self.data[dest]['list'][idx]
            if c is None: c = self.getcost(dest, nh, ts)
        return c

    def timestamp(self, dest):
        """Find timestamp for route to destination `dest`."""
        if isinstance(dest, Reference): dest = dest._deref
        ts = None
        c = None
        if self.hasroute(dest):
            idx = self._currindex(dest)
            c, ts, nh = self.data[dest]['list'][idx]
        return ts

    def getcost(self, dest, nexthop, ts):
        """Determine cost of path to `dest` using route through `nexthop`.

        :param dest: Destination of path.
        :param nexthop: Next hop in path.
        :param ts: Timestamp for when route was added to table.

        By default this method calculates the cost using the timestamp `ts`.
        That is, this method returns the "age" of the path.

        :note: Overload this method to change how cost of a path is computed.
        """
        cost = now() - ts
        return cost

    def _index(self, dest, maxcost=True, costfunc=None):
        """Internal method to find index of max [or min] cost path.

        :param maxcost: If true, return index of path with maximum cost;
                        otherwise return index of path with minimum cost.
        :param costfunc: Cost function used to automatically determine cost when
                         needed [default=`getcost`].
        :return: Index into `data` entry corresponding to `dest`.
        """
        if isinstance(dest, Reference): dest = dest._deref
        idx = None
        if costfunc is None: costfunc = self.getcost
        if dest in self.data:
            # sort routes to dest based on cost
            olist = [x for x in self.data[dest]['list'] ]
            for k in range(len(olist)):
                # calculate cost if needed
                c, ts, nh = olist[k]
                if (c is None): c = costfunc(dest, nh, ts)
                olist[k] = (c, ts, nh, k)
            olist.sort()
            # get appropriate index
            if maxcost: idx = olist[-1][-1]  # max cost path
            else:       idx = olist[0][-1]   # min cost path
        return idx

    def _currindex(self, dest):
        """Internal method to get/update current index for path to `dest`."""
        idx = None
        if dest in self.data:
            idx = self.data[dest]['index']
            if idx is None: idx = self._index(dest, self.maxweight)
            # set index in entry
            self.data[dest]['index'] = idx
        return idx

    def log(self, evt=None, p=None, *args, **kwargs):
        """Overloaded to check verbose level and set common annotations."""
        force = False
        if ('verbose' in kwargs): force = (kwargs['verbose']>ROUTE_TABLE_VERBOSE)
        if self.verbose>ROUTE_TABLE_VERBOSE or force:
            Traceable.log(self, evt, p, *args, **kwargs)
