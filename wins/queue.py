#!  /usr/bin/env python

"""
Implementation of FIFO and priority queue; contains `Queue` class.

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

:var QUEUE_LOW_VERBOSE:
    Constant enumeration to control low verbose threshold in this file.

    Setting the verbose level of a `Queue` above this threshold will cause the
    corresponding output in this file to be written (or logged).

:var QUEUE_HIGH_VERBOSE:
    Constant enumeration to control high verbose threshold in this file.
"""
__docformat__ = "restructuredtext en"

from SimPy.Simulation import Process, SimEvent, Store, Monitor, PriorityQ
from SimPy.Simulation import passivate, activate, waitevent
from SimPy.Simulation import hold, put, get, now
from wins.trace import Traceable
from wins import const

import heapq

QUEUE_LOW_VERBOSE = 50
QUEUE_HIGH_VERBOSE = 100

class _PriorityList(list):
    """Internal class for implementing `theBuffer` for `PriorityQueue`."""
    def remove(self, item):
        flatBuffer = [a for a,b in self]
        idx = flatBuffer.index(item)
        i, prio = self[idx]
        return list.remove(self, (i, prio) )
    def sort(self, cmp=None, *args, **kwargs):
        if cmp is None: cmp = self.__cmp
        return list.sort(self, cmp=cmp, *args, **kwargs)
    def __cmp(self, a, b):
        return cmp(b[1], a[1])

def _prioritySort(self, buff):
    """Internal method for keeping priority queue sorted."""
    global _PriorityList
    assert isinstance(buff, _PriorityList), "_prioritySort must get  _PriorityList!"
    buff.sort()
    return buff

class Queue(Traceable, Store):
    """Convenient API to SimPy Store.

    See the constructor `__init__()` for configuration options.

    Queues and Renege Clauses
    =========================
    This class provides an API to get/put to SimPy Stores with renege clauses
    (i.e. `remove()` and `insert()`). A valid renege argument to these methods
    must be one of the following:

        1. a timeout value,
        #. a SimEvent to wait on,
        #. a list or tuple of SimEvents to wait on,
        #. a tuple that is a valid renege clause, e.g. (hold, proc, tsec)

    The `acquired()` and `stored()` method can only be used with SimPy Processes
    when a renege clause is provided to the associated `remove()` or `insert()`.
    Whenver using an `FSM`, however, these methods are overloaded so that they
    will work regardless of the presence of a renege clause.

    :CVariables:
     * `name`: Name of Queue.
     * `tracename`: Name used in `Trace`.

       By default, `Queue` tracenames will be appended with the queue's `uid`.

     * `all`: Property to access private SimEvent that gets signalled when a
       get, put, or drop occurs.
     * `priority`: Property to access priority flag set during constructor.
     * `length`: Property to access number of elements buffered in Queue.
     * `traceid`: Overload property to return `tracename`.

    :IVariables:
     * `monitorQ`: Boolean flag; if true, signal events `enQ` and `deQ`.
     * `enQ`: SimEvent signalled when a new item is inserted in Queue.
     * `deQ`: SimEvent signalled when an item is removed from the Queue.
     * `__priorityQ`: Private boolean flag; if true, `Queue` operates in
       priority mode.
     * `__dummy`: Private SimEvent used so that `acquired()` and `stored()`
       methods can be used to check if `insert()` and `remove()` are successful.
    """
    name="queue"
    tracename="Q"
    def __init__(self, priority=False, monitorQ=False, \
                 unitName="packet", capacity="unbounded", \
                 initialBuffered=None, **kwargs):
        """Constructor.

        :param priority:  Boolean; if True, use priority queueing.
        :param monitorQ:  Boolean; if True, support `enQ`, `deQ`, and `drp`.
        :param unitName:  Description of units stored in `Queue`.
        :param capacity:  Capacity of underlying Store [default='unbounded'].
        :param initialBuffered: Initialize list of buffered objects.
        :param kwargs:    Keywords passed to `Traceable` constructors.
        """
        # check that initialBuffered is properly formatted
        if initialBuffered and priority:
            isPrio = all([isinstance(s,tuple) for s in initialBuffered] )
            if not isPrio: initialBuffered = [(p, 0) for p in initialBuffered]
        # call constructors
        Store.__init__(self, unitName=unitName, capacity=capacity, \
                       initialBuffered=initialBuffered, \
                       putQType=PriorityQ, getQType=PriorityQ)
        Traceable.__init__(self, **kwargs)
        # set other parameters
        self.tracename = self.tracename + "%d"%(self.uid)
        self.__priority = priority
        self.monitorQ = monitorQ
        self.enQ = SimEvent(name=self.name+".enQ")
        self.deQ = SimEvent(name=self.name+".deQ")
        self.drp = SimEvent(name=self.name+".drp")
        self.__dummy = SimEvent(name=self.name+".dummy")
        # set up Queue for priority
        if self.priority:
            self.theBuffer = _PriorityList(self.theBuffer)
        self.addSort(None) # setup addSort for priority queueing

    priority = property(fget=lambda self: self.__priority)
    length = property(fget=lambda self: self.nrBuffered)
    traceid = property(fget=lambda self: self.tracename)

    def insert(self, proc, S, prio=None, renege=None):
        """Insert a list of objects `S` into the `Queue`.

        :param proc:    SimPy Process that will execute insert.
        :param S:       List of objects to insert.
        :param prio:    Priority level of insert (only used when `priority` is
                        set to `True`).
        :param renege:  Renege clause or timeout.
        :return: Yield clause to block on.

        If `priority` is false, `prio` will be ignored. Otherwise, if `priority`
        is True and `prio` is not provided, then objects in `S` will be inserted
        with priority level 0.

        Also, when `priority` is True, `S` can contain 2-tuples of *(obj, prio)*
        that would allow each object in `S` to be inserted with its own priority.

        Normal usage with a `renege` clause is as follows:

                >>> yield queue.insert(proc, S, renege=<renege clause>)
                >>> if proc.stored(queue):
                ...     # do something
                ... else:
                ...     # do something else

        See notes above for more on `Queues and Renege Clauses`_.
        """
        assert isinstance(S, list) or isinstance(S, tuple), \
                "[QUEUE]: can only insert() list or tuple of objects!"
        assert isinstance(proc, Process), \
                "[QUEUE]: insert() requires valid Process!"
        # handle prio arg
        if self.priority:
            if prio is None: prio = 0
            isPrio = all([isinstance(s,tuple) for s in S] )
            if not isPrio: S = [(s, prio) for s in S]
        else:
            prio = 0
        # create put command
        pcmd = put, proc, self, S, prio
        # set renege clause
        if isinstance(renege, int) or isinstance(renege, float):
            pcmd = pcmd, (hold, proc, renege)
        elif isinstance(renege, SimEvent):
            pcmd = pcmd, (waitevent, proc, renege)
        elif isinstance(renege, tuple) or isinstance(renege, list):
            isEventList = all([isinstance(e, SimEvent) for e in renege] )
            if isEventList:
                pcmd = pcmd, (waitevent, proc, renege)
            else:
                pcmd = pcmd, renege  # assume renege is a valid tuple
        else:
            errmsg = "[QUEUE]: Found invalid renege clause!"
            assert (renege is None), errmsg
        return pcmd

    def remove(self, proc, fn=1, prio=None, priofilter=None, \
               renege=None, **kwargs):
        """Remove a list of objects from the `Queue`.

        :param proc:        `SimPy` `Process` that will block on insert.
        :param fn:          Positive integer or filter function (or "all").
        :param prio:        Priority level of remove.
        :param priofilter:  Priority-based filter function.
        :param renege:      Renege clause.
        :param kwargs:      Additional keywords passed to filter function.

        `proc` is the only mandatory argument.
        
        Normally, `fn` is used to indicate either the number of objects to
        remove from `Queue` or provide a filter function (i.e. the same as
        removing objects from a `SimPy` `Store`).
        
        If `priority` is True and `priofilter` is specified, a list of 2-tuples
        containing *(object, priority level)* is passed to the `priofilter`
        function. Normal filter functions (passed with `fn`) will receive a list
        of buffered objects without priority information.

        The `renege` clause can be one of the following:

            * a timeout value;
            * a (list of) SimEvent(s) to wait on;
            * or tuple representing a valid renege clause (e.g. hold, proc, 0).

        Normal usage with a `renege` clause is as follows:

                >>> yield queue.remove(proc, fn, renege=<renege clause>)
                >>> if proc.acquired(queue):
                ...     # do something
                ... else:
                ...     # do something else

        See notes above for more on `Queues and Renege Clauses`_.

        :note: When `fn` or `priofilter` are used to specify filter functions,
               all objects returned by the filter function will be returned to
               the requesting `Process`.
        """
        isallstr = lambda s: isinstance(s, str) and (s.lower()=="all")
        if isallstr(fn): fn = "all"
        if isallstr(priofilter): priofilter = "all"
        assert isinstance(proc, Process), \
                "[QUEUE]: remove() must be called with a valid Process!"
        assert isinstance(fn, int) or callable(fn) or (fn=="all"), \
                "[QUEUE]: remove() must be called with " + \
                "\'fn\' that is integer, callable, or \"all\"!"
        # set up filter parameters
        if fn=="all": fn = lambda buff: [a for a in buff]
        if priofilter=="all": priofilter = lambda buff: [a for a,p, in buff]
        # set up priority parameters
        if self.priority:
            callback, nget = None, None
            if callable(fn): callback = fn
            elif callable(priofilter): callback = priofilter
            elif isinstance(fn, int): nget = fn
            elif fn is None: nget = 1
            # call private __filter method
            if priofilter is not None:
                filt = lambda buff: self.__filter(buff, callback, nget, prio, True, **kwargs)
            else:
                filt = lambda buff: self.__filter(buff, callback, nget, prio, **kwargs)
            fn = filt
        else:
            prio = 0
            if fn is None: fn = 1
        # create get command
        gcmd = get, proc, self, fn
        if prio is not None: gcmd += (prio, )
        # set renege clause
        if isinstance(renege, int) or isinstance(renege, float):
            # renege is a timeout
            gcmd = gcmd, (hold, proc, renege)
        elif isinstance(renege, SimEvent):
            # renege is a single SimEvent
            gcmd = gcmd, (waitevent, proc, renege)
        elif isinstance(renege, tuple) and isinstance(renege, tuple):
            # renege is an event list or valid tuple
            isEventList = all([isinstance(e, SimEvent) for e in renege] )
            if isEventList:
                gcmd = gcmd, (waitevent, proc, renege)
            else:
                gcmd = gcmd, renege  # assume renege is a valid tuple
        else:
            errmsg = "[QUEUE]: Found invalid renege clause!"
            assert (renege is None), errmsg
        return gcmd

    def _buffer_object(self, x):
        """Internal method to extract actual object from a buffer entry."""
        if self.priority: return x[0]
        else:             return x

    def __filter(self, buff, callback, nget, prio, flag=False, **kwargs):
        """Private filter method to manage filtering in `remove()` with
        `priority`."""
        pbuff = [a for a,b in buff if b>=prio]
        if flag:
            rbuff = callback(buff, **kwargs)
        elif (callback is None) and (nget is not None):
            rbuff = pbuff
        else:
            rbuff = callback(pbuff, **kwargs)
        if len(rbuff)< nget: return []
        else: return rbuff[:nget]
        return rbuff[:nget]

    def __copy_theBuffer(self):
        """Private method to get objects in `Store.theBuffer`."""
        return [self._buffer_object(a) for a in self.theBuffer]

    def __copy_putQ(self, contents=False):
        """Private method to get any `Process` waiting to put objects in `Store`."""
        if contents:
            x = []
            for p in self.putQ:
                x += [self._buffer_object(a) for a in p._whatToPut]
        else:
            x = [p for p in self.putQ]
        return x

    def __copy_getQ(self):
        """Private method to get any `Process` waiting to get objects from `Store`."""
        return [p for p in self.getQ]

    def _get(self, arg):
        """Overload `SimPy` `Store._get()` method."""
        if self.verbose>QUEUE_HIGH_VERBOSE or self.monitorQ:
            # run _get
            prethebuff  = self.__copy_theBuffer()
            preputbuff  = self.__copy_putQ(True)
            rval = Store._get(self, arg)
            postthebuff = self.__copy_theBuffer()
            postputbuff = self.__copy_putQ(True)
            # perform checks
            obj = arg[1]
            getbuff = obj.got
            putbuff = [a for a in preputbuff if \
                        ((a in postthebuff) or (a in getbuff)) ]
            drpbuff = [a for a in preputbuff if \
                        ((a not in postputbuff) and (a not in postthebuff) ) ]
            #self.log("logget")
            self.__log_all({'get': getbuff, \
                            'put': putbuff, \
                            'drp': drpbuff} )
        else:
            rval = Store._get(self, arg)
        return rval

    def _put(self, arg):
        """Overload `SimPy` `Store._put()` method."""
        nevents = len(self.sim._timestamps)
        if self.verbose>QUEUE_HIGH_VERBOSE or self.monitorQ:
            # get objects to request
            reqput = []
            reqput += [self._buffer_object(p) for p in arg[0][3] ]
            # run _put
            prethebuff  = self.__copy_theBuffer()
            preget      = self.__copy_getQ()
            preputbuff  = self.__copy_putQ(True)
            reqput += [p for p in preputbuff]
            rval = Store._put(self, arg)
            postthebuff = self.__copy_theBuffer()
            postget     = self.__copy_getQ()
            postputbuff = self.__copy_putQ(True)
            # perform checks
            getbuff = []
            for proc in [p for p in preget if (p not in postget)]:
                getbuff += [a for a in proc.got]
            putbuff = [a for a in reqput if \
                        ((a in postthebuff) or (a in getbuff) ) ]
            drpbuff = [a for a in reqput if \
                        ((a not in putbuff) and (a not in postputbuff) ) ]
            #self.log("logput")
            self.__log_all({'get': getbuff, \
                            'put': putbuff, \
                            'drp': drpbuff} )
        else:
            rval = Store._put(self, arg)
        # fudge event queue to "de-prior" put event
        numnew = len(self.sim._timestamps) - nevents
        if (0<numnew):
            # process put
            p = heapq.nsmallest(1, [x for x in self.sim._timestamps if (x[0]==self.sim._t)])
            putevt = p[0]
            pt = putevt[0]
            psortpr = putevt[1]
            assert (pt==self.sim._t), "[QUEUE]: _put(), put(at, sim._t) = (%s, %s)"%(pt, self.sim._t)
            assert (psortpr==self.sim._sortpr)
            if 1<numnew:
                glargest = heapq.nlargest(numnew-1, [x for x in self.sim._timestamps if (x[0]==self.sim._t)])
                for k in range(numnew-1):
                    g, krev = glargest[k], (numnew-2-k)
                    # process put with get commands
                    getevt, gt, gsortpr = g, g[0], g[1]
                    assert (gt==self.sim._t), "[QUEUE]: _put(), get(at, sim._t) = (%s, %s)"%(gt, self.sim._t)
                    assert (abs(psortpr+krev+1)==abs(gsortpr)), "getpr = %s, putpr = %s"%(p, g)
                    # set new get pointer and make non-prior
                    getevt[0] = gt + const.EPSILON
                    getevt[1] = abs(gsortpr) + 1
            putevt[1] = abs(psortpr) - (numnew - 1)
        else:
            assert False, "[QUEUE]: _put found unexpected number of events!"
        return rval

    def __log_all(self, buf):
        """Private method to log all get, put, and drop events.

        :param buf: Dictionary of buffers corresponding to each `Queue` event.
        
        Use appropriate callback to execute logging of buffers in `buf`.
        """
        callback = {'put':self.log_insert, \
                    'get':self.log_remove, \
                    'drp':self.log_drop}
        for key,val in buf.items():
            if key in callback:
                for x in val: callback[key](x)

    def log_insert(self, p):
        """Log insert or enqueue event."""
        pargs = {'qlen': self.length}
        if self.verbose>QUEUE_HIGH_VERBOSE: self.log("+", p, **pargs)
        if self.monitorQ: self.enQ.signal(p)

    def log_remove(self, p):
        """Log remove or dequeue event."""
        pargs = {'qlen': self.length}
        if self.verbose>QUEUE_HIGH_VERBOSE: self.log("-", p, **pargs)
        if self.monitorQ: self.deQ.signal(p)

    def log_drop(self, p):
        """Log drop event."""
        pargs = {'qlen': self.length}
        if self.verbose>QUEUE_LOW_VERBOSE: self.log("drp", p, **pargs)
        if self.monitorQ: self.drp.signal(p)

    def addSort(self, sortFunc):
        """
        Overload `Store.addSort()` for priority queueing.
        
        :param sortFunc: Sort method to reorder objects.
        
        `sortFunc` should take two arguments, namely `Queue` and a list
        containing the objects currently buffered by the `Queue`. When operating
        in `priority` mode, `sortFunc` will receive a priority buffer. That is a
        list containing 2-tuples of (*object*, *prioritylevel*).
        """
        return Store.addSort(self, \
                lambda me, buff: me.__callsort(sortFunc, buff) )

    def __callsort(self, callback, buff):
        """Private method to call sorting methods set by `Store.addSort()`;
        enforce priority sorting before calling auxilary sort method."""
        global _prioritySort, _PriorityList
        if callable(callback) and self.priority:
            # sort buff -> call callback
            cbuff = callback(self, _prioritySort(self, buff) )
            isPrio = all([isinstance(s, tuple) for s in cbuff] )
            # make into _PriorityList again
            if isPrio:
                priobuff = _PriorityList(cbuff)
            else:
                priobuff = _PriorityList()
                noprio_buff = [a for a,b in buff]
                for p in cbuff:
                    if p in noprio_buff:
                        prio = buff[noprio_buff.index(p)][1]
                        priobuff.append((p,prio) )
            # re-sort for priority and return
            return _prioritySort(self, priobuff)
        elif callable(callback) and (not self.priority):
            # non-priority callback
            return callback(self, buff)
        elif self.priority:
            # sort for priority
            return _prioritySort(self, buff)
        else:
            # don't do anything
            return buff

    def __len__(self):
        """Allows Python len() function to be used on `Queue`."""
        return self.length
