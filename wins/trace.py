#!  /usr/bin/env python

"""
Implementation of trace functionality; contains `Trace` and `Traceable`.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-11-01 20:25:59 -0500 (Tue, 01 Nov 2011) $
* $LastChangedRevision: 5325 $

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

:var USE_ROOT_AS_ROOT:
    Enumeration for option in which the 'root' annotation of a Packet is used as
    its root for the purpose of tracing.

:var USE_PARENT_AS_ROOT:
    Enumeration for option in which the 'parent' annotation of a Packet is used
    as its root for the purpose of tracing; the parent must be `Traceable`.

:var USE_PAYLOAD_AS_ROOT:
    Enumeration for option in which the Packet payload is used as its root for
    the purpose of tracing.

:var PACKET_ROOT_MODE:
    Determines method by which trace will determine a Packet's root for the
    purposes of tracing.

:var BATCHSIZE:
    Determines the maximum number of events logged by each temporary trace file.
"""
__docformat__ = "restructuredtext en"

from SimPy.Simulation import now
from wins.base import Base, Reference
from wins.packet import ANNO
from wins.helper import isnum

from scapy.all import Packet, NoPayload, IP

import os
import sys
import shutil
from copy import deepcopy


USE_ROOT_AS_ROOT = 0
USE_PARENT_AS_ROOT = 1
USE_PAYLOAD_AS_ROOT = 2
PACKET_ROOT_MODE = 2

import tempfile
import weakref
TempfileRepo = []
BATCHSIZE = 1000000

class Trace(Base):
    """Maintains log of events.

    :cvar Global: Trace with global scope.
    :cvar knownevents: Known events
    :cvar events: Property to access `_events`
    :cvar nevents: Property to access `_nevents`

    :ivar _events: Internal list of events maintained by trace.
    :ivar _nevents: Internal counter to keep track of number of events logged.
    :ivar _maxlen: Internal dictionary
    """
    name = "trace"
    Global = None
    knownevents = {"rcv":"rcv", "snd":"snd", \
                   "enq":"+", "deq":"-", \
                   "drop":"drp", "forward":"fwd"}

    def __init__(self, **kwargs):
        """Constructor."""
        Base.__init__(self, **kwargs)
        self._events = []
        self._maxlen = None
        self._nevents = 0
        self._tempfiles = []

    events = property(fget=lambda self: self.get_events() )
    nevents = property(fget=lambda self: self._nevents)

    def log(self, **kwargs):
        """Log event as a dictionary.

        :param kwargs: Event is a dictionary of keyword arguments.
        """
        # check if new temporary file is needed
        if ((self._nevents%BATCHSIZE)==0):
            t = tempfile.NamedTemporaryFile()
            TempfileRepo.append(t)
            self._tempfiles.append(weakref.proxy(t))
            t = None
        # log event into current tempfile
        currtemp = self._tempfiles[-1]
        currtemp.write(str(kwargs)+"\n")
        self._nevents += 1

    def reset(self):
        """Reset list of `events`."""
        for f in self._tempfiles:
            f.close()
            m = None
            for g in TempfileRepo:
                if (g.name==f.name): m = g
            if m:
                TempfileRepo.remove(m)
        # reset remaining parameters
        m = None
        self._tempfiles = []
        self._events = []
        self._maxlen = None
        self._nevents = 0

    def sync(self):
        """Synchronize event list with any events batched and stored in
        temporary files.

        This method loads all events from temporary files into `events`.
        """
        needsync = (len(self._events) != self._nevents)
        if not needsync: return
        # read through tempfiles and copy to event list
        ncopy, self._events = 0, []
        for f in self._tempfiles:
            f.seek(0)
            s = f.readline()
            while s:
                try:
                    d = eval(s)
                    assert isinstance(d, dict)
                except:
                    d = None
                # add event
                if d:
                    self._events.append(d)
                    ncopy += 1
                s = f.readline()
        # check if all events were read
        errmsg = "[TRACE]: Invalid number of events found! " + \
                 "ncopy = %d, _nevents = %d"%(ncopy, self._nevents)
        assert (ncopy==self._nevents), errmsg
        return ncopy

    def get_events(self):
        """Get event list; calls `sync()`."""
        self.sync()
        return self._events

    def header(self):
        """Get header for trace."""
        s, maxlen = "", self._max_fieldlen()
        fields = ['ts', 'event', 'obj', 'uid', 'parent', 'cid', \
                  'node', 'nid', 'packet', 'pid', 'root', 'rid', \
                  'net', 'netid', 'plen']
        for f in fields:
            mlen = maxlen[f]
            if (f=="ts"):
                mlen = int(maxlen[f])
            exec("s += \"%%%ss\"%%(f)"%(mlen) )
        h = "\n"+"-"*len(s)
        return h+"\n"+s.upper()+h

    def format(self, evt):
        """Format event into readable string."""
        s, maxlen = "", self._max_fieldlen()
        fields = ['ts', 'event', 'obj', 'uid', 'parent', 'cid', \
                  'node', 'nid', 'packet', 'pid', 'root', 'rid', \
                  'net', 'netid', 'plen']
        event = deepcopy(evt)
        for f in fields:
            mlen = maxlen[f]
            fmt = "%%%ss"%(mlen)
            if (f=="ts"): fmt = "%%%.1ff"%(mlen)
            exec("s += \"%s\"%%(event[f])"%(fmt) )
            del event[f]
        # print remaining arguments
        if 'debug' in evt: del event['debug']
        if event: s += "    %s"%(event)
        if 'debug' in evt:
            tslen = int(maxlen['ts'])
            s += "\n" + " "*(tslen) + " debug: %s"%(evt['debug'])
        return s

    def output(self, stderr=False):
        """Print formatted version of trace.

        :param stderr: Boolean; if true, output to stderr, otherwise use stdout.
        """
        output = sys.stdout.write
        if stderr: output = sys.stderr.write
        self.sync()
        # format output
        header = self.header()
        output("%s\n"%(header) )
        for e in self._events:
            evt = self.format(e)
            output("%s\n"%(evt) )
        return self._nevents

    def write(self, outfile):
        """Write all trace events to a file.

        :param outfile: Output file to be written.
        :return: Number of trace entries written.

        :note: If outfile is the name of a file, this method will open the file
               for writing in the current directory.
        """
        # write out from self._events
        f, nout = outfile, 0
        if isinstance(outfile, str): f = file(outfile, "w")
        assert isinstance(f, file), "[TRACE]: Cannot write() to non-file!"
        for e in self.events:
            f.write(str(e) + "\n")
            nout += 1
        if isinstance(outfile, str): f.close()
        return nout

    def read(self, infile):
        """Read events from a trace file generated by `write()`.

        :param infile: Input file (or name of file) to be read from.
        :return: Number of trace entries read.

        :note: This method will overwrite all `events` by calling `reset()`.
        """
        f, nin = infile, 0
        if isinstance(infile, str): f = file(infile, "r")
        assert isinstance(f, file), "[TRACE]: Cannot read() from non-file!"
        self.reset()
        s = f.readline()
        while (len(s)>0):
            try:    d = eval(s)
            except: d = None
            if isinstance(d, dict):
                self.log(**d)
                nin += 1
            s = f.readline()
        if isinstance(infile, str): f.close()
        return nin

    def _max_fieldlen(self, force=False, pad=4):
        """Internal method to compute maximum length of fields in trace."""
        if (not force) and (self._maxlen is not None): return self._maxlen
        # iterate over events to find max field length
        self._maxlen = {'ts':16, 'event':10, 'obj': 5, 'uid':5, \
                        'parent':10, 'cid':5, 'node':10, 'nid':5, \
                        'packet':10, 'pid':5, 'root':10, 'rid':5, \
                        'net':10, 'netid':5, 'plen':5}
        for e in self.events:
            for key in e:
                l = len(str(e[key]) ) + pad
                if key not in self._maxlen: self._maxlen[key] = l
                if l>self._maxlen[key]: self._maxlen[key] = l
        self._maxlen['ts'] = 16.6
        return self._maxlen

    @staticmethod
    def _init_global(g=None):
        """Internal method to initialize `Global`."""
        if Trace.Global is None:
            if g is None: g = Trace()
            Trace.Global = g #Reference(g)

    @staticmethod
    def _set_global(tr):
        """Internal method to set `Global`."""
        if isinstance(tr, Reference): tr = tr._deref
        assert isinstance(tr, Trace), \
                "[TRACE]: Cannot set global trace to non-Trace object!"
        Trace.Global = tr #Reference(tr)

Trace._init_global()

class Traceable(Base):
    """Base class to provide convenient API for accessing a `Trace`.

    :cvar name: Name of traceable object.
    :cvar tracename: Name used in `Trace`.
    :cvar traceid: Formatted trace name with id information.
    :cvar trace: Property to access/modify `Trace` used for logging.

    :ivar __trace: Private pointer to `Trace` object[default=`Trace.Global`];
                   use property to access/modify this.
    """
    name = "traceable"
    tracename = "TRACEABLE"
    def __init__(self, tracename=None, trace=None, **kwargs):
        """Constructor.

        :param tracename: Name to be used in logging events.
        :param trace: `Trace` to be used [default=`Trace.Global`].
        """
        Base.__init__(self, **kwargs)
        if tracename is None: tracename = self.__class__.tracename
        self.tracename = tracename
        if trace is None: trace = Reference(Trace.Global)
        self.__trace = trace

    trace = property(fget=lambda self: self.__trace, \
                     fset=lambda self,t: self.set_trace(t) )
    traceid = property(fget=lambda self: "%s%s"%(self.tracename,self.uid) )

    def addchild(self, *args, **kwargs):
        """Overloaded to set `trace` if possible."""
        c = Base.addchild(self, *args, **kwargs)
        if hasattr(c, 'trace'): c.trace = self.trace
        return c

    def set_trace(self, t, recursive=True):
        """Set `trace` to log events.

        :param t: New `Trace`.
        :param recursive: Boolean; if true, call `set_trace()` on all children
                          that support it.
        """
        if isinstance(t, Reference): t = t._deref
        assert isinstance(t, Trace), \
                "[TRACEABLE]: Cannot set_trace() to non-Trace object!"
        self.__trace = t
        if not recursive: return
        # otherwise set children
        self.callchild('set_trace', t, recursive=recursive)

    def log(self, event=None, packet=None, **kwargs):
        """Log events in `trace`."""
        ts = now()
        if event is None: event = "-----"
        event = str(event).lower()
        if event in Trace.knownevents:
            event = Trace.knownevents[event]
        event = event.upper()
        obj = self.tracename
        uid = self.uid
        cid, parent  = "---", "-----"
        if isinstance(self.parent, Traceable):
            cid, parent = self.parent.uid, self.parent.tracename
        nid, node = "---", "---"
        if hasattr(self, 'container'):
            if hasattr(self.container,'tracename'):
                node = self.container.tracename
                nid = self.container.uid
        # get packet info
        pid, pkt, plen = "---", "-----", "---"
        rid, root = "---", "-----"
        netid, net = "---", "-----"
        if isinstance(packet, Reference): packet = packet._deref
        if isinstance(packet, Packet):
            pid, pkt, plen = packet._id, packet.tracename, len(packet)
            if PACKET_ROOT_MODE==USE_ROOT_AS_ROOT:
                if packet.hasanno('root'):
                    r = packet.getanno('root')
                    if isinstance(r, Packet):
                        rid, root = r._id, r.tracename
            elif PACKET_ROOT_MODE==USE_PARENT_AS_ROOT:
                if packet.hasanno('parent'):
                    r = packet.getanno('parent')
                    if isinstance(r, Traceable):
                        rid, root = r.uid, r.tracename
            elif PACKET_ROOT_MODE==USE_PAYLOAD_AS_ROOT:
                ispayload = lambda p: isinstance(p, Packet) and \
                        not(isinstance(p, NoPayload) or isinstance(p, ANNO) )
                if ispayload(packet.payload):
                    r = packet.payload
                    rid, root = r._id, r.tracename
            if packet.haslayer(IP):
                ip = packet[IP]
                netid, net = ip._id, ip.tracename
        # set remaining parameters
        largs = {'packet':str(pkt)}
        for s in ['ts','event','obj','uid','parent','cid','node','nid', \
                  'pid', 'root', 'rid', 'net', 'netid', 'plen']:
            if isnum(eval(s)):
                largs[s] = deepcopy(eval(s))
            else:
                largs[s] = str(eval(s))
        for k,v in kwargs.items():
            largs[k] = str(v)
        self.trace.log(**largs)
        """
        self.trace.log(ts=ts, event=event, obj=obj, uid=uid, \
                       parent=parent, cid=cid, node=node, nid=nid, \
                       packet=pkt, pid=pid, root=root, rid=rid, plen=plen, \
                       **kwargs)
        """

    def debug(self, s, *args, **kwargs):
        """Log debug string in `trace`."""
        kwargs['debug'] = str(s)
        self.log("debug", *args, **kwargs)
