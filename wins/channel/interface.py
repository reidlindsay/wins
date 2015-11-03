#!  /usr/bin/env python

"""
Interface for connecting nodes to `Channel`; contains `ChannelInterface` class.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-10-04 16:24:27 -0500 (Tue, 04 Oct 2011) $
* $LastChangedRevision: 5184 $

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

:var CHANNELIF_RX: Enumeration for receive state.

:var CHANNELIF_TX: Enumeration for transmit state.

:var CHANNELIF_VERBOSE:
    Constant enumeration to control verbose thresholds in this file.

    Setting the verbose level of a `ChannelInterface` above this threshold will
    cause the corresponding output in this file to be written (or logged).
"""
__docformat__ = "restructuredtext en"

from SimPy.Simulation import SimEvent, hold, now
from wins.helper  import time2usec, linear2db, db2linear, monitor_events
from wins.base    import Reference
from wins.element import Element
from wins.packet  import ANNO, Packet
from wins.fsm     import FSM
from wins         import const

CHANNELIF_RX = "receive"
CHANNELIF_TX = "transmit"
CHANNELIF_VERBOSE = 70

def strcollision(p):
    """Get string representation of 'cif-collision' list from packet `p`."""
    coll = None
    if ANNO.supports(p, 'cif-collision'):
        coll = []
        for c in p.getanno('cif-collision'):
            try:
                s = c.traceid
                assert ANNO.supported(c), "s = %s, coll = %s"%(s, coll)
            except:
                s = c
                assert not ANNO.supported(c), "s = %s, coll = %s"%(s, coll)
            coll.append(s)
    return coll

class ChannelInterface(Element):
    """Interface to connect to a `Channel`.

    This element has four `Port` objects:
        1. "TXU" - sends received packets to upstream element.
        #. "RXU" - receives traffic from an upstream element.
        #. "TXD" - sends traffic to `Channel`.
        #. "RXD" - receives traffic from the `Channel`.

    A `ChannelInterface` must be connected to a `Channel` using the
    `Channel.connect()` method.

    Packet Annotations and ChannelInterface
    =======================================
    The following annotations are marked or used by `ChannelInterface`.

    ============== =================================================
    Name            Description
    ============== =================================================
    cif-duration    Specifies duration of packet so that
                    `ChannelInterface` can simulate transmission and
                    reception of packet. **This must be set by an
                    upstream protocol.**
    -------------- -------------------------------------------------
    cif-collision   List of other Packet objects that have arrived
                    during the reception of a Packet. This collision
                    list should be used to resolve collisions.
    -------------- -------------------------------------------------
    cif-src         Reference to `ChannelInterface` set when sending
                    a packet to the `Channel`.
    -------------- -------------------------------------------------
    cif-dst         Reference to `ChannelInterface` set when 
                    receiving a packet from the `Channel`.
    -------------- -------------------------------------------------
    cif-drp         If marked, indicates that the packet was not
                    completely received because of transmission.
                    This only occurs in halfduplex operation.
    -------------- -------------------------------------------------
    cif-txts        Timestamp for when the packet was transmitted
                    from 'cif-src'.
    -------------- -------------------------------------------------
    cif-rxts        Timestamp for when the packet arrived at the
                    receiver 'cif-dst'.
    ============== =================================================

    :CVariables:
     * `ifstate`: Property to access current state of `ChannelInterface`.
     * `intransmit`: Property to check if in `CHANNELIF_TX` state.
     * `inreceive`:  Property to check if in `CHANNELIF_RX` state.

    :IVariables:
     * `halfduplex`: Boolean flag; if true, then use half-duplex operation,
       otherwise use full-duplex operation.
     * `txdata`: SimEvent signalled when Packet is sent to `Channel`.
     * `txdone`: SimEvent signalled when Packet duration is done simulating.
     * `rxdata`: SimEvent signalled when Packet starts being received.
     * `rxdone`: SimEvent signalled when Packet is done being received.
     * `rxbuffer`: Receive buffer of Packets that are actively being received.
    """
    name = "channel interface"
    tracename = "CIF"
    def __init__(self, **kwargs):
        """Constructor."""
        self.halduplex = None
        self.__ifstate = CHANNELIF_RX
        # set up events and buffer
        self.txdata = SimEvent(name="txdata")
        self.txdone = SimEvent(name="txdone")
        self.rxdata = SimEvent(name="rxdata")
        self.rxdone = SimEvent(name="rxdone")
        self.rxbuffer = []
        Element.__init__(self, **kwargs)
        # rename events
        self.txdata.name = "%s(%s).txdata"%(self.name, self.uid)
        self.txdone.name = "%s(%s).txdone"%(self.name, self.uid)
        self.rxdata.name = "%s(%s).rxdata"%(self.name, self.uid)
        self.rxdone.name = "%s(%s).rxdone"%(self.name, self.uid)
        # monitor events -> keep up to date
        #monitor_events(self.txdata, self.txdone)
        #monitor_events(self.rxdata, self.rxdone)

    ifstate = property(fget=lambda self: self.__ifstate)
    intransmit = property(fget=lambda self: (self.ifstate==CHANNELIF_TX) )
    inreceive  = property(fget=lambda self: (self.ifstate==CHANNELIF_RX) )

    def configure(self, halfduplex=True, **kwargs):
        """Set up ports and other parameters.

        :param halfduplex: Boolean; if true, operate in half-duplex mode;
                           otherwise use full-duplex operation.
        """
        # ports to upstream target
        self.addport("RXU"), self.addport("TXU")
        # ports to downstream target (i.e. a Channel)
        self.addport("TXD"), self.addport("RXD")
        # set up other parameters
        self.halfduplex = halfduplex
        self.newchild("txfsm", FSM, tracename=self.tracename+".TX")
        self.newchild("rxfsm", FSM, tracename=self.tracename+".RX")
        self.txfsm.goto(self.SEND)
        self.rxfsm.goto(self.RECV)

    def SEND(self, fsm):
        """Manage downstream (or outgoing traffic."""
        yield self.RXU.recv(fsm, 1)
        assert fsm.acquired(self.RXU) and (len(fsm.got)==1), \
                "[CHANNELIF]: SEND() error occurred during recv() from RXU!"
        # get packet and set annotations
        p, duration = fsm.got[0], 0
        errmsg = "[CHANNELIF]: Cannot send packet that does not support ANNO!"
        assert ANNO.supported(p), errmsg
        self.set_sendanno(p)
        if p.hasanno('cif-duration'): duration = p.getanno('cif-duration')
        # min-time is const.EPSILON
        if duration<const.EPSILON: duration = const.EPSILON
        p.setanno('cif-duration', duration)
        self.log_send(p)
        # send and simulate duration
        self.__ifstate = CHANNELIF_TX     # start TX
        self.drop("all")                # drop all packet in rxbuffer
        self.txdata.signal(p)
        yield self.TXD.send(fsm, [p])
        yield hold, fsm, duration       # simulate duration
        self.drop("all")                # drop all packet in rxbuffer
        self.__ifstate = CHANNELIF_RX     # resume RX
        self.txdone.signal(p)
        # continue in SEND
        yield fsm.goto(self.SEND)

    def RECV(self, fsm):
        """Manage upstream (or incoming) traffic.

        This method spawn worker processes to simulate the capture of each
        packet. These worker processes manage 'cif-drp' and 'cif-collision'
        annotations; along with simulating packet 'cif-duration'.
        """
        yield self.RXD.recv(fsm, 1)
        errmsg = "[CHANNELIF]: RECV() error occurred during recv() from RXD!"
        assert fsm.acquired(self.RXD) and (len(fsm.got)==1), errmsg
        # start capture thread
        for p in fsm.got:
            w = FSM()
            w.goto(self.CAPTURE, p)
            w.start(prior=True)
        # continue in RECV
        yield fsm.goto(self.RECV)
        assert False, "State transition failed!"

    def CAPTURE(self, fsm, p):
        """Simulate capture process for a Packet `p`."""
        self.log("CAPSTART", p)
        errmsg = "[CHANNELIF]: Cannot capture packet that doesn't support ANNO!"
        assert ANNO.supported(p), errmsg
        duration = const.EPSILON
        self.rxbuffer.append(p)                 # add to rxbuffer
        if self.intransmit: self.drop("all")    # drop all packet in rxbuffer
        # mark/use other annotations
        self.set_recvanno(p)
        if p.hasanno('cif-duration'): duration = p.getanno('cif-duration')
        assert not(duration<const.EPSILON), \
                "[CHANNELIF]: Invlaid duration in CAPTURE! (t=%s)"%(duration)
        # resume operation
        self.log("rxdata.sig", p)
        self.rxdata.signal(p)                   # signal rxdata
        yield hold, fsm, duration               # simulate duration
        if self.intransmit: self.drop("all")    # drop all packet in rxbuffer
        self.rxbuffer.remove(p)                 # remove from rxbuffer
        # drop or forward to upper layer
        if p.hasanno('cif-drp') and p.getanno('cif-drp'):
            self.log_drop( p, halfduplex="%s"%(self.halfduplex) )
            self.cleananno(p)
        else:
            pargs = {'cif-duration':time2usec(duration) }
            self.log_recv(p, **pargs)
            yield self.TXU.send(fsm, [p])
        # signal rxdone
        self.log("rxdone.sig", p)
        self.rxdone.signal(p)
        yield fsm.stop()
        assert False, "stop failed!"

    def set_sendanno(self, p):
        """Set annotations for outgoing traffic prior to sending downstream.

        :return: Modified packet `p`.

        By default, this method sets the 'cif-src' and 'cif-txts' annotations.
        Overload this method as necessary.
        """
        p.setanno('cif-src', self, ref=True)
        p.setanno('cif-txts', now() )
        # remove unwanted annotations in outgoing packets
        rannolist = ['cif-dst', 'cif-collision', 'cif-rxts', \
                     'cif-drp', 'cif-iheap']
        for a in rannolist:
            if p.hasanno(a): p.delanno(a)
        return p

    def set_recvanno(self, p):
        """Set annotation for incoming traffic at start of capture (i.e. right
        after packet has been inserted into `rxbuffer`).

        :return: Modified packet `p`.

        By default, this method initializes the 'cif-collision' annotation, sets the
        'cif-dst' annotation, and sets the 'cif-rxts' annotation. Overload this
        method as necessary.
        """
        assert (p in self.rxbuffer), "[CHANNELIF]: set_recvanno(p) " + \
                "could not find packet 'p' in rxbuffer!"
        p.setanno('cif-dst', self, ref=True)
        p.setanno('cif-collision', [], priv=True)
        # add p to collision list of all in rxbuffer/p
        idx = self.rxbuffer.index(p)
        range_not_idx = range(len(self.rxbuffer) )
        range_not_idx.remove(idx)
        for k in range_not_idx:
            c = self.rxbuffer[k]
            ### XXX ###
            #c.getanno('cif-collision').append(Reference(p) )
            ### XXX ###
            c.getanno('cif-collision').append(p)
        # add rxbuffer/p to collision list of p
        for k in range_not_idx:
            c = self.rxbuffer[k]
            ### XXX ###
            #p.getanno('cif-collision').append(Reference(c))
            ### XXX ###
            p.getanno('cif-collision').append(c)
        # set timestamp for arrival at cif-dst
        p.setanno('cif-rxts', now() )
        return p

    def cleananno(self, p):
        """Remove unwanted annotations from packet `p`.

        Removes any annotations that could cause cyclic references.
        """
        if ANNO.supports(p, 'cif-collision'):
            coll = strcollision(p)
            p.delanno('cif-collision')
            p.setanno('cif-collision', coll, priv=True)
        return p

    def drop(self, p):
        """Set 'cif-drp' annotation in packet p.

        :param p: Packet to mark; or if "all", then mark all packets in
                  `rxbuffer`.
        """
        if isinstance(p, Reference): p = p._deref
        if (p=="all"):
            for c in self.rxbuffer:
                self.drop(c)
        elif p in self.rxbuffer:
            p.setanno('cif-drp', True)
        else:
            raise RuntimeError, \
                  "[CHANNELIF]: drop() could not find packet in rxbuffer!"

    def interval_heap(self, p):
        """Create interval heap from collision list in packet `p`.

        :param p: Packet containing collision list in 'cif-collision' annotation.
        :return: Interval heap.

        This method uses the 'cif-rxts' and 'cif-duration' annotations of packet
        `p` and each packet in its collision list to create an interval heap. To
        do this, the method will create a list of partitions over the duration
        of packet `p` and sort colliding packets into the appropriate
        partitions. An interval heap looks like:

            [(t0,t1,[...]), (t1,t2,[...]), ...]

        This method sets the 'cif-iheap' annotation.
        """
        # get packet p parameters
        for a in ['cif-collision', 'cif-rxts']:
            errmsg = "[CHANNELIF]: interval_heap() requires '%s' annotation!"%(a)
            assert ANNO.supports(p, a), errmsg
        coll = p.getanno('cif-collision')
        duration = const.EPSILON
        if p.hasanno('cif-duration'): duration = p.getanno('cif-duration')
        ta = p.getanno('cif-rxts')
        tb = ta + duration
        # get times for all packets in collision list
        times = [ta, tb]
        for c in coll:
            errmsg = "[CHANNELIF]: interval_heap() requires 'cif-rxts' " + \
                     "annotation in collision list packet!"
            assert ANNO.supports(c, 'cif-rxts'), errmsg
            duration = const.EPSILON
            if c.hasanno('cif-duration'): duration = c.getanno('cif-duration')
            t0 = c.getanno('cif-rxts')  # start of packet c
            t1 = t0 + duration          # end of packet c
            # check if t0, t1 are in times
            t0intimes = any([(abs(t0-t)<2*const.EPSILON) for t in times])
            if (not t0intimes) and (ta<t0<tb): times.append(t0)
            t1intimes = any([(abs(t1-t)<2*const.EPSILON) for t in times])
            if (not t1intimes) and (ta<t1<tb): times.append(t1)
        # sort times and create interval heap
        times.sort()
        iheap = [(times[k], times[k+1], []) for k in range(len(times)-1) ]
        #print "%s: Interval heap for %s (%.8f, %.8f) @ %.8f"%(self.traceid, \
        #        p.traceid, ta,tb, now())
        for c in coll:
            errmsg = "[CHANNELIF]: interval_heap() requires 'cif-rxts' " + \
                     "annotation in collision list packet!"
            assert ANNO.supports(c, 'cif-rxts'), errmsg
            duration = const.EPSILON
            if c.hasanno('cif-duration'): duration = c.getanno('cif-duration')
            t0 = c.getanno('cif-rxts')  # start of packet c
            t1 = t0 + duration          # end of packet c
            # insert into interval heap
            #print "  + inserting %s, (%.8f, %.8f)"%(c.traceid,t0,t1)
            for k in range(len(iheap)):
                ia, ib = iheap[k][0], iheap[k][1]
                errmsg = "[CHANNELIF]: malformed interval in  " + \
                         "interval_heap()! (ia=%s, ib=%s)"%(ia, ib)
                assert (ia<ib), errmsg
                if (t0<ib) and (t1>ia):
                    iheap[k][2].append(Reference(c))
                    #print "    --> inserted  into (%.8f, %.8f)"%(ia,ib)
                else:
                    #print "    --> not added into (%.8f, %.8f)"%(ia,ib)
                    pass
        # set iheap annotation
        p.setanno('cif-iheap', iheap, priv=True)
        return iheap

    def sinr_heap(self, p, force=True):
        """Calculate signal-to-interference-and noise ratio (SINR) for each
        partition created by `interval_heap()`.

        :param p: Packet to inspect.
        :param force: If true, recalculate interval heap, else use existing
                      annotation if it exists.
        :return: SINR heap.

        SINR heap has looks like this:

            [(t0, t1, sinr0), (t1, t2, sinr1), ... ]

        Note: This method uses the 'rxpower' and 'noisepower' annotations.
        """
        # check packet
        errmsg = "[CHANNELIF]: sinr_heap() cannot process non-Packet!"
        assert ANNO.supported(p), errmsg
        for a in ['rxpower', 'noisepower']:
            errmsg = "[CHANNELIF]: sinr_heap() cannot find '%s' annotation!"%(a)
            assert ANNO.supports(p, a), errmsg
        # get parameters
        rxpower = p.getanno('rxpower')          # in dBm
        noisepower = p.getanno('noisepower')    # in dBm
        npow = db2linear(noisepower)
        # get interval heap
        if p.hasanno('cif-iheap') and not force:
            iheap = p.getanno('cif-iheap')
        else:
            iheap = self.interval_heap(p)
        # start creating sinr heap
        sinrheap = []
        #print "%s: SINR heap for %s @ %.8f"%(self.traceid, p.traceid, now())
        for ta,tb,coll in iheap:
            ipow = 0
            for c in coll:
                errmsg = "[CHANNELIF]: sinr_heap() cannot find 'rxpower' " + \
                         "annotation in collision list!"
                assert ANNO.supports(c, 'rxpower'), errmsg
                ipow += db2linear(c.getanno('rxpower') )
            sinr = rxpower - linear2db(ipow + npow)
            sinrheap.append((ta,tb,sinr) )
            #print "  --> (%.8f, %.8f): %.3f dB, coll = %s"%(ta,tb, sinr, [c.traceid for c in coll])
        return sinrheap

    def log_send(self, p, *args, **kwargs):
        """Convenience method for logging a send event for packet `p`."""
        if self.verbose>CHANNELIF_VERBOSE:
            kwargs.update(self.get_cif_anno(p))
            self.log("snd", p, *args, **kwargs)

    def log_recv(self, p, *args, **kwargs):
        """Convenience method for logging a receive event for packet `p`."""
        if self.verbose>CHANNELIF_VERBOSE:
            if p.hasanno('cif-src') and ('cif-src' not in kwargs):
                kwargs['cif-src'] = p.getanno('cif-src').traceid
            if p.hasanno('cif-dst') and ('cif-dst' not in kwargs):
                kwargs['cif-dst'] = p.getanno('cif-dst').traceid
            kwargs.update(self.get_cif_anno(p))
            self.log("rcv", p, *args, **kwargs)

    def log_drop(self, p, *args, **kwargs):
        """Convenience method for logging a drop event for packet `p`."""
        if self.verbose>CHANNELIF_VERBOSE:
            if 'halfduplex' not in kwargs: kwargs['halfduplex'] = self.halfduplex
            self.log("drp", p, *args, **kwargs)

    def log(self, event=None, p=None, *args, **kwargs):
        """Overloaded to check verbose level and set common annotations."""
        force = False
        if ('verbose' in kwargs): force = (kwargs['verbose']>CHANNELIF_VERBOSE)
        if self.verbose>CHANNELIF_VERBOSE or force:
            kwargs.update(self.get_cif_anno(p))
            Element.log(self, event, p, *args, **kwargs)

    def get_cif_anno(self, p):
        """Convenience method to extract annotations and convert to strings."""
        kwargs = {}
        if not isinstance(p, Packet): return kwargs
        if p.hasanno('cif-collision'):
            kwargs['cif-collision'] = strcollision(p)
        if p.hasanno('cif-duration'):
            kwargs['cif-duration'] = time2usec(p.getanno('cif-duration') )
        if p.hasanno('cif-src'):
            kwargs['cif-src'] = "%s"%(p.getanno('cif-src').traceid)
        if p.hasanno('cif-dst'):
            kwargs['cif-dst'] = "%s"%(p.getanno('cif-dst').traceid)
        return kwargs
