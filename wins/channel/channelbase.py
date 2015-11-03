#!  /usr/bin/env python

"""
Base classes needed for implementing channel to connect nodes together; contains
`Channel` and `ChannelModel` classes.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-12-18 16:54:42 -0600 (Sun, 18 Dec 2011) $
* $LastChangedRevision: 5380 $

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

:var CHANNEL_VERBOSE:
    Constant enumeration to control verbose thresholds in this file.
  
    Setting the verbose level of a `Channel` above this threshold will cause the
    corresponding output in this file to be written (or logged).
"""
__docformat__ = "restructuredtext en"

from SimPy.Simulation import now
from wins.base    import Reference
from wins.element import Element, Port
from wins.trace   import Traceable
from wins.packet  import Packet, ANNO
from wins.fsm     import FSM
from wins.channel.interface import ChannelInterface
from wins.helper  import time2usec
from wins.crc import crcupdate

import networkx as nx
from copy import deepcopy

CHANNEL_VERBOSE=80

class Channel(Element):
    """Maintain connection graph in which each node is a `ChannelInterface` and
    each link is a `ChannelModel`.

    The `model` and `bidirectional` parameters can be set through keywords
    passed to the constructor (i.e. passed to `configure()`). The basic
    operation of a `Channel` is as follows:

        1. Get packets transmitted from `ChannelInterface`.
        #. Create deepcopy of packets, and forward on all outgoing edges in
           `graph`, after applying annotations and appropriate `ChannelModel`.
        #. Check for 'cm-delay' annotation to apply propagation delay.
       

    Packet Annotations and Channel
    ==============================
    Channel will mark or use the following annotations when forwarding packets.

    ========== ==============================================================
    Name        Description
    ========== ==============================================================
    cif-src     Reference to `ChannelInterface` that sent packet.
    ---------- --------------------------------------------------------------
    cif-dst     Reference to `ChannelInterface` to which packet is forwarded.
    ---------- --------------------------------------------------------------
    cm-delay    Propagation delay annotation marked by `ChannelModel`.
    ========== ==============================================================

    If the 'cm-delay' annotation is set by a `ChannelModel`, then `Channel` will
    forward the packet after applying the appropriate delay.

    :CVariables:
     * `model`: Property to access default `ChannelModel` used to create new
       links between `ChannelInterface` objects in `graph`.
     * `bidirectional`: Immutable boolean flag; if true, indicates that the
       underlying graph is bidirectional (or undirected).
     * `graph`: Property to access underlying NetworkX graph.
     * `interfaces`: Property to access nodes in `graph` for `Channel`.
     * `listen`: Property to access private dictionary for processes listening
       for data from all `ChannelInterface` connected.
    """
    name = "channel"
    tracename = "CHAN"
    def __init__(self, **kwargs):
        """Constructor."""
        self.__model = None
        self.__bidirectional = None
        self.__graph = None
        self.__listen = {}
        Element.__init__(self, **kwargs)

    model = property(fget=lambda self: self.__model, \
                     fset=lambda self,m: self.set_model(m) )
    bidirectional = property(fget=lambda self: self.__bidirectional)
    graph = property(fget=lambda self: self.__graph)
    interfaces = property(fget=lambda self: self.graph.nodes() )
    listen = property(fget=lambda self: self.__listen)

    def configure(self, model=None, bidirectional=False, **kwargs):
        """Set up `graph` and other parameters

        :param bidirectional: Boolean; if true, use undirected graph, otherwise
                              use a directed graph.
        :param model: Default `ChannelModel` used to create links in `graph`.
        """
        self.model = model
        self.__bidirectional = bidirectional
        if self.bidirectional:
            self.__graph = nx.Graph()
        else:
            self.__graph = nx.DiGraph()

    def set_model(self, m):
        """Set default `ChannelModel` used to create links in `graph`."""
        if m is None: m = ChannelModel
        assert issubclass(m, ChannelModel), \
                "[CHANNEL]: model must be a ChannelModel!"
        self.__model = m

    def connect(self, cif):
        """Add `ChannelInterface` as a node in graph.

        :param cif: `ChannelInterface` to add.
        :return: Added `ChannelInterface`.

        :note: If `cif` is already in `graph`, this method will do nothing.
        """
        assert isinstance(cif, ChannelInterface), \
                "[CHANNEL]: Cannot add non-ChannelInterface!"
        if isinstance(cif, Reference): cif = cif._deref
        hascif = (cif in self.graph) and (self.haschild((cif, "FSM")) )
        if not hascif:
            self.graph.add_node(cif)
            # create input/output ports
            rxport = self.addport((cif,"RX") )  # data sent by cif to channel
            txport = self.addport((cif,"TX") )  # data sent from channel to cif
            # connect ports
            ciftx, cifrx = cif.TXD, cif.RXD
            ciftx.connect(rxport)
            txport.connect(cifrx)
            # create FSM to listen
            self.newchild((cif, "FSM"), FSM, tracename=cif.traceid+".LISTEN")
            self.listen[cif] = self.getchild((cif, "FSM") )
            self.listen[cif].goto(self.LISTEN, cif, rxport)
            self.listen[cif].start()
        return cif

    def disconnect(self, cif):
        """Disconnect a previously connected `ChannelInterface`.

        After removing the interface and its listener process, this method will
        delete the `ChannelInterface` from the `graph` and each `ChannelModel`
        for every outgoing edge from `cif`.

        :note: After disconnecting an interface, it will still be able to send
               traffic from its ports, but this traffic will not be received by
               the `Channel`.
        """
        if isinstance(cif, Reference): cif = cif._deref
        delcif = (cif in self.listen) and (cif in self.graph) \
                        and (self.haschild((cif, "FSM")) ) 
        if delcif:
            # stop listener and delete child FSM
            f = self.listen[cif]
            del self.listen[cif]
            self.delchild(f)
            assert isinstance(f, FSM), "[CHANNEL]: Found non-FSM listener!"
            f.halt()
            # disconnect ports
            cif.TXD.disconnect()                # downstream port from cif
            rxport = self.delport((cif,"RX") )  # data sent by cif to channel
            txport = self.delport((cif,"TX") )  # data sent from channel to cif
            # remove all outgoing edges in graph and corresponding node
            neighbors = [c for c in self.graph.neighbors(cif) ]
            for n in neighbors:
                self.del_edge(cif, n)
            self.graph.del_node(cif)

    def has_edge(self, a, b):
        """Check if a link exists between `ChannelInterface` `a` and `b`."""
        if isinstance(a, Reference): a = a._deref
        if isinstance(b, Reference): b = b._deref
        hasedge = self.graph.has_edge(a,b)
        haschild = (self.haschild((a,b)) is not None)
        return hasedge and haschild

    def get_edge(self, a, b):
        """Get edge between `ChannelInterface` a and b."""
        if isinstance(a, Reference): a = a._deref
        if isinstance(b, Reference): b = b._deref
        return self.graph[a][b]['model']

    def add_edge(self, a, b, model=None, **kwargs):
        """Add a link between `ChannelInterface` `a` and `b`.

        :param a: `ChannelInterface` on one end of link.
        :param b: `ChannelInterface` on other end of link.
        :param model: Optional `ChannelModel` subclass to create link with.
        :param kwargs: Keywords passed to channel model constructor.
        :return: Newly created link (i.e. channel model)

        If model is not specified as a parameter, the default `model` will be
        used to create the link.
        """
        u = self.connect(a)
        v = self.connect(b)
        if model is None: model = self.model
        assert issubclass(model, ChannelModel), \
                "[CHANNEL]: Cannot create graph edge without valid ChannelModel!"
        # remove existing edge
        if self.has_edge(u,v): self.del_edge(u,v)
        cm = self.newchild((u,v), model, **kwargs)
        self.graph.add_edge(u,v, model=cm)
        return cm

    def del_edge(self, a, b):
        """Delete edge between `ChannelInterface` `a` and `b`.

        :return: Deleted `ChannelModel` for removed link or `None` if not found.
        """
        cm = None
        if self.has_edge(a,b):
            u = self.connect(a)
            v = self.connect(b)
            cm = self.get_edge(u,v)
            # check if channel model is a child
            assert (self.haschild(cm) is not None), "[CHANNEL]: del_edge() " + \
                    "found non-child ChannelModel in graph!"
            self.graph.remove_edge(u,v) # remove edge from graph
            self.delchild(cm)           # remove channel model as child
        return cm

    def apply_filter(self, model, p, *args, **kwargs):
        """Apply filter for `model` to `p`.

        :param model: `ChannelModel` to do filtering.
        :param p: Packet to be passed to filter.
        :param args: Additional arguments passed to filter of `model`.
        :param kwargs: Additional keyword arguments passed to filter of `model`.
        :return: Reason for drop (string), or `None` if packet was not filtered
                 out by the `ChannelModel`.
        """
        errmsg = "[CHANNEL]: Got non-ChannelModel in apply_filter()!"
        assert isinstance(model, ChannelModel), errmsg
        return model.filter(p, *args, **kwargs)

    def apply_model(self, model, p, u, v):
        """Apply a `ChannelModel` to a given packet.

        :param model: `ChannelModel` to apply.
        :param p: Packet being modified.
        :param u: Source `ChannelInterface` that sent `p`.
        :param v: Destination `ChannelInterface`.
        :return: Modified packet from `ChannelModel.apply()`.
        """
        errmsg = "[CHANNEL]: Got non-ChannelModel in apply_model()!"
        assert isinstance(model, ChannelModel), errmsg
        r = model.apply(p, u, v)
        errmsg = "[CHANNEL]: ChannelModel.apply() returned invalid value!"
        assert (isinstance(r, Packet) or (r is None)), errmsg
        return r

    def LISTEN(self, fsm, cif, rxport):
        """Listen to traffic sent by `ChannelInterface` on `rxport`."""
        assert (self.hasport(rxport) is not None), \
                "[CHANNEL]: Cannot listen on invalid Port!"
        # get data from rxport
        yield rxport.recv(fsm, 1)
        assert fsm.acquired(rxport) and (len(fsm.got)==1), \
                "[CHANNEL]: Error occurred during LISTEN on %s!"%(rxport)
        # forward on all outgoing edge
        p, u = fsm.got[0], cif
        neighbors = [str(c.traceid) for c in self.graph.neighbors(u) ]
        errmsg = "[CHANNEL]: Cannot forward non-Packet!"
        assert isinstance(p, Packet), errmsg
        for v in self.graph[u]:
            cmodel = self.graph[u][v]['model']
            errmsg = "[CHANNEL]: Cannot find corresponding TX Port!"
            assert (self.hasport((v, "TX") ) is not None)
            # if channel model for link is not active -> ignore packet
            if not cmodel.active: continue
            # filter packets using channel model
            drop = self.apply_filter(cmodel, p, u, v)
            if drop:
                self.log("drp",p,drop=drop,src=u.traceid,dst=v.traceid)
                continue            # continue with next link
            # copy and mark annotations
            c, txport = p.copy(), self.getport((v, "TX") )
            c.setanno('cif-src', u, ref=True)
            c.setanno('cif-dst', v, ref=True)
            # apply channel model
            propdelay = None
            r = self.apply_model(cmodel, c, u, v)
            if cmodel.verbose>CHANNEL_VERBOSE:
                cmodel.log_forward(c, src=u.traceid, dst=v.traceid)
            if ANNO.supports(r, 'cm-delay'):
                propdelay = r.getanno('cm-delay')
            # forward to destination
            if (r is None):
                self.log("drp",c,action="dropped by %s"%(cmodel.traceid) )
            elif (propdelay>0):
                # apply propagation delay
                f = FSM()
                f.goto(self.FORWARD, txport, [r])
                f.start(delay=propdelay)
            else:
                # send immediately (no delay)
                yield txport.send(fsm, [r])
        if self.verbose>CHANNEL_VERBOSE:
            self.log("fwd", p, dest="%s"%(neighbors) )
        # continue in LISTEN
        yield fsm.goto(self.LISTEN, cif, rxport)

    def FORWARD(self, fsm, txport, S):
        """FORWARD state; internal state to forward packets for channel."""
        errmsg = "[CHANNEL]: Cannot FORWARD() to non-Port!"
        assert isinstance(txport, Port), errmsg
        yield txport.send(fsm, S)

    def log(self, evt=None, p=None, *args, **kwargs):
        """Overloaded to check verbose level and set common annotations."""
        force = False
        if ('verbose' in kwargs): force = (kwargs['verbose']>CHANNEL_VERBOSE)
        if self.verbose>CHANNEL_VERBOSE or force:
            Element.log(self, evt, p, *args, **kwargs)

class ChannelModel(Traceable):
    """Implement channel model between two `ChannelInterface`.

    Packet Annotations and ChannelModel
    ===================================
    `ChannelModel` may set any number of annotations to be used in determining a
    variety of parameters, including: packet error rate, path loss, signal to
    noise ratio, propagation delay, etc. The only annotation used by `Channel`
    is the 'cm-delay' annotation (see `Packet Annotations and Channel`).

    `Channel` calls the `apply()` method before forwarding the packet from the
    source `ChannelInterface` to the destination.
    """
    name = "channel model"
    tracename = "LNK"
    def __init__(self, **kwargs):
        """Constructor."""
        self.__active = True
        Traceable.__init__(self, **kwargs)
        self.tracename = "%s%s"%(self.tracename, self.uid)

    active = property(fget=lambda self: self.__active)

    def filter(self, p, src, dst):
        """Filter out packets that the `ChannelModel` will drop.

        :param p: Packet to filter.
        :param src: Source `ChannelInterface` that sent `p`.
        :param dst: Destination `ChannelInterface`.
        :return: Reason for drop (string), or `None` if packet was not filtered
                 out by the `ChannelModel`.

        By default this method does nothing (i.e. all packets pass through
        filter). *Overload this in derived class to change filter.*
        """
        return None

    def apply(self, p, src, dst):
        """Apply `ChannelModel` to specified packet.

        :param p: Packet to be forwarded from `src`.
        :param src: Source `ChannelInterface` that sent `p`.
        :param dst: Destination `ChannelInterface`.

        :return: Modified packet with appropriate annotations to be forwarded to
                 destintation `ChannelInterface` by `Channel`.

        By default this method returns `p` after setting the 'cm-delay'
        annotation to zero (i.e. denoting zero propagation delay).

        :note: If this method returns `None`, then the packet will be dropped by
               the parent `Channel` in `Channel.apply_model()`.
        """
        p.setanno('cm-delay', 0)
        return p

    def enable(self):
        """Set active flag."""
        self.__active = True

    def disable(self):
        """Clear active flag."""
        self.__active = False

    def log_forward(self, p, *args, **kwargs):
        """Convenience method for logging a forward event for packet `p`."""
        if p.hasanno('cm-delay'):
            delay = p.getanno('cm-delay')
            kwargs['cm-delay'] = time2usec(delay)
        self.log("fwd", p, *args, **kwargs)
