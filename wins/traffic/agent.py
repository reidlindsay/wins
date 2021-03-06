#!  /usr/bin/env python

"""
Agent for generating and receiving traffic; contains `AGT` and `Agent`.

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

:var AGENT_ID_INDEX: Boolean; if true, use packet id to index receive buffers.
"""
__docformat__ = "restructuredtext en"

from SimPy.Simulation import SimEvent, waitevent, hold, now
from wins import Element, Base, FSM, NET, MAC
from wins import const, strcollision
from scapy.all import Packet, IntField

import struct
import numpy as np
import sys

AGENT_VERBOSE = -1  # always log to trace
AGENT_ID_INDEX = 1

class AGT(Packet):
    """Packet header for traffic generated by Agent."""
    fields_desc = [IntField("id", 0), IntField("seqno", 0)]

class Agent(Element):
    """Act as a source and sink.

    This module generates (and receives) traffic for (and from) other Agents.
    There are two ports for traffic: the egress port 'TX' and ingress port 'RX'.
    This traffic module can operate in one of three modes: 'cbr', 'poisson', and
    the default 'backlog'.
    
    The traffic source generates fixed length packets, but varies the
    interpacket wait time based on the `mode` of operation and the average
    `delay`. The 'backlog' mode only works if the parent container of the
    `Agent` is a node containing a valid 'mac' object.

    If the parent container of `Agent` is a node containing a valid 'net'
    object, this module will use the source and destination addresses of the
    received packets for indexing into the `recvbuffer`.

    :ivar nsent: Number of packets sent by `Agent`.
    :ivar nrcvd: Number of packets received by `Agent`.
    :ivar delay: Average delay (in sec).
    :ivar plen:  Fixed packet length (in bytes).
    :ivar recvbuffer: Receive buffer maintained to keep track of duplicate
                      packets; buffer size is `RecvBufferSize`.

    :cvar DefaultPacketLength: Default packet length (in bytes).
    :cvar MaxSeqNo: Maximum sequence number.
    :cvar TrafficModes: List of supported traffic generating modes.
    :cvar DefaultTrafficMode: Default traffic generating mode.
    :cvar RecvBufferSize: Buffer size of `recvbuffer`.
    :cvar mode: Property access/modify current mode of operation
                (see valid `TrafficModes`).
    :cvar mac: Property to access `MAC` module of parent node.
    :cvar addr: Property to access address of `NET` module of parent node.
    :cvar seqno: Property to access internal sequence number.
    :cvar broadcast: Property to access broadcast address of `NET` module of
                     parent node.
    """
    name = "agent"
    tracename = "AGT"
    DefaultPacketLength = 1024
    MaxSeqNo = 65536
    TrafficModes = ["backlog", "cbr", "poisson"]
    DefaultTrafficMode = "backlog"
    RecvBufferSize = 16
    def __init__(self, *args, **kwargs):
        """Constructor."""
        self._seqno, self._mode = 0, None
        self._payload = None
        # set parameters
        self.delay = None
        self.nsent, self.nrcvd = 0, 0
        self.dest, self.plen = None, None
        self.recv, self.send = SimEvent(), SimEvent()
        self.recvbuffer = {}
        Element.__init__(self, *args, **kwargs)

    mode = property(fget=lambda self: self.get_mode(), \
                    fset=lambda self, m: self.set_mode(m))
    mac = property(fget=lambda self: self.get_mac())
    addr = property(fget=lambda self: self.get_address())
    seqno = property(fget=lambda self: self._seqno)
    broadcast = property(fget=lambda self: self.get_address(broadcast=True))

    def configure(self, dest=None, plen=None, delay=None, mode=None, **kwargs):
        """Set up parameters and start agent if needed.

        :param dest: If `None`, the agent will not send any data; otherwise it
                     will send data to the desired `dest`.
        :param plen: Packet length for data to send.
        :param delay: Average delay to wait before sending a new packet.
        :param mode: Operation mode ("cbr", "poisson", or "backlog")
                     [default="backlog"].
        """
        # set up parameters
        if plen<1: plen = self.DefaultPacketLength
        self._seqno = 0
        self.dest = dest
        self.delay = delay
        self.mode = mode
        # initialize random data in payload buffer
        self.plen = plen
        x = np.random.randint(0, 128, self.plen*2)
        args = tuple(x.tolist())
        self._payload = struct.pack("b"*len(x), *args)
        # set up ports and start FSM's
        self.addport("TX"), self.addport("RX")
        f = self.newchild("txfsm", FSM, tracename=self.tracename+".TX")
        g = self.newchild("rxfsm", FSM, tracename=self.tracename+".RX")
        f.goto(self.INIT)
        g.goto(self.RECV)

    def connect(self, net):
        """Connect to a network layer protocol."""
        self.TX.connect(net.RXU)
        net.TXU.connect(self.RX)

    def INIT(self, fsm):
        """INIT state; *overload to implement initialization before `SEND`*."""
        yield fsm.goto(self.SEND)

    def SEND(self, fsm):
        """SEND state; keep backlog of data in MAC RXU queue."""
        yield hold, fsm, const.EPSILON         # yield to other threads
        if (not self.dest): yield fsm.stop()    # HALT and sleep
        # check parameters
        if (self.mode=="backlog") and self.mac:
            macQ = self.mac.RXU
            macQ.monitorQ = True
        # create packet to send downstream
        idx = self.seqno%self.plen
        pay = self._payload[idx:idx+self.plen]
        p = AGT(id=self.uid, seqno=self.seqno)/pay
        p.setanno('net-dst', self.dest)
        # set agent annotations
        txts = now()
        p.setanno('agt-txts', txts)
        p.setanno('agt-plen', self.plen)
        # additional send processing
        pkt = self.senddata(p)
        # send packet
        self.nsent += 1
        self.log("snd", pkt, seqno=self.seqno, nsent=self.nsent)
        self.nextseqno(update=True)       # update counter
        self.send.signal(pkt)
        yield self.TX.send(fsm, [pkt])
        # pause agent
        yield fsm.goto(self.PAUSE)

    def PAUSE(self, fsm):
        """PAUSE state; pause agent based on operation mode."""
        yield hold, fsm, 0  # yield to other threads
        if (not self.dest): yield fsm.stop()    # HALT and sleep
        # wait for MAC to empty in 'backlog' mode
        p, moredata = None, False        # need to send more data
        if (self.mode=="backlog") and self.mac:
            macQ = self.mac.RXU
            # check MAC input queue
            while (not moredata):
                yield waitevent, fsm, macQ.deQ
                p = macQ.deQ.signalparam
                if (macQ.length<1):
                    moredata = True
        # wait for prescribed delay
        d = None
        if self.delay:
            # pause based on operation mode
            # -> 'poisson' uses exponentail backoff
            # -> otherwise use fixed delay
            if (self.mode=="poisson"): d = np.random.exponential(self.delay)
            else:                      d = self.delay
        # enforce minimum wait time
        d = max(d, const.EPSILON)
        self.log("WAIT", p, delay=d, avgdelay=self.delay)
        yield hold, fsm, d
        # goto CONTINUE
        yield fsm.goto(self.CONTINUE)

    def CONTINUE(self, fsm):
        """CONTINUE state; add processing after PAUSE before next SEND."""
        yield fsm.goto(self.SEND)

    def RECV(self, fsm):
        """RECV state; Receive traffic."""
        # get packet from lower layer
        yield self.RX.recv(fsm, 1)
        p = fsm.got[0]
        self.recv.signal(p)
        # check packet parameters
        assert isinstance(p, AGT)
        if self.addr and (not AGENT_ID_INDEX):
            # use addressing for indexing into recvbuffer
            assert p.hasanno('net-dst') and p.hasanno('net-src')
            src, dst = p.getanno('net-src'), p.getanno('net-dst')
            isforme = (dst==self.addr) or (dst==self.broadcast)
            errmsg = "[AGT]: Got packet not intended for me!"
            assert (isforme), errmsg
            index = src
        else:
            # use AGT.id for indexing into recvbuffer
            index = int(p.id)
        # set agent annotations
        rxts = now()
        latency = rxts - p.getanno('agt-txts')
        p.setanno('agt-rxts', rxts)
        p.setanno('agt-latency', latency)
        # check for duplicates -> drop DUP
        key = ID, seqno = int(p.id), int(p.seqno)
        if index not in self.recvbuffer: self.recvbuffer[index] = []
        if key in self.recvbuffer[index]:
            self.recvdup(p)
            self.log("dup", p)
        else:
            # log receive
            self.nrcvd += 1
            self.recvbuffer[index].append(key)
            pkt = self.recvdata(p)
            self.log("rcv", pkt)
        # trim receive buffer
        while (len(self.recvbuffer[index])>self.RecvBufferSize):
            x = self.recvbuffer[index].pop(0)
        # conintue in RECV
        yield fsm.goto(self.RECV)

    def senddata(self, p):
        """Process received data packets as desired.

        :return: Data packet to send.

        By default, this method does nothing; *overload as needed*.
        """
        return p

    def recvdup(self, p):
        """Process duplicate packets as desired.

        By default, this method does nothing; *overload as needed*.
        """
        return p

    def recvdata(self, p):
        """Process received data packets as desired.

        :return: Processed packet.

        By default, this method does nothing; *overload as needed*.
        """
        return p

    def nextseqno(self, update=False):
        """Increment sequence number."""
        newseqno = (self.seqno+1)%self.MaxSeqNo
        if update: self._seqno = newseqno
        return newseqno

    def get_mode(self):
        """Get current traffic generating mode."""
        mode = self.DefaultTrafficMode
        if self._mode in self.TrafficModes: mode = self._mode
        return mode

    def set_mode(self, mode):
        """Set current traffic generating mode."""
        ismode = isinstance(mode, str) and (mode.lower() in self.TrafficModes)
        if ismode: mode = mode.lower()
        else:      mode = self.DefaultTrafficMode
        self._mode = mode

    def get_mac(self):
        """Return MAC of parent node."""
        mac = None
        if isinstance(self.parent, Base):
            if self.parent.haschild('mac'):
                mac = self.parent.mac
        if not isinstance(mac, MAC):
            mac = None
        return mac

    def get_address(self, broadcast=False):
        """Return network address of parent node."""
        addr, net = None, None
        if isinstance(self.parent, Base):
            if self.parent.haschild('net'):
                net = self.parent.net
        if isinstance(net, NET):
            addr = net.address
            if broadcast: addr = net.broadcast
        return addr

    def log(self, event=None, p=None, *args, **kwargs):
        """Overloaded to check verbose level and set common annotations."""
        force = False
        if ('verbose' in kwargs): force = (kwargs['verbose']>AGENT_VERBOSE)
        if self.verbose>AGENT_VERBOSE or force:
            kwargs.update(self.get_agt_anno(p))
            Element.log(self, event, p, *args, **kwargs)

    def get_agt_anno(self, p):
        """Get common annotations for logging data."""
        kwargs = {}
        kwargs['mode'] = self.mode
        if self.addr: kwargs['addr'] = self.addr
        if not isinstance(p, AGT): return kwargs
        kwargs['id'] = p.id
        kwargs['seqno'] = p.seqno
        kwargs['nsent'] = self.nsent
        kwargs['nrcvd'] = self.nrcvd
        if p.hasanno('cif-collision'):
            kwargs['cif-collision'] = strcollision(p)
            assert (kwargs['cif-collision'] is not None)
        if p.hasanno('dot11n-sinr'):
            kwargs['dot11n-sinr'] = "%.4f dB"%(p.getanno('dot11n-sinr') )
        if p.hasanno('dot11n-channel-fading'):
            kwargs['dot11n-channel-fading'] = "%.4f dB"%(p.getanno('dot11n-channel-fading') )
        if p.hasanno('phy-sinr'):
            kwargs['sinr'] = "%.4f dB"%(p.getanno('phy-sinr') )
        if p.hasanno('phy-rate'):
            kwargs['phy-rate'] = p.getanno('phy-rate')
        if p.hasanno('net-dst'):
            kwargs['net-dst'] = p.getanno('net-dst')
        if p.hasanno('net-src'):
            kwargs['net-src'] = p.getanno('net-src')
        if p.hasanno('net-root'):
            kwargs['net-root'] = p.getanno('net-root')
        if p.hasanno('mac-root'):
            kwargs['mac-root'] = p.getanno('mac-root')
        if p.hasanno('mac-txts'):
            kwargs['mac-txts'] = p.getanno('mac-txts')
        if p.hasanno('mac-rxts'):
            kwargs['mac-rxts'] = p.getanno('mac-rxts')
        if p.hasanno('agt-txts'):
            kwargs['agt-txts'] = p.getanno('agt-txts')
        if p.hasanno('agt-rxts'):
            kwargs['agt-rxts'] = p.getanno('agt-rxts')
        if p.hasanno('agt-plen'):
            kwargs['agt-plen'] = p.getanno('agt-plen')
        if p.hasanno('agt-latency'):
            kwargs['agt-latency'] = p.getanno('agt-latency')
        return kwargs
