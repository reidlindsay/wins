#!  /usr/bin/env python

"""
Master Controller and Agent for DSR simulation.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-12-20 23:57:48 -0600 (Tue, 20 Dec 2011) $
* $LastChangedRevision: 5393 $

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

from SimPy.Simulation import *
from scapy.all import *
from wins import *
from wins.ieee80211 import *
from copy import copy, deepcopy

from wins.backend import RNG_init
from wins.traffic import Agent
from wins.mac import RBAR, ARF
from wins.net import DSR

import sys
from optparse import OptionParser

import numpy as np
import struct
import gc
import time

class MyAgent(Agent):
    name = "my agent"
    def __init__(self, *args, **kwargs):
        self.newdata = SimEvent()
        self.sndnext = SimEvent()
        self.rcvdata = SimEvent()
        Agent.__init__(self, *args, **kwargs)
    def configure(self, **kwargs):
        Agent.configure(self, **kwargs)
    def INIT(self, fsm):
        """INIT; wait in CONTINUE before starting."""
        yield fsm.goto(self.CONTINUE)
    def CONTINUE(self, fsm):
        """CONTINUE; wait on `sndnext` before returning to CONTINUE."""
        yield waitevent, fsm, self.sndnext
        yield fsm.goto(self.SEND)
    def senddata(self, p):
        """Signal new data packet before sending."""
        self.newdata.signal(self.dest)
        return Agent.senddata(self, p)
    def recvdata(self, p):
        """Signal received data packet."""
        self.rcvdata.signal()
        return Agent.recvdata(self, p)

class MasterControl(Element):
    name = "master control"
    tracename = "MC"
    DLIMIT = 10
    TLIMIT = 10
    TIMEOUT = 35
    def __init__(self, src, dst, netaddr, chan, border, **kwargs):
        self.src = src
        self.dst = dst
        self.netaddr = netaddr
        self.chan = chan
        self.border = border
        # set up counters
        self.nsent, self.nrcvd = None, None
        self.ntopo = 0     # number of topologies simulated
        Element.__init__(self, **kwargs)
    def configure(self, **kwargs):
        f = self.newchild('fsm', FSM)
        f.goto(self.INIT)
    def set_layout(self):
        """Set random layout for non-src/non-dst nodes."""
        nodes = [n for n in self.netaddr.values() \
                    if ((n is not self.src) and (n is not self.dst))]
        xmin, xmax, ymin, ymax = self.border[0:4]
        self.src.motion.position = (xmin, ymin)
        #self.log("POS", addr=self.src.net.address, pos=self.src.motion.position)
        for n in nodes:
            x = np.random.uniform(xmin, xmax)
            y = np.random.uniform(ymin, ymax)
            n.motion.position = (x,y)
            #self.log("POS", addr=n.net.address, pos=n.motion.position)
        self.dst.motion.position = (xmax, ymax)
        #self.log("POS", addr=self.dst.net.address, pos=self.dst.motion.position)

    def INIT(self, fsm):
        """INIT state."""
        ##fsm.verbose = 150
        # remove all routes to destination
        dstaddr = self.dst.net.address
        self.src.net.delroute(dstaddr, reset=True)
        # randomize layout
        self.set_layout()
        # connect all nodes
        nodelist = self.netaddr.values()
        self.connect(nodelist)
        self.ntopo += 1
        self.nsent = self.src.agent.nsent
        self.nrcvd = self.dst.agent.nrcvd
        yield fsm.goto(self.SNDNEXT)

    def SNDNEXT(self, fsm):
        """SNDNEXT; signal for next data packet to be sent."""
        # wait for delay
        delay = self.src.agent.delay
        yield hold, fsm, delay
        # send next packet from agent
        nrcvd = self.dst.agent.nrcvd - self.nrcvd
        nsent = self.src.agent.nsent - self.nsent
        if nsent<self.DLIMIT:
            self.log("SND", nsent=nsent, nrcvd=nrcvd)
            self.src.agent.sndnext.signal()
            yield fsm.goto(self.WAIT)
        else:
            yield fsm.goto(self.DONE)

    def WAIT(self, fsm):
        """WAIT; wait on newdata from src."""
        newdata = self.src.agent.newdata
        yield waitevent, fsm, newdata
        dstaddr = newdata.signalparam
        yield fsm.goto(self.DATA, dstaddr)

    def DATA(self, fsm, dstaddr):
        """DATA; set routing for data to dstaddr."""
        srcaddr = self.src.net.address
        if self.src.net.hasroute(dstaddr):
            # connect only nodes on route
            self.disconnect()
            srcroute = self.src.net.srcroute(dstaddr)
            path = [srcaddr] + srcroute + [dstaddr]
            nodes = []
            for p in path: nodes.append(self.netaddr[p])
            self.connect(nodes)
            yield fsm.goto(self.MON, nodes, dstaddr)
        else:
            # connect everything
            nodelist = self.netaddr.values()
            self.connect(nodelist)
            yield fsm.goto(self.RREQ)

    def RREQ(self, fsm):
        """RREQ; wait for outcome of RREQ."""
        drprreq = self.src.net.drprreq
        finrreq = self.src.net.finrreq
        yield waitevent, fsm, (drprreq, finrreq)
        if drprreq in fsm.eventsFired:
            yield fsm.goto(self.DONE)
        elif finrreq in fsm.eventsFired:
            nodelist = self.netaddr.values()
            dstaddr = finrreq.signalparam
            yield fsm.goto(self.MON, nodelist, dstaddr)

    def MON(self, fsm, nodes, dstaddr):
        """MON; monitor nodes along an active route."""
        flags = [n.net.sndfail for n in nodes]
        rcvdata = self.dst.agent.rcvdata
        yield waitevent, fsm, (flags+[rcvdata])
        if rcvdata in fsm.eventsFired:
            # success -> send next packet
            yield fsm.goto(self.SNDNEXT)
        # link is broken -> reconnect all
        nodelist = self.netaddr.values()
        self.connect(nodelist)
        # wait for timeout and then continue
        yield hold, fsm, self.TIMEOUT
        yield fsm.goto(self.SNDNEXT)

    def DONE(self, fsm):
        """DONE; log data for current topology."""
        nrcvd = self.dst.agent.nrcvd - self.nrcvd
        nsent = self.src.agent.nsent - self.nsent
        assert (nsent>0)
        pdr = 1.0*nrcvd/nsent
        nodelist = self.netaddr.values()
        numnodes = len(nodelist)
        self.log("PDR", pdr=pdr, nsent=nsent, topo=self.ntopo, numnodes=numnodes)
        # check if finished
        if self.ntopo<self.TLIMIT:
            yield fsm.goto(self.INIT)
        else:
            self.src.agent.dest = None
            # stop simulation
            stopSimulation()

    def disconnect(self):
        """Disconnect everything."""
        nodelist = self.netaddr.values()
        self.connect(nodelist, disable=True)

    def connect(self, nodelist, disable=False):
        """Connect selected nodes."""
        for n in nodelist:
            for m in nodelist:
                if (n is not m):
                    cm = self.chan.get_edge(n.cif,m.cif)
                    if disable: cm.disable()
                    else:       cm.enable()
        #self.log("CONNECT", nodes=[n.net.address for n in nodelist if (not disable)])
