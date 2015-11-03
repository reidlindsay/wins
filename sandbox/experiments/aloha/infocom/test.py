#!  /usr/bin/env python

"""
Simulate ALOHA over a network of nodes.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-10-19 17:04:02 -0500 (Wed, 19 Oct 2011) $
* $LastChangedRevision: 5220 $

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
"""
__docformat__ = "restructuredtext en"

from SimPy.Simulation import *
from scapy.all import *
from wins import *
from wins.ieee80211 import *
from copy import copy, deepcopy

from wins.backend import RNG_init

import sys
from optparse import OptionParser

import numpy as np
import struct
import gc
import time

RNG_INIT = 1

class AGT(Packet):
    fields_desc = [IntField("id", 0), \
                   IntField("seqno", 0)]
    def __init__(self, *args, **kwargs):
        Packet.__init__(self, *args, **kwargs)

class Agent(Element):
    """Agent to act as source and sink.

    This module creates a Poisson traffic generator.

    :ivar dest: IP address of destination.
    :ivar plen: Packet length (in bytes).
    :ivar rate: Send rate (in packets-per-second).
    :ivar nsent: Number of packets sent by Agent.
    :ivar nrcvd: Number of packets received by Agent.
    """
    name="agent"
    tracename="AGT"
    PACKETLEN = 1024
    def __init__(self, **kwargs):
        self.nsent, self.nrcvd = 0, 0
        self.netaddr, self.dest, self.addr = None, None, None
        self._plen = None
        self.duration = None
        self.payload = None
        self.seqno, self.maxseqno = 0, 4096
        Element.__init__(self, **kwargs)

    plen = property(fget=lambda self: self._plen)

    def configure(self, rate=None, dest=None, datarate=None, \
                        plen=None, **kwargs):
        """Configure Agent.

        :param rate: Rate for packet generation (in transmission per packet time).
        :param dest: Network address for destination node.
        :param datarate: Rate index for PHY.
        """
        if plen is None: plen=self.PACKETLEN
        if rate is None: rate=1.0
        self.rate = rate
        self.dest = dest
        self.datarate = datarate
        self._plen = plen - len(IP()/AGT())
        # set up ports and start FSM's
        self.addport("TX"), self.addport("RX")
        f = self.newchild("txfsm", FSM, tracename=self.tracename+".TX")
        g = self.newchild("rxfsm", FSM, tracename=self.tracename+".RX")
        f.goto(self.SEND)
        g.goto(self.RECV)
    def connect(self, net):
        """Connect to a network layer protocol."""
        self.TX.connect(net.RXU)
        net.TXU.connect(self.RX)
    def SEND(self, fsm):
        """SEND state; generate Poisson traffic."""
        if not (self.rate>0): yield fsm.stop()
        if self.dest is None: yield fsm.stop()
        # set parameters
        if self.duration is None:
            phy = self.parent.phy
            plen = len(IP()/AGT()) + self.plen
            self.duration = phy.calcduration(plen, self.datarate)
        if self.payload is None:
            x = np.arange(self.plen)%128
            args = tuple(x.tolist())
            self.payload = struct.pack("b"*self.plen, *args)
        # find holding time
        avgdelay = self.duration/self.rate
        delay = np.random.exponential(avgdelay)
        self.log("wait", delay=delay, duration="%.2f usec"%(self.duration*1e6))
        yield hold, fsm, delay
        # create new packet
        idx = self.seqno%len(self.payload)
        payload = self.payload[idx:] + self.payload[:idx]
        p = AGT(id=self.uid, seqno=self.seqno)/payload
        assert len(payload) == self.plen
        p.setanno('phy-rate', self.datarate)
        p.setanno('net-dst',  self.dest)
        # send packet
        self.nsent += 1
        kwargs = {'seqno':self.seqno, 'nsent':self.nsent}
        self.log("snd", p, **kwargs)
        yield self.TX.send(fsm, [p])
        # update counters
        self.seqno = (self.seqno+1)%self.maxseqno
        # continue in SEND
        yield fsm.goto(self.SEND)
    def RECV(self, fsm):
        """RECV state; Receive traffic."""
        # set parameters
        if self.addr is None:
            self.addr = self.parent.net.address
        # get packet
        yield self.RX.recv(fsm, 1)
        p = fsm.got[0]
        assert isinstance(p, AGT)
        assert (p.getanno('net-dst')==self.addr)
        # log receive
        self.nrcvd += 1
        kwargs = {'id':p.id, 'seqno':p.seqno, 'nrcvd':self.nrcvd}
        if p.hasanno('cif-collision'):
            kwargs['cif-collision'] = [c for c in p.getanno('cif-collision')]
        if p.hasanno('dot11n-sinr'):
            kwargs['sinr'] = p.getanno('dot11n-sinr')
        self.log("rcv", p, **kwargs)
        # conintue in RECV
        yield fsm.goto(self.RECV)

class Node(Element):
    name = "node"
    tracename = "NODE"
    def __init__(self, **kwargs):
        Element.__init__(self, **kwargs)
    def configure(self, rate=None, pos=None, datarate=None, \
                        dest=None, useshared=False, cfocorrection=True, **kwargs):
        cif = self.newchild('cif', Dot11NRadio)
        phy = self.newchild('phy', Dot11NPHY, radio=cif, cfocorrection=cfocorrection)
        mac = self.newchild('mac', Aloha, phy=phy)
        net = self.newchild('net', Routing)
        arp = self.newchild('arp', ARP, useshared=useshared)
        agt = self.newchild('agent', Agent, rate=rate, dest=dest, datarate=datarate)
        mobi = self.newchild('motion', Motion, pos=pos)
        # connect ports
        agt.connect(net)
        arp.connect(net, mac)
        mac.connect(phy)
        phy.connect(cif)

def read_topo(options, topofile):
    f = file(topofile, 'r')
    s = f.readline()
    topo = {'border':None, 'layout':None}
    done = not (s)
    while not done:
        # covert s to dict (check for border and layout)
        try:
            d = eval(s)
            assert isinstance(d, dict)
            assert ('border' in d) and ('layout' in d)
        except:
            d = None
        # add dict to datain
        if d: topo = d
        # get next input
        s = f.readline()
        done = not(s)
    return topo

def run_experiment(options):
    # record start time
    starttime = time.time()
    # initialize RNG
    if RNG_INIT: RNG_init()
    # experiment parameters
    ntx, nrx = 1, 1
    datarate = 0
    load = options.load     # total transmissions per packet time
    numnodes = options.numnodes
    if options.plen>0: Agent.PACKETLEN = options.plen
    if 0<=options.mcs<8*ntx: datarate = options.mcs
    # set other parameters
    Dot11NPHY.usewaveform = options.usewaveform
    Dot11NRadio.Ntx, Dot11NRadio.Nrx = ntx, nrx
    Dot11NRadio.fomax = options.fomax
    cfocorrection = True
    if options.disable_cfo_correction: cfocorrection = False
    # set up parameters
    initialize()
    stoptime = 2.0
    if options.stop>0: stoptime = options.stop
    stoptime *= 1.05     # allow events around stoptime to finish
    verbose = 20
    useshared = True
    kwargs = {'verbose':verbose, 'useshared':useshared}
    nodeargs = {'datarate':datarate, 'cfocorrection':cfocorrection}
    nodeargs.update(kwargs)
    # set up channel
    ch = Channel(model=Dot11NChannel, **kwargs)
    modeltype  = options.tgnmodel       # default -> LOS Channel
    alpha      = options.alpha
    usedoppler = options.usedoppler
    usefading  = options.usefading
    chargs = {'modeltype':modeltype, 'n':alpha, \
              'usedoppler':usedoppler, 'usefading':usefading}
    # calculate maximum network diameter
    node = Node(rate=0)
    phy = node.phy
    snr, plen = 0.0, Agent.PACKETLEN
    per = phy.calcper(plen, datarate, snr)
    per_floor = 0.01
    while per>per_floor:
        snr += 3
        per = phy.calcper(plen, datarate, snr)
    snr += 3
    Ptx = DOT11N_MAXPOWER
    Lrx, Ltx = Dot11NRadio.Lrx, Dot11NRadio.Ltx
    Grx, Gtx = Dot11NRadio.Grx, Dot11NRadio.Gtx
    No = Dot11NRadio.thermalnoise(DOT11N_BANDWIDTH) + DOT11N_NOISEFIGURE
    Prx = Ptx + Gtx + Grx - Ltx - Lrx
    refdist = Dot11NChannel.bpdist       # use breakpoint distance as refdist
    refloss = Propagation.freespace(refdist, n=2.0, fc=Dot11NChannel.fc)
    # use reference pathloss model (with BPdist as RefDist)
    PL = Prx - No - refloss    # pathloss normalized to refloss
    dnetwork = refdist*(db2linear(PL-snr)**(1.0/alpha))
    # set (or create) topology
    if options.usetopo:
        topofile = options.usetopo
        topo = read_topo(options, topofile)
        border = topo['border']
        layout = topo['layout']
        xmin, xmax, ymin, ymax = border[0:4]
    else:
        # create new uniform layout
        assert (options.xmin<=options.xmax)
        assert (options.ymin<=options.ymax)
        xmin, xmax = options.xmin, options.xmax
        ymin, ymax = options.ymin, options.ymax
        # use dnetwork as network diameter
        beta = 0.95
        xmax = xmin + beta*dnetwork/np.sqrt(2)
        ymax = ymin + beta*dnetwork/np.sqrt(2)
        border = (xmin, xmax, ymin, ymax)
        xpos = xmin + np.random.uniform(0,1,numnodes)*(xmax-xmin)
        ypos = ymin + np.random.uniform(0,1,numnodes)*(ymax-ymin)
        layout = [(xpos[k],ypos[k]) for k in range(numnodes)]
    # verify layout parameters
    assert (len(layout)==numnodes)
    d2 = (xmax-xmin)**2 + (ymax-ymin)**2
    assert (dnetwork**2>d2), "boundary = %s, dnetwork = %s"%(border, dnetwork)
    topo   = {'border':border, 'layout':layout}
    # set up receviers
    rlist, nlinks = [], int(numnodes/2)
    for k in range(nlinks):
        pos  = layout[nlinks+k]
        rnode = Node(rate=0.0,  pos=pos, **nodeargs)
        rlist.append(rnode)
    # set up senders
    slist, rate = [], load/nlinks      # load per sender
    for k in range(nlinks):
        dest = rlist[k].net.address
        pos  = layout[k]
        snode = Node(rate=rate, pos=pos, dest=dest, **nodeargs)
        slist.append(snode)
    # connect all nodes and add routes
    for n in (slist+rlist):
        for m in (slist+rlist):
            if (n is not m):
                ch.add_edge(n.cif, m.cif, modeltype=modeltype, n=alpha)
                n.net.addroute(m.net.address)
    # start monitor
    if options.monitor:
        mon = Monitor(period=stoptime/1e4)
        mon.start()

    # run simulation - log relavent information
    ch.log("load", G="%.5g"%(load) )
    simulate(until=stoptime)
    ch.log("stoptime", stoptime="%.6f"%(stoptime))
    n = gc.collect()
    ch.log("GC", collected=n)
    totaltime = time.time() - starttime
    t = time.gmtime(totaltime)
    ch.log("runtime", runtime="%02d:%02d:%02d (h/m/s)"%(t.tm_hour, t.tm_min, t.tm_sec) )

    # print output
    sys.stdout.flush()
    if options.verbose: ch.trace.output()

    # write tracefile
    if options.output is not None: ch.trace.write(options.output)

    # write topofile
    if options.savetopo:
        f = file(options.savetopo, 'w')
        f.write("%s\n"%(topo) )
        f.close()

def main():
    usage = "%prog [OPTIONS]"
    parser = OptionParser(usage=usage)
    # simulation parameters
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true",  \
            default=False, help="Output formatted trace to stdout")
    parser.add_option("-o", "--output", dest="output", \
            default=None, help="Name of output file for trace")
    parser.add_option("-s", "--stop", dest="stop", \
            type="float", default=2.0, \
            help="Run simulation until stop time [default=%default]")
    parser.add_option("-m", "--monitor", dest="monitor", action="store_true",  \
            default=False, help="Enable simulation montior")
    # protocol parameters
    parser.add_option("-G", "--load", dest="load", type="float", \
            default=1.0, help="Set network traffic load in total " + \
                              "transmissions per packet time [default=%default]")
    parser.add_option("-l", "--packet-length", dest="plen", type="int", \
            default=1024, help="Set packet size in bytes [default=%default]")
    parser.add_option("-n", "--num-nodes", dest="numnodes", type="int", \
            default=50, help="Set number of nodes [default=%default]")
    parser.add_option("", "--mcs", dest="mcs", type="int", \
            default=0, help="Set rate index for MCS [default=%default]")
    parser.add_option("", "--fomax", dest="fomax", \
            type="float", default=0.0, \
            help="Specify maximum frequency offset in ppm [default=%default]")
    parser.add_option("", "--use-waveform", dest="usewaveform", action="store_true", \
            default=False, help="Enable waveform-level simulation [default=%default]")
    parser.add_option("", "--disable-cfo-correction", \
            dest="disable_cfo_correction", action="store_true", \
            default=False, help="Disable CFO correction in waveform-level simulation [default=%default]")
    # channel parameters
    parser.add_option("", "--tgn-model", dest="tgnmodel", \
            default=None, help="Specify TGn model.")
    parser.add_option("", "--alpha", dest="alpha", type="float", \
            default=2.0, help="Specify pathloss exponent [default=%default]")
    parser.add_option("", "--use-doppler", dest="usedoppler", action="store_true",  \
            default=False, help="Enable doppler filter for fading in TGn channel model.")
    parser.add_option("", "--disable-fading", dest="usefading", action="store_false",  \
            default=True, help="Normalize channel and remove impact of fading on pathloss in TGn channel model.")
    # topology/layout parameters
    parser.add_option("", "--xmin", dest="xmin", type="float", \
            default=0.0, help="Set x-axis left boundary [default=%default]")
    parser.add_option("", "--xmax", dest="xmax", type="float", \
            default=500.0, help="Set x-axis right boundary [default=%default]")
    parser.add_option("", "--ymin", dest="ymin", type="float", \
            default=0.0, help="Set y-axis lower boundary [default=%default]")
    parser.add_option("", "--ymax", dest="ymax", type="float", \
            default=500.0, help="Set y-axis upper boundary [default=%default]")
    parser.add_option("", "--use-topo", dest="usetopo", \
            default=None, help="Specify topology file instead of generating random topology.")
    parser.add_option("", "--save-topo", dest="savetopo", \
            help="Save topology to file.")
    (options, args) = parser.parse_args()

    if len(args)>0:
        print "Invalid number of arguments."
        parser.print_help()
        raise SystemExit

    run_experiment(options)

if __name__ == '__main__':
    main()
