#!  /usr/bin/env python

"""
Simulate point-to-point links.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-10-18 12:12:13 -0500 (Tue, 18 Oct 2011) $
* $LastChangedRevision: 5213 $

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
#from wins.traffic import Agent

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
        self.dest = None
        self._plen = None
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
        # adjust packet length to account for header overhead
        self._plen = plen - len(Ether()/IP()/AGT())
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
        if self.dest is None:
            self.dest = self.parent.net.broadcast
        # set parameters
        if self.payload is None:
            x = np.arange(self.plen)%128
            args = tuple(x.tolist())
            self.payload = struct.pack("b"*self.plen, *args)
        # hold for constant wait time (1.0/rate)
        delay = 1.0/self.rate
        self.log("wait", delay=delay)
        yield hold, fsm, delay
        # create new packet
        idx = self.seqno%len(self.payload)
        payload = self.payload[idx:] + self.payload[:idx]
        p = AGT(id=self.uid, seqno=self.seqno)/payload
        assert len(payload) == self.plen
        p.setanno('phy-rate', self.datarate)    # set MCS
        p.setanno('net-dst',  self.dest)        # set destination addr
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
        addr  = self.parent.net.address
        bcast = self.parent.net.broadcast
        # get packet
        yield self.RX.recv(fsm, 1)
        p = fsm.got[0]
        assert isinstance(p, AGT)
        netdst = p.getanno('net-dst')
        assert ((netdst==addr) or (netdst==bcast))
        # log receive
        self.nrcvd += 1
        kwargs = {'id':p.id, 'seqno':p.seqno, 'nrcvd':self.nrcvd}
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
        agt = self.newchild('agent', Agent, rate=rate, \
                                     dest=dest, datarate=datarate)
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
        # add dict to topo
        if d: topo.update(d)
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
    if options.plen>0: Agent.PACKETLEN = options.plen
    if 0<=options.mcs<8*ntx: datarate = options.mcs
    # set other parameters
    Dot11NPHY.usewaveform = options.usewaveform
    Dot11NRadio.Ntx, Dot11NRadio.Nrx = ntx, nrx
    Dot11NRadio.fomax = options.fomax
    cfocorrection = True
    if options.disable_cfo_correction: cfocorrection = False
    # set simulation parameters
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
    modeltype  = options.tgnmodel
    alpha      = options.alpha
    usedoppler = options.usedoppler
    usefading  = options.usefading
    chargs = {'modeltype':modeltype, 'n':alpha, \
              'usedoppler':usedoppler, 'usefading':usefading}
    # set (or create) topology
    if options.usetopo:
        topofile = options.usetopo
        topo = read_topo(options, topofile)
        border = topo['border']
        layout = topo['layout']
        numnodes = len(layout)
    else:
        # create line topology based on snr parameters
        snrmin, snrmax = options.snrmin, options.snrmax
        dsnr = options.snrstep
        Ptx = DOT11N_MAXPOWER
        Lrx, Ltx = Dot11NRadio.Lrx, Dot11NRadio.Ltx
        Grx, Gtx = Dot11NRadio.Grx, Dot11NRadio.Gtx
        No = Dot11NRadio.thermalnoise(DOT11N_BANDWIDTH) + DOT11N_NOISEFIGURE
        Prx = Ptx + Gtx + Grx - Ltx - Lrx
        refdist = Dot11NChannel.bpdist       # use breakpoint distance as refdist
        refloss = Propagation.freespace(refdist, n=2.0, fc=Dot11NChannel.fc)
        # use reference pathloss model (with BPdist as RefDist)
        PL = Prx - No - refloss    # pathloss normalized to refloss
        snr = numpy.arange(snrmin, snrmax+const.EPSILON, dsnr)
        snrlist, pos = snr.tolist(), [(0,0)]    # initial point is sender
        snrlist.reverse()
        for s in snrlist:
            x = refdist*(db2linear(PL-s)**(1.0/alpha))
            pos.append((x,0))
        numnodes = len(snrlist) + 1
        # create topology layout
        xmin = min([x for (x,y) in pos])
        xmax = max([x for (x,y) in pos])
        ymin = min([y for (x,y) in pos])
        ymax = max([y for (x,y) in pos])
        border = (xmin, xmax, ymin, ymax)
        layout = pos
    assert (len(layout)==numnodes)
    topo   = {'border':border, 'layout':layout}
    # set up sender
    dest = NET.broadcast
    rate, pos = 1.0, layout[0]
    sender = Node(rate=rate, pos=pos, dest=dest, **nodeargs)
    # set up remaining nodes as receviers
    nlist = [sender]
    for k in range(1, numnodes):
        rate, pos  = 0.0, layout[k]
        rnode = Node(rate=rate, pos=pos, **nodeargs)
        nlist.append(rnode)
    # connect sender to all receivers and add routes
    for n in nlist[1:]:
        ch.add_edge(sender.cif, n.cif, **chargs)
        sender.net.addroute(n.net.address)      # add route sender->receiver

    """
    # print channel model (to verify flags are being set correctly)
    p = AGT()/"helloworld"
    cm = ch.get_edge(sender.cif, nlist[1].cif)
    cm.apply(p, sender.cif, nlist[1].cif)
    print "Channel(%s, %s) = %s"%(sender.traceid, nlist[1].traceid, \
                                  ch.get_edge(sender.cif, nlist[1].cif).chan.model())
    """
    

    # run simulation - log relavent information
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
    # protocol parameters
    parser.add_option("", "--mcs", dest="mcs", type="int", \
            default=0, help="Set rate index for MCS [default=%default]")
    parser.add_option("-l", "--packet-length", dest="plen", type="int", \
            default=1500, help="Set packet size in bytes [default=%default]")
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
    parser.add_option("", "--use-topo", dest="usetopo", \
            default=None, help="Specify topology file instead of generating random topology.")
    parser.add_option("", "--save-topo", dest="savetopo", \
            help="Save topology to file.")
    # experiment range parameters
    parser.add_option("", "--snrmin", dest="snrmin", type="float", \
            default=0, help="Set minimum SNR value for simulation [default=%default]")
    parser.add_option("", "--snrmax", dest="snrmax", type="float", \
            default=30, help="Set maximum SNR value for simulation [default=%default]")
    parser.add_option("", "--snrstep", dest="snrstep", type="float", \
            default=5, help="Set SNR step size for simulation [default=%default]")
    (options, args) = parser.parse_args()

    if len(args)>0:
        print "Invalid number of arguments."
        parser.print_help()
        raise SystemExit

    run_experiment(options)

if __name__ == '__main__':
    main()
