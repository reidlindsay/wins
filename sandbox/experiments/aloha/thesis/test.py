#!  /usr/bin/env python

"""
Simulate ALOHA network.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-10-21 15:20:36 -0500 (Fri, 21 Oct 2011) $
* $LastChangedRevision: 5245 $

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

import sys
from optparse import OptionParser

import numpy as np
import struct
import gc
import time

RNG_INIT = 1
EXIT_WITH_TRACE = 0

class MyAgent(Agent):
    name = "my agent"
    def __init__(self, *args, **kwargs):
        Agent.__init__(self, *args, **kwargs)
    def configure(self, datarate=None, **kwargs):
        self.datarate = datarate
        Agent.configure(self, **kwargs)
        self.txfsm.goto(self.PAUSE)     # start in wait state
    def senddata(self, p):
        if not ANNO.supported(p): return p
        if (self.datarate is not None):
            p.setanno('phy-rate', self.datarate)
        return p
    def get_agt_anno(self, p):
        kwargs = Agent.get_agt_anno(self, p)
        if not isinstance(p, AGT): return kwargs
        if p.hasanno('dot11n-sinr'):
            kwargs['sinr'] = p.getanno('dot11n-sinr')
        return kwargs

class Node(Element):
    name = "node"
    tracename = "NODE"
    def __init__(self, **kwargs):
        Element.__init__(self, **kwargs)
    def configure(self, pos=None,                       # motion    \
                        useshared=True,                 # arp       \
                        cfocorrection=True,             # phy       \
                        datarate=None,                  # myagent   \
                        dest=None, plen=None, delay=None, mode=None, # agent \
                        **kwargs):
        cif = self.newchild('cif', Dot11NRadio)
        phy = self.newchild('phy', Dot11NPHY, radio=cif, cfocorrection=cfocorrection)
        mac = self.newchild('mac', Aloha, phy=phy)
        net = self.newchild('net', Routing)
        arp = self.newchild('arp', ARP, useshared=useshared)
        agt = self.newchild('agent', MyAgent, dest=dest, plen=plen, \
                                     delay=delay, mode=mode, datarate=datarate)
        mobi = self.newchild('motion', Motion, pos=pos)
        # connect ports
        agt.connect(net)
        arp.connect(net, mac)
        mac.connect(phy)
        phy.connect(cif)

USE_TRAFFIC_MONITORS = 0

class TrafficMonitor(Element):
    """Monitor traffic to see if sufficient number of packets has been
    received; allows simulation to exit early."""
    name = "traffic monitor"
    tracename = "TMON"
    Limit = 2000
    def __init__(self, *args, **kwargs):
        """Constructor."""
        self.src, self.dst = None, None
        self.limit = None
        Element.__init__(self, *args, **kwargs)

    def configure(self, src=None, dst=None, limit=None, **kwargs):
        """Set up source, destination, and receive limit."""
        self.src, self.dst = src, src
        if (limit is None): limit = self.Limit
        self.limit = limit
        self.nrcvd = 0
        f = self.newchild("ctrl", FSM, tracename=self.tracename+".MON")
        f.goto(self.MON)

    def MON(self, fsm):
        """Monitor `dst` to see if agent has received enough packets."""
        if not isinstance(self.src, Node): yield fsm.stop()
        if not isinstance(self.dst, Node): yield fsm.stop()
        if not self.src.haschild('agent'): yield fsm.stop()
        if not self.dst.haschild('agent'): yield fsm.stop()
        while (self.nrcvd<self.limit):
            yield waitevent, fsm, self.dst.agent.recv
            self.nrcvd = self.dst.agent.nrcvd
        # shutdown source and sleep
        self.src.agent.dest = None
        self.log("SHUTDOWN", rcvd=self.nrcvd, limit=self.limit)

def read_topo(options, topofile):
    """Read topology layout from file."""
    f = file(topofile, 'r')
    s = f.readline()
    topo = {'border':None, 'layout':None}
    done = not (s)
    while not done:
        # convert s to dict (check for border and layout)
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
    f.close()
    return topo

def get_topology(options):
    """Get/create new topology."""
    # create star topology: nodes on disc from snr and channel parameters
    snr = options.snr
    numnodes = options.numnodes + 1
    # channel parameters
    alpha      = options.alpha
    modeltype  = options.tgnmodel       # default -> LOS Channel
    # calculate SNR values
    Ptx = DOT11N_MAXPOWER
    Lrx, Ltx = Dot11NRadio.Lrx, Dot11NRadio.Ltx
    Grx, Gtx = Dot11NRadio.Grx, Dot11NRadio.Gtx
    No = Dot11NRadio.thermalnoise(DOT11N_BANDWIDTH) + DOT11N_NOISEFIGURE
    Psignal = Ptx + Gtx + Grx - Ltx - Lrx
    # use breakpoint distance as refdist
    refdist = Dot11NChannel(modeltype=modeltype, n=alpha).bpdist
    refloss = Propagation.freespace(refdist, n=2.0, fc=Dot11NChannel.fc)
    # calculate pathloss = Psignal - Prx
    PL = Psignal - (snr + No)
    n = 2.0
    if PL>refloss: n = alpha
    netradius = refdist*(db2linear((1.0*PL-refloss)/n))
    assert (netradius>0)
    # place nodes on disc of radius netradius
    pos = [(0,0)] # initial position is receiver
    for k in range(1, numnodes):
        r = np.random.uniform(0, netradius)
        th = np.random.uniform(0, 2*np.pi)
        x, y = r*np.cos(th), r*np.sin(th)
        pos.append((x,y))
    # create topology layout
    xmin = min([x for (x,y) in pos])
    xmax = max([x for (x,y) in pos])
    ymin = min([y for (x,y) in pos])
    ymax = max([y for (x,y) in pos])
    border = (xmin, xmax, ymin, ymax)
    layout = pos
    # verify layout parameters
    assert (len(layout)>=numnodes)
    topo = {'border':border, 'layout':layout}
    return topo

def run_experiment(options):
    # record start time
    starttime = time.time()
    # initialize RNG
    if RNG_INIT: RNG_init()

    # set SIMULATION parameters
    mon = Element(tracename="MON")
    verbose = options.verbose
    stoptime = 2.0
    if not (options.stop<0): stoptime = options.stop
    stoptime *= 1.05     # allow events around stoptime to finish
    simargs = {'verbose':verbose}

    # set EXPERIMENT parameters
    ntx, nrx = options.ntx, options.nrx
    numnodes = options.numnodes + 1
    mcs = 0
    if options.mcs: mcs = options.mcs
    assert (0<=mcs<8*ntx)
    assert (0<=mcs<8*nrx)
    workload = options.workload         # transmissions per frame-time
    dummy = Node()   # dummy node for PHY

    # set CHANNEL parameters
    alpha      = options.alpha
    modeltype  = options.tgnmodel       # default -> LOS Channel
    usedoppler = options.usedoppler
    usefading  = options.usefading
    envspeed   = options.envspeed
    chargs = {'modeltype':modeltype, 'n':alpha, \
              'usedoppler':usedoppler, 'usefading':usefading, \
              'environmentspeed': envspeed}
    chargs.update(simargs)

    # set AGENT parameters
    mode  = options.agent_mode
    if mode is None: mode = "poisson"
    plen  = Agent.DefaultPacketLength
    if options.plen>0: plen = options.plen
    duration = dummy.phy.calcduration(plen, mcs)
    load  = workload/duration
    rate  = load/options.numnodes       # transmission rate in frames/second
    delay = None
    if (rate>0): delay = 1.0/rate
    # set agent delay if not already specified
    if delay is None:
        cm = Dot11NChannel(**chargs)
        chan = Dot11N_Channel(cm.modelnum, nrx, ntx, cm.flags)
        delay = 2*chan.coherencetime()
        if rate is None: rate = 1.0/delay
    agtargs = {'plen': plen, 'mode':mode, 'delay':delay, 'datarate':mcs}

    # set NET parameters
    netargs = {}

    # set other protocol parameters (MAC, ARP, etc.)
    useshared = True
    arpargs = {'useshared':useshared}
    macargs = {}

    # set phy parameters
    Dot11NPHY.usewaveform = options.usewaveform
    Dot11NRadio.Ntx, Dot11NRadio.Nrx = ntx, nrx
    Dot11NRadio.fomax = options.fomax
    cfocorrection = True
    if options.disable_cfo_correction: cfocorrection = False
    phyargs = {'cfocorrection':cfocorrection}

    # set node parameters
    nodeargs = {}
    nodeargs.update(agtargs)
    nodeargs.update(netargs)
    nodeargs.update(arpargs)
    nodeargs.update(macargs)
    nodeargs.update(phyargs)
    nodeargs.update(simargs)

    ############################
    # Set Up Simulation
    ############################
    initialize()

    # create channel
    bidirectional = options.bidirectional
    ch = Channel(model=Dot11NChannel, bidirectional=bidirectional, **simargs)

    # get topology
    topo = get_topology(options)
    border = topo['border']
    layout = topo['layout']
    assert (numnodes==len(layout) )

    # create nodes in star topology
    nodelist = []
    for k in range(numnodes):
        pos = layout[k]
        n = Node(pos=pos, **nodeargs)
        nodelist.append(n)
        n.motion.log("pos", pos=["%.3f"%(p) for p in n.motion.position] )
    # set up senders to receiver
    receiver = nodelist[0]
    for k in range(1, numnodes):
        s = nodelist[k]
        s.agent.dest = receiver.net.address

    # connect all senders to receiver via channel
    for n in nodelist[1:]:
        assert (n is not receiver)
        ch.add_edge(n.cif, receiver.cif, **chargs)
        n.net.addroute(receiver.net.address)  # add route to receiver

    # create monitor
    if options.monitor:
        mon = Monitor(period=stoptime/1e4)
        mon.start()

    # create traffic monitors
    if options.limit and USE_TRAFFIC_MONITORS:
        tlist = []
        for n in nodelist[1:]:
            t = TrafficMonitor(src=n, dst=receiver, limit=options.limit)
            tlist.append(t)

    ############################
    # Run Simulation
    ############################
    mon.log("model", **chargs)
    mon.log("workload", G="%.5g"%(workload), T=time2usec(duration))
    mon.log("stoptime", stoptime="%.6f"%(stoptime))
    simerror = None
    if EXIT_WITH_TRACE:
        try:
            simulate(until=stoptime)
        except Exception, e:
            mon.log("SIMERR", error=str(e))
            simerror = e
    else:
        simulate(until=stoptime)
    # log remaining trace information
    n = gc.collect()
    mon.log("GC", collected=n)
    totaltime = time.time() - starttime
    t = time.gmtime(totaltime)
    mon.log("runtime", runtime="%02d:%02d:%02d (h/m/s)"%(t.tm_hour, t.tm_min, t.tm_sec) )

    ############################
    # Teardown/Cleanup
    ############################

    # print output
    sys.stdout.flush()
    if options.trace: ch.trace.output()

    # write tracefile
    if options.output is not None: ch.trace.write(options.output)

    # write topofile
    if options.savetopo:
        f = file(options.savetopo, 'w')
        f.write("%s\n"%(topo) )
        f.close()

    # if Exception occurred during simulation ...
    if simerror: raise simerror

def main():
    usage = "%prog [OPTIONS]"
    parser = OptionParser(usage=usage)
    # simulation parameters
    parser.add_option("-v", "--verbose", dest="verbose", type="int", \
            default=DOT11N_VERBOSE+1, help="Set verbose level [default=%default].")
    parser.add_option("-t", "--trace", dest="trace", action="store_true",  \
            default=False, help="Output formatted trace to stdout")
    parser.add_option("-o", "--output", dest="output", \
            default=None, help="Name of output file for trace")
    parser.add_option("-s", "--stop", dest="stop", \
            type="float", default=2.0, \
            help="Run simulation until stop time [default=%default]")
    parser.add_option("-m", "--monitor", dest="monitor", action="store_true",  \
            default=False, help="Enable simulation montior")
    parser.add_option("", "--limit", dest="limit", \
            type="int", default=TrafficMonitor.Limit, \
            help="Specify minimum number of samples needed before stopping [default=%default].")
    # experiment parameters
    parser.add_option("", "--snr", dest="snr", type="float", \
            default=20, help="Minimum SNR for links in network [default=%default]")
    parser.add_option("-G", "--workload", dest="workload", type="float", \
            default=1.0, help="Set total network traffic load in transmissions/frame-time [default=%default]")
    parser.add_option("-n", "--num-nodes", dest="numnodes", type="int", \
            default=100, help="Set number of nodes [default=%default]")
    # agent parameters
    parser.add_option("-r", "--rate", dest="rate", type="float", \
            default=None, help="Packets/second generated by a source [default=%default]")
    parser.add_option("-l", "--packet-length", dest="plen", type="int", \
            default=1024, help="Set packet size in bytes [default=%default]")
    parser.add_option("", "--agent-mode", dest="agent_mode", \
            default=None, help="Specify traffic mode [options=%s]."%(Agent.TrafficModes))
    # net parameters
    parser.add_option("", "--rreqrate", dest="rreqrate", type="int", \
            default=None, help="Set rate index for RREQ in DSR [default=%default]")
    parser.add_option("", "--datarate", dest="datarate", type="int", \
            default=None, help="Set rate index for non-RREQ packets in DSR [default=%default]")
    parser.add_option("", "--use-route", dest="useroute", \
            default=None, help="Specify routing file to initialize route tables.")
    parser.add_option("", "--save-route", dest="saveroute", \
            default=None, help="Save route tables to file.")
    # mac parameters
    # phy parameters
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
    parser.add_option("", "--ntx", dest="ntx", type="int", \
            default=1, help="Set number of transmit antennas [default=%default]")
    parser.add_option("", "--nrx", dest="nrx", type="int", \
            default=1, help="Set number of receive antennas [default=%default]")
    # channel parameters
    parser.add_option("", "--tgn-model", dest="tgnmodel", \
            default=None, help="Specify TGn model.")
    parser.add_option("", "--alpha", dest="alpha", type="float", \
            default=2.0, help="Specify pathloss exponent [default=%default]")
    parser.add_option("", "--use-doppler", dest="usedoppler", action="store_true",  \
            default=False, help="Enable doppler filter for fading in TGn channel model.")
    parser.add_option("", "--disable-fading", dest="usefading", action="store_false",  \
            default=True, help="Normalize channel and remove impact of fading on pathloss in TGn channel model.")
    parser.add_option("-E", "--environment-speed", dest="envspeed", type="float", \
            default=1.2, help="Environmental speed in (km/hr) [default=%default]")
    parser.add_option("", "--bidirectional-channel", dest="bidirectional", action="store_true", \
            default=False, help="Use bidirectional links in channel [default=%default]")
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
            default=None, help="Save topology to file.")
    (options, args) = parser.parse_args()

    if len(args)>0:
        print "Invalid number of arguments."
        parser.print_help()
        raise SystemExit

    run_experiment(options)

if __name__ == '__main__':
    main()
