#!  /usr/bin/env python

"""
Simulate interference scenario with interferers synchronized to sender.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-11-01 19:59:54 -0500 (Tue, 01 Nov 2011) $
* $LastChangedRevision: 5321 $

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
    def configure(self, datarate=None, pausefirst=True, **kwargs):
        self.datarate = datarate
        Agent.configure(self, **kwargs)
        if pausefirst:
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

class SyncAgent(MyAgent):
    name = "sync agent"
    def __init__(self, *args, **kwargs):
        self.trigger = None
        MyAgent.__init__(self, *args, **kwargs)
    def configure(self, trigger=None, maxdelay=None, **kwargs):
        MyAgent.configure(self, **kwargs)
        self.trigger = trigger
        self.maxdelay = maxdelay
        assert isinstance(trigger, SimEvent)
        assert (self.mode=="poisson")
        self.txfsm.goto(self.PAUSE)
    def PAUSE(self, fsm):
        """PAUSE state; wait on trigger before pausing."""
        yield hold, fsm, 0  # yield to other threads
        if (not self.dest): yield fsm.stop()    # HALT and sleep
        # wait for trigger
        avgdelay = self.delay
        maxdelay = self.maxdelay
        #self.log("SYNC", avgdelay=avgdelay, maxdelay=maxdelay)
        yield waitevent, fsm, self.trigger
        # calculate prescribed delay
        # -> 'poisson' mode uses exponentail backoff
        assert (self.delay>0)
        d = np.random.exponential(self.delay)
        # wait and send if delay d < maxdelay
        offset = 0  #DOT11N_TSHORT # STF duration
        yield hold, fsm, offset
        if (d<self.maxdelay-offset):
            self.log("WAIT", delay=d, avgdelay=avgdelay, maxdelay=maxdelay)
            yield hold, fsm, d
            yield fsm.goto(self.SEND)
        else:
            # otherwise -> continue in PAUSE
            yield fsm.goto(self.PAUSE)

class Node(Element):
    name = "node"
    tracename = "NODE"
    def __init__(self, **kwargs):
        Element.__init__(self, **kwargs)
    def configure(self, pos=None,                       # motion    \
                        useshared=True,                 # arp       \
                        cfocorrection=True,             # phy       \
                        datarate=None, pausefirst=True, # myagent   \
                        dest=None, plen=None, delay=None, mode=None, # agent \
                        trigger=None, maxdelay=None,   # syncagent \
                        **kwargs):
        cif = self.newchild('cif', Dot11NRadio)
        phy = self.newchild('phy', Dot11NPHY, radio=cif, cfocorrection=cfocorrection)
        mac = self.newchild('mac', Aloha, phy=phy)
        net = self.newchild('net', Routing)
        arp = self.newchild('arp', ARP, useshared=useshared)
        if isinstance(trigger, SimEvent):
            agt = self.newchild('agent', SyncAgent, dest=dest, plen=plen, \
                                         delay=delay, mode=mode, \
                                         trigger=trigger, maxdelay=maxdelay, \
                                         datarate=datarate, pausefirst=pausefirst)
        else:
            agt = self.newchild('agent', MyAgent, dest=dest, plen=plen, \
                                         delay=delay, mode=mode, \
                                         datarate=datarate, pausefirst=pausefirst)
        mobi = self.newchild('motion', Motion, pos=pos)
        # connect ports
        agt.connect(net)
        arp.connect(net, mac)
        mac.connect(phy)
        phy.connect(cif)

def get_topology(options):
    """Get/create new topology."""
    # create interference topology: fixed SR link + interferers on a line
    snr = options.snr
    numnodes = options.numnodes
    sinrmin = options.sinrmin
    sinrmax = options.sinrmax
    sinr = np.linspace(sinrmax, sinrmin, numnodes)
    snr_lin  = db2linear(snr)
    sinr_lin = db2linear(sinr)
    assert all(sinr<snr)
    sir_lin  = 1.0/(1.0/sinr_lin - 1.0/snr_lin)
    sir = linear2db(sir_lin)
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
    # calculate pathloss for direct link = Psignal - Prx
    Prx = snr + No
    PL  = Psignal - Prx
    n = 2.0
    if PL>refloss: n = alpha
    srdist = refdist*(db2linear((1.0*PL-refloss)/n))
    assert (srdist>0)
    # place nodes on line using SIR
    pos = [(-srdist, 0), (0,0)]     # initialize with send and receiver
    for s in sir:
        Pi = Prx - s
        pl = Psignal - Pi
        n = 2.0
        if pl>refloss: n = alpha
        x = refdist*(db2linear((1.0*pl-refloss)/n))
        pos.append((x,0))
    assert (numnodes+2==len(pos))
    # create topology layout
    xmin = min([x for (x,y) in pos])
    xmax = max([x for (x,y) in pos])
    ymin = min([y for (x,y) in pos])
    ymax = max([y for (x,y) in pos])
    border = (xmin, xmax, ymin, ymax)
    layout = pos
    # verify layout parameters
    assert (len(layout)>=numnodes)
    topo = {'border':border, 'layout':layout, 'sir':[None, None]+sir.tolist()}
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
    numnodes = options.numnodes
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

    # set AGENT parameters for sender
    smode = "cbr"
    plen  = Agent.DefaultPacketLength
    if options.plen>0: plen = options.plen
    # calculate length/duration of PHY packet
    phylen = len(Ether()/IP()/AGT()) + plen + len(CRC32())
    frametime = dummy.phy.calcduration(phylen, mcs)
    srate = 1.0/3                     # sender rate = 1 packet/N frametimes
    sdelay = frametime/srate
    sndargs = {'mode': smode, 'delay': sdelay}

    # set AGENT parameters for interferers
    imode  = options.agent_mode
    if imode is None: imode = "poisson"
    # calculat rates and delays
    load  = workload/frametime
    irate = load/numnodes               # interferer rate in frames/second
    idelay = 1.0/irate
    agtargs = {'plen': plen, 'mode':imode, 'delay':idelay, 'datarate':mcs}

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
    sir    = topo['sir']
    assert (numnodes+2==len(layout) )

    # create receiver
    pos = layout[0]
    receiver = Node(pos=pos, **nodeargs)
    receiver.motion.log("pos", pos=["%.3f"%(p) for p in receiver.motion.position] )
    # create sender
    pos = layout[1]
    sargs = nodeargs.copy()
    sargs.update(sndargs)
    sender = Node(pos=pos, **sargs)
    sender.motion.log("pos", pos=["%.3f"%(p) for p in sender.motion.position] )
    # create nodes in interference topology
    ilist = []
    trigger = sender.agent.send
    maxdelay = frametime
    for k in range(numnodes):
        pos = layout[2+k]
        n = Node(pos=pos, trigger=trigger, maxdelay=maxdelay, **nodeargs)
        ilist.append(n)
        n.motion.log("pos", pos=["%.3f"%(p) for p in n.motion.position] , \
                            sir=sir[2+k])
    # set up sender and interferer agents
    sender.agent.dest = receiver.net.address
    sender.net.addroute(receiver.net.address)   # add route to receiver
    for n in ilist:
        n.agent.dest = sender.net.address   # throw away addresses
        n.net.addroute(sender.net.address)  # add sender as neighbor

    # connect all senders/interferers to receiver via channel
    ch.add_edge(sender.cif, receiver.cif, **chargs)
    for n in ilist:
        assert (n is not receiver)
        ch.add_edge(n.cif, receiver.cif, **chargs)

    # create monitor
    if options.monitor:
        mon = Monitor(period=stoptime/1e4)
        mon.start()

    ############################
    # Run Simulation
    ############################
    mon.log("model", **chargs)
    mon.log("workload", G="%.5g"%(workload), T=time2usec(frametime))
    mon.log("stoptime", stoptime="%.6f"%(stoptime))
    mon.log("sender",   net=sender.net.address, mac=sender.mac.address, \
                        nodeid=sender.uid)
    mon.log("receiver", net=receiver.net.address, mac=receiver.mac.address, \
                        nodeid=receiver.uid)
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
    mon.log("stoptime", stoptime="%.6f"%(now()) )

    ############################
    # Teardown/Cleanup
    ############################

    # print output
    sys.stdout.flush()
    if options.trace: ch.trace.output()

    # write tracefile
    if options.output is not None: ch.trace.write(options.output)

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
    # experiment parameters
    parser.add_option("", "--snr", dest="snr", type="float", \
            default=10, help="Set SNR of direct link [default=%default]")
    parser.add_option("", "--sinrmax", dest="sinrmax", type="float", \
            default=10, help="Set maximum SINR value for simulation [default=%default]")
    parser.add_option("", "--sinrmin", dest="sinrmin", type="float", \
            default=10, help="Set minimum SINR value for simulation [default=%default]")
    parser.add_option("-n", "--num-nodes", dest="numnodes", type="int", \
            default=10, help="Set number of interference nodes [default=%default]")
    parser.add_option("-G", "--workload", dest="workload", type="float", \
            default=1.0, help="Set transmission rate of interferers in transmissions/frame-time [default=%default]")
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
