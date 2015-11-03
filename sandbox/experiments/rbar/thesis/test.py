#!  /usr/bin/env python

"""
Simulate point-to-point link using RBAR.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-10-24 05:45:57 -0500 (Mon, 24 Oct 2011) $
* $LastChangedRevision: 5305 $

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

import sys
from optparse import OptionParser

import numpy as np
import struct
import gc
import time

RNG_INIT = 1
EXIT_WITH_TRACE = 0

RBAR_THRESHOLD        = {0: 4.2, 1: 7.2, 2:10.1, 3:13.8, 4:16.9, 5:21.7, 6:22.9, 7:24.7}
# SINR thresholds tuned for CM-D and PER<10%
RBAR_CMD_THRESHOLD_10 = {0:10.0, 1:13.0, 2:18.0, 3:20.5, 4:25.7, 5:28.5, 6:31.6, 7:35.6}
# SINR thresholds tuned for CM-D and PER<5%
RBAR_CMD_THRESHOLD_5  = {0:11.7, 1:14.8, 2:19.8, 3:22.0, 4:26.7, 5:29.2, 6:32.6, 7:36.3}
# SINR thresholds tuned for CM-D and PER<1%
RBAR_CMD_THRESHOLD_1  = {0:14.2, 1:16.8, 2:22.6, 3:24.5, 4:29.5, 5:31.8, 6:36.1, 7:38.1}

class Node(Element):
    name = "node"
    tracename = "NODE"
    def __init__(self, **kwargs):
        Element.__init__(self, **kwargs)
    def configure(self, pos=None,                       # motion    \
                        useshared=True,                 # arp       \
                        cfocorrection=True,             # phy       \
                        base=None, thresh=None,         # mac       \
                        dest=None, plen=None, delay=None, mode=None, # agent \
                        **kwargs):
        cif = self.newchild('cif', Dot11NRadio)
        phy = self.newchild('phy', Dot11NPHY, radio=cif, cfocorrection=cfocorrection)
        mac = self.newchild('mac', RBAR, base=base, thresh=thresh, phy=phy)
        net = self.newchild('net', Routing)
        arp = self.newchild('arp', ARP, useshared=useshared)
        agt = self.newchild('agent', Agent, dest=dest, plen=plen, \
                                     delay=delay, mode=mode)
        mobi = self.newchild('motion', Motion, pos=pos)
        # connect ports
        agt.connect(net)
        arp.connect(net, mac)
        mac.connect(phy)
        phy.connect(cif)

class MotionControl(Element):
    """Monitor/control motion of nodes."""
    name = "motion control"
    tracename = "MC"
    Limit = 100
    def __init__(self, *args, **kwargs):
        """Constructor."""
        self.src, self.dst = None, None
        self.nrcvd, self.idx = 0, None
        Element.__init__(self, *args, **kwargs)

    def configure(self, src=None, dst=None, limit=None, \
                        snr=None, pos=None, **kwargs):
        """Configure parameters."""
        if limit is None: limit = self.Limit
        self.src, self.dst = src, dst
        self.snr, self.pos = snr, pos
        self.limit = limit
        assert (len(snr)==len(pos))
        # initialize counters and parameters
        self.idx = None
        self.nrcvd, self.nsent = 0, 0
        # start motion controller
        f = self.newchild("ctrl", FSM, tracename=self.tracename+".CTRL")
        f.goto(self.MOVE)

    def MOVE(self, fsm):
        """Initialize position of destination."""
        # check for valid dst
        if not isinstance(self.dst, Node):  yield fsm.stop()
        if not self.dst.haschild('motion'): yield fsm.stop()
        # get index into position list
        if self.idx is None: self.idx = 0
        else:                self.idx += 1
        if not (self.idx<len(self.pos)): yield fsm.goto(self.SHUTDOWN)
        # set new position
        pos = self.pos[self.idx]
        snr = "%.4f dB"%(self.snr[self.idx])
        self.dst.motion.position = pos
        # log trace information
        self.log("AVGSNR", pos=pos, snr=snr)
        yield fsm.goto(self.CTRL)

    def CTRL(self, fsm):
        """Control motion of receiver."""
        # check for valid src and dst
        if not isinstance(self.src, Node): yield fsm.stop()
        if not isinstance(self.dst, Node): yield fsm.stop()
        if not self.src.haschild('agent'): yield fsm.stop()
        if not self.dst.haschild('agent'): yield fsm.stop()
        # wait for packet to be received by Agent
        self.nrcvd, self.nsent = 0, 0
        sndagt, rcvagt = self.src.agent, self.dst.agent
        done = False
        while (not done):
            yield waitevent, fsm, (rcvagt.recv, sndagt.send)
            if (rcvagt.recv in fsm.eventsFired):   self.nrcvd +=1
            elif (sndagt.send in fsm.eventsFired): self.nsent +=1
            # check if enough packet have been sent/received
            done  = (self.nrcvd>self.limit)
            done |= (self.nrcvd<1) and (self.nsent>self.limit*2)
            done |= (self.nrcvd<self.limit) and (self.nsent>self.limit*25)
        # advance position
        yield fsm.goto(self.MOVE)

    def SHUTDOWN(self, fsm):
        """Shutdown sender when all data has been collected."""
        # check for valid src
        if not isinstance(self.src, Node): yield fsm.stop()
        if not self.src.haschild('agent'): yield fsm.stop()
        # shutdown agent
        self.src.agent.dest = None
        self.log("SHUTDOWN", nrcvd=self.nrcvd, nsent=self.nsent)
        yield fsm.stop()        # HALT

def get_topology(options):
    """Get/create topology."""
    # create new line topology from snr and channel parameters
    snrmin, snrmax = options.snrmin, options.snrmax
    dsnr = options.snrstep
    assert (snrmin<=snrmax)
    assert (dsnr>0)
    snr = numpy.arange(int((1.0*snrmax-snrmin)/dsnr)+1)*dsnr + snrmin
    snrlist = snr.tolist()
    snrlist.reverse()   # go from max to min snr
    snr = np.array(snrlist)
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
    # use reference pathloss model (with BPdist as RefDist)
    pos = [(0,0)]       # initial position is sender
    for pl in PL.tolist():
        n = 2.0
        if pl>refloss: n = alpha
        x = refdist*(db2linear((1.0*pl-refloss)/n))
        pos.append((x,0))
    numnodes = len(snrlist) + 1
    # create topology layout
    xmin = min([x for (x,y) in pos])
    xmax = max([x for (x,y) in pos])
    ymin = min([y for (x,y) in pos])
    ymax = max([y for (x,y) in pos])
    border = (xmin, xmax, ymin, ymax)
    layout = pos
    # verify layout parameters
    assert (len(layout)>=numnodes)
    topo = {'border':border, 'layout':layout, 'snr': snrlist}
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
    plen  = Agent.DefaultPacketLength
    if options.plen>0: plen = options.plen

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
    if mode is None: mode = "backlog"
    rate  = options.rate     # transmission rate in packets/second
    delay = None
    if (rate>0): delay = 1.0/rate
    # use no delay if rate is not specified
    agtargs = {'plen': plen, 'mode':mode, 'delay':delay}

    # set NET parameters
    netargs = {}

    # set other protocol parameters (MAC, ARP, etc.)
    useshared = True
    arpargs = {'useshared':useshared}
    thresh = None
    if options.conservative: thresh = RBAR_CMD_THRESHOLD_5   # 5 % PER
    macargs = {'thresh':thresh}

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
    border  = topo['border']
    layout  = topo['layout']
    snrlist = topo['snr']
    assert (len(layout)>1)
    assert (len(layout)>len(snrlist))

    # create sender
    pos = layout[0]
    sender = Node(pos=pos, **nodeargs)
    sender.motion.log("pos", pos=["%.3f"%(p) for p in sender.motion.position] )
    # create receiver
    pos = layout[1]
    receiver = Node(pos=pos, **nodeargs)
    receiver.motion.log("pos", snr="%.3f dB"%(snrlist[0]), \
            pos=["%.3f"%(p) for p in receiver.motion.position] )
    # set up sender for unicast to receiver
    sender.agent.dest = receiver.net.address
    sender.net.addroute(receiver.net.address)   # add route to receiver

    # connect sender and receiver via channel
    ch.add_edge(sender.cif, receiver.cif, **chargs)
    ch.add_edge(receiver.cif, sender.cif, **chargs)

    # set up motion controller
    npos = len(snrlist)
    mc = MotionControl(src=sender, dst=receiver, limit=options.limit, \
                       snr=snrlist, pos=layout[1:], **simargs)

    # create monitor
    if options.monitor:
        mon = Monitor(period=stoptime/1e4)
        mon.start()

    ############################
    # Run Simulation
    ############################
    sender.mac.log("THRESH", thresholds="%s"%(sender.mac.thresh))
    mon.log("model", **chargs)
    mon.log("limit", limit=mc.limit)
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
            type="int", default=MotionControl.Limit, \
            help="Specify minimum number of samples needed before moving [default=%default].")
    # experiment parameters
    parser.add_option("", "--conservative", dest="conservative", action="store_true",  \
            default=False, help="Enable conservative RBAR thresholds.")
    parser.add_option("", "--snrmin", dest="snrmin", type="float", \
            default=0, help="Set minimum SNR value for simulation [default=%default]")
    parser.add_option("", "--snrmax", dest="snrmax", type="float", \
            default=30, help="Set maximum SNR value for simulation [default=%default]")
    parser.add_option("", "--snrstep", dest="snrstep", type="float", \
            default=5, help="Set SNR step size for simulation [default=%default]")
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
