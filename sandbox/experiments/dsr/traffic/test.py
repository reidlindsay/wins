#!  /usr/bin/env python

"""
Simulate DSR network.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-12-08 20:24:22 -0600 (Thu, 08 Dec 2011) $
* $LastChangedRevision: 5362 $

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
import threading

RNG_INIT = 1
EXIT_WITH_TRACE = 0

RBAR_THRESHOLD        = {0: 4.2, 1: 7.2, 2:10.1, 3:13.8, 4:16.9, 5:21.7, 6:22.9, 7:24.7}
# SINR thresholds tuned for CM-D and PER<10%
RBAR_CMD_THRESHOLD_10 = {0:10.0, 1:13.0, 2:18.0, 3:20.5, 4:25.7, 5:28.5, 6:31.6, 7:35.6}
# SINR thresholds tuned for CM-D and PER<5%
RBAR_CMD_THRESHOLD_5  = {0:11.7, 1:14.8, 2:19.8, 3:22.0, 4:26.7, 5:29.2, 6:32.6, 7:36.3}
# SINR thresholds tuned for CM-D and PER<1%
RBAR_CMD_THRESHOLD_1  = {0:14.2, 1:16.8, 2:22.6, 3:24.5, 4:29.5, 5:31.8, 6:36.1, 7:38.1}
# SINR thresholds tuned for CFO and PER<5%
RBAR_CFO_THRESHOLD_5  = {0: 7.9, 1:11.4, 2:12.4, 3:15.6, 4:17.8, 5:22.6, 6:23.7, 7:25.2}

class MyAgent(Agent):
    name = "my agent"
    def __init__(self, *args, **kwargs):
        self.probe = None       # probe routes?
        self.destinations = []  # possible destinations
        self.probed = []        # destinations probed
        Agent.__init__(self, *args, **kwargs)
    ndest = property(fget=lambda self: len(self.destinations))
    nprobed = property(fget=lambda self: len(self.probed))
    def configure(self, probe=False, **kwargs):
        self.probe = probe
        Agent.configure(self, **kwargs)
    def INIT(self, fsm):
        """INIT state."""
        if self.dest: self.nextdest(reset=True)
        yield fsm.goto(self.PAUSE)
    def nextdest(self, reset=False):
        """Set next destination."""
        assert self.destinations
        dest = self.destinations[np.random.randint(0,self.ndest)]
        if self.probe:
            if reset: self.probed = []
            # probe next destination
            if (self.nprobed<self.ndest):
                dest = self.destinations[self.nprobed]
            else:
                dest = None     # halt since all dest were probed
                assert all([(a in self.probed) for a in self.destinations])
            self.probed.append(dest)
            self.log("probe", dest=dest)
        # set destination
        self.dest = dest
    def senddata(self, p):
        """Get destination for next packet."""
        if self.dest: self.nextdest()
        return Agent.senddata(self, p)

class Node(Element):
    name = "node"
    tracename = "NODE"
    def __init__(self, **kwargs):
        Element.__init__(self, **kwargs)
    def configure(self, pos=None,                       # motion    \
                        useshared=True,                 # arp       \
                        cfocorrection=True,             # phy       \
                        base=None, thresh=None,         # mac       \
                        rreqrate=None, datarate=None,   # net       \
                        dest=None, plen=None, delay=None, mode=None, probe=False, # agent \
                        userbar=False, usecsma=False,    # node setup \
                        **kwargs):
        cif = self.newchild('cif', Dot11NRadio)
        phy = self.newchild('phy', Dot11NPHY, radio=cif, cfocorrection=cfocorrection)
        if userbar:
            mac = self.newchild('mac', RBAR, base=base, thresh=thresh, phy=phy)
        else:
            mac = self.newchild('mac', DCF, usecsma=usecsma, phy=phy)
        net = self.newchild('net', DSR, rreqrate=rreqrate, datarate=datarate)
        arp = self.newchild('arp', ARP, useshared=useshared)
        agt = self.newchild('agent', MyAgent, dest=dest, plen=plen, \
                                     delay=delay, mode=mode, probe=probe)
        mobi = self.newchild('motion', Motion, pos=pos)
        # connect ports
        agt.connect(net)
        arp.connect(net, mac)
        mac.connect(phy)
        phy.connect(cif)

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

def read_route(options, routefile):
    """Read routing tables from file."""
    f = file(routefile, 'r')
    s = f.readline()
    routedata = {}
    done = not (s)
    while not done:
        # convert s to dict
        try:
            d = eval(s)
            assert isinstance(d, dict)
            for x,y in d.items():
                # maps src x -> routing table y
                assert isinstance(y, dict)
                for a,b in y.items():
                    # maps dst a -> info b (for route table y)
                    assert ('index' in b)
                    assert ('list' in b)
        except:
            d = None
        # add dict to routedata
        if d: routedata.update(d)
        # get next input
        s = f.readline()
        done = not(s)
    f.close()
    return routedata

def get_topology(options):
    """Get/create topology."""
    # load topology from file
    numnodes = options.numnodes
    if options.usetopo:
        topofile = options.usetopo
        topo = read_topo(options, topofile)
        border = topo['border']
        layout = topo['layout']
        xmin, xmax, ymin, ymax = border[:4]
        assert (len(layout)>=numnodes)
        return topo
    # create new topology
    assert (options.xmin<=options.xmax)
    assert (options.ymin<=options.ymax)
    xmin, xmax = options.xmin, options.xmax
    ymin, ymax = options.ymin, options.ymax
    border = (xmin, xmax, ymin, ymax)
    # use uniform distribution for layout
    xpos = np.random.uniform(xmin, xmax, numnodes)
    ypos = np.random.uniform(ymin, ymax, numnodes)
    layout = [(xpos[k],ypos[k]) for k in range(numnodes)]
    # verify layout parameters
    assert (len(layout)>=numnodes)
    topo = {'border':border, 'layout':layout}
    return topo

def set_routing(options, nodelist):
    """Set routing tables if needed."""
    if (not options.useroute) or options.probe: return
    routefile = options.useroute
    rdata = read_route(options, routefile)
    for n in nodelist:
        addr = n.net.address
        if addr not in rdata: continue
        for dst, data in rdata[addr].items():
            paths = data['list']
            for c,ts,nh in paths:
                n.net.addroute(dst, nexthop=nh, cost=c)
    return rdata

class StopSim(threading.Thread):
    """Thread for stopping simulation early."""
    PERIOD = 0.10
    def __init__(self, *args, **kwargs):
        self._done = 0
        threading.Thread.__init__(self, None, self.stopsim, None, args, kwargs)
        self.setDaemon(1)
        # start FSM to actually stop simulation
        f = FSM.launch(self.STOP)
    def stopsim(self, maxruntime):
        stoptime = maxruntime
        # convert string to float "
        if isinstance(maxruntime, str):
            x = maxruntime.split(":")
            t = 0
            for k in range(len(x)):
                t = t*60 + float(x[k])
            stoptime = t
        # wait until stoptime
        if (stoptime>0):
            time.sleep(stoptime)
        # set flag to stop simulation
        self._done = 1
    def STOP(self, fsm):
        # pause periodically until done
        while not self._done:
            yield hold, fsm, self.PERIOD
        # stop simulation
        stopSimulation()

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
    nconnect = options.nconnect
    assert (nconnect>0)
    assert (numnodes>=2*nconnect)
    plen  = Agent.DefaultPacketLength
    if options.plen>0: plen = options.plen

    # set CHANNEL parameters
    alpha      = None # options.alpha   # use default pathloss exponent
    modeltype  = options.tgnmodel       # default -> LOS Channel
    envspeed   = options.envspeed
    if not (envspeed>0):
        usedoppler = False
        usefading  = False
        envspeed = 0
        vo = None
    else:
        usedoppler = options.usedoppler
        usefading  = options.usefading
        vo = envspeed
    chargs = {'modeltype':modeltype, 'n':alpha, \
              'usedoppler':usedoppler, 'usefading':usefading, \
              'environmentspeed': vo}
    chargs.update(simargs)

    # set AGENT parameters
    mode  = options.agent_mode
    probe = options.probe
    if mode is None: mode = "poisson"
    rate  = options.rate     # transmission rate in packets/second
    #load  = options.rate     # transmission rate in packets/second
    #rate  = load/nconnect
    delay = None
    if (rate>0): delay = 1.0/rate
    # use no delay if rate is not specified
    agtargs = {'plen': plen, 'mode':mode, 'delay':delay, 'probe': probe}

    # set DSR parameters
    rreqrate, datarate = None, None
    if 0<=options.rreqrate<8*ntx: rreqrate=options.rreqrate
    if 0<=options.datarate<8*ntx: datarate=options.datarate
    netargs = {'rreqrate':rreqrate, 'datarate':datarate}

    # set other protocol parameters (MAC, ARP, etc.)
    useshared = True
    arpargs = {'useshared':useshared}
    userbar = options.userbar
    usecsma = False
    thresh = None
    if options.conservative:
        #thresh = RBAR_CMD_THRESHOLD_5   # 5 % PER
        thresh = RBAR_CFO_THRESHOLD_5   # 5 % PER
    macargs = {'userbar':userbar, 'usecsma':usecsma, 'thresh':thresh}

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
    assert (len(layout)>=numnodes)

    # create nodes
    nodelist = []
    for k in range(numnodes):
        pos = layout[k]
        n = Node(pos=pos, **nodeargs)
        nodelist.append(n)
        n.motion.log("pos", pos=["%.3f"%(p) for p in n.motion.position] )
    # initialize source/destination pairs
    assert (nconnect<len(nodelist))
    for k in range(nconnect):
        src = nodelist[k]               # first N are sources
        dst = nodelist[-k-1]            # last  N are destinations
        # set destination list and initial address
        src.agent.destinations = [n.net.address for n in nodelist if (n is not src)]
        src.agent.dest = dst.net.address

    # set routing tables
    set_routing(options, nodelist)

    # connect all nodes via channel
    for n in nodelist:
        for m in nodelist:
            if (n is not m):
                ch.add_edge(n.cif, m.cif, **chargs)

    # create monitor
    if options.monitor:
        mon = Monitor(period=stoptime/1e4)
        mon.start()

    # create StopSim
    if options.max_runtime:
        stopper = StopSim(options.max_runtime)
        stopper.start()

    ############################
    # Run Simulation
    ############################
    mon.log("thresh", thresholds="%s"%(thresh))
    if options.usetopo:
        mon.log("topo", topofile=options.usetopo)
    chargs['environmentspeed'] = envspeed
    chargs['bidirectional'] = bidirectional
    mon.log("model", **chargs)
    mon.log("rate", rate=rate, rreqrate=rreqrate, datarate=datarate)
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

    # write routefile
    if options.saveroute:
        # write data
        f = file(options.saveroute, 'w')
        for n in nodelist:
            addr = n.net.address
            rdata = {addr: n.net.table.data.copy()}
            f.write("%s\n"%(rdata))
        f.close()

    # if Exception occurred during simulation ...
    if simerror: raise simerror

def main():
    usage = "%prog [OPTIONS]"
    parser = OptionParser(usage=usage)
    # simulation parameters
    parser.add_option("-v", "--verbose", dest="verbose", type="int", \
            default=ROUTING_VERBOSE+1, help="Set verbose level [default=%default].")
    parser.add_option("-t", "--trace", dest="trace", action="store_true",  \
            default=False, help="Output formatted trace to stdout")
    parser.add_option("-o", "--output", dest="output", \
            default=None, help="Name of output file for trace")
    parser.add_option("-s", "--stop", dest="stop", \
            type="float", default=2.0, \
            help="Run simulation until stop time [default=%default]")
    parser.add_option("-m", "--monitor", dest="monitor", action="store_true",  \
            default=False, help="Enable simulation montior")
    parser.add_option("", "--max-runtime", dest="max_runtime", \
            default=None, help="Specify maximum runtime (e.g. HH:MM:SS)")
    # experiment parameters
    parser.add_option("-n", "--num-nodes", dest="numnodes", type="int", \
            default=50, help="Set number of nodes [default=%default]")
    parser.add_option("-c", "--num-connections", dest="nconnect", type="int", \
            default=1, help="Set number of active connections [default=%default]")
    # agent parameters
    parser.add_option("-r", "--rate", dest="rate", type="float", \
            default=None, help="Packets/second generated by a source [default=%default]")
    parser.add_option("-l", "--packet-length", dest="plen", type="int", \
            default=1024, help="Set packet size in bytes [default=%default]")
    parser.add_option("", "--agent-mode", dest="agent_mode", \
            default=None, help="Specify traffic mode [options=%s]."%(Agent.TrafficModes))
    parser.add_option("", "--probe", dest="probe", action="store_true",  \
            default=False, help="Probe routing instead of normal operation.")
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
    parser.add_option("", "--use-rbar", dest="userbar", action="store_true",  \
            default=False, help="Use RBAR instead of DCF for MAC.")
    parser.add_option("", "--conservative", dest="conservative", action="store_true",  \
            default=False, help="Enable conservative RBAR thresholds.")
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
