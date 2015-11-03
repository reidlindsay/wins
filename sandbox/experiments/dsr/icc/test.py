#!  /usr/bin/env python

"""
Simulate DSR over a network of nodes.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-10-26 21:51:40 -0500 (Wed, 26 Oct 2011) $
* $LastChangedRevision: 5314 $

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
from wins.backend import *

from wins.mac import RBAR, ARF
from wins.net import DSR

from wins.traffic import Agent

import sys
from optparse import OptionParser

import numpy as np
import struct
import gc
import time

RNG_INIT = 1
EXIT_WITH_TRACE = 1

class Node(Element):
    name = "node"
    tracename = "NODE"
    def __init__(self, **kwargs):
        Element.__init__(self, **kwargs)
    def configure(self, pos=None,                       # motion    \
                        useshared=False,                # arp       \
                        cfocorrection=True,             # phy       \
                        usecsma=False,                  # mac       \
                        rreqrate=None, datarate=None,   # net       \
                        dest=None, plen=None, delay=None, mode=None, # agent \
                        **kwargs):
        cif = self.newchild('cif', Dot11NRadio)
        phy = self.newchild('phy', Dot11NPHY, radio=cif, cfocorrection=cfocorrection)
        mac = self.newchild('mac', DCF, usecsma=usecsma, phy=phy)
        net = self.newchild('net', DSR, rreqrate=rreqrate, datarate=datarate)
        arp = self.newchild('arp', ARP, useshared=useshared)
        agt = self.newchild('agent', Agent, dest=dest, plen=plen, \
                                            delay=delay, mode=mode)
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

def get_topology(options, numnodes):
    """Get/create topology."""
    # load topology from file
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
    if not options.useroute: return
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
    ntx, nrx = 1, 1
    numnodes = options.numnodes
    nconnect = options.nconnect
    assert (nconnect>0)
    assert (numnodes>=2*nconnect)

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
    plen  = Agent.DefaultPacketLength
    rate  = options.rate     # transmission rate in packets/second
    delay = None
    if mode is None: mode = "cbr"
    if options.plen>0: plen = options.plen
    if (rate>0): delay = 1.0/rate
    # set agent delay if not already specified
    if delay is None:
        cm = Dot11NChannel(**chargs)
        chan = Dot11N_Channel(cm.modelnum, nrx, ntx, cm.flags)
        delay = 2*chan.coherencetime()
        if rate is None: rate = 1.0/delay
    agtargs = {'plen': plen, 'mode':mode, 'delay':delay}

    # set DSR parameters
    rreqrate, datarate = None, None
    if 0<=options.rreqrate<8*ntx: rreqrate=options.rreqrate
    if 0<=options.datarate<8*ntx: datarate=options.datarate
    netargs = {'rreqrate':rreqrate, 'datarate':datarate}

    # set other protocol parameters (MAC, ARP, etc.)
    useshared = True
    arpargs = {'useshared':useshared}
    usecsma   = False
    macargs = {'usecsma':usecsma}

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
    topo = get_topology(options, numnodes)
    border = topo['border']
    layout = topo['layout']

    # create nodes
    nodelist = []
    for k in range(numnodes):
        pos = layout[k]
        n = Node(pos=pos, **nodeargs)
        nodelist.append(n)
        n.motion.log("pos", pos=["%.3f"%(p) for p in n.motion.position] )
    # connect source/destination pairs
    assert (nconnect<len(nodelist))
    for k in range(nconnect):
        src = nodelist[k]               # first N are sources
        dst = nodelist[-k-1]            # last  N are destinations
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

    ############################
    # Run Simulation
    ############################
    if options.usetopo:
        mon.log("topo", topofile=options.usetopo)
    mon.log("model", **chargs)
    mon.log("rate", rate="%.5g"%(rate) )
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
    mon.log("stoptime", stoptime="%.6f"%(stoptime))
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
    # net parameters
    parser.add_option("", "--rreqrate", dest="rreqrate", type="int", \
            default=None, help="Set rate index for RREQ in DSR [default=%default]")
    parser.add_option("", "--datarate", dest="datarate", type="int", \
            default=None, help="Set rate index for non-RREQ packets in DSR [default=%default]")
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
    # routing parameters
    parser.add_option("", "--use-route", dest="useroute", \
            default=None, help="Specify routing file to initialize route tables.")
    parser.add_option("", "--save-route", dest="saveroute", \
            default=None, help="Save route tables to file.")
    (options, args) = parser.parse_args()

    if len(args)>0:
        print "Invalid number of arguments."
        parser.print_help()
        raise SystemExit

    run_experiment(options)

if __name__ == '__main__':
    main()
