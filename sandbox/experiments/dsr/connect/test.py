#!  /usr/bin/env python

"""
Simulate DSR network.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-12-18 16:57:02 -0600 (Sun, 18 Dec 2011) $
* $LastChangedRevision: 5381 $

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

from mastercontrol import MasterControl, MyAgent

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
    plen  = Agent.DefaultPacketLength
    if options.plen>0: plen = options.plen

    # set CHANNEL parameters
    alpha      = options.alpha          # use default pathloss exponent
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
    if mode is None: mode = "cbr"
    rate  = options.rate     # transmission rate in packets/second
    if rate is None: rate = 1.0
    delay = None
    if (rate>0): delay = 1.0/rate
    # use no delay if rate is not specified
    agtargs = {'plen': plen, 'mode':mode, 'delay':delay}

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

    # get topology border
    numnodes = options.numnodes
    assert (options.xmin<=options.xmax)
    assert (options.ymin<=options.ymax)
    xmin, xmax = options.xmin, options.xmax
    ymin, ymax = options.ymin, options.ymax
    border = (xmin, xmax, ymin, ymax)

    # create src/dst
    xmin, xmax, ymin, ymax = border[0:4]
    src = Node(pos=(xmin,ymin), **nodeargs)
    dst = Node(pos=(xmax,ymax), **nodeargs)
    # initialize traffic
    src.agent.dest = dst.net.address
    # make address book
    netaddr = {src.net.address:src, dst.net.address:dst}

    # create remaining nodes
    for k in range(numnodes):
        n = Node(**nodeargs)
        netaddr[n.net.address] = n
    nodelist = netaddr.values()

    # connect all nodes via channel
    for n in nodelist:
        for m in nodelist:
            if (n is not m):
                ch.add_edge(n.cif, m.cif, **chargs)

    # create MasterControl
    tlimit = options.numtopo
    dlimit = options.numdata
    if tlimit is not None: MasterControl.TLIMIT = tlimit
    if dlimit is not None: MasterControl.DLIMIT = dlimit
    mc = MasterControl(src, dst, netaddr, ch, border)

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
    parser.add_option("", "--num-topo", dest="numtopo", type="int", \
            default=None, help="Set number of random topologies [default=%default]")
    parser.add_option("", "--num-data", dest="numdata", type="int", \
            default=None, help="Set number of data per topology [default=%default]")
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
            default=None, help="Specify pathloss exponent [default=%default]")
    parser.add_option("", "--use-doppler", dest="usedoppler", action="store_true",  \
            default=False, help="Enable doppler filter for fading in TGn channel model.")
    parser.add_option("", "--disable-fading", dest="usefading", action="store_false",  \
            default=True, help="Normalize channel and remove impact of fading on pathloss in TGn channel model.")
    parser.add_option("-E", "--environment-speed", dest="envspeed", type="float", \
            default=1.2, help="Environmental speed in (km/hr) [default=%default]")
    parser.add_option("", "--bidirectional-channel", dest="bidirectional", action="store_true", \
            default=True, help="Use bidirectional links in channel [default=%default]")
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
