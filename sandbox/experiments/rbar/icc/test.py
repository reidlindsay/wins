#!  /usr/bin/env python

"""
Simulate RBAR over a single link.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-10-21 15:51:37 -0500 (Fri, 21 Oct 2011) $
* $LastChangedRevision: 5246 $

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

from wins.mac import RBAR, ARF

import sys
from optparse import OptionParser

import numpy as np
import struct
import gc
import time

RNG_INIT = 1

class AGT(Packet):
    fields_desc = [IntField("id", 0), IntField("seqno", 0)]

class Agent(Element):
    """Act as a source and sink.
    
    This module keeps the MAC RXU queue backlogged with data to send.
    """
    name = "agent"
    tracename = "AGT"
    DefaultPacketLength = 1500
    MaxSeqNo = 4096
    TrafficModes = ["backlog", "cbr", "poisson"]
    def __init__(self, *args, **kwargs):
        """Constructor."""
        self.nsent, self.nrcvd = 0, 0
        self._seqno = 0
        self.payload = None
        self.dest, self.plen = None, None
        self.recv, self.send = SimEvent(), SimEvent()
        self.cmodel, self.delay = None, None
        self.mode = None
        Element.__init__(self, *args, **kwargs)

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
        ismode = isinstance(mode, str) and (mode.lower() in self.TrafficModes)
        if ismode: mode = mode.lower()
        else:      mode = "backlog"
        self._seqno = 0
        self.dest = dest
        self.delay = delay
        self.mode = mode
        # initialize payload buffer
        self.plen = plen
        x = np.arange(self.plen*2)%128
        args = tuple(x.tolist())
        self.payload = struct.pack("b"*len(x), *args)
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
        """SEND state; keep backlog of data in MAC RXU queue."""
        yield hold, fsm, const.EPSILON         # yield to other threads
        if (not self.dest): yield fsm.stop()    # HALT and sleep
        # check parameters
        errmsg = "[AGT]: No valid network address!"
        assert self.addr, errmsg
        errmsg = "[AGT]: No valid MAC!"
        assert self.mac,  errmsg
        macQ = self.mac.RXU
        macQ.monitorQ = True
        # send packet downstream
        idx = self.seqno%self.plen
        pay = self.payload[idx:idx+self.plen]
        p = AGT(id=self.uid, seqno=self.seqno)/pay
        p.setanno('net-dst', self.dest)
        # send packet
        self.nsent += 1
        self.log("snd", p, seqno=self.seqno, nsent=self.nsent)
        self.nextseqno(update=True)       # update counter
        self.send.signal(p)
        yield self.TX.send(fsm, [p])
        # pause agent
        yield fsm.goto(self.PAUSE)

    def PAUSE(self, fsm):
        """PAUSE state; pause agent based on operation mode."""
        yield hold, fsm, 0  # yield to other threads
        # wait for MAC to empty in 'backlog' mode
        if (self.mode=="backlog"):
            # check MAC input queue
            p, moredata = None, False        # need to send more data
            while (not moredata):
                yield waitevent, fsm, macQ.deQ
                p = macQ.deQ.signalparam
                if (macQ.length<1):
                    moredata = True
        # if delay not specified -> compute delay based on channel
        if (self.delay is None) and self.cmodel.chan:
            self.delay = 2*self.cmodel.chan.coherencetime()
        # wait for prescribed delay
        if self.delay:
            # pause based on operation mode
            # -> 'poisson' and 'backlog' use exponentail backoff
            if (self.mode=="cbr"): d = self.delay
            else:                  d = np.random.exponential(self.delay)
            self.log("WAIT", p, delay=d, avgdelay=self.delay)
            yield hold, fsm, d
        # continue in SEND
        yield fsm.goto(self.SEND)

    def RECV(self, fsm):
        """RECV state; Receive traffic."""
        # get packet from lower layer
        yield self.RX.recv(fsm, 1)
        p = fsm.got[0]
        self.recv.signal(p)
        # check packet parameters
        assert isinstance(p, AGT) and p.hasanno('net-dst')
        dst = p.getanno('net-dst')
        isforme = (dst==self.addr) or (dst==self.broadcast)
        errmsg = "[AGT]: Got packet not intended for me!"
        assert (isforme), errmsg
        # log receive
        self.nrcvd += 1
        kwargs = self.get_agt_anno(p)
        self.log("rcv", p, **kwargs)
        # conintue in RECV
        yield fsm.goto(self.RECV)

    def get_address(self, broadcast=False):
        """Return network address of parent node."""
        addr = None
        if isinstance(self.parent, Base):
            if self.parent.haschild('net'):
                addr = self.parent.net.address
                if broadcast: addr = self.parent.net.broadcast
        return addr

    def get_mac(self):
        """Return MAC of parent node."""
        mac = None
        if isinstance(self.parent, Base):
            if self.parent.haschild('mac'):
                mac = self.parent.mac
        return mac

    def nextseqno(self, update=False):
        """Increment sequence number."""
        newseqno = (self.seqno+1)%self.MaxSeqNo
        if update: self._seqno = newseqno
        return newseqno

    def get_agt_anno(self, p):
        """Get common annotations for logging data."""
        kwargs = {}
        if not isinstance(p, AGT): return kwargs
        kwargs['id'] = p.id
        kwargs['seqno'] = p.seqno
        kwargs['nrcvd'] = self.nrcvd
        if p.hasanno('cif-collision'):
            kwargs['cif-collision'] = [c for c in p.getanno('cif-collision')]
        if p.hasanno('dot11n-sinr'):
            kwargs['sinr'] = "%.4f dB"%(p.getanno('dot11n-sinr') )
        return kwargs

class Node(Element):
    name = "node"
    tracename = "NODE"
    def __init__(self, **kwargs):
        Element.__init__(self, **kwargs)
    def configure(self, rate=None, pos=None, dest=None, useshared=False, \
                        plen=None, cfocorrection=True, **kwargs):
        cif = self.newchild('cif', Dot11NRadio)
        phy = self.newchild('phy', Dot11NPHY, radio=cif, cfocorrection=cfocorrection)
        mac = self.newchild('mac', RBAR, phy=phy)
        net = self.newchild('net', Routing)
        arp = self.newchild('arp', ARP, useshared=useshared)
        agt = self.newchild('agent', Agent, dest=dest, plen=plen)
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
    Limit = 200
    def __init__(self, *args, **kwargs):
        """Constructor."""
        self.src, self.dst = None, None
        self.nrcvd, self.idx = 0, None
        Element.__init__(self, *args, **kwargs)

    def configure(self, src=None, dst=None, snr=None, pos=None, **kwargs):
        """Configure parameters."""
        assert isinstance(src, Node)    # sender is valid Node
        assert isinstance(dst, Node)    # receiver is valid Node
        self.src, self.dst = src, dst
        self.snr, self.pos = snr, pos
        self.nrcvd, self.nsent = 0, 0
        self.idx = None
        f = self.newchild("ctrl", FSM, tracename=self.tracename+".CTRL")
        f.goto(self.MOVE)

    def MOVE(self, fsm):
        """Initialize position of destination."""
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
        # wait for packet to be received by Agent
        self.nrcvd, self.nsent = 0, 0
        agt = self.dst.agent
        while (self.nrcvd < self.Limit) and (self.nsent<self.Limit*20):
            yield waitevent, fsm, (agt.recv, agt.send)
            if   (agt.recv in fsm.eventsFired): self.nrcvd += 1
            elif (agt.send in fsm.eventsFired): self.nsent += 1
        # advance position
        yield fsm.goto(self.MOVE)

    def SHUTDOWN(self, fsm):
        """Shutdown sender when all data has been collected."""
        self.src.agent.dest = None
        yield fsm.stop()        # HALT

def run_experiment(options):
    # record start time
    starttime = time.time()
    # initialize RNG
    if RNG_INIT: RNG_init()
    # experiment parameters
    ntx, nrx = 1, 1
    if options.plen>0: plen = options.plen
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
    verbose = 12
    useshared = True
    kwargs   = {'verbose':verbose, 'useshared':useshared}
    nodeargs = {'cfocorrection':cfocorrection, 'plen':plen}
    nodeargs.update(kwargs)
    # set up channel
    ch = Channel(model=Dot11NChannel, **kwargs)
    modeltype  = options.tgnmodel       # default -> LOS Channel
    alpha      = options.alpha
    usedoppler = options.usedoppler
    usefading  = options.usefading
    chargs = {'modeltype':modeltype, 'n':alpha, \
              'usedoppler':usedoppler, 'usefading':usefading}
    # calculate positions of receiver based on snr parameters
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
    snr = snrmin + numpy.arange(int((snrmax+const.EPSILON+dsnr-snrmin)/dsnr))*dsnr
    snrlist, pos = snr.tolist(), [(0,0)]    # initial point is sender
    snrlist.reverse()
    for s in snrlist:
        x = refdist*(db2linear(PL-s)**(1.0/alpha))
        pos.append((x,0))
    numpos = len(snrlist) + 1
    # create topology layout
    xmin = min([x for (x,y) in pos])
    xmax = max([x for (x,y) in pos])
    ymin = min([y for (x,y) in pos])
    ymax = max([y for (x,y) in pos])
    border = (xmin, xmax, ymin, ymax)
    layout = pos
    assert (len(layout)==numpos)
    topo = {'border':border, 'layout':layout}
    # set up receiver
    pos = layout[1]
    receiver = Node(pos=pos, **nodeargs)
    # set up sender
    pos  = layout[0]
    dest = receiver.net.address
    sender = Node(pos=pos, dest=dest, **nodeargs)
    # connect to channel and set up sender
    cm = ch.add_edge(sender.cif, receiver.cif, **chargs)
    ch.add_edge(receiver.cif, sender.cif, **chargs)
    sender.agent.cmodel = cm
    sender.net.addroute(dest)
    # motion control
    mc = MotionControl(src=sender, dst=receiver, \
                       snr=snrlist, pos=layout[1:], **kwargs)
    # start monitor
    if options.monitor:
        mon = Monitor(period=stoptime/1e4)
        mon.start()

    # run simulation - log relavent information
    ch.log("MODEL", **chargs)
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
    parser.add_option("", "--mcs", dest="mcs", type="int", \
            default=0, help="Set rate index for MCS [default=%default]")
    parser.add_option("-l", "--packet-length", dest="plen", type="int", \
            default=1024, help="Set packet size in bytes [default=%default]")
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
