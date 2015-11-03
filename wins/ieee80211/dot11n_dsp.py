#!  /usr/bin/env python

"""
Implements baseband digital signal processing for waveform-level simulation of
IEEE 802.11n PHY protocol; contains `Dot11N_DSP` class.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-10-19 13:31:24 -0500 (Wed, 19 Oct 2011) $
* $LastChangedRevision: 5214 $

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

:var DOT11N_DSP_VERBOSE:
    Constant enumeration to control verbose thresholds in this file.
  
    Setting the verbose level of a `Dot11N_DSP` above this threshold will cause
    the corresponding output in this class to be written (or logged).
"""
__docformat__ = "restructuredtext en"

from wins.element import Element
from wins.helper import time2usec, linear2db, db2linear, monitor_events
from wins.packet import ANNO
from wins.crc    import CRC32
from wins import const

from wins.ieee80211.dot11n_support import *

from wins.backend import Dot11N_Transmitter, Dot11N_Receiver, Waveform
from wins.backend import hamming_distance

from numpy import sqrt

DOT11N_DSP_VERBOSE = 67

class Dot11N_DSP(Element):
    """Provides API for baseband DSP needed for waveform-level processing of
    IEEE 802.11n packets.

    Dot11N_DSP and Annotations
    =========================
    This class uses/sets the following annotations:

    ==================== ===========================================================
    Name                  Description
    ==================== ===========================================================
    dot11n-start          Start time for packet as determined during detection.
    -------------------- -----------------------------------------------------------
    dot11n-txwaveform     Transmitted waveform.
    -------------------- -----------------------------------------------------------
    cif-duration          Duration of packet is marked by `encode()` if it is
                          not already set.
    -------------------- -----------------------------------------------------------
    collision             List of packets that arrived at the same time as the
                          current packet (i.e. collided with the marked packet).
    -------------------- -----------------------------------------------------------
    phy-rate              Rate index used to specify coding and modulation scheme.
    ==================== ===========================================================

    :IVariables:
     * `transmitter`: Internal transmitter for IEEE 802.11n PHY.
     * `receiver`: Internal receiver for IEEE 802.11n PHY.
     * `rxinput`: Internal receive buffer to maintain input samples to receiver.

    :note: This module should be used with `Dot11NPHY` as its parent and *only*
           when `Dot11NPHY.usewaveform` is set.
    """
    name = "IEEE 802.11n DSP"
    tracename = "DSP"
    MAXINPUT  = DOT11N_MAXINPUT
    MINPAD    = DOT11N_WAVEFORM_PAD
    def __init__(self, **kwargs):
        """Constructor."""
        self.transmitter = Dot11N_Transmitter()
        self.receiver    = Dot11N_Receiver()
        self.rxinput = None
        Element.__init__(self, **kwargs)

    cfocorrection = property(fset=lambda self, x: self.enable_cfo_correction(x))

    def configure(self, cfocorrection=True, **kwargs):
        """Configure parameters and initialize internal variables.

        :param cfocorrection: Boolean flag; if true, enable CFO correction in
                              `receiver` [default=True].
        """
        self.rxinput = None
        self.cfocorrection = cfocorrection

    def enable_cfo_correction(self, enable=True):
        """Enable/Disable CFO correction in `receiver`."""
        assert isinstance(self.receiver, Dot11N_Receiver)
        if enable: self.receiver.enable_cfo_correction()
        else:      self.receiver.disable_cfo_correction()

    def encode(self, p, rate=0, ntx=1):
        """Encode waveform and set corresponding annotation.

        :param p: Packet to encode (i.e. IEEE 802.11n payload).
        :param rate: Rate index.
        :param ntx: Number of transmit antennas.

        This method will check the validity of `rate` and `ntx`. In addition, it
        will also verify that the encoded waveform matches the 'cif-duration'
        annotation (if it is found).
        """
        # verify rate and ntx parameters
        errmsg = "[DOT11N_DSP]: Invalid number of antennas!"
        assert (0<ntx<=4), errmsg
        errmsg = "[DOT11N_DSP]: Invalid rate parameter, %d "%(rate) + \
                 "not in [%d,%d)!"%(0, 8*ntx)
        assert (0<=rate<8*ntx), errmsg
        # set parameters and encode
        m = str(p)
        self.transmitter.set_ntx(ntx)
        self.transmitter.set_rate(rate)
        txwaveform = self.transmitter.encode(m)
        self.currwaveform = txwaveform      # keep local copy of txwaveform
        p.setanno('dot11n-txwaveform', txwaveform, ref=True)
        tduration = txwaveform.cols()/DOT11N_BANDWIDTH
        # check duration annotation
        if p.hasanno('cif-duration'):
            duration = p.getanno('cif-duration')
            errmsg = "[DOT11N_DSP]: Waveform duration does not " + \
                     "match calculated value! " + \
                     "(%.4f usec != %.4f usec)"%(duration*1e6, tduration*1e6)
            assert (abs(tduration-duration)<const.EPSILON), errmsg
        else:
            p.setanno('cif-duration', tduration)
        return p

    def calcrxwaveform(self, p, bandwidth=None):
        """Pre-compute receive waveform and set 'dot11n-rxwaveform'.

        Use 'dot11n-channel' and 'dot11n-txwaveform' annotations to pre-compute
        the received waveform (modulo intereference and noise).
        """
        if bandwidth is None: bandwidth = DOT11N_BANDWIDTH
        # check for annotations
        for a in ['rxpower', 'dot11n-txwaveform', 'dot11n-channel']:
            errmsg = "[DOT11N_DSP]: cannot find '%s' annotation!"%(a)
            assert ANNO.supports(p,a), errmsg
        A = db2linear(p.getanno('rxpower'))
        H = p.getanno('dot11n-channel')
        x = p.getanno('dot11n-txwaveform')
        # check that channel dimensions match Ntx
        nrx = H.rows()
        ntx, hcols = x.rows(), H.cols()
        assert (hcols == ntx), "[DOT11N]: Got channel with invalid dimensions!"
        z = sqrt(A)*H.conv(x)
        # apply local frequency offset (and doppler) if annotations are found
        cfo = 0
        if p.hasanno('tx-cfo') and p.hasanno('rx-cfo'):
            cfo += p.getanno('tx-cfo') - p.getanno('rx-cfo')
        if p.hasanno('doppler'):
            cfo += p.getanno('doppler')
        z.apply_offset(cfo, bandwidth)
        # set rxwaveform annotation
        p.setanno('dot11n-rxwaveform', z)
        return p

    def addrxinput(self, p, **kwargs):
        """Update local waveform with new packet.

        :param p: Packet to add to `rxinput`.
        :param kwargs: Keyword arguments passed to `calcrxwaveform()`.
        """
        # check if waveform has been calculated or already been added to input
        if p.hasanno('dot11n-rxadded'): return
        if not p.hasanno('dot11n-rxwaveform'): self.calcrxwaveform(p, **kwargs)
        # check other annotations
        for a in ['cif-rxts', 'dot11n-rxwaveform', 'noisepower']:
            errmsg = "[DOT11N_DSP]: cannot find '%s' annotation!"%(a)
            assert ANNO.supports(p,a), errmsg
        # get parameters
        y = p.getanno('dot11n-rxwaveform')
        npow = db2linear(p.getanno('noisepower'))
        nrx  = y.rows()
        yts  = p.getanno('cif-rxts')
        errmsg = "[DOT11N]: Cannot add waveform with %d "%(y.cols() ) + \
                 "(> %d) samples in addrxinput()!"%(self.MAXINPUT)
        assert (y.cols()<self.MAXINPUT), errmsg
        # apply pad to front and back of new input
        pad  = Waveform.zeros(nrx,2*self.MINPAD)
        y    = pad.concat_horizontal(y.concat_horizontal(pad) )
        yts -= pad.cols()/DOT11N_BANDWIDTH
        # initialize rxinput if it is not already set
        if self.rxinput is None:
            #self.log("rxadded", p)
            noise =  Waveform.randn(nrx, y.cols(), 0, npow)
            self.rxinput = (yts, y+noise)
            p.setanno('dot11n-rxadded', True)
            return p
        # update rxinput
        (rxts, w) = self.rxinput
        rxend = rxts + w.cols()/DOT11N_BANDWIDTH
        assert (y.rows() == nrx), "[DOT11N]: Invalid row size in addrxinput()!"
        assert (w.rows() == nrx), "[DOT11N]: Invalid row size in addrxinput()!"
        # get indices for begin/end of rxinput and new input
        rxa, rxb = 0, w.cols()-1
        ya = rxa + int((yts-rxts)*DOT11N_BANDWIDTH)
        yb = ya + y.cols() - 1
        istart, istop = min(rxa, ya), max(rxb, yb)
        errmsg = "[DOT11N]: Cannot use addrxinput() to add " + \
                "waveform that starts prior to current input stream!"
        assert (istart==rxa), errmsg
        # new waveform starts after current input ends -> reinitialize input
        if ya>rxb:
            self.rxinput = None
            return self.addrxinput(p)
        # add pads to y as needed
        if ya>istart:
            pad = Waveform.zeros(nrx, ya - istart)
            y = pad.concat_horizontal(y)
            self.log("-DSP.PAD.IN", p, npad=ya-istart)
        if yb<istop:
            pad = Waveform.zeros(nrx, istop - yb)
            y = y.concat_horizontal(pad)
            self.log("+DSP.PAD.IN", p, npad=istop-yb)
        # add pads to (end of) rxinput as needed
        if rxb<istop:
            noise = Waveform.randn(nrx, istop-rxb, 0, npow)
            w = w.concat_horizontal(noise)
            self.log("+DSP.PAD.RX", p, npad=istop-rxb)
        # combine waveforms
        errmsg = "[DOT11N]: Incompatible waveform dimensions, " + \
                 "y ~ %d x %d, "%(y.rows(),y.cols()) + \
                 "w ~ %d x %d, "%(w.rows(),w.cols()) + \
                 "y: [%d, %d], w: [%d, %d]"%(ya,yb,rxa,rxb)
        assert (y.cols() == w.cols() ), errmsg
        assert (y.rows() == w.rows() ), errmsg
        assert (y.cols()>0), errmsg
        assert (y.rows()>0), errmsg
        z = y + w
        # trim input if necessary
        if z.cols()> self.MAXINPUT:
            trimidx = z.cols() - self.MAXINPUT
            z = z.get_cols(trimidx, z.cols()-1)
            rxts += trimidx/DOT11N_BANDWIDTH
        self.rxinput = (rxts, z)
        p.setanno('dot11n-rxadded', True)
        #self.log("rxadded", p)

    def delrxwaveform(self, p):
        """Remove waveform annotations from packet `p`."""
        wanno = ['dot11n-txwaveform', 'dot11n-rxwaveform', \
                 'dot11n-rxnoise', 'dot11n-rxadded', 'dot11n-channel']
        for a in [w for w in wanno if ANNO.supports(p,w)]:
            if p.hasanno(a): p.delanno(a)
        return p

    def getrxinput(self, p):
        """Get input to receiver for specified packet.

        :param p: Packet to grab from input stream.

        Get the samples corresponding to the input for a given packet. This uses
        the timestamp annotation 'cif-rxts' and 'cif-duration' annotation to
        find the appropriate input samples from `rxinput`.
        """
        #self.stdout("%s: called getrxinput() on %s @ %.8f\n"%(self.traceid, p.traceid, now()))
        errmsg = "[DOT11N]: in getrxinput(), no input stream found!"
        assert (self.rxinput is not None), errmsg
        # check annotations
        for a in ['cif-rxts', 'cif-duration', 'dot11n-rxadded']:
            errmsg = "[DOT11N]: Cannot find '%s' annotation!"%(a)
            assert ANNO.supports(p, a), errmsg
        assert p.getanno('dot11n-rxadded')
        # calculate start/end of waveform
        tstart = p.getanno('cif-rxts')
        duration = p.getanno('cif-duration')
        tend = tstart + duration
        # pad returned waveform
        minpad = self.MINPAD
        tstart -= minpad/DOT11N_BANDWIDTH
        tend   += minpad/DOT11N_BANDWIDTH
        # extract from input stream
        rxts, w = self.rxinput
        rxend = rxts + w.cols()/DOT11N_BANDWIDTH
        errmsg = "[DOT11N]: in getrxinput(), packet starts before input " + \
                 "stream! tpacket, tbuffer = (%.6f, %.6f)"%(tstart, rxts)
        assert (tstart>=rxts), errmsg
        # calculare indices
        istart = int((tstart-rxts)*DOT11N_BANDWIDTH)
        istop  = int((tend-rxts)*DOT11N_BANDWIDTH)-1
        assert (istart>-1), "[DOT11N]: In getrxinput(), requested " + \
                "packet starting before start of input stream!"
        assert (istop<=w.cols()), "[DOT11N]: In getrxinput(), requested " + \
                "packet beyond end of input stream!"
        # get rxinput
        y = w.get_cols(istart, istop)
        return y

    def decode_header(self, p, checkcrc=False):
        """Waveform-level decoding of header.

        :param p: Packet to decode.
        :param checkcrc: If true, check header CRC.

        Set `checkcrc` to false when using `decode_header()` to just do
        detection. By default, the PHY does detection & decoding in the same
        functions, so this method is allows for some slight variation in usage.
        """
        # check annotations
        assert isinstance(p, Dot11N)
        for a in ['cif-rxts', 'dot11n-rxadded']:
            errmsg = "[DOT11N]: Cannot find '%s' annotation!"%(a)
            assert ANNO.supports(p, a), errmsg
        errmsg = "[DOT11N]: Packet was not added to input!"
        assert p.hasanno('dot11n-rxadded'), errmsg
        # get rx input stream
        y = self.getrxinput(p)
        detect = self.receiver.decode_header(y, checkcrc)   # check header crc?
        prepad = int(self.MINPAD + DOT11N_TSHORT*DOT11N_BANDWIDTH)
        startidx = self.receiver.start_index() - prepad
        crcok = detect
        # check header parameters
        plen = len(p.payload)
        # get header parameters from receiver after decoding
        rate, length = self.receiver.mcs(), self.receiver.length()
        okrate = (0<=rate<len(DOT11N_DATARATE) )
        oklen  = (length == plen)
        okpar  = crcok
        # calculate start time
        tstart = p.getanno('cif-rxts')
        tstart  += startidx/DOT11N_BANDWIDTH
        p.setanno('dot11n-start-index', startidx)
        p.setanno('dot11n-cfo', self.receiver.get_cfo() )
        # return decoding parameters
        if detect: header = "detect success"
        else:      header = "detect failed"
        if checkcrc:
            if crcok: header = "success"
            else:     header = "decoding failed"
        param = {'status':header, 'detect':detect, 'crcok':crcok, \
                 'rate':rate, 'length':length, 'start':tstart}
        return param

    def decode_data(self, p):
        """Waveform-level decoding of packet payload.

        :param p: Packet to decode.
        :return: Decoded data string.
        """
        # check annotations
        assert isinstance(p, Dot11N)
        for a in ['cif-rxts', 'dot11n-header', 'dot11n-detect', 'dot11n-rxadded']:
            errmsg = "[DOT11N]: Cannot find '%s' annotation!"%(a)
            assert ANNO.supports(p, a), errmsg
        errmsg = "[DOT11N]: Packet was not added to input!"
        assert p.hasanno('dot11n-rxadded'), errmsg
        errmsg = "[DOT11N]: Packet was not detected!"
        assert p.getanno('dot11n-detect'), errmsg
        header = p.getanno('dot11n-header')
        errmsg = "[DOT11N]: Header decoding failed before decode_data()!"
        assert (header=="success"), errmsg
        # get parameters to check startidx
        Ngi      = int(DOT11N_TGI*DOT11N_BANDWIDTH)
        prepad   = int(self.MINPAD + DOT11N_TSHORT*DOT11N_BANDWIDTH)
        startidx = self.receiver.start_index() - prepad
        ## startidx = -8     # overwrite startidx
        ## self.receiver.set_start_index(prepad+startidx)
        # decode payload -> check if startidx is valid
        y = self.getrxinput(p)
        data = self.receiver.decode_data(y)
        # check if decoding error occurred
        #  -> do contents of data match and does CRC pass?
        crcchk, error = "FAIL", True
        if (data == str(p.payload)):
            if not CRC32.haserror(data):
                crcchk = "OK"
                error = False
        # hamming distance between decoded data and packet p
        hdist = hamming_distance(data, str(p.payload) )
        self.log("CRC", p.payload, crcchk=crcchk, error=error, hdist=hdist)
        # clean up waveform annotations
        self.delrxwaveform(p)
        # return decoding parameters
        param = {'data':data, 'error':error, 'nerror':hdist}
        # return decoded data
        return param

    def log(self, event=None, p=None, *args, **kwargs):
        """Overloaded to check verbose level and set common annotations."""
        force = False
        if ('verbose' in kwargs): force = (kwargs['verbose']>DOT11N_DSP_VERBOSE)
        if (self.verbose>DOT11N_DSP_VERBOSE) or force:
            Element.log(self, event, p, *args, **kwargs)
