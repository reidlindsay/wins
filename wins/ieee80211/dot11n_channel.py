#!  /usr/bin/env python

"""
IEEE 802.11n channel models; contains `Dot11NChannel` class.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-11-28 20:13:29 -0600 (Mon, 28 Nov 2011) $
* $LastChangedRevision: 5347 $

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

from SimPy.Simulation import now
from wins.channel.breakpoint import Breakpoint
from wins.channel.radio import Radio

from wins.ieee80211.dot11n_support import DOT11N_CARRIER
from wins.backend.dot11n_backend import *

DOT11N_CHANNEL_VERBOSE = 75

class Dot11NChannel(Breakpoint):
    """IEEE 802.11n channel models for propagation over wireless channel.

    This channel model uses an internal `Dot11N_Channel` from `wins.backend` to
    run the internal workings of the IEEE 802.11n channel model.

    Dot11NChannel and Annotations
    =========================
    This class uses/sets the following annotations:

    ================ ==========================================================
    Name              Description
    ================ ==========================================================
    dot11n-channel    Copy of internal Dot11N_Channel applied to packet.
    ---------------- ----------------------------------------------------------
    dot11n-rmsdelay   RMS delay spread of channel (in ns).
    ---------------- ----------------------------------------------------------
    pathloss          Pathloss (in dB). Note that pathloss has inverse 
                      semantics since it subtracted from the transmitted power.
                      So a higher pathloss implies less 'rxpower'.
    ---------------- ----------------------------------------------------------
    fading            Fading contribution from channel (in dB). If `usefading` 
                      is enabled, the fading component will be incorporated    
                      into the pathloss annotation: PL = PL (prop) - X (fade).
                      Fading has normal semantics since it is considered an
                      additive component of receive power. A higher fade
                      implies more 'rxpower'. However, since this value is
                      incorporated into 'pathloss' it does not need to be
                      directly considered to compute instantaneous 'rxpower'.
    ================ ==========================================================

    :CVariables:
     * `modeltype`: IEEE 802.11n channel model type (A-F) [default=None].
     * `modelnum`: enumeration of IEEE 802.11n channel model type (A-F)
                   [default=DOT11N_TGN_MODEL_A].
     * `chan`: Internal `Dot11N_Channel` that runs (taken from `wins.backend`).
     * `flags`: Property to access flags of internal `Dot11N_Channel`.
     * `environmentspeed`: Property to access environmental speed.
     * `DefaultCutoff`: Default value for SNR threshold used to filter packets.

    :IVariables:
     * `usedoppler`: Boolean flag; if true use doppler filter for internal
                     computation, which can require more computation [default=False].
     * `usefading`: Boolean flag; if false normalize channel to remove impact of
                    fading on pathloss [default=True].
     * `cutoff`: Minimum SNR threshold used to filter out packets to reduce
                 computational overhead.
    """
    name = "802.11n propagation"
    tracename = "80211NCHAN"
    fc = DOT11N_CARRIER
    n  = 3.5
    DefaultCutoff = 0
    def __init__(self, modeltype=None, modelnum=DOT11N_LOS_MODEL, \
                       modelflags=None, usedoppler=False, usefading=True, \
                       environmentspeed=None, cutoff=None, **kwargs):
        """Constructor.

        :param modeltype: IEEE 802.11n channel model type (A-F) [default=None].
        :param modelnum: IEEE 802.11n channel model enumeration
                         [default=DOT11N_LOS_MODEL].
        :param modelflags: IEEE 802.11n channel model flags (e.g.
                           `DOT11N_USE_LOS`, `DOT11N_USE_NLOS`,
                           `DOT11N_LOS_RICE`, `DOT11_NLOS_FADING`, etc.).
        :param usedoppler: Enable doppler filtering of NLOS component in TGn
                           channel model [default=False].
        :param usefading: Boolean flag; if false, normalize channel to remove
                          impact of fading on pathloss [default=True].
        :param environmentspeed: Environmental speed (in km/hr).
        :param cutoff: Minimum SNR threshold used to filter packets.

        Create internal `Dot11N_Channel` from `wins.backend` to run internal
        channel model.
        """
        self.usefading  = usefading     # enable fading (otherwise normalize)
        self.usedoppler = usedoppler    # enable doppler filter
        self._lastupdate = None         # internal parameter for last update time
        self.__modelnum = None
        self.__chan = None              # internal channel implementation
        self._envspeed = environmentspeed   # environment speed for internal channel
        self.cutoff = cutoff            # SNR threshold for filtering packets
        # set flags
        if modelflags is None:
            modelflags = DOT11N_TGN_DEFAULT
        self.__flags = modelflags
        # set modeltype
        if modeltype is None:
            self.__modelnum = modelnum
        else:
            self.set_modeltype(modeltype)
        # set cutoff
        if cutoff is None:
            self.cutoff = self.DefaultCutoff
        # set remaining parameters
        m = self.modelnum
        if  (m==DOT11N_TGN_MODEL_D):
            bpdist = 10
        elif (m==DOT11N_TGN_MODEL_E):
            bpdist = 20
        elif (m==DOT11N_TGN_MODEL_F):
            bpdist = 30
        else:
            bpdist = 5      # Used for TGN Models A, B, and C
        kwargs['bpdist'] = bpdist
        Breakpoint.__init__(self, **kwargs)

    modeltype = property(fget=lambda self: self.get_modeltype() )
    modelnum  = property(fget=lambda self: self.get_modelnum() )
    chan      = property(fget=lambda self: self.get_chan(), \
                         fset=lambda self,*args: self.set_chan(*args) )
    flags     = property(fget=lambda self: self.get_flags(), \
                         fset=lambda self,*args: self.set_flags(*args) )
    environmentspeed = property(fget=lambda self: self.get_environmentspeed(), \
                                fset=lambda self, *args: self.set_environmentspeed(*args))

    def filter(self, p, u, v):
        """Filter out packets that should be dropped by the channel model.

        :param p: Packet to filter.
        :param src: Source `ChannelInterface` that sent `p`.
        :param dst: Destination `ChannelInterface`.
        :return: Reason for drop (string), or `None` if packet was not filtered
                 out by the channel model.

        Filter out packets based on pathloss of `Breakpoint` model and expected
        receive power. This should be conservative so as to not adversely impact
        the interference modeling in the receiver.
        """
        # drop packets based
        assert isinstance(u, Radio), "[DOT11N_CHANNEL]: Endpoint is not a Radio!"
        assert isinstance(v, Radio), "[DOT11N_CHANNEL]: Endpoint is not a Radio!"
        dist = self.distance(u,v)
        pl = self.pathloss(dist, bpdist=self.bpdist, fc=self.fc, n=self.n)
        # update channel
        tupdate = self.update(u,v)
        # compute channel power and pathloss
        Pdb = self.chan.powerdb()
        Pnorm = self.chan.normpowerdb()
        if self.usefading: Xdb = Pdb - Pnorm
        else:              Xdb = 0
        pathloss = pl - Xdb
        # get SNR from destination interface
        snr = v.snr(p, pathloss=pathloss)
        # filter packets with insufficient SNR
        drop = None
        if (self.cutoff is not None):
            if (snr<self.cutoff): drop = "below cutoff SNR (%.3f<%.3f)"%(snr,self.cutoff)
        return drop

    def update(self, u, v):
        """Update channel between `u` and `v` as needed.

        :param u: Transmitting `ChannelInterface`.
        :param v: Receiving `ChannelInterface`.
        """
        assert isinstance(u, Radio), "[DOT11N_CHANNEL]: Endpoint is not a Radio!"
        assert isinstance(v, Radio), "[DOT11N_CHANNEL]: Endpoint is not a Radio!"
        nrx, ntx = v.Nrx, u.Ntx
        if self.chan is None:
            h = Dot11N_Channel(self.modelnum, nrx, ntx, self.flags)
            self.set_chan(h)
            self.set_environmentspeed(self._envspeed)
        # update channel using time passed since last update
        errmsg = "[DOT11N_CHANNEL]: Channel not correctly set!"
        assert isinstance(self.chan, Dot11N_Channel), errmsg
        tnow, tupdate = now(), -1
        if self._lastupdate is not None:
            tupdate = tnow - self._lastupdate
        self.chan.update(tupdate)
        self._lastupdate = tnow
        return tupdate

    def apply(self, p, u, v):
        """Apply breakpoint pathloss model.

        :param p: Packet to modify.
        :param u: Transmitting `ChannelInterface`.
        :param v: Receiving `ChannelInterface`.
        :return: Modified packet.

        This method sets the 'dot11n-channel' annotation after calling
        `Breakpoint.apply()`.
        """
        pkt = Breakpoint.apply(self, p, u, v)
        assert(pkt.hasanno('pathloss') )
        pl = pkt.getanno('pathloss')
        # update internal channel between u and v
        tupdate = self.update(u, v)
        # normalize channel power and update pathloss annotation
        Pdb = self.chan.powerdb()
        self.chan.normalize()
        Pnorm = self.chan.normpowerdb()
        # if fading -> add fade contribution from channel (Xdb) to pathloss (pl)
        # otherwise -> ignore fade contribution and just use normalized channel
        if self.usefading: Xdb = Pdb - self.chan.powerdb()
        else:              Xdb = 0
        pkt.setanno('pathloss', pl - Xdb)
        pkt.setanno('fading', Xdb)
        self.log("PWR", power=self.chan.powerdb(), pathloss=pl-Xdb)
        # copy channel and set annotation
        ch = self.chan + 0
        pkt.setanno('dot11n-channel', ch)
        pkt.setanno('dot11n-channel-fading', Xdb)
        pkt.setanno('dot11n-channel-power', Pdb-Pnorm)
        pkt.setanno('dot11n-rmsdelay', ch.rmsdelay() )
        pkt.setanno('dot11n-tgn-model', ch.model())
        return pkt

    def get_modeltype(self):
        """Get IEEE 802.11n channel model type."""
        m = self.modelnum
        tgnenum = {'A': DOT11N_TGN_MODEL_A, \
                   'B': DOT11N_TGN_MODEL_B, \
                   'C': DOT11N_TGN_MODEL_C, \
                   'D': DOT11N_TGN_MODEL_D, \
                   'E': DOT11N_TGN_MODEL_E, \
                   'F': DOT11N_TGN_MODEL_F}
        t = "non-std"
        if m in tgnenum.values():
            for (k,v) in tgnenum.items():
                if (v==m): t = k
        return t

    def set_modeltype(self, mtype):
        """Set IEEE 802.11n channel model type from a string (or int)."""
        if not isinstance(mtype, str):
            self.__modelnum = mtype
            return
        # process as a string
        m = mtype.upper()[0]
        tgnenum = {'A': DOT11N_TGN_MODEL_A, \
                   'B': DOT11N_TGN_MODEL_B, \
                   'C': DOT11N_TGN_MODEL_C, \
                   'D': DOT11N_TGN_MODEL_D, \
                   'E': DOT11N_TGN_MODEL_E, \
                   'F': DOT11N_TGN_MODEL_F}
        v = DOT11N_LOS_MODEL    # use LOS model by default
        if m in tgnenum:
            v = tgnenum[m]
        self.__modelnum = v

    def get_modelnum(self):
        """Get IEEE 802.11n channel model enumeration."""
        return self.__modelnum

    def get_chan(self):
        """Get internal `Dot11N_Channel`."""
        return self.__chan

    def set_chan(self, ch):
        """Set internal `Dot11N_Channel`."""
        assert isinstance(ch, Dot11N_Channel), \
                "[DOT11NCHANNEL]: Error setting internal Dot11N_Channel!"
        self.__chan = ch

    def get_flags(self):
        """Get flags for internal `Dot11N_Channel`."""
        f = self.__flags
        if f is not None:
            if self.usedoppler: f |= DOT11N_NLOS_DOPPLER
            else:               f &= ~DOT11N_NLOS_DOPPLER
        return f

    def set_flags(self, flags):
        """Set flags for internal `Dot11N_Channel`."""
        self.__flags = flags
        # set flags for internal channel
        if isinstance(self.chan, Dot11N_Channel):
            self.chan.set_flags(self.flags)

    def get_environmentspeed(self):
        """Get environmental speed in km/hr from internal `Dot11N_Channel`."""
        if isinstance(self.chan, Dot11N_Channel):
            self._envspeed = self.chan.environmentspeed()
        return self._envspeed

    def set_environmentspeed(self, v=None):
        """Set environmental speed in km/hr [default=1.2 km/hr]."""
        if v is None: v = 1.2
        self._envspeed = v
        if isinstance(self.chan, Dot11N_Channel):
            self.chan.set_environmentspeed(v)

    def log(self, event=None, p=None, *args, **kwargs):
        """Overloaded to check verbose level and set common annotations."""
        force = False
        if ('verbose' in kwargs): force = (kwargs['verbose']>DOT11N_CHANNEL_VERBOSE)
        if self.verbose>DOT11N_CHANNEL_VERBOSE or force:
            Breakpoint.log(self, event, p, *args, **kwargs)
