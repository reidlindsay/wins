#!  /usr/bin/env python

"""
Contains class for creating and modifying cyclic-redundancy check (CRC);
contains `CRC32` class and helper functions.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-10-19 17:02:26 -0500 (Wed, 19 Oct 2011) $
* $LastChangedRevision: 5219 $

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

from wins.packet import Packet, ANNO, HiddenField
from scapy.all import SignedIntField, ActionField, crc32

class CRC32(Packet):
    """
    32-bit cyclic redundancy check (or CRC-32).

    This class maintains a 4-byte CRC value that can be appended at the end of a
    Packet; however it should not be directly instantiated. Use the
    accessor/modifier methods `crcupdate()` and `crcremove()` instead.

    This class contains two fields, which can be accessed by Packets containing
    a `CRC32` layer:

        * 'crcvalue' - `SignedIntField` for the CRC-32 value.
        * 'crcerror' - `HiddenField` wrapping a `SignedIntField` that contains an
          error condition (e.g. 0 indicates no error).
    """
    name = "CRC-32"
    fields_desc = [SignedIntField('crcvalue', 0), \
                   ActionField(HiddenField(SignedIntField('crcerror', 0) ), \
                       '_setcrcanno') ]

    def __init__(self, *args, **kwargs):
        """Constructor.

        A `CRC32` should not be directly instantiated. Use the accessor/modifier
        methods `crcupdate()` and `crcremove()`.
        """
        Packet.__init__(self, *args, **kwargs)
        self.hide_defaults()
        # rewrite crcerror so that annotation is set
        c = self.crcerror
        self.crcerror = c

    def _setcrcanno(self, val, *args, **kwargs):
        """Internal method to set 'crcerror' annotation from packet field and
        increment 'crcvalue' if error is set."""
        self.setanno('crcerror', val)
        if (val!=0): self.crcvalue += 1     # if error, then corrupt crcvalue

    @staticmethod
    def calculate(p):
        """Calculate CRC-32 of a string or Packet.

        :param p: Packet or string used to calculate CRC-32.

        :return: Computed CRC-32.
        """
        # get string without CRC
        s = CRC32._noncrcstring(p)
        # compute CRC-32 from string
        crc = crc32(s)
        return crc

    @staticmethod
    def supported(p):
        """Check if packet `p` has a `CRC32` layer."""
        hascrc = isinstance(p, Packet) and p.haslayer(CRC32)
        return hascrc

    @staticmethod
    def _noncrcstring(p, crc=False):
        """Internal method to extract a string from a packet without the CRC.

        :param crc: Boolean flag; if True, include CRC in returned string
                    [default=False].
        """
        # error messages
        errornonstr = "[CRC32]: _noncrcstring() could not find string!"
        # convert to string
        s = p
        if isinstance(p, Packet):
            # detach anno
            anno = None
            if p.haslayer(ANNO):
                anno = p[ANNO]
                anno.underlayer.remove_payload()
            # copy packet and convert to string
            pkt = p.copy()
            scrc = ""
            if p.haslayer(CRC32):
                if crc: scrc = str(pkt[CRC32])
                del pkt[CRC32]
            s = str(pkt) + str(scrc)
            # reattach anno
            if anno: p.add_payload(anno)
        assert isinstance(s, str), errornonstr
        return s

    @staticmethod
    def haserror(p, force=False):
        """Check if CRC of a packet `p` indicates that an error has occurred.

        :param p: `Packet` (or string) to be checked.
        :param force: Boolean; if True, packet `p` is converted to a string and
                      the computed CRC-32 is compared against the existing one.
        :return: Boolean; True if `CRC32` has an error; False otherwise.

        Three normal cases of operation exist for this method:

            1. *p is a `Packet` with a `CRC32` layer.* This method will check
               the `CRC32` layer for an error.

            #. *p is a `Packet` without a `CRC32` layer.* This method will check
               for a 'crcerror' annotation and report this value.

            #. *p is a String.* The last four bytes of the string are assumed to
               be a 32-bit CRC. This method will compute the CRC from the rest
               of the string and compare it to the reported CRC value.

        Otherwise this method will raise a RuntimeError.
        """
        # error messages
        errorshort = "[CRC32]: Packet less than 32 bits long!"
        errornocrc = "[CRC32]: haserror() could not find CRC!"
        # check for error
        error = 1
        hascrc   = CRC32.supported(p)
        hasanno  = isinstance(p, Packet) and p.hasanno('crcerror')
        isstring = isinstance(p, str)
        if force or isstring:
            # convert to string and check
            s = CRC32._noncrcstring(p, crc=True)
            v = CRC32._noncrcstring(p)
            assert len(s)>len(CRC32()), errorshort
            acrc = CRC32.calculate(s[:-4])  # computed crc
            rcrc = CRC32(s[-4:]).crcvalue   # reported crc
            error = acrc - rcrc
        elif hascrc:
            # check CRC32 error field
            rcrc = p[CRC32].crcvalue
            acrc = CRC32.calculate(p)
            error = acrc - rcrc
            p[CRC32].crcerror = (error != 0) # kludge to set annotation
            # FIXME: This might be an easier way to do it. Just look at the
            #        error field instead of recomputing the CRC.
            #error = p[CRC32].crcerror
        elif hasanno:
            # get crcerror from ANNO
            error = p.getanno('crcerror')
        else:
            raise RuntimeError, errornocrc
        # (error = 0) -> no error!
        return (error != 0)

def crcupdate(p, *args, **kwargs):
    """Add/update `CRC32` for Packet `p`.

    :param p: Packet to modify.
    :param args: Arguments passed to `CRC32` constructor if needed.
    :param kwargs: Keywords passed to `CRC32` constructor if needed.

    If `p` already has a `CRC32`, this method will overwrite the current
    `CRC32`; otherwise it will append a new one.

    :returns: Modified `Packet` with updated `CRC32`.
    """
    # error messages
    errornonpacket = "[CRC32]: crcupdate() must have Packet!"
    errornonanno = "[CRC32]: CRC-32 is not at end of Packet (before ANNO)!"
    # continue processing packet
    assert isinstance(p, Packet), errornonpacket
    if not p.haslayer(CRC32):
        p.add_payload(CRC32(*args, **kwargs) )
    crcvalue = CRC32.calculate(p)
    # update CRC value
    if p[CRC32].payload:
        assert isinstance(p[CRC32].payload, ANNO), errornonanno
    p[CRC32].crcvalue = crcvalue
    p[CRC32].crcerror = 0
    return p

def crcremove(p):
    """Remove `CRC32` from Packet `p`.

    :returns: Modified `Packet` with `CRC32` removed.

    :note: If `p` does not have a `CRC32`, this method does nothing.
    """
    if CRC32.supported(p):
        crc = p[CRC32]
        crcerror = crc.crcerror
        crc.crcerror = crcerror
        pay = crc.payload
        crc.remove_payload()
        del p[CRC32]
        p.add_payload(pay)
        p.setanno('crcerror', crcerror)
    return p
