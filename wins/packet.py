#!  /usr/bin/env python

"""
Contains classes for creating and modifying packets; contains `Packet` and
`ANNO`.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-09-19 21:56:57 -0500 (Mon, 19 Sep 2011) $
* $LastChangedRevision: 5134 $

:author: Ketan Mandke <mandke@ece.utexas.edu>

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

Packet objects from Scapy now implicitly support an `ANNO` layer. They also have
an id parameter (i.e. Packet._id), which is copied by `copy()` or `deepcopy()`
operations.

:var PACKET_USE_SEQUENTIAL_IDS:
    Switch which will enable sequential IDs when turned on.

:var MAX_PACKET_ID:
    Maximum packet ID; sequence numbers wrap after this limit.
"""
__docformat__ = "restructuredtext en"

from wins.base import Reference
from scapy.all import Packet, Field
from scapy.all import ConditionalField, StrField
from wins import const

from scapy.all import NoPayload

import sys
from copy import deepcopy
import re

##  FIXME: Using sequential IDs?
#   Using sequential IDs is more complex and requires some more complicated
#   modifications to the Packet class. Although this is nice, it might be better
#   just to have unique IDs that are not sequential.
PACKET_USE_SEQUENTIAL_IDS = 1
MAX_PACKET_ID = 65536

class HiddenField:
    """Null or zero length field that maintains a value that will not be added
    to the `Packet` when it is converted to a string.
    """
    fld = None
    def __init__(self, fld):
        self.fld = fld

    def i2m(self, *args, **kwargs):
        return ""

    def addfield(self, pkt, s, val):
        return s+self.i2m(pkt, val)

    def __getattr__(self, attr):
        return getattr(self.fld, attr)

class FuncField(StrField):
    """Allows a function to be assigned to a Packet Field.

    This class is used by `ANNO` to create callable Packet fields. The following
    code shows how to create a custom Packet with a `FuncField`.

        >>> class MYPKT(Packet):
        ...     fields_desc = [FuncField("foo")]
        ...     def __init__(self, *args, **kwargs):
        ...         Packet.__init__(self, *args, **kwargs)
        ...         # create callable method foo that calls private __foo
        ...         self.foo = lambda *args, **kwargs: self.__foo(*args, **kwargs)
        ...     def __foo(self, *args, **kwargs):
        ...         print "foo called"
        ...
        >>> p = MYPKT()
        >>> p.foo()
        foo called

    This example created a new Packet class `MYPKT` and new instance of the
    class `p`. Then called the method `p.foo()`.
    """
    def __init__(self, funcname):
        StrField.__init__(self, funcname, "")

    def i2m(self, pkt, x):
        return ""

    def i2repr(self, pkt, x):
        f = getattr(pkt,self.name)
        return repr(f() )


class ANNO(Packet):
    """Packet annotations implicit to every Scapy Packet object.

    This class implements annotations through a series of `FuncField`, namely
    `setanno`, `getanno`, `delanno`, `hasanno`, `mergeanno`, and `writeanno`.
    These methods access/modify the underlying dictionary object used to
    implement `ANNO`.

    The following code gives an example of how to generate a Packet and
    access/modify its annotations:

        >>> p = Packet()
        >>> p.hasanno("x")
        False
        >>> p.setanno("x", 100)
        >>> p
        <Packet |<ANNO anno={'x': 100} |>>
        >>> p.hasanno("x")
        True
        >>> p.delanno("x")
        <Packet |<ANNO anno={} |>>
        >>> p.hasanno("x")
        False

    Some common annotations used by other classes are listed in the table below.
    These can be used in a variety of ways and are not set by default. Note that
    these special annotation are kept by `Reference` automatically.

    ================ ===========================================================
    Annotation Name  Description/Usage
    ============================================================================
    `parent`         Denotes the object that created the Packet. When setting
                     this annotation, use the ``ref=True`` option with
                     `setanno()`. This will maintain a `Reference` to the parent
                     object rather than a direct pointer. This prevents
                     unwanted copies from being created.
    ---------------- -----------------------------------------------------------
    `root`           Maintain a pointer to the *underlayer* that might be
                     considered the root of a packet. For instance, an IP Packet
                     might encapsulate a TCP packet. So we may want to mark that
                     the *root* of the IP packet is the TCP packet. It may or
                     may not be desirable to use the ``ref=True`` option when
                     setting this annotation with `setanno()`.
    ================ ===========================================================

    :IVariables:
     * `getanno`: `FuncField` that can be called to retrieve an annotation; see
       `__getanno` for more.

     * `setanno`: `FuncField` that can be called to set an annotation; see
       `__setanno` for more.

     * `hasanno`: `FuncField` that can be called to check if an annotation is
       present; see `__hasanno` for more.

     * `delanno`: `FuncField` that can be called to remove an annotation; see
       `__delanno` for more.

     * `mergeanno`: `FuncField` that can be called to merge the current set of
       annotations with another packet's annotation or with a dictionary of
       annotations; see `__mergeanno` for more.

     * `writeanno`: `FuncField` that can be called to overwrite the current set
       of annotations with another packet's annotation or with a dictionary of
       annotations; see `__writeanno` for more.

     * `_anno`: Private dictionary object to implement `ANNO` functionality.

    :cvar fields_desc: Internal list of Packet fields used to implement
                       interface to `ANNO`.

    :cvar __specialanno: Private list of special annotations maintained by
                         `Reference` in `ANNO`.
    """
    fields_desc = [ConditionalField(FuncField('getanno'), lambda *args: False), \
                   ConditionalField(FuncField('setanno'), lambda *args: False), \
                   ConditionalField(FuncField('hasanno'), lambda *args: False), \
                   ConditionalField(FuncField('delanno'), lambda *args: False), \
                   ConditionalField(FuncField('mergeanno'), lambda *args: False), \
                   ConditionalField(FuncField('writeanno'), lambda *args: False), \
                   ConditionalField(FuncField('setprivanno'), lambda *args: False), \
                   FuncField('anno'), \
                   FuncField('privanno') ]
    __specialanno = ['root', 'parent']

    def __init__(self, *args, **kwargs):
        """Constructor."""
        self._anno = {}
        self.__priv = {}
        Packet.__init__(self, *args, **kwargs)
        self.hide_defaults()
        self._setup()

    def _setup(self):
        """Sets `FuncField` fields of packet so that they are callable and point
        to the appropriate internal method."""
        self.getanno = lambda *args, **kwargs: self.__getanno(*args, **kwargs)
        self.setanno = lambda *args, **kwargs: self.__setanno(*args, **kwargs)
        self.hasanno = lambda *args, **kwargs: self.__hasanno(*args, **kwargs)
        self.delanno = lambda *args, **kwargs: self.__delanno(*args, **kwargs)
        self.mergeanno = lambda *args, **kwargs: self.__mergeanno(*args, **kwargs)
        self.writeanno = lambda *args, **kwargs: self.__writeanno(*args, **kwargs)
        self.setprivanno = lambda *args, **kwargs: self.__setprivanno(*args, **kwargs)
        self.anno = lambda *args, **kwargs: self._anno
        self.privanno = lambda *args, **kwargs: self.__priv

    def __getanno(self, key):
        """Internal method for retrieving an annotation `key`; *do not call
        directly*, instead use `getanno`."""
        if key in self.__priv:
            if isinstance(self.__priv[key], Reference): return self.__priv[key]._deref
            else: return self.__priv[key]
        if isinstance(self._anno[key], Reference): return self._anno[key]._deref
        else: return self._anno[key]

    def __setanno(self, key, val, ref=False, priv=False):
        """Internal method for setting an annotation `key`; *do not call
        directly*, instead use `setanno`.

        :Parameters:
         * `key`: Name of annotation to be set.
         * `val`: New value of annotation.
         * `ref`: Optional argument; when True, a `Reference` to `val` is stored
           rather than a direct pointer.

        When `ref` is True a `Reference` to `val` is stored in the internal
        dictionary. This means that when a deepcopy of the Packet occurs, the
        new Packet will contain a `Reference` to the original object.
        """
        if priv: self.__setprivanno(key, val, ref)
        # keep special annotations by reference
        if key in ANNO.__specialanno: ref = True
        # store annotation
        if ref: self._anno[key] = Reference(val)
        else: self._anno[key] = val

    def __setprivanno(self, key, val, ref=False):
        """Internal method for setting private annotation `key`; *do not call
        directly*, instead use `setprivanno`.

        :Parameters:
         * `key`: Name of annotation to be set.
         * `val`: New value of annotation.
         * `ref`: Optional argument; when True, a `Reference` to `val` is stored
           rather than a direct pointer.

        When `ref` is True a `Reference` to `val` is stored in the internal
        dictionary. This means that when a deepcopy of the Packet occurs, the
        new Packet will contain a `Reference` to the original object.
        """
        # keep special annotations by reference
        if key in ANNO.__specialanno: ref = True
        # store annotation
        if ref: self.__priv[key] = Reference(val)
        else: self.__priv[key] = val

    def __hasanno(self, key):
        """Internal method to check if annotation `key` is present. *do not call
        directly*, instead use `setanno`."""
        if key in self.__priv: return True
        return key in self._anno

    def __delanno(self, key):
        """Internal method for removing an annotation `key`; *do not call
        directly* (instead use `delanno`)."""
        if key in self.__priv:
            del self.__priv[key]
            return  # exit
        del self._anno[key]

    # merging annotation will override existing fields
    def __mergeanno(self, p):
        """Internal method for merging an annotation from a Packet or dictionary
        with the current set of annotations; *do not call directly*, instead use
        `mergeanno`.

        :param p: Packet or dictionary of new annotations.

        If `p` is a Packet, this method will extract the new set of annotations
        from the packet using the `anno` method. The new set of annotations will
        be merged with the current set.

        *Annotations from `p` will overwrite any existing annotations of the
        same name*. So if an annotation 'x' existed in `p` and in the current
        set of annotations, the annotation value `p.getanno('x')` would
        overwrite the current annotation value.
        """
        a = p
        priv = None
        if isinstance(p, Packet): a = p.anno()
        if isinstance(p, Packet): priv = p.privanno()
        self._anno.update(a)
        if priv: self.__priv.update(priv)

    def __writeanno(self, a):
        """Internal method for overwriting the entire set of annotations; *do
        not call directly*, instead use `writeanno`."""
        newanno = a
        newpriv = None
        if isinstance(a, Packet): newanno = a.anno()
        if isinstance(a, Packet): newpriv = a.priv()
        if isinstance(newanno, dict): self._anno = newanno
        if isinstance(newpriv, dict): self.__priv = newpriv
        else:                         self.__priv = {}

    def copy(self, *args, **kwargs):
        """Overwrite Packet.copy from Scapy Packet class."""
        clone = Packet.copy(self, *args, **kwargs)
        clone._setup()
        newanno, keepkeys = {}, []
        for k in self._anno:
            keep = True
            keep &= not isinstance(self.getanno(k), list)   # do not copy lists
            keep &= not isinstance(self.getanno(k), dict)   # do not copy dicts
            keep &= not (k in self.privanno())
            keep &= not (k in const.ANNO_BLACKLIST)
            if keep: keepkeys.append(k)
        for k in keepkeys:
            if isinstance(self._anno[k], Reference):
                newanno[k] = Reference(self._anno[k])
            else:
                try: newanno[k] = deepcopy(self._anno[k])    # deepcopy item
                except: sys.stderr.write("%s caused deepcopy error\n"%(k))
        clone.writeanno(newanno)
        return clone

    @staticmethod
    def supported(p):
        """Convenience method to check if object `p` supports `Packet`
        annotations as implemented by `ANNO` class."""
        x = p
        if isinstance(p, Reference): x = p._deref
        return isinstance(x, Packet)

    @staticmethod
    def supports(p, aname):
        """Convenience method to check if object `p` has an `ANNO` sublayer and
        has annotation `aname`.

        This is the same as calling `ANNO.supported()` and `hasanno`.
        """
        x = p
        if isinstance(p, Reference): x = p._deref
        return isinstance(x, Packet) and x.hasanno(aname)


__orig_add_payload = Packet.add_payload

def __anno_add_payload(self, payload, *args, **kwargs):
    """Overwrite Packet.add_payload() for `ANNO` support."""
    global __orig_add_payload
    # add annotation for the first time
    if isinstance(payload, ANNO) and (not self.haslayer(ANNO) ):
        return __orig_add_payload(self, payload, *args, **kwargs)
    # extract current ANNO
    if self.haslayer(ANNO):
        oldanno = self[ANNO]
        del self[ANNO]
    else:
        oldanno = ANNO()
    # merge/add to payload ANNO
    if payload.haslayer(ANNO):
        payload[ANNO].mergeanno(oldanno)
    else:
        payload.add_payload(oldanno)
    # add payload
    return __orig_add_payload(self, payload, *args, **kwargs)

Packet.add_payload = \
        lambda self, payload, *args, **kwargs: __anno_add_payload(self, payload, *args, **kwargs)

__orig_init_fields = Packet.init_fields

def __anno_init_fields(self, *args, **kwargs):
    """Overwrite Packet.init_fields() for `ANNO` and Packet id support."""
    global __orig_init_fields
    __orig_init_fields(self, *args, **kwargs)
    # update Packet counter and set id
    if PACKET_USE_SEQUENTIAL_IDS:
        cls = self.__class__
        clsname = cls.__name__
        cntname = "_%s__count"%(clsname)
        # More complex approach - create unique counter for every class
        if not hasattr(cls, cntname):
            setattr(cls, cntname, 0)
        c = getattr(cls, cntname)
        if not hasattr(self, '_id'):
            self._id = c
            setattr(cls, cntname, (c+1)%MAX_PACKET_ID)
    else:
        # Straight forward approach - does not keep sequential numbering!
        if not hasattr(Packet, '_count'): Packet._count = 0
        self._id = Packet._count
        Packet._count = (Packet._count+1)%(MAX_PACKET_ID)

Packet.init_fields = lambda self, *args, **kwargs: __anno_init_fields(self, *args, **kwargs)

__orig_copy = Packet.copy

def __packet_id_copy(self, *args, **kwargs):
    """Overwrite Packet.copy() for Packet id support."""
    global __orig_copy
    if PACKET_USE_SEQUENTIAL_IDS:
        cls = self.__class__
        clsname = cls.__name__
        cntname = "_%s__count"%(clsname)
        if hasattr(cls, cntname):
            c = getattr(cls, cntname)
            setattr(cls, cntname, (c-1)%MAX_PACKET_ID)
    # make copy and return
    clone = __orig_copy(self, *args, **kwargs)
    if (hasattr(clone, '_id') and hasattr(self,'_id')):
        clone._id = self._id
    return clone

Packet.copy = lambda self, *args, **kwargs: __packet_id_copy(self, *args, **kwargs)
Packet.__deepcopy__ = lambda self, memo: self.copy()

__orig_iter = Packet.__iter__

def __packet_iter(self, *args, **kwargs):
    """Overwrite Packet.__iter__() for Packet id support."""
    global __orig_iter
    if PACKET_USE_SEQUENTIAL_IDS:
        cls = self.__class__
        clsname = cls.__name__
        cntname = "_%s__count"%(clsname)
        if hasattr(cls, cntname):
            c = getattr(cls, cntname)
            setattr(cls, cntname, (c-1)%MAX_PACKET_ID)
    return __orig_iter(self, *args, **kwargs)

Packet.__iter__ = lambda self, *args, **kwargs: __packet_iter(self, *args, **kwargs)

__orig_getfield_and_val = Packet.getfield_and_val

def __packet_getfield_and_val(self, *args, **kwargs):
    """Overwrite Packet.getfield_and_val() to support ANNO."""
    try:
        return __orig_getfield_and_val(self, *args, **kwargs)
    except AttributeError:
        if not self.haslayer(ANNO): self.add_payload(ANNO())
        return __orig_getfield_and_val(self, *args, **kwargs)

Packet.getfield_and_val = lambda self, *args, **kwargs: __packet_getfield_and_val(self, *args, **kwargs)

Packet.flatten = property(fget=lambda self: re.compile("\x1b\[\d{1,2}m").sub("",repr(self)) )
Packet.traceid = property(fget=lambda self: "%s(%s)"%(self.name, self._id) )
Packet.tracename = property(fget=lambda self: self.get_tracename())

Packet.get_tracename = lambda self: self.name
