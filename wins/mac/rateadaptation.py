#!  /usr/bin/env python

"""
Rate Adaptation module; contains `RateAdapt` and `RAI` classes.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-08-25 16:51:48 -0500 (Thu, 25 Aug 2011) $
* $LastChangedRevision: 5115 $

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

from wins.element import Element
from wins.packet import Packet

from scapy.all import IntField

class RateAdapt(Element):
    """Maintains table for current rate being used to send data.

    Indices into the table should be treated like keys in a dictionary. Although
    more complex indices can be used (e.g. tuples, hashes), the most common
    indices will likely be simple strings (e.g. addresses, readable strings).

    Subclasses of this module should overload `get_rate` and `set_rate` to
    enable more complex adaptation policies.

    :CVariables:
     * `db`: Database (i.e. table) containing (key,rate) pairings.

    :IVariables:
     * `base`: Base rate used as default value for (key,rate) pairs.
    """
    name = "rate adapation"
    tracename = "RATE"
    def __init__(self, **kwargs):
        """Constructor."""
        self.base = None
        self.__db = {}
        Element.__init__(self, **kwargs)

    db = property(fget=lambda self: self.__db)

    def configure(self, base=None, **kwargs):
        """Initialize parameters for `RateAdapt`.
        
        :param base: Base rate used as default value for (key,rate) pairs.
        """
        self.base = base

    def get_rate(self, key):
        """Use `key` to determine the appropriate rate.

        :param key: Key to new or existing entry.

        If `key` is not found, this method will add a new entry to the table for
        `key` using the `base` rate and return that rate.
        """
        # check if key is in db
        if key not in self.db:
            self.set_rate(key, None)
        assert (key in self.db)
        # get rate
        rate = self.db[key]
        if rate is None: rate = self.base
        return rate

    def set_rate(self, key, rate=None):
        """Set rate entry for `key`.

        :param key: Key to new or exisitng entry.
        :param rate: New rate for entry [default=None].

        If `rate` is not specified (or `None`), this method will enter `None`
        into the table, which signifies that the `base` rate will be used.
        """
        self.db[key] = rate

    def __getitem__(self, key):
        """Dictionary semantics for `get_rate()`."""
        return self.get_rate(key)

    def __setitem__(self, key, value):
        """Dictionary semantics for `set_rate()`."""
        return self.set_rate(key, value)

class RAI(Packet):
    """Rate Adaptation Information packet.

    This class can be used to create a `Packet` object to maintain information
    relevant to rate adaptation.
    
    This base class only contains two fields: (i) rate and (ii) length.
    """
    name = "RA-INFO"
    fields_desc = [IntField("rate", 0), \
                   IntField("length", 0)]
    def __init__(self, *args, **kwargs):
        """Constructor."""
        Packet.__init__(self, *args, **kwargs)
