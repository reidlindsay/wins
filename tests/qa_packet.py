#!  /usr/bin/env python

"""
Runs Packet tests.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2010-10-19 16:34:02 -0500 (Tue, 19 Oct 2010) $
* $LastChangedRevision: 4811 $

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

import unittest

from scapy.all import *
from wins import *
from copy import copy, deepcopy

class TestPacket(unittest.TestCase):

    def setUp(self):
        Trace.Global.reset()

    def not_test_001(self):
        p = Packet()
        h = Raw("helloworld")
        h.setanno("x", 100)
        r = TCP()/h
        s = TCP()/h
        s.setanno("x", 50)
        r.setanno("foo", s, ref=True)
        foo = r.getanno("foo")
        self.assertFalse(isinstance(foo, Reference), "setanno() error!")
        self.assertTrue(foo==s,"getanno() error!")
        self.assertTrue(foo is s,"getanno() error!")
        t = copy(r)
        u = deepcopy(s)
        self.assertEqual(t._id, r._id, "packet copy() error!")
        self.assertEqual(u._id, s._id, "packet deepcopy() error!")
        self.assertFalse(t is r, "packet copy() error!")
        self.assertTrue(t == r, "packet copy() error!")
        self.assertFalse(u == s, "packet deepcopy() error!")

    def not_test_002(self):
        p = Packet()/"helloworld"
        q = Packet()/"helloworld"
        p.setanno("x", q, ref=True)
        crcupdate(p)
        s = deepcopy(p)

    def test_crc(self):
        p = Packet()/"helloworld"
        crcupdate(p)
        p[CRC32].crcerror = 1

