#!  /usr/bin/env python

"""
Runs basic tests.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-05-19 23:11:00 -0500 (Thu, 19 May 2011) $
* $LastChangedRevision: 4986 $

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

from wins import *
from copy import copy, deepcopy

class TestBase(unittest.TestCase):

    def setUp(self):
        pass

    def test_base(self):
        class Derived(Base):
            name = "Derived"
            def __init__(self, *args, **kwargs):
                Base.__init__(self, *args, **kwargs)

        a = Base()
        b = Derived()
        self.assertEqual(a.uid, 0, "ID counter error; %s != %s"%(a.uid, 0))
        self.assertEqual(b.uid, 0, "ID counter error; %s != %s"%(b.uid, 0))

        c = a.newchild('x', Base)
        self.assertTrue(a.x in a.children.values(), "addchild() error!")
        self.assertTrue(a.haschild('x'), "haschild(nickname) error!")
        self.assertTrue(a.haschild(c), "haschild(child) error!")

        f = a.addchild("x", b)
        self.assertTrue(c.parent is None, "delchild() error!")
        self.assertTrue(a.x is f, "addchild() error!")

    def test_reference(self):
        a = Base()
        r = Reference(a)
        s = Reference(r)
        s = copy(r)
        t = deepcopy(r)
        self.assertTrue(a is r._deref, "Reference copy error!")
        self.assertTrue(r._deref is s._deref, "Reference copy error!")
        self.assertTrue(r._deref is t._deref, "Reference deepcopy error!")
        self.assertFalse(r._deref is deepcopy(a), "Reference deepcopy error!")
        self.assertTrue(isinstance(r._deref, Base), "Reference failed type check!")
