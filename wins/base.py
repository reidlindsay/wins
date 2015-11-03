#!  /usr/bin/env python

"""
Implement common functionality for all `wins` objects ; contains
`Base` and `Reference` classes.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-09-23 12:36:44 -0500 (Fri, 23 Sep 2011) $
* $LastChangedRevision: 5140 $

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

import sys
import copy

class Base(object):
    """Implement basic functionality needed by objects.

    :cvar name: Name of class objects.
    :cvar verbose: Property to access/modify verbose level.
    :cvar children: Property to access children dictionary.

    :ivar parent: Pointer to parent.
    :ivar uid: Unique identification number.
    """
    name = "base"
    def __init__(self, name=None, parent=None, verbose=None, uid=None, **kwargs):
        """Constructor.

        :param name: Name of object.
        :param parent: Pointer to parent.
        :param verbose: Verbose level of output.
        :param uid: Unique identification number.
        """
        if name is None: name = self.__class__.name
        self.name = name
        self.__parent = None
        self.parent = parent
        self.__verbose = verbose
        self.__children = {}
        if uid is None: uid = self._next_uid()
        self.uid = uid

    verbose = property(fget=lambda self: self.__verbose, \
                       fset=lambda self, v: self.set_verbose(v) )
    children = property(fget=lambda self: self.__children)
    container = property(fget=lambda self: self.get_container() )
    parent = property(fget=lambda self: self.get_parent(), \
                      fset=lambda self, p: self.set_parent(p) )

    def get_parent(self):
        """Return pointer to parent."""
        if isinstance(self.__parent, Reference):
            return self.__parent._deref
        else:
            return self.__parent

    def set_parent(self, p):
        """Keep weak reference to parent."""
        if p: self.__parent = Reference(p)
        else: self.__parent = p

    def haschild(self, c):
        """Check if `c` is a child.

        :param c: Child or nickname of child to check for.
        :return: Nickname of child; or 'None' if not found.
        """
        nickname = None
        if c in self.children.keys():
            nickname = c
        elif c in self.children.values():
            idx = self.children.values().index(c)
            nickname = self.children.keys()[idx]
        return nickname

    def delchild(self, c):
        """Remove child from children list.

        :param c: Child or nickname of child to delete.
        :return: Removed child.

        :note: This method may raise a KeyError if no matching child exists.
        """
        key = c
        if self.haschild(c):
            key = self.haschild(c)
            child = self.children[key]
        del self.children[key]
        # disconnect child's parent pointer also
        if hasattr(child, 'parent'): child.parent = None
        return child

    def addchild(self, nickname, child):
        """Add new `child` to set of children.
        
        :param nickname: Nickname of child; used when referring to child.
        :param child: New child to add.

        :note: If the child nickname is already in use, this method will
               overwrite the existing child and call `delchild` as needed.
        """
        if self.haschild(nickname): self.delchild(nickname)
        self.children[nickname] = child
        if hasattr(child, 'parent'):  child.parent = self
        if hasattr(child, 'set_verbose'): child.set_verbose(self.verbose)
        return child

    def getchild(self, c):
        """Get child `c`.

        :param c: Child or nickname of child.
        :return: Matching child or `None` if not found.
        """
        child, cname = None, self.haschild(c)
        if cname is not None:
            child = self.children[cname]
        return child

    def newchild(self, nickname, childclass, *args, **kwargs):
        """Factory method for constructing and adding a child.

        :param nickname: Nickname of new child.
        :param childclass: Class of child to be constructed.
        :param args: Arguments to pass to constructor.
        :param kwargs: Keyword arguments to pass to constructor.

        :returns: newly constructed child.

        Child is added using `addchild()`.

        :note: This method overloads the 'parent' pointer and passes the
               'verbose' level as keyword arguments.
        """
        kwargs['parent'] = self
        if 'verbose' not in kwargs: kwargs['verbose'] = self.verbose
        c = childclass(*args, **kwargs)
        return self.addchild(nickname, c)

    def callchild(self, func, *args, **kwargs):
        """Call a function on children that support it.

        :param func: String; name of function to be called.
        :param args: Arguments passed to function.
        :param kwargs: Keywords passed to function.
        """
        assert isinstance(func, str), \
                "[Base.callchild]: Function name must be a string!"
        for c in self.children.values():
            if hasattr(c, func):
                f = getattr(c, func)
                f(*args, **kwargs)

    def get_container(self):
        """Retrieve top level parent."""
        p = self
        while hasattr(p, 'parent') and (p.parent is not None):
            p = p.parent
        # FIXME: Should it return self as container?
        #if p is self: p = None
        return p

    def set_verbose(self, v, recursive=True):
        """Set verbosity of output.

        :param v: New verbose level.
        :param recursive: Boolean; if true, call `set_verbose()` on all children
                          that support it.
        """
        self.__verbose = v
        if not recursive: return
        # otherwise set children
        self.callchild('set_verbose', v)

    def __getattr__(self, key):
        """Kludge to allow dot access to children."""
        if key in self.children:
            return self.children[key]
        else:
            return object.__getattribute__(self, key)

    @staticmethod
    def stderr(s):
        """Convenience method to write string to stderr."""
        sys.stderr.write(s)

    @staticmethod
    def stdout(s):
        """Convenience method to write string to stdout."""
        sys.stdout.write(s)

    @classmethod
    def _next_uid(cls):
        """Internal method to generate unique IDs."""
        clsname = cls.__name__
        c = 0
        if hasattr(cls, "_%s__count"%(clsname)):
            c = getattr(cls, "_%s__count"%(clsname))
        setattr(cls, "_%s__count"%(clsname), c+1)
        return c

    def __str__(self):
        """Convert to string."""
        return "%s(%s)"%(self.name, self.uid)

    def __repr__(self):
        """Get string representation."""
        return str(self)

import weakref
class Reference(object):
    """
    Weakly reference another object.

    :cvar __rsvd: Reserved keywords that cannot be used by objects being
                  referenced.
    :cvar _deref: Property to access underlying object.

    :ivar _ref: Internal reference for underlying object
    """
    __rsvd = ['_ref', '_deref', '_dereference', '__copy__', '__deepcopy__']
    def __init__(self, x=None):
        """Constructor."""
        if isinstance(x, Reference): x = x._deref
        if x is None: x = self.__class__    # default value
        self._ref = weakref.ref(x)
    _deref = property(fget=lambda self: self._dereference() )
    def _dereference(self):
        """Return weakly referenced object."""
        x = self._ref
        while isinstance(x, Reference): x = x._deref
        return x()
    def __copy__(self, *args):
        """Return new `Reference` to same underlying object."""
        return Reference(self)
    def __deepcopy__(self, memo):
        """Return new `Reference` to same underlying object."""
        r = Reference(self)
        memo[id(self)] = r
        return r
    def __getattr__(self, key):
        """Getter for referenced item."""
        if key in Reference.__rsvd: return object.__getattribute__(self, key)
        else: return getattr(self._deref, key)
    def __getattribute__(self, key):
        """Getter for referenced item."""
        if key in Reference.__rsvd: return object.__getattribute__(self, key)
        else: return self._deref.__getattribute__(key)
    def __setattr__(self, key, val):
        """Setter for attributes of referenced item."""
        if key in Reference.__rsvd: return object.__setattr__(self, key, val)
        else: return setattr(self._deref, key, val)
    def __str__(self):
        """Dereference before converting to string."""
        return "Ref[%s]"%(str(self._deref))
    def __repr__(self):
        """Dereference before converting to representation."""
        return repr(self._deref)
