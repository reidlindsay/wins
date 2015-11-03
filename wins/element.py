#!  /usr/bin/env python

"""
Basic building blocks for all protocol objects; contains `Element` and `Port`
classes.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-02-08 13:52:17 -0600 (Tue, 08 Feb 2011) $
* $LastChangedRevision: 4908 $

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

from wins.trace import Traceable
from wins.queue import Queue

import re
import string

class Port(Queue):
    """Interface for enabling connections between `Element` objects.

    Port Targets
    ============
    Each `Port` has a `target`. By default, each Port points to itself, thus a
    `send()` call on an unconnected Port will result in queueing of objects at
    the Port itself. Multiple Ports may be connected together in a chain. The
    target of any Port in this chain would be the Port at the tail of this
    chain. Ports can be connected using the member method `Port.connect()`.

    :CVariables:
     * `name`: Name of Port.
     * `tracename`: Name for logging in `Trace`.
     * `output`: Property to access `Port` to which this port is connected.
     * `target`: Property to access the target for this port (see above).

    :IVariables:
     * `__output`: Private pointer to connected `Port`.
    """
    name="port"
    tracename="P"
    def __init__(self, **kwargs):
        Queue.__init__(self, **kwargs)
        self.__output = None

    output = property(fget=lambda self: self.__output)
    target = property(fget=lambda self: self._get_target())

    def connect(self, port):
        """Connect to another `Port`."""
        assert isinstance(port, Port), \
                "[PORT]: Can only connect() to another Port!"
        self.__output = port

    def disconnect(self):
        """Disconnect port from its output."""
        self.__output = None

    def _get_target(self):
        """Internal method to determine the target of the Port."""
        p = self
        while isinstance(p.output, Port): p = p.output
        assert isinstance(p, Port), \
                "[PORT]: Error! Port has non-Port target!"
        return p

    def send(self, *args, **kwargs):
        """Send items to `target`.

        :param args: Arguments passed to `insert()` on `target`.
        :param kwargs: Keywords passed to `insert()` on `target`.
        :return: Yield clause from `insert()` call.

        This method is effectively a wrapper for `insert()` on the `target`.
        """
        t = self.target
        return t.insert(*args, **kwargs)

    def recv(self, *args, **kwargs):
        """Receive items from port.

        :param args: Arguments passed to `remove()`.
        :param kwargs: Keywords passed to `remove()`.
        :return: Yield clause from `remove()` call.

        This method is effectively a wrapper for `remove()`.
        """
        return self.remove(*args, **kwargs)

class Element(Traceable):
    """Basic building block for all protocol objects.

    When an `Element` is instantiated, the `configure()` and `start()` method
    may be called (see `__init___()` for more). Subclasses should overload these
    methods as needed, since these methods do nothing by default.

    :cvar started: Property to access internal start flag; if true, it indicates
                   that `start()` has already been called. By default, the
                   `start()` method will not do anything if `started` is true.
                   This parameter should be used to prevent multiple calls to
                   `start()` which may result in errors.
    """
    def __init__(self, start=True, **kwargs):
        """Constructor.

        :param start: Boolean; if true, call `start()` at the end of constructor.
        :param kwargs: Keywords passed to `Traceable` constructor and
                       `configure()` method.

        After this method invokes the `Traceable` constructor, it calls the
        `configure()` method, and then optionally calls the `start()` method.
        All keywords will be passed to the `Traceable` constructor as well as
        the `configure()` method.
        """
        self.__started = False
        Traceable.__init__(self, **kwargs)
        self.configure(**kwargs)
        if start: self.start()

    started = property(fget=lambda self: self.__started)

    def configure(self, **kwargs):
        """Configure `Element` and create components; **overload this method**.

        :param kwargs: Optional keyword arguments passed from constructor.

        By default this method does nothing and should be overloaded.
        """
        pass

    def start(self, **kwargs):
        """Start execution; **overload this method as needed**.

        :param kwargs: Optional keyword arguments.

        By default, this method will call 'start' on all child objects. If
        `start()` has already been called, this method will do nothing.
        """
        if (self.started): return
        self.__started = True
        self.callchild('start', **kwargs)

    def connect(self, p, *args, **kwargs):
        """Convenience method to create connections with the ports of a `p`.

        By default this method does nothing and should be overloaded. Subclasses
        must specify how to connect to each `Port` on `p` and vice versa.

        This is meant to be a convenience method that Elements can use to
        specify how they should connect with other elements.
        """
        pass

    def hasport(self, p):
        """Check if `p` is a child `Port`.

        :param p: `Port` or nickname for `Port` in question.
        :return: Nickname of child `Port`; or 'None' if not found.
        """
        cname = self.haschild(p)
        c = self.getchild(cname)
        if not isinstance(c, Port): cname = None
        return cname

    def addport(self, nickname, *args, **kwargs):
        """Create a new `Port`.

        :param nickname: Nickname of new `Port`.
        :param args: Additional arguments passed to `Port` constructor.
        :param kwargs: Additional keywords passed to `Port` constructor.

        This is essentially a wrapper for `newchild()` with some mangled
        parameters. All arguments and keywords are passed to the `Port`
        constructor.

        Any child `Port` may be accessed using the same approach as a normal
        child (i.e. it may be accessed as a member of the `Base` object).
        """
        if 'tracename' not in kwargs: kwargs['tracename'] = self.tracename+".P"
        port = self.newchild(nickname, Port, *args, **kwargs)
        if isinstance(nickname, str):
            portnick = re.sub("[%s]"%(string.whitespace), "", nickname)
            port.tracename = self.tracename + "." + portnick
        assert isinstance(port, Port), \
                "[ELEMENT]: addport() failed to create valid Port!"
        return port

    def delport(self, p):
        """Remove port `p` by disconnecting and deleting from child list.

        :param p: `Port` or nickname of Port to delete.
        :return: Removed `Port`.

        `disconnect()` Port and call `delchild()`.
        """
        c = self.delchild(p)
        c.disconnect()
        return c

    def getport(self, p):
        """Get `Port` p.

        :param p: `Port` or nickname of port.
        :return: Matching port or `None` if not found.
        """
        port, pname = None, self.hasport(p)
        if pname is not None:
            port = self.getchild(pname)
        return port

    def ports(self):
        """Get dictionary containing all Ports."""
        ports = {}
        for c in self.children:
            child = self.children[c]
            if isinstance(child, Port):
                ports[c] = child
        return ports
