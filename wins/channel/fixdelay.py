#!  /usr/bin/env python

"""
Wrapper to be used to apply a fixed delay to any `ChannelModel` by modifying the
'cm-delay' annotation.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-09-06 16:44:24 -0500 (Tue, 06 Sep 2011) $
* $LastChangedRevision: 5121 $

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

from wins.channel.channelbase import ChannelModel
from wins.packet import ANNO

def _apply_fixed_delay(cm, *args, **kwargs):
    """Internal method to implement `FixedDelay` `ChannelModel` wrapper."""
    assert isinstance(cm, ChannelModel)
    assert hasattr(cm, '__subclass__')
    assert hasattr(cm, '__fixed_delay__')
    subclass = cm.__subclass__
    delay = cm.__fixed_delay__
    r = subclass.apply(cm, *args, **kwargs)
    if ANNO.supports(r, 'cm-delay'):
        r.delanno('cm-delay')
    if ANNO.supported(r) and (delay>0):
        r.setanno('cm-delay', delay)
    return r

def FixedDelay(model, delay=0):
    """Create a new `ChannelModel` with overloaded `ChannelModel.apply()` method
    to modify 'cm-delay' annotation.

    :param model: `ChannelModel` class.
    :param delay: Fixed delay value [default = 0].
    """
    assert issubclass(model, ChannelModel)
    clsname = model.__name__
    newname = "_FixedDelay_%s"%(clsname)
    tracename = "D%s"%(model.tracename)
    createclass = "class %s(model):\n\tpass"%(newname)
    exec(createclass)
    globals()[newname] = eval(newname)
    exec("%s.__subclass__ = model"%(newname) )
    exec("%s.apply = _apply_fixed_delay"%(newname) )
    exec("%s.__fixed_delay__ = delay"%(newname) )
    return eval(newname)
