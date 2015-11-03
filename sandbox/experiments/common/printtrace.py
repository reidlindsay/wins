#!  /usr/bin/env python

"""
Read in and print formatted output of a WiNS trace.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-06-17 17:21:51 -0500 (Fri, 17 Jun 2011) $
* $LastChangedRevision: 5005 $

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

from wins import *
from wins.ieee80211 import *

from optparse import OptionParser
import sys

def read_trace(options, tracefile):
    # load trace from file
    tr = Trace()
    tr.read(tracefile)
    # return trace
    return tr

def output_trace():
    usage = "%prog TRACEFILE"
    parser = OptionParser(usage=usage)
    (options, args) = parser.parse_args()

    if len(args)<>1:
        print "Incorrect number of arguments."
        parser.print_help()
        raise SystemExit

    tracefile = args[0]
    # read in trace from file
    tr = read_trace(options, tracefile)
    # output formatted trace
    tr.output()

if __name__ == '__main__':
    output_trace()
