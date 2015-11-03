#!  /usr/bin/env python

"""
Calculate theoretical throughput of pure ALOHA vs. normalized workload (G).

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-10-20 21:39:01 -0500 (Thu, 20 Oct 2011) $
* $LastChangedRevision: 5236 $

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

from optparse import OptionParser
import numpy as np
import sys

def calctput():
    usage = "%prog [OPTIONS] \n" + \
            "      Writes ALOHA throughput data to standard output."
    parser = OptionParser(usage=usage)
    parser.add_option("", "--dg", dest="dg", type="float", \
            default=0.05, help="Set step size for load (G) [default=%default]")
    parser.add_option("-G", "--max", dest="max", type="float", \
            default=2.0, help="Set max load (G) [default=%default]")
    (options, args) = parser.parse_args()

    if len(args)>0:
        print "Invalid number of arguments."
        parser.print_help()
        raise SystemExit

    # set parameters
    gmax = options.max
    dg   = options.dg

    # calculate throughput
    ndata = int(gmax/dg) + 1
    G = np.arange(ndata)*dg
    tput = G*np.exp(-2*G)

    # format data and output to stdout
    assert (len(G) == len(tput))
    assert (len(G) == ndata)
    data = []
    for k in range(ndata):
        dp = {'x': G[k], 'y': tput[k], 'ndata': 0}
        data.append(dp)
    param = {'xlabel': "G", \
             'ylabel': "${\\rm T}_{put}$", \
             'title': "Throughput vs. Workload"}
    param['format'] = 'k:'
    param['label']  = "${\\rm S}_{pure}$"
    parsed_data = {'parameters': param, 'data': data}
    sys.stdout.write("%s\n"%(parsed_data) )

if __name__ == '__main__':
    calctput()
