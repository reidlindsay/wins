#!  /usr/bin/env python

"""
Combine PDR/FDR data to create PDR/FDR vs. Nc; output to stdout.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-10-28 17:34:53 -0500 (Fri, 28 Oct 2011) $
* $LastChangedRevision: 5317 $

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
from optparse import OptionParser
from copy import copy

from numpy import *

def read_data(options, fin=sys.stdin):
    """Read data from file (or standard in by default).

    :return: list containing {parameters, dict} dictionary definitions.
    """
    datain = []
    infile = fin
    if isinstance(fin, str): infile = open(fin, 'r')
    # read from data file
    s = infile.readline()
    done = not (s)
    sys.stderr.write("Reading from standard in ...\n")
    while not done:
        # covert s to dict
        try:
            d = eval(s)
            assert isinstance(d, dict)
        except:
            d = None
        # add dict to datain
        if d: datain.append(d)
        # get next input
        s = infile.readline()
        done = not(s)
    # close file (if needed)
    if isinstance(fin, str): infile.close()
    return datain

def hasfield(d, *args):
    ok = all([(a in d) for a in args])
    return ok

USE_NC_MAX = 1

def combine_data(dataset, options):
    """Modify data and return modified data set."""
    param, data = dataset['parameters'].copy(), dataset['data']
    hasfield(param, 'ncollision', 'label')
    ncoll = param['ncollision']
    label = param['label']
    if not isinstance(ncoll, int):
        # ncoll is 'all' or 'max'
        isall = (ncoll.lower()=='all')
        ismax = (ncoll.lower()=='max')
        assert (isall or ismax)
        if isall: return {}
        idx = label.find(">")
        if (idx<0): return {}
        ncoll = int(label[idx+1:].replace(")", "")) + 1
        # or return nothing
        if not USE_NC_MAX: return {}
    # combined data points
    z, nz, zv = 0, 0, 0
    for dp in data:
        # get y-data values (ignore x-data)
        y = dp['y']
        ny, yv = 0, 0
        if 'yvar' in dp:  yv = dp['yvar']
        if 'ndata' in dp: ny = dp['ndata']
        # combine data
        c = (y+z)/2.0
        cv = ((yv + y**2) + (zv + z**2))/2.0 - z**2
        if (nz+ny>0):
            c = (1.0*y*ny + z*nz)/(ny+nz)
            ysqr = ny*(yv + y**2)   # (sum of Y)^2
            zsqr = nz*(zv + z**2)   # (sum of Z)^2
            cv = (1.0*ysqr + zsqr)/(ny+nz) - c**2
        z, nz, zv = c, (ny+nz), cv
    # create combined data point
    cp = {'x': ncoll, 'y':z, 'ndata': nz, 'xvar': 0, 'yvar': zv}
    # modify parameters
    del param['ncollision']
    ishydra = (label.find("Hydra")>-1) and (label.find("PHY")>-1)
    ismodel = (label.find("WiNS")>-1)  and (label.find("Model")>-1)
    if   ishydra: param['label'] = "Hydra PHY"
    elif ismodel: param['label'] = "WiNS Model"
    if options.title:
        param['title'] = options.title
    param['xlabel'] = "$N_{coll}$"
    # return data set
    cdata = {'parameters': param, 'data': [cp]}
    return cdata

def main():
    usage = "%prog [OPTIONS] DATAFILE1 [DATAFILE2 ...] > MODIFIEDFILE \n" + \
            "      Writes modified data to standard output."
    parser = OptionParser(usage=usage)
    parser.add_option("", "--title", dest="title", \
                      default=None, help="Set title for combined data set.")
    (options, args) = parser.parse_args()

    if len(args)<1:
        print "Insufficient number of arguments."
        parser.print_help()
        raise SystemExit

    # process args and options
    datafile = args[0:]
    numdata  = len(datafile)

    # combine all datasets
    alldata = []
    for dfile in datafile:
        data = read_data(options, dfile)
        alldata += data

    # merge data as needed
    while (len(alldata)>0):
        # get data set containing parameters and data points
        dataset = alldata.pop(0)
        # merge data points
        cdata = combine_data(dataset, options)
        if cdata:
            sys.stdout.write("%s\n"%(cdata))

if __name__ == '__main__':
    main()
