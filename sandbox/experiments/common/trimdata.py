#!  /usr/bin/env python

"""
Trim data based on program options.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-10-21 20:11:44 -0500 (Fri, 21 Oct 2011) $
* $LastChangedRevision: 5252 $

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

def trimdata(dataset, options):
    """Trim data points and return modified data set."""
    param, data = dataset['parameters'].copy(), dataset['data']
    newdata = []
    for dp in data:
        x, y = dp['x'], dp['y']
        ndata, xvar, yvar = None, None, None
        if 'xvar' in dp:  xvar = dp['xvar']
        if 'yvar' in dp:  yvar = dp['yvar']
        if 'ndata' in dp: ndata = dp['ndata']
        # evaluate if data point should be trimmed
        trim = False
        if options.trim:
            trim = eval(options.trim)
        # keep/discard data points
        if not trim:
            newdata.append(dp)
    pdata = {'parameters': param, 'data': newdata}
    return pdata

def main():
    usage = "%prog [OPTIONS] DATAFILE1 [DATAFILE2 ...] > MODIFIEDFILE \n" + \
            "      Writes trimmed data sets to standard output."
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--trim", dest="trim", \
                      default=None, help="Specify condition used to trim out data points.")
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
        # modify data set
        newdata = trimdata(dataset, options)
        sys.stdout.write("%s\n"%(newdata))

if __name__ == '__main__':
    main()
