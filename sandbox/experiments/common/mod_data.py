#!  /usr/bin/env python

"""
Modify data and parameters in a data file.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-10-17 20:06:45 -0500 (Mon, 17 Oct 2011) $
* $LastChangedRevision: 5212 $

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

def mod_data(dataset, options):
    """Modify data and return modified data set."""
    param, data = dataset['parameters'].copy(), dataset['data']
    for dp in data:
        x, y = dp['x'], dp['y']
        ndata, xvar, yvar = None, None, None
        if 'xvar' in dp:  xvar = dp['xvar']
        if 'yvar' in dp:  yvar = dp['yvar']
        if 'ndata' in dp: ndata = dp['ndata']
        # evaluate new values as requested
        if options.setndata: dp['ndata'] = eval(options.setndata)
        if options.setx:     dp['x'] = eval(options.setx)
        if options.sety:     dp['y'] = eval(options.sety)
        if options.setxvar:  dp['xvar'] = eval(options.setxvar)
        if options.setyvar:  dp['yvar'] = eval(options.setyvar)
    if options.setparam: param[options.setparam] = options.pvalue
    newdata = {'parameters': param, 'data': data}
    return newdata

def main():
    usage = "%prog [OPTIONS] DATAFILE1 [DATAFILE2 ...] > MODIFIEDFILE \n" + \
            "      Writes modified data to standard output."
    parser = OptionParser(usage=usage)
    parser.add_option("", "--set-ndata", dest="setndata", \
                      default=None, help="Set ndata-value in data.")
    parser.add_option("", "--set-x", dest="setx", \
                      default=None, help="Set x-value in data.")
    parser.add_option("", "--set-y", dest="sety", \
                      default=None, help="Set y-value in data.")
    parser.add_option("", "--set-xvar", dest="setxvar", \
                      default=None, help="Set xvar-value in data.")
    parser.add_option("", "--set-yvar", dest="setyvar", \
                      default=None, help="Set yvar-value in data.")
    parser.add_option("", "--set-param", dest="setparam", \
                      default=None, help="Set parameter in data set to pvalue.")
    parser.add_option("", "--pvalue", dest="pvalue", \
                      default=None, help="Parameter value used with --set-param.")
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
        newdata = mod_data(dataset, options)
        sys.stdout.write("%s\n"%(newdata))

if __name__ == '__main__':
    main()
