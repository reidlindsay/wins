#!  /usr/bin/env python

"""
Set confidence interval for binomial proportion estimator. Works with 'y' values
in the [0,1] range.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-08-06 15:40:35 -0500 (Sat, 06 Aug 2011) $
* $LastChangedRevision: 5103 $

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

def calc_error(dataset):
    """Modify data and return modified data set."""
    param, data = copy(dataset['parameters']), copy(dataset['data'])
    if 'label' in param:
        lbl = param['label']
        isrealdata = (lbl.lower().find("actual")>-1) or (lbl.lower().find("hydra")>-1)
        if isrealdata:
            param['use-yerror'] = True
        else:
            return copy(dataset)
    for dp in data:
        x, y, ndata = dp['x'], dp['y'], dp['ndata']
        yvar = (y*(1-y)/ndata)
        dp['yvar'] = yvar
    newdata = {'parameters': param, 'data': data}
    return newdata

def main():
    usage = "%prog [OPTIONS] DATAFILE1 [DATAFILE2 ...] > MODIFIEDFILE \n" + \
            "      Writes modified data to standard output."
    parser = OptionParser(usage=usage)
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
        # calculate error values for data
        newdata = calc_error(dataset)
        sys.stdout.write("%s\n"%(newdata))

if __name__ == '__main__':
    main()
