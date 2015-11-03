#!  /usr/bin/env python

"""
Merge x-y data files and output result to stdout.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-11-13 13:36:44 -0600 (Sun, 13 Nov 2011) $
* $LastChangedRevision: 5333 $

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
import sys
from copy import copy

def read_data(options, dfile):
    """Read data from input file.

    :return: list containing {parameters, dict} dictionary definitions.
    """
    f = file(dfile, 'r')
    data = []
    # read data from file
    s = f.readline()
    done = not (s)
    while not done:
        # convert s to dict
        try:
            d = eval(s)
            assert isinstance(d, dict)
            assert ('parameters' in d)
            assert ('data' in d)
        except:
            d = None
        # add dict to data
        if d: data.append(d)
        # get next input
        s = f.readline()
        done = not (s)
    return data

def merge_data(data1, data2=[], epsilon=1e-12):
    """Merge [two] list(s) of datapoints consisting of (x,y,ndata) information.

    :epsilon: Define epsilon size for each data bin.

    :return: List of merged data.
    """
    dp = []
    # check validity of data
    ok = True
    isnum = lambda x: isinstance(x,float) or isinstance(x,int)
    for c in ['x', 'y', 'ndata']:
        ok &= all([((c in d) and isnum(d[c])) for d in data1])
        ok &= all([((c in d) and isnum(d[c])) for d in data2])
    errmsg = "Invalid data point found!"
    assert ok, errmsg
    # merge data points
    newdata = []
    olddata = data1 + data2
    done = (len(olddata)==0)
    while (len(olddata)>0):
        # get next data point
        d1 = olddata.pop(0)
        keepdata = []
        for d2 in olddata:
            # extract values {x, y, ndata [,xvar, yvar]}
            x1, y1, n1 = d1['x'], d1['y'], d1['ndata']
            x2, y2, n2 = d2['x'], d2['y'], d2['ndata']
            xv1, xv2, yv1, yv2 = 0, 0, 0, 0
            if ('xvar' in d1): xv1 = d1['xvar']
            if ('yvar' in d1): yv1 = d1['yvar']
            if ('xvar' in d2): xv2 = d2['xvar']
            if ('yvar' in d2): yv2 = d2['yvar']
            # merge data point d2?
            if (abs(x1-x2)<epsilon):
                # combine assuming no ndata info
                x = (x1 + x2)/2.0
                y = (y1 + y2)/2.0
                # Var[Z] = E[Z^2] - (E[Z])^2
                xv = ((xv1 + x1**2) + (xv2 + x2**2))/2.0 - x**2
                yv = ((yv1 + y1**2) + (yv2 + y2**2))/2.0 - y**2
                # combine using # of data points
                if (n1+n2>0):
                    x = (1.0*x1*n1 + x2*n2)/(n1+n2)
                    y = (1.0*y1*n1 + y2*n2)/(n1+n2)
                    x1sqr = n1*(xv1 + x1**2)  # (sum of X1)^2
                    x2sqr = n2*(xv2 + x2**2)  # (sum of X2)^2
                    y1sqr = n1*(yv1 + y1**2)  # (sum of Y1)^2
                    y2sqr = n2*(yv2 + y2**2)  # (sum of Y2)^2
                    xv = (1.0*x1sqr + x2sqr)/(n1+n2) - x**2
                    yv = (1.0*y1sqr + y2sqr)/(n1+n2) - y**2
                # update data point d1 with merged data
                d1 = {'x':x, 'y':y, 'ndata': (n1+n2), \
                      'xvar':max(0,xv), 'yvar':max(0,yv)}
            else:
                # keep data point d2 for next pass
                keepdata.append(copy(d2))
        # add merged data point d1 to new data
        newdata.append(copy(d1))
        olddata = keepdata
    # return merged data
    return newdata

def is_param_same(param1, param2, ignore=[]):
    """Check if two parameter dictionaries have the same parameters.

    :ignore: List of parameters to ignore.
    """
    # get copy of parameters
    p1 = copy(param1)
    p2 = copy(param2)
    # remove parameters to ignore
    for p in ignore:
        if (p in p1): del p1[p]
        if (p in p2): del p2[p]
    # compare parameters
    issame = (p1 == p2)
    return issame

def merge_multi():
    usage = "%prog [OPTIONS] DATAFILE1 [DATAFILE2 ...] > MERGEDFILE \n" + \
            "      Writes merged data to standard output."
    parser = OptionParser(usage=usage)
    parser.add_option("-e", "--epsilon", dest="epsilon", \
                      type="float", default=1e-12, \
                      help="Epsilon value for comparing floats [default = %default]")
    parser.add_option("-i", "--ignore", dest="ignore", \
                      type="string", default=None, \
                      help="Specify list of parameters to ignore (e.g. \"p1 p2 p3\")")
    parser.add_option("", "--ignore-all", dest="ignoreall", \
                      action="store_true", default=False, \
                      help="Ignore all parameters and merge all data")
    (options, args) = parser.parse_args()

    if len(args)<1:
        print "Insufficient number of arguments."
        parser.print_help()
        raise SystemExit

    # process args and options
    datafile = args[0:]
    numdata  = len(datafile)
    epsilon  = options.epsilon
    ignore   = []
    if options.ignore:
        ignore = options.ignore.split()

    # combine all datasets
    alldata = []
    for dfile in datafile:
        data = read_data(options, dfile)
        alldata += data

    # merge data as needed
    while (len(alldata)>0):
        # get data set containing parameters and data points
        data1 = alldata.pop(0)
        p1, d1 = data1['parameters'], data1['data']
        keepdata = []
        for data2 in alldata:
            p2, d2 = data2['parameters'], data2['data']
            issame = is_param_same(p1, p2, ignore)
            if issame or options.ignoreall:
                # merge with data points d1 & d2
                d1 = merge_data(d1, d2, epsilon)
            else:
                # keep data set data2 for next pass
                keepdata.append(data2)
        alldata = keepdata
        # merge with data points within d1
        d1 = merge_data(d1, epsilon=epsilon)
        # print new data to stdout
        newdata = {'parameters': copy(p1), 'data': copy(d1)}
        sys.stdout.write("%s\n"%(newdata))

if __name__=='__main__':
    merge_multi()
