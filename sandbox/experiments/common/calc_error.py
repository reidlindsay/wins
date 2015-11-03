#!  /usr/bin/env python

"""
Calculate absolute-error distance (L1-norm distance) and output result to stdout.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-11-22 17:23:33 -0600 (Tue, 22 Nov 2011) $
* $LastChangedRevision: 5344 $

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

import numpy as np

from pylab import *

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
            assert ('data' in d)
            assert ('parameters' in d)
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

def resample(xydata, newx):
    """Resample xydata and return newy."""
    assert (xydata) # cannot process empty set
    newy = []
    oldx, oldy = [a for (a,b) in xydata], [d for (c,d) in xydata]
    # set up parameters
    idx = 0
    xcurr, xnext = None, oldx[idx]
    ycurr, ynext = None, oldy[idx]
    # resample with new x
    for x in newx:
        yval = None
        # advance xnext until end is reached
        if (xnext is not None) and (x>xnext):
            while ( (idx<len(oldx)-1) and (x>xnext) ):
                xcurr, ycurr = xnext, ynext
                idx += 1
                xnext, ynext = oldx[idx], oldy[idx]
            if not (idx<len(oldx)):
                # passed end of data
                xnext, ynext = None, oldy[-1]
        # get new y value
        if xnext is None:
            yval = ynext
        elif (x<xnext):
            # in current bin
            if xcurr is None:
                yval = ynext  # before start of xydata -> use first value
            else:
                assert (xcurr<xnext)
                m = 1.0*(ynext-ycurr)/(xnext - xcurr)
                yval = ycurr + m*(x - xcurr)
        else:
            # assume x == xnext
            yval = ynext
        # add new value of y
        assert (yval is not None)
        newy.append(yval)
    # return resampled y
    return np.array(newy)

def main():
    usage = "%prog [OPTIONS] DATAFILE1 [DATAFILE2 ...] > MODIFIEDFILE \n" + \
            "      Writes modified data to standard output."
    parser = OptionParser(usage=usage)
    parser.add_option("", "--dx", dest="dx", type="float", \
            default=None, help="Specify minimum x tick.")
    parser.add_option("", "--no-combine", dest="nocombine", \
                      action="store_true", default=False, \
                      help="Do not sum error")
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
    numdata = len(alldata)
    assert (numdata>=2), "Can only calculate benchmark between two data sets!"

    dataset1 = alldata[0]
    dataset2 = alldata[1]
    # extract x,y data and sort data
    xy1 = [(dp['x'],dp['y']) for dp in dataset1['data']]
    xy2 = [(dp['x'],dp['y']) for dp in dataset2['data']]
    xy1.sort(), xy2.sort()
    # get labels
    label1, label2 = None, None
    if ('label' in dataset1['parameters']): label1 = dataset1['parameters']['label']
    if ('label' in dataset2['parameters']): label2 = dataset2['parameters']['label']
    # combine x-values
    xvals = [a for a,b in xy1] + [c for c,d in xy2]
    xvals.sort()
    xmin = min(xvals)
    xmax = max(xvals)
    # find minimum dx
    dx = options.dx
    if options.dx is None:
        xv0 = np.array(xvals[0:-1])
        xv1 = np.array(xvals[1:])
        dx = min(0.1, min([a for a in (xv1-xv0) if (abs(a)>0)]) )
    #newx = np.arange(int((xmax-xmin)/dx))*dx + xmin
    newx = np.linspace(xmin, xmax, int((xmax-xmin)/dx)+1)
    # resample xy1 and xy2
    newy1 = resample(xy1, newx)
    newy2 = resample(xy2, newx)
    # calculate benchmark
    abserror = abs(newy1 - newy2)
    benchmark = dx*abserror.sum()
    # write to stdout
    sys.stdout.write("label1 = %s\n"%(label1))
    sys.stdout.write("label2 = %s\n"%(label2))
    if options.nocombine:
        sys.stdout.write("----------\n")
        assert (len(newx)<=len(abserror))
        xstr = "["
        for x in newx: xstr += "%.2g, "%(x)
        xstr += "]"
        estr = "["
        for e in abserror: estr += "%.4g, "%(e)
        estr += "]"
        sys.stdout.write("x     = %s\n"%(xstr))
        sys.stdout.write("error = %s\n"%(estr))
    else:
        sys.stdout.write("error  = %.8g\n"%(benchmark))

    if 1:
        plot([a for a,b in xy1], [b for a,b in xy1], 'bo')
        hold(True), grid(True)
        plot(newx, newy1, 'bx')
        plot([a for a,b in xy2], [b for a,b in xy2], 'ro')
        plot(newx, newy2, 'rx')
        show()

if __name__ == '__main__':
    main()
