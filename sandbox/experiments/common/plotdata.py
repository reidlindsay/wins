#!  /usr/bin/env python

"""
Plot data from a formatted data file.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2012-02-23 16:08:43 -0600 (Thu, 23 Feb 2012) $
* $LastChangedRevision: 5419 $

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

from pylab import *
from matplotlib import pyplot

from numpy import sqrt

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
    sys.stderr.write("Reading from %s ...\n"%(infile.name))
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

from scipy.special import erfinv

def calc_zvalue(c):
    # calculate z-value for 1-a/2 confidence interval
    a = 1 - c/100.0     # error level
    z = -sqrt(2)*erfinv(1-2*(1-a/2) )
    return z

def plot_data(data, param, options):
    # set parameters
    use_xerror, use_yerror = False, False
    if 'use-xerror' in param:
        use_xerror = param['use-xerror']
    if 'use-yerror' in param:
        use_yerror = param['use-yerror']
    # create x/y arrays
    xdata, ydata = [], []
    xerrdata, yerrdata = [], []
    z = calc_zvalue(options.confidence_level)
    errfunc = lambda x: z*sqrt(x)       # x = variance
    if options.usestddev: errfunc = sqrt
    if options.usevar:    errfunc = lambda x: x
    # extract x/y data
    for d in data:
        x, y, ndata = d['x'], d['y'], d['ndata']
        xdev, ydev = 0, 0
        if options.plot_xerror and ('xvar' in d): xdev = errfunc(d['xvar'])
        if options.plot_yerror and ('yvar' in d): ydev = errfunc(d['yvar'])
        xerr, yerr = xdev, ydev
        # add new data point
        if (ndata is not None):
            assert (y is not None)
            assert (x is not None)
            xdata.append(x)
            ydata.append(y)
            xerrdata.append(xerr)
            yerrdata.append(yerr)
    assert (len(xdata)==len(ydata))
    # order data
    z = [(xdata[k], ydata[k], xerrdata[k], yerrdata[k]) for k in range(len(xdata))]
    z.sort()
    xdata, ydata = [t[0] for t in z], [t[1] for t in z]
    xerrdata, yerrdata = [t[2] for t in z], [t[3] for t in z]
    # set parameters
    kwargs = {}
    if 'format' in param:
        kwargs['fmt'] = param['format']
    label = None
    if 'label' in param:
        label = param['label']
    if 'markerfacecolor' in param:
        kwargs['markerfacecolor']=param['markerfacecolor']
        if 'format' in param:
            kwargs['markeredgecolor']=param['format'][0]
    if 'linewidth' in param:
        kwargs['linewidth']=param['linewidth']
    kwargs['capsize'] = 3.0*options.fscale
    # edit arguments for bar plot
    if options.barplot:
        if 'fmt' in kwargs:
            fmt = kwargs['fmt']
            kwargs['facecolor'] = fmt[0]
            del kwargs['fmt']
        if 'markerfacecolor' in kwargs:
            del kwargs['markerfacecolor']
    # plot data
    pline = None
    if (xdata):
        if options.plot_xerror and use_xerror: kwargs['xerr'] = xerrdata
        if options.plot_yerror and use_yerror: kwargs['yerr'] = yerrdata
        if options.barplot:
            p = pyplot.bar(xdata, ydata, **kwargs)
        else:
            p = pyplot.errorbar(xdata, ydata, **kwargs)
        if options.semilogx: pyplot.semilogx()
        if options.semilogy: pyplot.semilogy()
        pline = p[0]
        if ('estyle' in param) and (options.plot_yerror or options.plot_xerror):
            p[2][0].set_linestyle(param['estyle'])
    # return line handle and label for plot
    return (pline, label)


def plot_multi():
    usage = "%prog [OPTIONS] DATAFILE1 [DATAFILE2 ...]"
    parser = OptionParser(usage=usage)
    parser.add_option("", "--ignore-title", dest="ignore_title", \
                      action="store_true", default=False, \
                      help="Ignore title when plotting data (use first title from first data set).")
    parser.add_option("", "--plot-xerror", dest="plot_xerror", \
                      action="store_true", default=False, \
                      help="Use x-error bars when plotting data.")
    parser.add_option("", "--plot-yerror", dest="plot_yerror", \
                      action="store_true", default=False, \
                      help="Use y-error bars when plotting data.")
    parser.add_option("", "--use-stddev", dest="usestddev", \
                      action="store_true", default=False, \
                      help="Use standard deviation (instead of variance) when plotting error bars.")
    parser.add_option("", "--use-var", dest="usevar", \
                      action="store_true", default=False, \
                      help="Use variance for plotting error bars.")
    parser.add_option("", "--save-fig", dest="savefig", \
                      help="Save figure to file.")
    parser.add_option("-s", "--fscale", dest="fscale", type="float", \
            default=1.0, help="Scale figure text/markers/lines when saving figure [default=%default]")
    parser.add_option("-w", "--wscale", dest="wscale", type="float", \
            default=1.0, help="Scale figure width [default=%default]")
    parser.add_option("-t", "--hscale", dest="hscale", type="float", \
            default=1.0, help="Scale figure height [default=%default]")
    parser.add_option("", "--auto-scale", dest="autoscale", \
                      action="store_true", default=False, \
                      help="Automatically determine width/height scaling.")
    parser.add_option("-c", "", dest="confidence_level", type="int", \
            default=95, help="Confidence level (percentage) [default=%default]")
    parser.add_option("-l", "--loc", dest="loc", type="int", \
            default=None, help="Specify location for legend.")
    parser.add_option("", "--xmin", dest="xmin", type="float", \
            default=None, help="Set minimum value on x-axis [default=%default]")
    parser.add_option("", "--xmax", dest="xmax", type="float", \
            default=None, help="Set maximum value on x-axis [default=%default]")
    parser.add_option("", "--ymin", dest="ymin", type="float", \
            default=None, help="Set minimum value on y-axis [default=%default]")
    parser.add_option("", "--ymax", dest="ymax", type="float", \
            default=None, help="Set maximum value on y-axis [default=%default]")
    parser.add_option("", "--semilogx", dest="semilogx", \
                      action="store_true", default=False, \
                      help="Plot x-axis using a log scale")
    parser.add_option("", "--semilogy", dest="semilogy", \
                      action="store_true", default=False, \
                      help="Plot y-axis using a log scale")
    parser.add_option("", "--barplot", dest="barplot", \
                      action="store_true", default=False, \
                      help="Make bar plot")
    (options, args) = parser.parse_args()

    if len(args)<1:
        print "Insufficient number of arguments."
        parser.print_help()
        raise SystemExit

    # process args and options
    datafile = args[0:]

    # combine all datasets
    datain = []
    for dfile in datafile:
        data = read_data(options, dfile)
        datain += data
    numdata = len(datain)

    # set parameters for saving figure to file
    if options.savefig:
        if (options.savefig.lower().find('svg')>-1):
            rcParams.update({'backend': 'svg'})
        elif (options.savefig.lower().find('pdf')>-1):
            rcParams.update({'backend': 'pdf'})
        else:
            print "Cannot create non-pdf or non-svg file %s"%(options.savefig)
            raise SystemExit
        fscale = options.fscale
        rcParams.update({'font.size': 12.0*fscale, \
                         'legend.fontsize': 'medium', \
                         'lines.linewidth': 1.0*fscale, \
                         'lines.markeredgewidth': 0.5*fscale, \
                         'lines.markersize': 6*fscale*1.00, \
                         'axes.linewidth': 1.0*fscale, \
                         'patch.linewidth': 1.0*fscale, \
                         'grid.linewidth': 0.5*fscale, \
                         'xtick.major.size' : 4*fscale, \
                         'xtick.minor.size' : 2*fscale, \
                         'xtick.major.pad'  : 4*fscale, \
                         'xtick.minor.pad'  : 4*fscale, \
                         'ytick.major.size' : 4*fscale, \
                         'ytick.minor.size' : 2*fscale, \
                         'ytick.major.pad'  : 4*fscale, \
                         'ytick.minor.pad'  : 4*fscale, \
                         })
    # plot data
    lgd = {'pline':[], 'label':[]}
    xname, yname, tname = None, None, None
    hold(True), grid(True)      # set up plot
    for k in range(numdata):
        # plot data
        p = datain[k]
        data, param = p['data'], p['parameters']
        pline, label = plot_data(data, param, options)
        # assign legend data
        lgd['pline'].append(pline)
        lgd['label'].append(label)
        # get title, xlabel, ylabel parameters
        if xname: assert (xname == param['xlabel'])
        else:     xname = param['xlabel']
        if yname: assert (yname == param['ylabel'])
        else:     yname = param['ylabel']
        if tname: assert ((tname == param['title']) or options.ignore_title)
        else:     tname = param['title']
    # set other plot parameters
    if (numdata>0):
        # check legend data
        assert (len(lgd['pline'])==len(lgd['label']))
        keep = {'pline':[], 'label':[]}
        for k in range(len(lgd['pline'])):
            pline = lgd['pline'][k]
            label = lgd['label'][k]
            if (pline is not None) and (label is not None):
                keep['pline'].append(pline)
                keep['label'].append(label)
        if tname: title(tname)
        if xname: xlabel(xname)
        if yname: ylabel(yname)
        lwargs = {}
        if options.loc is not None: lwargs['loc'] = options.loc
        if keep['pline']: legend(keep['pline'], keep['label'], **lwargs)
    # fix dimensions (if needed)
    fig = gcf()
    if options.autoscale:
        options.hscale = 1.0
        options.wscale = 1.0
        if (numdata>2): options.wscale = 1.67
    if (options.wscale>0):
        w = fig.get_figwidth()
        fig.set_figwidth(w*options.wscale)
    if (options.hscale>0):
        h = fig.get_figheight()
        fig.set_figheight(h*options.hscale)

    # adjust axis if xmin/xmax/ymin/ymax is set
    xmin, xmax = options.xmin, options.xmax
    ymin, ymax = options.ymin, options.ymax
    z = list(axis())
    if xmin is not None: z[0] = xmin
    if xmax is not None: z[1] = xmax
    if ymin is not None: z[2] = ymin
    if ymax is not None: z[3] = ymax
    axis(z)

    # save to figure and show plot
    if options.savefig:
        savefig(options.savefig,transparent=True)
    # show plot
    show()

if __name__ == '__main__':
    plot_multi()
