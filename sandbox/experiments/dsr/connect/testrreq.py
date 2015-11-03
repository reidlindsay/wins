#!  /usr/bin/env python

from wins import *
from wins.ieee80211 import *

from pylab import *
import numpy as np
import networkx as nx

from numpy.random import uniform, randint

import sys
import time
from optparse import OptionParser

def db(x):
    return 10*np.log10(x)

def invdb(x):
    return 10**(x/10.0)

#DATA_MODEL =  {'cma':{'base':  3.9, 'alt':  3.9},
#               'cfo':{'base':  3.9, 'alt':  3.9}}
#DATA_HYDRA =  {'cma':{'base':  4.2, 'alt':  4.2},
#               'cfo':{'base':  7.9, 'alt':  7.9}}
#RREQ_MODEL =  {'cma':{'base':  2.9, 'alt': 12.3},
#               'cfo':{'base':  2.9, 'alt': 12.3}}
#RREQ_HYDRA =  {'cma':{'base':  3.2, 'alt': 12.3},
#               'cfo':{'base':  7.3, 'alt': 14.4}}

DATA_MODEL =  {'cma':{'base':  3.9, 'alt':  3.9},
               'cfo':{'base':  3.9, 'alt':  3.9}}
DATA_HYDRA =  {'cma':{'base':  4.0, 'alt':  4.0},
               'cfo':{'base':  8.1, 'alt':  8.1}}
RREQ_MODEL =  {'cma':{'base':  3.0, 'alt': 12.5},
               'cfo':{'base':  3.0, 'alt': 12.5}}
RREQ_HYDRA =  {'cma':{'base':  3.2, 'alt': 13.5},
               'cfo':{'base':  7.3, 'alt': 14.0}}

DATA0_MODEL = {'cma':{'base':  2.5, 'alt':  2.5},
               'cfo':{'base':  2.5, 'alt':  2.5}}
DATA0_HYDRA = {'cma':{'base':  1.7, 'alt':  1.7},
               'cfo':{'base':  4.6, 'alt':  4.6}}
RREQ0_MODEL = {'cma':{'base':  2.0, 'alt': 10.9},
               'cfo':{'base':  2.0, 'alt': 10.9}}
RREQ0_HYDRA = {'cma':{'base':  1.3, 'alt':  9.3},
               'cfo':{'base':  4.2, 'alt': 10.1}}

SNRPENALTY = 0

def calcdist(a, b):
    return np.linalg.norm(a-b)

def findroute(g, a, b, hops):
    assert (g.has_node(a) and g.has_node(b))
    xg, rg = g.copy(), nx.Graph()
    rg.add_node(a)
    slist = [(0,a)]
    done = False
    while not done:
        # decrement counters
        for k in range(len(slist)):
            t, n = slist[k]
            slist[k] = (t-1,n)
        # exit if empty
        if not slist:
            done = True
            continue
        # check if head is expired
        t, sender = slist[0]
        if (t>0): continue
        slist.pop(0)
        assert rg.has_node(sender)
        # process neighbors
        neighbors = xg.neighbors(sender)
        for n in neighbors:
            if (n==b): done = True          # target reached
            if rg.has_node(n): continue     # node already added
            assert (n not in [w for v,w in slist])
            rg.add_edge(sender, n)
            t = randint(0, DOT11_CWMIN)
            slist.append((t,n))
        # remove sender from original graph
        xg.remove_node(sender)
        slist.sort()
    # check for route from a -> b
    path = False
    if rg.has_node(b):
        path = nx.shortest_path(rg,a,b)
        if path and (len(path)>hops+1): path = False
    return path

def addedge(dist, dmin, dmax=None, ch=None):
    dmax = None
    if (ch is None): alpha = 2.0
    else:            alpha = ch.n
    if dmax is None: return (dist<dmin)
    n = alpha
    num = 1.0*n*db(dist/dmax)
    den = 1.0*n*db(dmin/dmax)
    p = max(0, min(1,num/den) )
    return (uniform(0,1)<p)

def main():
    usage = "%prog [OPTIONS]"
    parser = OptionParser(usage=usage)
    parser.add_option("", "--alt", dest="alt", action="store_true", \
            default=False, help="Use alternate DSR.")
    parser.add_option("", "--cfo", dest="cfo", action="store_true", \
            default=False, help="Use CFO.")
    parser.add_option("", "--model", dest="model", action="store_true", \
            default=False, help="Use PHY model.")
    parser.add_option("", "--ntopo", dest="ntopo", type="int", \
            default=30, help="Specify number of topologies [default=%default].")
    parser.add_option("-f", "--fmt", dest="fmt", \
            default=None, help="Format of data parsed from tracefiles [default=%default]")
    (options, args) = parser.parse_args()

    if len(args)>0:
        print "Invalid number of arguments."
        parser.print_help()
        raise SystemExit

    # process options
    mode = 'base'
    usemodel = 0
    chmode = 'cma'
    if options.alt: mode = 'alt'
    if options.model: usemodel = 1
    if options.cfo:  chmode = 'cfo'
    ntopo = options.ntopo

    # parameters
    border = (0, 150, 0, 150)
    xmin, xmax, ymin, ymax = border[0:4]
    numnodes = [2, 5, 10, 15, 20, 25, 30]
    #numnodes = [1, 2, 5, 10, 15, 20, 25, 30, 35, 40]
    MAXRREQ = 8

    # channel parameters
    alpha = None
    modeltype = "A"
    # calculate SNR values
    Ptx = DOT11N_MAXPOWER
    Lrx, Ltx = Dot11NRadio.Lrx, Dot11NRadio.Ltx
    Grx, Gtx = Dot11NRadio.Grx, Dot11NRadio.Gtx
    No = Dot11NRadio.thermalnoise(DOT11N_BANDWIDTH) + DOT11N_NOISEFIGURE
    Psignal = Ptx + Gtx + Grx - Ltx - Lrx
    # use breakpoint distance as refdist
    ch = Dot11NChannel(modeltype=modeltype, n=alpha)
    alpha = ch.n
    refdist = ch.bpdist
    refloss = Propagation.freespace(refdist, n=2.0, fc=ch.fc)

    # get RREQ and DATA SNR
    if usemodel:
        snrdata  = DATA_MODEL[chmode][mode]
        snrrreq  = RREQ_MODEL[chmode][mode]
        snrdata0 = DATA0_MODEL[chmode][mode]
        snrrreq0 = RREQ0_MODEL[chmode][mode]
    else:
        snrdata  = DATA_HYDRA[chmode][mode]
        snrrreq  = RREQ_HYDRA[chmode][mode]
        snrdata0 = DATA0_HYDRA[chmode][mode]
        snrrreq0 = RREQ0_HYDRA[chmode][mode]

    # add fixed penalty
    snrdata  += SNRPENALTY
    snrrreq  += SNRPENALTY
    snrdata0 += SNRPENALTY
    snrrreq0 += SNRPENALTY

    # calculate critical distance for DATA and RREQ
    pldata  = Psignal - (snrdata + No)
    plrreq  = Psignal - (snrrreq + No)
    pl0data = Psignal - (snrdata0 + No)
    pl0rreq = Psignal - (snrrreq0 + No)
    if (pldata>refloss): n = alpha
    else:                n = 2.0
    ddata = refdist*(invdb((1.0*pldata-refloss)/n))
    if (plrreq>refloss): n = alpha
    else:                n = 2.0
    drreq = refdist*(invdb((1.0*plrreq-refloss)/n))
    if (pl0data>refloss): n = alpha
    else:                 n = 2.0
    d0data = refdist*(invdb((1.0*pl0data-refloss)/n))
    if (pl0rreq>refloss): n = alpha
    else:                 n = 2.0
    d0rreq = refdist*(invdb((1.0*pl0rreq-refloss)/n))

    sys.stderr.write("alpha    = %.1f\n"%(alpha))
    sys.stderr.write("usemodel = %s\n"%(usemodel))
    sys.stderr.write("snrrreq  = %s\n"%(snrrreq))
    sys.stderr.write("%s.%s.drreq = %s\n"%(chmode,mode,drreq))
    
    # keep track of results
    data = []

    spos = np.array([xmin,ymin])
    dpos = np.array([xmax,ymax])
    srdist = np.linalg.norm(spos-dpos)
    assert (srdist>drreq)
    for NN in numnodes:
        for NT in range(ntopo):
            # make new topology
            xpos = uniform(xmin,xmax, NN)
            ypos = uniform(ymin,ymax, NN)
            pos = [np.array(xpos[k],ypos[k]) for k in range(NN)]
            # check for route
            for r in range(MAXRREQ):
                # connect topology
                gdata = nx.Graph()
                grreq = nx.Graph()
                gdata.add_node('src'), gdata.add_node('dst')
                grreq.add_node('src'), grreq.add_node('dst')
                for a in range(NN):
                    apos = pos[a]
                    sdist = calcdist(spos,apos)
                    ddist = calcdist(dpos,apos)
                    if addedge(sdist,ddata,d0data,ch=ch): gdata.add_edge('src', a)
                    if addedge(ddist,ddata,d0data,ch=ch): gdata.add_edge('dst', a)
                    if addedge(sdist,drreq,d0rreq,ch=ch): grreq.add_edge('src', a)
                    if addedge(ddist,drreq,d0rreq,ch=ch): grreq.add_edge('dst', a)
                    for b in range(NN):
                        bpos = pos[b]
                        if (a==b): continue
                        abdist = calcdist(apos,bpos)
                        if addedge(abdist,ddata,d0data,ch=ch): gdata.add_edge(a,b)
                        if addedge(abdist,drreq,d0rreq,ch=ch): grreq.add_edge(a,b)
                # do route discovery
                x = findroute(grreq, 'src', 'dst', 2**r)
                pathok = False
                if x:
                    pathok = True
                    for k in range(len(x)-1):
                        if not gdata.has_edge(x[k],x[k+1]): pathok = False
                if pathok: break    # route found
                #sys.stderr.write("no path: n = %d, r = %d\n"%(NN, r))
            # check if src and dst are connected
            y = 0.0
            #if nx.shortest_path(g,'src','dst'): y = 1.0
            if pathok: y = 1.0
            dp = {'x':NN, 'y':y, 'ndata': 1}
            data.append(dp)
            # print status
            if ((NT%50)==0):
                sys.stderr.write("finished topo %3s, NN = %3s ...\n"%(NT,NN))

    # set param
    param = {'xlabel': "Number of Nodes", \
             'ylabel': "% connected", \
             'title': ""}

    fmt = options.fmt
    if fmt is None:
        fmt = ""
        if usemodel: fmt += "r"
        else:        fmt += "b"
        if (chmode=="cma"): fmt += "o"
        else:               fmt += "^"
        if usemodel: fmt += "--"
        else:        fmt += "-"
    param['format'] = fmt

    mlabel = "(%s.%s)"%(chmode.upper(), mode.upper())
    if usemodel: label = "WiNS Model %s"%(mlabel)
    else:        label = "Hydra PHY %s"%(mlabel)
    param['label'] = label

    param['use-yerror'] = True
    if usemodel: param['markerfacecolor'] = 'None'

    # output results
    pdata = {'data': data, 'parameters': param}
    sys.stdout.write("%s\n"%(pdata))

if __name__=='__main__':
    main()
