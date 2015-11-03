#!  /usr/bin/env python

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
    param['use-yerror'] = True
    for dp in data:
        x, y, ndata = dp['x'], dp['y'], dp['ndata']
        if (x>0):
            p = y/x  # y = p*x
            pvar = (p*(1-p)/ndata)
            yvar = (x**2)*pvar
        else:
            yvar = 0
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
