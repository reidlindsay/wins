#!  /usr/bin/env python

"""
Helper methods to read/write data from/to files as dictionaries.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-03-24 11:29:29 -0500 (Thu, 24 Mar 2011) $
* $LastChangedRevision: 4942 $

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

def readparam(infile, importstring="", overwrite=True):
    """Read parameters from a file.

    :param infile: File object to read from.
    :param importstring: Specify a string containg import statements that should
                         be executed in order for the parameter file `infile` to
                         be interpreted correctly.
    :param overwrite: Boolean; if true, all dictionary objects in the parameter
                      file `infile` will be merged together into one dictionary;
                      otherwise all dictionaries in `infile` will be appended
                      together to form a list.

    The parmeter file `infile` **MUST** contain only dictionaries of objects.

    :note: If `infile` is the name of the parameter file, this method will
           attempt to open the specified file with read-only permissions.
    """
    exec(importstring)
    try:    f = file(infile, 'r')
    except: f = infile
    assert isinstance(f, file), "[readfile]: No file object found!"
    if overwrite: opts = {}
    else:         opts = []
    s = f.readline()
    while s:
        if (len(s)>1) and (s[-2]=='\\'):
            s = s[:-2] + f.readline()
            continue
        try: d = eval(s)
        except: d = None
        if (d is None):
            try:
                # check if line is a comment
                exec(s)
                s = f.readline()
            except:
                # check if more needs to be read
                u = f.readline()
                if u: s = s[:-1] + u
                else: break
            continue
        elif isinstance(d, dict):
            if overwrite: opts.update(d)
            else:         opts.append(d)
        s = f.readline()
    return opts
