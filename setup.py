#!  /usr/bin/env python

"""
Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-10-08 10:39:45 -0500 (Sat, 08 Oct 2011) $
* $LastChangedRevision: 5189 $

:author: Ketan Mandke <kmandke@mail.utexas.edu>

:copyright:
    Copyright 2009-2010 The University of Texas at Austin

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

from distutils.core import setup
from extras import CheckCommand, CleanCommand, DocsCommand
from distutils.core import Extension
import os

PROJECT_ROOT = ''
if 'PROJECT_ROOT' in os.environ: PROJECT_ROOT = os.environ['PROJECT_ROOT']

#ex_module = Extension("wins._example", sources=["ext/example.i","ext/example.cc"])

crc_module = Extension('wins.backend._crc', \
        sources = ['backend/crc.i', 'backend/crc.cc'], \
        swig_opts = ['-c++', '-I`pwd`/backend', '-Wall'], \
        library_dirs = ["%s/lib"%(PROJECT_ROOT)], \
        libraries=['python2.6', 'itpp'], \
        extra_compile_args=['-Wall'], \
        include_dirs=['/usr/include/python2.6', 'backend', "%s/include"%(PROJECT_ROOT)] \
        )

itpp_module = Extension('wins.backend._itpp', \
        sources = ['backend/itpp.i'],
        swig_opts = ['-c++', '-I`pwd`/backend'],
        library_dirs = ["%s/lib"%(PROJECT_ROOT)],
        libraries=['python2.6', 'wins-80211n'],
        include_dirs=['/usr/include/python2.6', 'backend', "%s/include"%(PROJECT_ROOT)] \
        )

waveform_module = Extension('wins.backend._waveform', \
        sources = ['backend/waveform.i'], \
        swig_opts = ['-c++', '-I`pwd`/backend'], \
        library_dirs = ["%s/lib"%(PROJECT_ROOT)], \
        libraries=['python2.6', 'wins-80211n'], \
        include_dirs=['/usr/include/python2.6', 'backend', "%s/include"%(PROJECT_ROOT), "%s/include/wins-80211n"%(PROJECT_ROOT)] \
        )

channel_module = Extension('wins.backend._channel', \
        sources = ['backend/channel.i'], \
        swig_opts = ['-c++', '-I`pwd`/backend'], \
        library_dirs = ["%s/lib"%(PROJECT_ROOT)], \
        libraries=['python2.6', 'wins-80211n'], \
        include_dirs=['/usr/include/python2.6', 'backend', "%s/include"%(PROJECT_ROOT), "%s/include/wins-80211n"%(PROJECT_ROOT)] \
        )

dot11n_module = Extension('wins.backend._dot11n_backend', \
        sources = ['backend/dot11n_backend.i', \
                   'backend/dot11n_transmitter.cc', \
                   'backend/dot11n_receiver.cc', \
                   'backend/dot11n_channel.cc'], \
        swig_opts = ['-c++', '-I`pwd`/backend'], \
        library_dirs = ["%s/lib"%(PROJECT_ROOT)], \
        libraries = ['python2.6', 'itpp', 'wins-80211n'], \
        include_dirs=['/usr/include/python2.6', 'backend', "%s/include"%(PROJECT_ROOT), "%s/include/wins-80211n"%(PROJECT_ROOT)] \
        )

setup(name="wins",
      version="0.1",
      description="Wireless Network Simulator",
      author="Ketan Mandke",
      author_email="kmandke@mail.utexas.edu",
      url="http://netlab.ece.utexas.edu/mandke/wins/",
      packages=['wins', \
                'wins.mobile', \
                'wins.channel', \
                'wins.digital', \
                'wins.protocol', \
                'wins.util', \
                'wins.ieee80211', \
                'wins.mac', \
                'wins.net', \
                'wins.traffic', \
                'wins.backend'],
      ext_modules = [crc_module, \
                     itpp_module, \
                     waveform_module, \
                     dot11n_module],
      package_dir = {'wins':'wins', 'wins.backend':'backend'},
      cmdclass = { 'check':CheckCommand, 'clean':CleanCommand,
                   'docs':DocsCommand } )
