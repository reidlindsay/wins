#!  /usr/bin/env python

"""
Miscellaneous helpful functions and classes.

Revision Info
=============
* $LastChangedBy: mandke $
* $LastChangedDate: 2011-02-08 13:52:17 -0600 (Tue, 08 Feb 2011) $
* $LastChangedRevision: 4908 $

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

from distutils.core import Command
from unittest import TextTestRunner, TestLoader
from glob import glob
from os.path import splitext, basename, join as pjoin, walk
import os
import sys


def _check_import(module):
    """Check if importing module works correctly; useful for checking
    dependencies."""
    try:
        exec("import %s"%(module) )
    except Exception, e:
        sys.stderr.write("[WINS]: Error importing %s!\n"%(module))
        sys.stderr.write("      : %s\n"%(e) )
        raise SystemExit


class CheckCommand(Command):
    user_options = [("test=", 't', "list of additional tests")]
    description = "run all UnitTest modules in \'tests\' subdirectory"

    def initialize_options(self):
        self._dir = os.getcwd()
        blist = glob(pjoin(self._dir, 'build', 'lib*'))
        self._build = blist[0]
        sys.path.insert(0, self._build) # look in build dir first for wins
        self.test = None

    def finalize_options(self):
        pass

    def run(self):
        global _check_import
        _check_import("wins")
        '''
        Finds all the test modules in tests/, and runs them.
        '''
        verbosity = 1
        # extract test names from 'test' option
        if (self.test is None) or len(self.test)<1: testname="*"
        else: testname = self.test
        # run specified tests
        for t in glob(pjoin(self._dir, 'tests', "qa_%s.py"%(testname))):
            testfile = '.'.join(['tests', splitext(basename(t))[0]])
            test = TestLoader().loadTestsFromNames([testfile])
            t = TextTestRunner(verbosity = verbosity)
            sys.stderr.write("\nRunning tests for %s "%str(testfile))
            r = t.run(test)
            if not r.wasSuccessful():
                sys.stderr.write("\n[WINS]: ERROR OCCURRED! Stopping testing ...\n")
                break
 

class CleanCommand(Command):
    user_options = [ ]
    description = "clean up output of 'build' command and *.pyc"

    def initialize_options(self):
        self._clean_me = [ ]
        for root, dirs, files in os.walk('.'):
            for f in files:
                if f.endswith('.pyc'):
                    self._clean_me.append(pjoin(root, f))

    def finalize_options(self):
        pass

    def run(self):
        for clean_me in self._clean_me:
            try:
                os.unlink(clean_me)
            except:
                pass
 

class DocsCommand(Command):
    user_options = [ ]
    description = "generate documentation"

    def initialize_options(self):
        self._dir = os.getcwd()
        blist = glob(pjoin(self._dir, 'build', 'lib*'))
        self._build = blist[0]
        sys.path.insert(0, self._build) # look in build dir first for wins

    def finalize_options(self):
        pass

    def run(self):
        os.chdir('docs')
        global _check_import
        _check_import("wins")
        try:
            os.system("epydoc --config winsdocs.cfg")
            #os.system("epydoc --config winsdocs.cfg > winsdocs.log 2>&1")
        except:
            pass
