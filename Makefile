# Makefile for WiNS project.  Mostly a wrapper for python setuptools.
#
# Revision Info
# =============
# * $LastChangedBy: mandke $
# * $LastChangedDate: 2011-04-08 18:23:02 -0500 (Fri, 08 Apr 2011) $
# * $LastChangedRevision: 4948 $
# 
# :author: Ketan Mandke <mandke@ece.utexas.edu>
# 
# :copyright:
#     Copyright 2009-2010 The University of Texas at Austin
# 
#     Licensed under the Apache License, Version 2.0 (the "License");
#     you may not use this file except in compliance with the License.
#     You may obtain a copy of the License at
# 
#        http://www.apache.org/licenses/LICENSE-2.0
# 
#     Unless required by applicable law or agreed to in writing, software
#     distributed under the License is distributed on an "AS IS" BASIS,
#     WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#     See the License for the specific language governing permissions and
#     limitations under the License.

INPUT=setup
ALLPY=extras.py wins/*.py wins/*/*.py
#ALLEXT=ext/*.cc ext/*.h ext/*.i
ALLBCK=backend/*.cc backend/*.h backend/*.i backend/__init__.py
INPUT_DEPS=$(ALLPY) $(ALLEXT) $(ALLBCK) $(INPUT).py

PYTHON=`which python`
PY_SETUP=$(PYTHON) $(INPUT).py

all: $(INPUT).cfg

$(INPUT).cfg: $(INPUT_DEPS)
	@$(PY_SETUP) build_ext
	@$(PY_SETUP) build
	@touch $(INPUT).cfg

build: $(INPUT).cfg

build_%:
	@$(PY_SETUP) $@

force: force_ext
	@$(PY_SETUP) build --force
	@touch $(INPUT).cfg

force_ext:
	@$(PY_SETUP) build_ext --force

install:
	@$(PY_SETUP) install

docs: docs.log

docs.log: build
	$(PY_SETUP) docs

pdf: docs
	if [ -d docs/latex ]; then \
		cd docs/latex; \
		pdflatex refman.tex; \
		pdflatex refman.tex; \
		cd ../..; \
	fi

check: build
	@$(PY_SETUP) check

check-%: build
	@$(PY_SETUP) check --test=$(subst check-,,$@)

clean: clean-ext clean-wins

clean-wins:
	$(PY_SETUP) clean
	rm -rf build
	@touch -m -t `date +20010101%H%M` $(INPUT).cfg

clean-docs:
	rm -rf docs/html docs/latex

clean-ext:
	rm -f backend/*_wrap.cpp
	rm -f backend/[a-z]*.py
	rm -f ext/*_wrap.cpp
	rm -f ext/*_wrap.c
	rm -f ext/[a-z]*.py

clean-all: clean clean-docs
	@touch -m -t `date +20010101%H%M` $(INPUT).cfg

help:
	$(PY_SETUP) --help

help-%:
	$(PY_SETUP) $(subst help-,,$@) --help

help-commands:
	$(PY_SETUP) --help-commands

