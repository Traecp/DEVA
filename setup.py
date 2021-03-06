#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'Tra NGUYEN THANH <thanh-tra.nguyen@esrf.fr>'
__version__ = '3.2'
__adv__ = 'setup.py'


import os, sys, glob
from os.path import join, expanduser, isdir
tmp_dir = expanduser("~")
tmp_dir = join(tmp_dir, "DEVA")
if not isdir(tmp_dir):
	cmd = "mkdir %s"%tmp_dir
	os.system(cmd)
cmd1 = "cp Instrument_configuration.DEVA %s/."%tmp_dir
cmd2 = "chmod 777 %s"%join(tmp_dir, "Instrument_configuration.DEVA")
print cmd1
print cmd2
os.system(cmd1)
os.system(cmd2)
#sys.argv.append('bdist_wininst')
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

# Package name
name = 'DEVA'

# Packages (subdirectories in lib/)
packages = [name, name+'.utilities', name+'.xpad']

# Scripts (in scripts/)
scripts = ['Deva.py']

command_options = {}


setup(name=name,
	version = __version__,
	description='The DEVA software: D2AM Edf images Visualisation and Analysis',
	author='Tra NGUYEN THANH',
	author_email='thanh-tra.nguyen@esrf.fr',
	maintainer="Tra NGUYEN",
	maintainer_email='thanh-tra.nguyen@esrf.fr',
	url='http://www.esrf.fr/',
	license="CeCILL-C FREE SOFTWARE LICENSE",
	packages=packages,
	scripts=['scripts/'+script for script in scripts],
	zip_safe = False,
	command_options=command_options
      )
