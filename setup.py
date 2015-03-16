#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'Tra NGUYEN THANH <thanh-tra.nguyen@esrf.fr>'
__version__ = '2.0.4'
__adv__ = 'setup.py'


import os, sys, glob
#sys.argv.append('bdist_wininst')
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

# Package name
name = 'DEVA'

# Packages (subdirectories in lib/)
#packages = [name, name+'.pyFAI', name+'.utilities', name+'.xpad', name+'.xrayutilities', name+'.xrayutilities/io']
packages = [name, name+'.utilities', name+'.xpad']

# Scripts (in scripts/)
scripts = ['Deva.py', 'xpad3_geometry.py']

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
