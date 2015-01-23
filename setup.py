#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'Tra NGUYEN THANH <thanh-tra.nguyen@esrf.fr>'
__version__ = '1.3.8'
__adv__ = 'setup.py'


import os, sys, glob
#sys.argv.append('bdist_wininst')
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

from distutils.core import  Extension
from Cython.Distutils import build_ext

# for numpy
from numpy.distutils.misc_util import get_numpy_include_dirs

if sys.platform in ["linux2", "posix"]:
    openmp = '-fopenmp'
    spliter = "/"
elif sys.platform in ["win32", "nt"]:
    openmp = '/openmp'
    spliter = "\\"
src = {}
ext_dir = "DEVA/pyFAI/extensions"
cython_files = [os.path.splitext(i)[0] for i in glob.glob(os.path.join(ext_dir,"*.pyx"))]
#print cython_files
if build_ext:
	for ext in cython_files:
		ext0 = ext.split(spliter)[-1]
		src[ext0] = os.path.join(".", ext + ".pyx")
		
else:
	for ext in cython_files:
		ext0 = ext.split(spliter)[-1]
		src[ext0] = os.path.join(".", ext + ".c")
#print src
hist_ext = Extension("histogram",
                    include_dirs=get_numpy_include_dirs(),
                    sources=[src['histogram']],
                    extra_compile_args=['-fopenmp'],
                    extra_link_args=['-fopenmp'])


#halfsplit_ext = Extension("halfSplitPixel",
                    #include_dirs=get_numpy_include_dirs(),
                    #sources=[src['halfSplitPixel']],
                    #extra_compile_args=['-fopenmp'],
                    #extra_link_args=['-fopenmp'])


split_ext = Extension("splitPixel",
                    include_dirs=get_numpy_include_dirs(),
                    sources=[src['splitPixel']],
                    extra_compile_args=['-fopenmp'],
                    extra_link_args=['-fopenmp'])

splitPixelFull_ext = Extension("splitPixelFull",
                    include_dirs=get_numpy_include_dirs(),
                    sources=[src['splitPixelFull']])
                    
splitBBoxCSR_ext = Extension("splitBBoxCSR",
                    include_dirs=get_numpy_include_dirs(),
                    sources=[src['splitBBoxCSR']])

splitPixelFullCSR_ext = Extension("splitPixelFullCSR",
                    include_dirs=get_numpy_include_dirs(),
                    sources=[src['splitPixelFullCSR']])

relabel_ext = Extension("relabel",
                        include_dirs=get_numpy_include_dirs(),
                        sources=[src['relabel']])

bilinear_ext = Extension("bilinear",
                        include_dirs=get_numpy_include_dirs(),
                        sources=[src['bilinear']])
#rebin_ext = Extension("fastrebin",
                        #include_dirs=get_numpy_include_dirs(),
                        #sources=[src["slist"]])

splitBBox_dic = dict(name="splitBBox",
                    include_dirs=get_numpy_include_dirs(),
                    sources=[src['splitBBox']],)

paraSplitBBox_dic = dict(name="paraSplitBBox",
                    include_dirs=get_numpy_include_dirs(),
                    sources=[src['paraSplitBBox']],
                    extra_compile_args=[openmp],
                    extra_link_args=[openmp])
splitBBoxLUT_dic = dict(name="splitBBoxLUT",
                    include_dirs=get_numpy_include_dirs(),
                    sources=[src['splitBBoxLUT']],
                    extra_compile_args=[openmp],
                    extra_link_args=[openmp]
                    )
marchingsquares_dict = dict(name="marchingsquares",
                    include_dirs=get_numpy_include_dirs(),
                    sources=[src['marchingsquares']],
                    #extra_compile_args=[openmp],
#                    extra_link_args=[openmp]
                    )

sparse_csr_dict = dict(name="sparse_csr",
                    include_dirs=get_numpy_include_dirs(),
                    sources=[src['sparse_csr']],
                    #extra_compile_args=[openmp],
#                    extra_link_args=[openmp]
                    )
_convolution_dict = dict(name="_convolution",
                    include_dirs=get_numpy_include_dirs(),
                    sources=[src['_convolution']],
                    extra_compile_args=[openmp],
                    extra_link_args=[openmp]
                    )
morphology_dict = dict(name="morphology",
                    include_dirs=get_numpy_include_dirs(),
                    sources=[src['morphology']],
#                    extra_compile_args=[openmp],
#                    extra_link_args=[openmp]
                    )
                    
ext_modules = [Extension(**splitBBox_dic),
				   Extension(**paraSplitBBox_dic),
                   Extension(**splitBBoxLUT_dic),
                   Extension(**marchingsquares_dict),
                   Extension(**sparse_csr_dict),
                   Extension(**_convolution_dict),
                   Extension(**morphology_dict),
                   hist_ext, split_ext, relabel_ext, bilinear_ext, splitPixelFull_ext, splitBBoxCSR_ext, splitPixelFullCSR_ext
                   ]
# Package name
name = 'DEVA'

# Packages (subdirectories in lib/)
packages = [name, name+'.pyFAI', name+'.utilities', name+'.xpad', name+'.xrayutilities', name+'.xrayutilities/io']

# Scripts (in scripts/)
scripts = ['Deva.py', 'xpad3_geometry.py']

cmdclass = {'build_ext': build_ext}
command_options = {}


setup(name=name,
	  version = __version__,
      description='The DEVA software: D2AM Edf images Visualisation and Analysis',
      author='Tra NGUYEN THANH',
      author_email='thanh-tra.nguyen@esrf.fr',
      maintainer="Tra NGUYEN",
      maintainer_email='thanh-tra.nguyen@esrf.fr',
      url='http://esrf.fr/',
      license="CeCILL-C FREE SOFTWARE LICENSE",
      packages=packages,
      package_data={
        "DEVA": ["xrayutilities/*.conf"]
        },
      scripts=['scripts/'+script for script in scripts],
      # Data - setuptools specific
      ext_package = name+'/pyFAI',
      ext_modules=ext_modules,
      cmdclass=cmdclass,
      zip_safe = False,
      command_options=command_options
      )
