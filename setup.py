#!/usr/bin/env python

from distutils.core import setup
from Cython.Build import cythonize

setup(
	name = 'flux',
	ext_modules = cythonize("flux.pyx"),
)
