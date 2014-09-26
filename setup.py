#!/usr/bin/env python

from distutils.core import setup
from Cython.Build import cythonize

setup(
	name = 'fields',
	ext_modules = cythonize("fields.pyx"),
)
