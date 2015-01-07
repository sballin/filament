#!/usr/bin/env python

from distutils.core import setup
from Cython.Build import cythonize

extra_compile_args = ["-O3"]

setup(
    name = 'fields',
    ext_modules = cythonize("fields.pyx"),
)
