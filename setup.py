#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
BEM-2D
A 2D boundary element method code

"""
from distutils.core import setup

setup(name='BEM2D',
      version='1.0.0',
      description='A 2D boundary element method code',
      author='Lehigh Biofluids Group',
      author_email='kmoored@lehigh.edu',
      url='http://www.lehigh.edu/',
      packages=['bem2d', 'data_IO_class', 'FSI_bem2d', 'FSIClass', 'functions_general', 'functions_graphics', 'functions_influence', 'parameter_classes', 'PyFEA', 'rigid_bem2d', 'SolidClass', 'swimmer_class', 'swimmer_subclasses', 'terminal_output'],
     )