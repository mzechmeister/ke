# python setup.py develop --user

from setuptools import setup, Extension

setup(name='ke_dbl_i',
      ext_modules=[Extension('lib.ke_dbl_i', ['lib/ke_dbl_i.c'])]
     )
