# python setup.py develop --user

from setuptools import setup, Extension

setup(name='ke',
      ext_modules=[Extension('lib.ke', ['lib/ke.c']),
                   Extension('lib.ke_newton', ['lib/ke_newton.c'])]
     )

