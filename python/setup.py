# python setup.py develop --user

from setuptools import setup, Extension

setup(name='ke',
      ext_modules=[Extension('lib.ke', ['lib/ke.c'], extra_compile_args = ['-O2']),
                   Extension('lib.ke_newton', ['lib/ke_newton.c'], extra_compile_args = ['-O2'])]
     )

