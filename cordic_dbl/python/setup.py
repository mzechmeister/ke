# python setup.py develop --user

from setuptools import setup, Extension

setup(name='ke_dbl_i',
      py_modules = ['ke_dbl_i'],
      ext_modules=[Extension('lib.ke_dbl_i', ['lib/ke_dbl_i.c'], extra_compile_args=['-Wno-parentheses'])]
     )

setup(name='ke_dbl_f',
      py_modules = ['ke_dbl_f'],
      ext_modules=[Extension('lib.ke_dbl_f', ['lib/ke_dbl_f.c'])]
     )
