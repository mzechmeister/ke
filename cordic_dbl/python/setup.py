# python setup.py develop --user

from setuptools import setup, Extension

# with options -pie and -lm the shared library is also executable

setup(name='ke_dbl_i',
      py_modules=['ke_dbl_i'],
      ext_modules=[Extension('lib.ke_dbl_i', ['lib/ke_dbl_i.c'],
                                extra_compile_args=['-Wno-parentheses', '-pie'],
                                extra_link_args=['-lm', '-pie', '-Wl,-E'])]
     )

setup(name='ke_dbl_f',
      py_modules=['ke_dbl_f'],
      ext_modules=[Extension('lib.ke_dbl_f', ['lib/ke_dbl_f.c'],
                                extra_compile_args=['-pie'],
                                extra_link_args=['-lm', '-pie', '-Wl,-E'])]
     )
