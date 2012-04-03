#!/usr/bin/env python

"""
setup.py file for SWIG helloworld
"""

from distutils.core import setup, Extension


helloworld_module = Extension('_helloworld',
                           sources=['helloworld_wrap.cxx'],
                           libraries=["helloworld"],
                           runtime_library_dirs = ["."],
                           library_dirs=['.'],
                           extra_compile_args=['-g'],
                           extra_link_args=['-g'],
                           undef_macros=["NDEBUG"]
                           )

setup (name = 'helloworld',
       version = '0.1',
       author      = "SWIG helloworld",
       description = """SWIG helloworld""",
       ext_modules = [helloworld_module],
       py_modules = ["helloworld"],
       )
