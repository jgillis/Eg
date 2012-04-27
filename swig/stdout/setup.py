#!/usr/bin/env python

"""
setup.py file for SWIG helloworld
"""

from distutils.core import setup, Extension


greetcpp_module = Extension('_greetcpp',
                           sources=['greetcpp_wrap.cxx'],
                           libraries=["greetcpp"],
                           runtime_library_dirs = ["."],
                           library_dirs=['.']
                           )

setup (name = 'greetcpp',
       version = '0.1',
       author      = "SWIG greetcpp",
       description = """SWIG greetcpp""",
       ext_modules = [greetcpp_module],
       py_modules = ["greetcpp"],
       )
