import os
from os.path import split, join
from glob import glob

from itertools import chain

from setuptools import setup
from Cython.Build import cythonize
from setuptools.extension import Extension
import numpy


version = "0.1"
license='Apache 2.0',
packages = ["volestipy"]
description="volestipy: wrapper for the VolEsti library to sample from convex sets and compute volume."
author = "Pedro Zuidberg Dos Martires"
author_email="pedro.zuidbergdosmartires@cs.kuleuven.be"
name = 'volestipy'
scripts = []
url = ""
download_url = ""
classifiers = []
zip_safe=False


source_directory_list = ['volestipy', join('volestipy','src')]


compiler_args = [
    "-std=c++11",
    "-O3",
    "-lm",
    "-ldl",
    "-DBOOST_NO_AUTO_PTR",
]


extra_volesti_include_dirs = [
    join("volestipy","include"),

    join("..","external"),
    join("..","external","minimum_ellipsoid"),
    join("..","external","LPsolve_src","run_headers"),
    join("..","external","boost"),

    join("..","include","generators"),
    join("..","include","volume"),
    join("..","include"),
    join("..","include","convex_bodies"),
    join("..","include","annealing"),
    join("..","include","samplers"),
]


src_files = ["volestipy/volestipy.pyx","volestipy/src/bindings.cpp"]
extra_include_dirs = [numpy.get_include()]


library_includes = ["lpsolve55"]
ext_module = Extension(
    "volestipy",
    language = "c++",
    sources = src_files,
    include_dirs = extra_include_dirs + extra_volesti_include_dirs,
    libraries = library_includes,
    extra_compile_args = compiler_args,
    # extra_link_args = link_args,
)
ext_modules = cythonize([ext_module], gdb_debug=False)

setup(
    version = version,
    author = author,
    author_email = author_email,
    name = name,
    packages = packages,
    ext_modules = ext_modules,
    zip_safe=zip_safe,
)
