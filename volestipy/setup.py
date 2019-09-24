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
lpsolve_base = join('external','lp_solve_5.5')
# lpsolve_base = join('..','external','LPSolve_src')



compiler_args = [
    "-std=c++11",
    "-O3",
    "-DYY_NEVER_INTERACTIVE",
    "-DPARSER_LP",
    "-DINVERSE_ACTIVE=INVERSE_LUSOL",
    "-DRoleIsExternalInvEngine",
    "-lm",
    "-ldl",
    "-DBOOST_NO_AUTO_PTR",
]
link_args = ['-O3']



extra_include_dirs = [join(lpsolve_base, d) for d in ['.', 'shared', 'bfp', join('bfp','bfp_LUSOL'), join('bfp','bfp_LUSOL/LUSOL'), 'colamd']]

extra_library_dirs = []
extra_volesti_include_dirs = [
    join("volestipy","include"),

    join("..","external","Eigen"),
    join("..","external"),
    join("..","external","minimum_ellipsoid"),
    join("..","external","boost"),

    join("..","include","generators"),
    join("..","include","volume"),
    join("..","include"),
    join("..","include","convex_bodies"),
    join("..","include","annealing"),
    join("..","include","samplers"),
]

library_includes = []



print(extra_volesti_include_dirs)

src_files = ["volestipy/volestipy.pyx","volestipy/src/bindings.cpp"]
extra_sources = {'volestipy' : [join(lpsolve_base, f) for f in
                                          ['lp_MDO.c', join('shared','commonlib.c'), join('shared','mmio.c'), join('shared','myblas.c'),
                                           'ini.c', 'fortify.c', 'colamd/colamd.c', 'lp_rlp.c', 'lp_crash.c',
                                           join('bfp','bfp_LUSOL','lp_LUSOL.c'), join('bfp','bfp_LUSOL','LUSOL','lusol.c'), 'lp_Hash.c',
                                           'lp_lib.c', 'lp_wlp.c', 'lp_matrix.c', 'lp_mipbb.c', 'lp_MPS.c', 'lp_params.c',
                                           'lp_presolve.c', 'lp_price.c', 'lp_pricePSE.c', 'lp_report.c', 'lp_scale.c',
                                           'lp_simplex.c', 'lp_SOS.c', 'lp_utils.c', 'yacc_read.c']]}
extra_include_dirs.append(numpy.get_include())


ld_library_path = os.getenv("LD_LIBRARY_PATH")
if ld_library_path is not None:
    lib_paths = ld_library_path.split(":")
    lib_paths = [l for l in lib_paths if len(l)!=0]
else:
    lib_paths = []
include_path = os.getenv("INCLUDE_PATH")
if include_path is not None:
    include_paths = [p.strip() for p in include_path.split(":") if len(p.strip()) > 0]
else:
    include_paths = []


def get_include_dirs():
    return extra_include_dirs + extra_volesti_include_dirs + include_paths

def get_library_dirs():
    return extra_library_dirs + lib_paths

def get_libraries():
    return library_includes

def get_source_files(m):
    return src_files + extra_sources[m]


modname = "volestipy"
ext_module = Extension(
    modname,
    language = "c++",
    sources = get_source_files(modname),
    include_dirs = get_include_dirs(),
    library_dirs = get_library_dirs(),
    libraries = get_libraries(),
    extra_compile_args = compiler_args,
    extra_link_args = link_args,
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
