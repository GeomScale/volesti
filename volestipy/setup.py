from setuptools import setup
from setuptools import Extension
from Cython.Build import cythonize
from os.path import join
import numpy


version = "0.2.0"
license='LGPL3',
packages = ["volestipy"]
description="volestipy: wrapper for the VolEsti library to sample from convex sets and compute volume."
author = "Haris Zafeiropoulos"
author_email="haris-zaf@hcmr.gr"
name = 'volestipy'
# scripts = []
# url = ""
# download_url = ""
# classifiers = []
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
# inside the volestipy directory, there is an extra volestipy folder containing a "include" folder; we add those as included dirs; the bindings.h file is located there
	join("volestipy","include"),

# the volesti code uses some external classes. these are located on the "external" directory and we need to add them as well
	join("..","external")
	join("..","external","Eigen")
	join("..","external","boost")
	join("..","external","boost","random")
	join("..","external","")


# finally, we also move back and include and add the directories on the "include" directory (generatorors, random_walks, sampling etc)
	join("..","include"),
	join("..","include", "misc"),
	join("..","include", "random_walks"),
	join("..","include", "volume"),
	join("..","include","generators"),	
]






