ecos-matlab
===========

Matlab interface for ECOS, the Embedded Conic Solver. For detailed documentation please visit [this page](https://www.embotech.com/ECOS).

Repository Content
====

You will find the following directories in this repository:

* `bin`: directory with the script for making a MEX binary (makemex.m) and other helper files for ecos-matlab.
* `conelp`: Matlab implementation of ECOS with different linear system solver options.
* `src`: mex interface C file
* `test`: testing code - run batchtest.m to run different tests

Using ECOS in MATLAB
====

ECOS can be used in three different ways within MATLAB:

* [native MEX interface](https://www.embotech.com/ECOS/Matlab-Interface/Matlab-Native)
* [through CVX](https://www.embotech.com/ECOS/Matlab-Interface/CVX)
* [through YALMIP](https://www.embotech.com/ECOS/Matlab-Interface/Yalmip)

In either case, a compiled mex binary is required, which can be downloaded from [embotech](https://embotech.com/ECOS/Download). The bottom part of the latter page also explains how to build it from source.
