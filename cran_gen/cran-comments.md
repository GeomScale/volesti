## Downstream dependencies

There are currently no downstream dependencies for this package

## check results

1. We ran devtools::check() and devtools::build_win() locally.  
2. We used package `r-hub` for other platforms.  

0 errors | 0 warnings in all platforms (Build fails only in `Oracle Solaris 10, x86, 32 bit, R-patched (experimental)`)

###  There were no ERRORs or WARNINGs. 

There was 3 NOTES in total:  


- FROM: all Platforms,  

* checking compiled code ... NOTE  
File ‘volesti/libs/volesti.so’:  
  Found ‘rand’, possibly from ‘rand’ (C)  
    Object: ‘vol_R.o’  

  Library lpSolveAPI uses rand() and srand() in lp_utils.c. We replace both functions with GetRNGstate(); PutRNGstate(); unif_rand(); from R’s internal random number generation routines as it is proposed in `Writing R Extensions`. Moreover if you run in folder `/src`:  
$ grep -r 'rand()'
You just get:  
`utils.c:  range *= (LPSREAL) unif_rand();`  
which is our replacement. If you replace `rand()` with `srand` in grep search you get a null result.  
This NOTE appears because of our functions in `/src/include/samplers` where word `rand` appears a lot of times, for example `rand_point_generator()`.  

--------------------------------------------

- FROM: all Platforms,  

* checking compiled code ... NOTE  
Found the following apparent object files/libraries:  
  src/lp_solve/colamd.o src/lp_solve/commonlib.o src/lp_solve/ini.o  
  src/lp_solve/liblp_solve.a src/lp_solve/lp_Hash.o  
  src/lp_solve/lp_LUSOL.o src/lp_solve/lp_MDO.o src/lp_solve/lp_MPS.o  
  src/lp_solve/lp_SOS.o src/lp_solve/lp_crash.o src/lp_solve/lp_lib.o  
  src/lp_solve/lp_matrix.o src/lp_solve/lp_mipbb.o  
  src/lp_solve/lp_params.o src/lp_solve/lp_presolve.o  
  src/lp_solve/lp_price.o src/lp_solve/lp_pricePSE.o  
  src/lp_solve/lp_report.o src/lp_solve/lp_rlp.o  
  src/lp_solve/lp_scale.o src/lp_solve/lp_simplex.o  
  src/lp_solve/lp_utils.o src/lp_solve/lp_wlp.o src/lp_solve/lusol.o  
  src/lp_solve/mmio.o src/lp_solve/myblas.o src/lp_solve/yacc_read.o  
Object files/libraries should not be included in a source package.  

We compile the C code from the lpsolveAPI R package using the same Makefile in folder `/src/lp_solve` and the same commands in the Makevars as the authors do.

--------------------------------------------

- FROM: `devtools::build_win()`, `Fedora Linux, R-devel, clang, gfortran`, `Windows Server 2008 R2 SP1, R-devel, 32/64 bit`,  

* checking CRAN incoming feasibility ... NOTE  
Maintainer: 'Fisikopoulos Vissarion <vissarion.fisikopoulos@gmail.com>'  

New submission  

Possibly mis-spelled words in DESCRIPTION:  
  Minkowski (9:17)  
  Polytopes (4:52)  
  Volesti (7:71)  
  polytopes (8:42, 10:86)  
  volesti (7:50)  
  zonotopes (9:44)  

These words except of `volesti` (which is the name of our package) describe geometrical concepts.  

--------------------------------------------

- FROM: `Windows Server 2008 R2 SP1, R-devel, 32/64 bit`,  

* checking sizes of PDF files under 'inst/doc' ... NOTE
Unable to find GhostScript executable to run checks on size reduction

We do not use GhostScript.


###  External libraries

We use two external libraries:
- The [lpSolveAPI](https://cran.r-project.org/web/packages/lpSolveAPI/index.html) package (LGPL-2). In order to build [lpSolveAPI](https://cran.r-project.org/web/packages/lpSolveAPI/index.html) we use the same pattern: folder src/lp_solve contains the .c files and folder external/LPsolve_src/include contains the header files that used to build .o files. Folder LPsolve_src/run_headers contains the header files that are suggested by [lpSolve](http://lpsolve.sourceforge.net/5.5/Build.htm) in order to call the main functions of the package.  
- We use a part of code from [BNMin1](https://github.com/bnikolic/oof/tree/master/bnmin1) library (GPL-2).

