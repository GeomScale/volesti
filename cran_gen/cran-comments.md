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

  Library lpSolveAPI uses rand() and srand() in lp_utils.c. We replace both functions with GetRNGstate(); PutRNGstate(); double unif_rand(); from R’s internal random number generation routines as it is proposed in `Writing R Extensions`. Moreover if you run in folder `/src`:  
$ grep -r 'rand()'
You just get:  
`utils.c:  range *= (LPSREAL) unif_rand();`  
which is our replacement. If you replace `rand()` with `srand` in grep search you get a null result.  
This NOTE appears because of our functions in `/src/include/samplers` where word `rand` appears a lot of times, for example `rand_point_generator`.  

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
- The [lpSolveAPI](https://cran.r-project.org/web/packages/lpSolveAPI/index.html) package (LGPL-2). In order to build [lpSolveAPI](https://cran.r-project.org/web/packages/lpSolveAPI/index.html) we use the same pattern: folder src/lp_solve containes the .c files and folder external/LPsolve_src/include containes the header files that used to build .o files. Folder LPsolve_src/run_headers containes the header files that are suggested by [lpSolve](http://lpsolve.sourceforge.net/5.5/Build.htm) in order to call the main functions of the package.  
- We use a part of code from [BNMin1](https://github.com/bnikolic/oof/tree/master/bnmin1) library (GPL-2).

