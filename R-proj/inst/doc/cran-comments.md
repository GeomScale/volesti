## R CMD check results

###  There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* hecking if this is a source package ... NOTE
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

  This comes from lpSolveAPI which is build exactly the same way as in lpSolveAPI R package

* checking installed package size ... NOTE
  installed size is 28.2Mb
  sub-directories of 1Mb or more:
    extdata   1.9Mb
    libs     26.0Mb

  This is because of the externals and extdata


* checking compiled code ... NOTE
File ‘volesti/libs/volesti.so’:
  Found ‘rand’, possibly from ‘rand’ (C)
    Object: ‘vol_R.o’

  This comes from external source in folder /src/external/minimum_ellipsoid which is part of [BNMin1](https://github.com/bnikolic/oof/tree/master/bnmin1) library


###  External libraries

We use two external libraries:
- The [lpSolveAPI](https://cran.r-project.org/web/packages/lpSolveAPI/index.html) package (LGPL-2). In order to build [lpSolveAPI](https://cran.r-project.org/web/packages/lpSolveAPI/index.html) we use the same pattern: folder src/lp_solve containes the .c files and folder external/LPsolve_src/include containes the header files that used to build .o files. Folder LPsolve_src/run_headers containes the header files that are suggested by [lpSolve](http://lpsolve.sourceforge.net/5.5/Build.htm) in order to call the main functions of the package.  
- We use a part of code from [BNMin1](https://github.com/bnikolic/oof/tree/master/bnmin1) library (GPL-2).

