Linear Complementarity Problem solver which implements lexicographic Lemke's algorithm. 


INSTALLATION:
To install the package within Matlab, just add the path to the main folder "lcp". 
In case that you want to recompile the source code, go to the "lcp" directory and type 
"lcp_compile" at the Matlab prompt.

FEATURES:
 - The algorithm avoids cycling and degeneracy problems via lexicographic 
   perturbation technique.
 - Solid linear algebra - based on BLAS and LAPACK 
                        - recursive basis update based on LUMOD package
 - Robust features - automatic data scaling
                   - repetitive basis refactorization
 - Tailored to Matlab & Simulink specifications and Real-Time Workshop.


USAGE:
For a demonstration how to use the code, type "lcp_test" for a simple example.


LICENSE:
The software is covered by General Public License (GPL) and needs linking with external 
libraries BLAS and LAPACK (also covered with GPL).
	  http://www.gnu.org/licenses/gpl.html

LUMOD files are covered by BSD License for SOL Numerical Software	
          http://www.opensource.org/licenses/bsd-license.php

OPENWATCOM files are covered by Open Watcom Public License
	   http://www.openwatcom.org/index.php/Open_Watcom_Public_License


AUTHORS:
 Copyright (C) 2006 by Colin N. Jones, Automatic Control Laboratory, ETH Zurich,
 cjones@control.ee.ethz.ch, colin.jones@epfl.ch
 
 Revised 2010-2013 by Martin Herceg, Automatic Control Laboratory, ETH Zurich,  
 herceg@control.ee.ethz.ch  
