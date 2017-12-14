function [z, w, basis, exitflag, pivots, time] = lcp(M, q, options)
%
%  LCP: Matlab interface to lexicographic Lemke algorithm 
%  =======================================================
%  
%  
%  SYNTAX
%  ------
%     
%      [z, w, basis, exitflag, pivots, time] = lcp(M, q, options)
%    
%  
%  DESCRIPTION
%  -----------
%     The lexicographic Lemke algorithm solves linear-complementarity problem (LCP) of the form 
%                                            find   w, z            (1) 
%                                            s.t.   w - Mz = q      (2) 
%                                                   w >= 0          (3) 
%                                                   z >= 0          (4) 
%                                                    T                  
%                                                   w z = 0         (5) 
%                                                                       
%    where w, z  are the unknown variables, and M, q  are given. For a feasible solution, matrix M 
%  should be positive-semidefinite. The complentarity condition w'z=0 ensures that either one, or both
%  elements of vectors w, z  are zero.
%  
%  INPUT
%  -----
%     
%    
%      M                       positive semi-definite square patrix (for a        
%                              feasible solution)                                 
%                              Class: double                                      
%      q                       right hand side vector with dimension equal M      
%                              Class: double                                      
%      options                 options structure                                  
%                              Class: struct                                      
%      options.zerotol         Less than this treshold the value is considered as 
%                              zero.                                              
%                              Class: double                                      
%                              Default: 1e-10                                     
%      options.lextol          Lexicographic tolerance - a small treshold from    
%                              which values are considered as equal.              
%                              Class: double                                      
%                              Default: 1e-10                                     
%      options.maxpiv          Maximum number of pivots to be performed.          
%                              Class: double                                      
%                              Default: 1e6                                       
%      options.nstepf          If routine is 0, then every 50 pivot steps the     
%                              basis is refactorized to avoid numerical problems  
%                              for LUMOD. For other routines the factorization is 
%                              performed at each step.                            
%                              Class: double                                      
%                              Default: 50                                        
%      options.clock           Show the information about the computational time. 
%                                                                           
%                              Class: double                                      
%                              Allowed values:                                    
%                                                                           
%                                 0                                                
%                                 1                                                
%                                                                           
%                              Default: 0                                         
%      options.verbose         Verbose output. Show progress of pivoting          
%                              algorithm including entering, leaving variables,   
%                              actual basis and basis solution.                   
%                              Class: double                                      
%                              Allowed values:                                    
%                                                                           
%                                0                                                
%                                1                                                
%                                                                           
%                              Default: 0                                         
%      options.routine         Routine which should be used to obtain a basis     
%                              solution.                                          
%                              Class: double                                      
%                              Allowed values:                                    
%                                                                           
%                                 0  Corresponds to LUmod package that performs    
%                                  factorization in the form LA = U. Depending on  
%                                  the change in A  factors L, U  are updated.     
%                                  This is the fastest method.                     
%                                 1  Corresponds to DGESV simple driver from       
%                                  LAPACK package which solves the system AX = B   
%                                  by factorizing A and overwriting B with the     
%                                  solution X. Since the factorization is          
%                                  performed at each pivot step, this method tends 
%                                  to be much slower than method 0.                
%                                 2  Corresponds to DGELS simple driver which      
%                                  solves overdetermined or underdetermined real   
%                                  linear systems min ||b - Ax||_2 involving an    
%                                  M-by- N  matrix A, or its transpose, using a QR 
%                                  or LQ factorization of A. Since the             
%                                  factorization is performed at each pivot step,  
%                                  this method tends to be much slower than method 0.                                              
%                                                                           
%                              Default: 0                                         
%      options.timelimit       Time limit in seconds. If this limit is exceeded,  
%                              the pivoting algorithm is terminated and current   
%                              basis is returned.                                 
%                              Class: double                                      
%                              Default: 3600                                      
%      options.normalize       Input matrices M, q  get scaled by D_1  and        
%                              D_2when invoking this option: M_n = D_1MD_2, q_n = 
%                              D_1q, and solution is recovered as z = D_2z_n, w = 
%                              Mz+q.                                              
%                              Class: double                                      
%                              Allowed values:                                    
%                                                                           
%                                0                                                
%                                1                                                
% 
%                              Default: 1                               
%      options.normalizethres  Scaling of matrices M and q is done only if
%                              normalize option is on and if 1-norm of
%                              matrix (maximum absolute column sum)
%                              M is greater than this value. This
%                              option ensures that the normalization is
%                              done only for badly scaled problems with
%                              large outliers and large norm. 
%                          
%  
%  
%  OUTPUT
%  ------
%     
%    
%      z        Solution vector to the linear complementarity      
%               problem                                            
%               Class: double                                      
%      w        Complementary vector to z.                         
%               Class: double                                      
%      basis    Index set describing the feasible basis.           
%               Class: double                                      
%      exitflag Information on how the algorithm terminated.       
%               Class: double                                      
%               Allowed values:                                    
%                                                                  
%                 1  feasible solution                             
%                 -1  infeasible                                   
%                 -2  unbounded                                    
%                 -3  preterminated (due to time limit or maximum  
%                  pivot limit )                                   
%                 -4  other (numerical) error                      
%                                                                  
%      pivots   Total number of pivots performed by the algorithm. 
%                                                                  
%               Class: double                                      
%      time     Time needed for algorithm to terminate (in         
%               seconds).                                          
%               Class: double                                      
%                 
%  
%  
%  EXAMPLE(s)
%  ----------
%  
%  
%  Example 1
%  =========
%    % A simple example
%      M = [1 -1 -1 -1; -1 1 -1 -1; 1 1 2 0; 1 1 0 2]
%      q = [3, 5, -9, -5]
%      [z,w,basis,exitflag,pivots] = lcp([1 -1 -1 -1;-1 1 -1 -1; 1 1 2 0;1 1 0 2],[3;5;-9;-5])
%    
%  
%  AUTHOR(s)
%  ---------
%     
%    
%   (c) 2006  Colin Neil Jones: ETH Zurich
%   mailto:cjones@control.ee.ethz.ch 
%     
%    
%   (c) 2010  Martin Herceg: ETH Zurich
%   mailto:herceg@control.ee.ethz.ch 
%  
%  
%  LICENSE
%  -------
%    
%    This program is free software; you can redistribute it and/or modify it under the terms of the GNU
%  General Public License as published by the Free Software Foundation; either version 2.1 of the
%  License, or (at your option) any later version.
%    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
%  even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
%  General Public License for more details.
%    You should have received a copy of the GNU General Public License along with this library; if not,
%  write to the  Free Software Foundation, Inc.,  59 Temple Place, Suite 330,  Boston, MA 02111-1307
%  USA
 
 
 
 
 
 
end
