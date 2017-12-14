function [A,D1,D2] = mpt_scale_matrix(A)
%
%  MPT_SCALE_MATRIX: Scales matrix row-wise and column-wise 
%  =========================================================
%  
%  
%  SYNTAX
%  ------
%     
%      [An,D1,D2] = mpt_scale_matrix(A)
%    
%  
%  DESCRIPTION
%  -----------
%     Scales matrix A  by finding diagonal matrices D_1  and D_2in A_n = D_1AD_2 
%  such that infinity norm of each row and column approaches 1. The problem is
%  given as 
%                                                     
%                                 min    ||D AD || )  
%                                           1  2      
%                                                     
%                                 s.t.    A  = D AD   
%                                          n    1  2  
%     Scaling matrix is used in solving linear equations of the type Ax=bfor badly
%  scaled matrix A  as follows: 
%                         Ax= b                                             (1) 
%                                                                               
%                    D Ax = D b        multiply from left by   D            (2) 
%                     1      1                                  1               
%                     -1                                 -1                     
%               D AD D  x = D b              insert   D D                   (3) 
%                1  2 2      1                         2 2                      
%                                                          -1                   
%                  D AD y = D b          substitute   y = D  x              (4) 
%                   1  2     1                             2                    
%                                                                               
%                      A y = b     substitute   A  = D AD ,  b  = D b       (5) 
%                       n     n                  n    1  2    n    1            
%                                                                           (6) 
%     First solve A_ny = b_n, then obtain x = D_2y.
%  
%  INPUT
%  -----
%     
%        
%          A Input matrix do be scaled. The matrix    
%            can be also rectangular.                 
%            Class: double                            
%              
%  
%  
%  OUTPUT
%  ------
%     
%        
%          An Scaled matrix A  such that A_n=D_1AD_2.  
%             Class: double                            
%          D1 Diagonal matrix D_1  such that           
%             A_n=D_1AD_2.                             
%             Class: double                            
%          D2 Diagonal matrix D_2  such that           
%             A_n=D_1AD_2.                             
%             Class: double                            
%               
%  
%  
%  SEE ALSO
%  --------
%     mptopt,  mpt_solve
%  

%  AUTHOR(s)
%  ---------
%     
%    
%   (c) 2010-2013  Martin Herceg: ETH Zurich
%   mailto:herceg@control.ee.ethz.ch 
%  
%  

%  LICENSE
%  -------
%    
%    This program is free software; you can redistribute it and/or modify it under
%  the terms of the GNU General Public License as published by the Free Software
%  Foundation; either version 2.1 of the License, or (at your option) any later
%  version.
%    This program is distributed in the hope that it will be useful, but WITHOUT
%  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
%  FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%    You should have received a copy of the GNU General Public License along with
%  this library; if not, write to the  Free Software Foundation, Inc.,  59 Temple
%  Place, Suite 330,  Boston, MA 02111-1307 USA
%  -------------------------------------------------------------------------------
%    
%      This document was translated from LaTeX by HeVeA (1).
%  ---------------------------------------
%    
%    
%   (1) http://hevea.inria.fr/index.html
 
 
global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

% check if A has any row or column equal zero
if ~all(any(A,1)) || ~all(any(A,2))
   error('Given matrix contains a row or column with all elements equal zero.');
end

validate_realmatrix(A);

% get dimensions
[m,n] = size(A);

% initialize 
D1 = eye(m);
D2 = eye(n);
row_norm = 1;
column_norm = 1;
k=0;

% tolerance
tol = MPTOPTIONS.rel_tol;

while (row_norm>tol) && (column_norm>tol)
    nA = max(abs(A), [], 2);
    r = sqrt(nA);
    rn = 1 - nA;
    DR = diag(1./r);

    nA = max(abs(A), [], 1)';
    c = sqrt(nA);
    cn = 1-nA;
    DC = diag(1./c);
    
    % matrix updates
    A = DR * A * DC;
    D1 = D1*DR;
    D2 = D2*DC;
    
    row_norm = norm(rn, Inf);
    column_norm = norm(cn, Inf);
    
    k = k+1;
end
