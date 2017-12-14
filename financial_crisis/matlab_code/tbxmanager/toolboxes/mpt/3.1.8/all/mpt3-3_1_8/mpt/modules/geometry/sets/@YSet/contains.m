function ts = contains(obj,x)
%
%  CONTAINS: Test if the point is contained inside convex set YSet. 
%  =================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      ts = S.contains(x)
%      ts = contains(S, x)
%    
%  
%  DESCRIPTION
%  -----------
%     Returns true if x in S  and false otherwise.
%  
%  INPUT
%  -----
%     
%        
%          S A convex set described as YSet object.   
%            Class: YSet                              
%          x A point given as vector. Note that for   
%            YSet with symmetric matrix variable, the 
%            point x must be given as vector with     
%            symmetric terms.                         
%            Class: double                            
%              
%  
%  
%  OUTPUT
%  ------
%     
%        
%          ts True if the point x is contained inside  
%             YSet                                     
%             Class: logical                           
%             Allowed values:                          
%                                                      
%               true                                   
%               false                                  
%                                                      
%               
%  
%  
%  SEE ALSO
%  --------
%     YSet
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

narginchk(2, 2);
validate_realmatrix(x);

if numel(obj)==0
	ts = [];
	return
end

ny = numel(obj);
nx = size(x, 2);
ts = false(ny, nx);
if obj(1).Dim ~= size(x, 1)
	error('The point must have %d rows.', obj(1).Dim);
end

for i = 1:ny
	for j = 1:nx
		% check residuals
		assign(obj(i).vars, x(:, j));
		residual = checkset(obj(i).constraints);
		ts(i, j) = all(residual>-MPTOPTIONS.abs_tol);
	end
end

end
