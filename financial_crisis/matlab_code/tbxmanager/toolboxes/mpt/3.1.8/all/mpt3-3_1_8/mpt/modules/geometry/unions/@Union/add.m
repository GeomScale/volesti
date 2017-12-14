function obj = add(obj, C)
%
%  ADD: Insert set to Union object. 
%  =================================
%  
%  
%  SYNTAX
%  ------
%     
%      U = add(U,Set)
%      U.add(Set)
%    
%  
%  DESCRIPTION
%  -----------
%     Insert the Set derived from from the ConvexSet inside the Union object. The
%  Set can be also an array. Each element of the array is stored under
%  Union.Setproperty as a cell array because objects derived from the ConvexSet
%  class cannot be concatenated. If the Set is empty, it is not added to the union.
%  
%  
%  INPUT
%  -----
%     
%        
%          U   The object of the Union class.           
%              Class: Union                             
%          Set A convex set derived from ConvexSet      
%              class.                                   
%              Class: ConvexSet                         
%                
%  
%  
%  OUTPUT
%  ------
%     
%        
%          U Union of the sets.                       
%            Class: Union                             
%              
%  
%  
%  SEE ALSO
%  --------
%     Union,  remove
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
 
 
if numel(obj)>1
    for i=1:numel(obj)
        obj(i) = obj(i).add(C);
    end
    return;
end

if ~isa(C,'ConvexSet')
    error('The argument must be derived from "ConvexSet" class.');
end

% remove empty sets
c = isEmptySet(C);
C(c) = [];

% check function hadles
fnames = obj.listFunctions;
nf = numel(fnames);
N = length(C);

% if union contains some function handles, the new added sets must also
% have the same function handles
if numel(fnames)>0 && N>0
	for i = 1:N
		if any(~C(i).hasFunction(fnames))
            error('The set %i to be added holds different function names than the union.',i);
		elseif length(C(i).Functions) ~= nf
			error('All sets to be added must have associated %d number of functions.',nf);
		end
	end
end

if N>0
	obj.Set(obj.Num+1:obj.Num+N,1) = num2cell(C(:));
	%obj.Set(obj.Num+1:obj.Num+N,1) = C(:);
end
   
end
