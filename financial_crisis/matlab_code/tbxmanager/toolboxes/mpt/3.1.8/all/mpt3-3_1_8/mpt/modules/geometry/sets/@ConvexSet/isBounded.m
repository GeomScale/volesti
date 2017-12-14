function ts = isBounded(obj)
%
%  ISBOUNDED: Test if a convex set is bounded. 
%  ============================================
%  
%  
%  SYNTAX
%  ------
%     
%      ts = S.isBounded
%      ts = isBounded(S)
%    
%  
%  DESCRIPTION
%  -----------
%     Return true if the convex set S  is bounded and false otherwise. Empty set is
%  considered as bounded.
%  
%  INPUT
%  -----
%     
%        
%          S Any set derived from ConvexSet class,    
%            e.g. YSet or Polyhedron.                 
%            Class: ConvexSet                         
%              
%  
%  
%  OUTPUT
%  ------
%     
%        
%          ts True if the convex set P  is bounded and 
%             false otherwise.                         
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
%     isEmptySet
%  

%  AUTHOR(s)
%  ---------
%     
%    
%   (c) 2010-2013  Colin Neil Jones: EPF Lausanne
%   mailto:colin.jones@epfl.ch 
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
 
 
no = numel(obj);
if no>1
    ts = true(size(obj));
    for i=1:no
        ts(i) = obj(i).isBounded;
    end
    return
end

% empty set -> bounded
if obj.isEmptySet,
    ts = true;
    return;
end

% Compute the support in all +- elementary directions.
% Bounded iff all are bounded
ts = true;
I = eye(obj.Dim);
for i=1:obj.Dim
    if isinf(obj.support( I(i,:)'))
        ts = false;
        return;
    end
    if isinf(obj.support(-I(i,:)')),
        ts = false;
        return;
    end
end

end
