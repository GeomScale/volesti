function ts = isEmptySet(obj)
%
%  ISEMPTYSET: Test if a convex set is empty. 
%  ===========================================
%  
%  
%  SYNTAX
%  ------
%     
%      ts = S.isEmptySet
%      ts = isEmptySet(S)
%    
%  
%  DESCRIPTION
%  -----------
%     Return true if the convex set S  is empty and false otherwise.
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
%          ts True if the convex set P  is empty and   
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
%     isBounded
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
 
 
global MPTOPTIONS

if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

no = numel(obj);
if no>1
    ts = false(size(obj));
    for i=1:no
        ts(i) = obj(i).isEmptySet;
    end
    return
elseif no<1
    ts = true;
    return
end

% Try to compute support - empty if this is infeasible
ret = obj.extreme(ones(obj.Dim,1));
ts = false;
if ret.exitflag == MPTOPTIONS.INFEASIBLE
    ts = true; 
end
end
