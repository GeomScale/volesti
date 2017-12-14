function ts = isFullDim(obj)
%
%  ISFULLDIM: Test if the union is built from full-dimensional polyhedra. 
%  =======================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      ts = U.isFullDim
%      ts = isFullDim(U)
%    
%  
%  DESCRIPTION
%  -----------
%     Return true if all polyhedra in the union U are full-dimensional and false
%  otherwise. Once this method has been called, the information about the
%  full-dimensionality can be retrieved from U.Internal.FullDim property.
%  
%  INPUT
%  -----
%     
%        
%          U Union of polyhedra in the same           
%            dimension.                               
%            Class: PolyUnion                         
%              
%  
%  
%  OUTPUT
%  ------
%     
%        
%          ts True if all polyhedra in the union are   
%             full-dimensional and false otherwise.    
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
%     isConvex,  isOverlapping,  isConnected,  isBounded
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
    ts = -ones(size(obj));
    for i=1:numel(obj)
        ts(i) = obj(i).isFullDim;
    end
    return;
end

% empty obj
if obj.Num==0
    ts = false;
    return;
end

% check full-dimensionality
if isempty(obj.Internal.FullDim)            
    ts = all(obj.Set.isFullDim);
    obj.Internal.FullDim = ts;
else
   ts = obj.Internal.FullDim; 
end



end
