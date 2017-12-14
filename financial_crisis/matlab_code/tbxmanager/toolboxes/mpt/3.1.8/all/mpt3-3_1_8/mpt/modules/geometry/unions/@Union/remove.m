function obj = remove(obj, index)
%
%  REMOVE: Remove set from Union object. 
%  ======================================
%  
%  
%  SYNTAX
%  ------
%     
%      U = remove(U,index)
%      U.remove(index)
%    
%  
%  DESCRIPTION
%  -----------
%     Removes the Set from the union based on the index. The elements of the Union
%  object are stored under U.Set property as a cell array. The index must
%  correspond to elements in this array to remove.
%  
%  INPUT
%  -----
%     
%        
%          U     The object of the Union class.           
%                Class: Union                             
%          index An index set derived that specifies      
%                which sets to remove from the union.     
%                Class: double                            
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
%     Union,  add
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
 
 
validate_indexset(index);

% deal with arrays
if numel(obj)>1
    for i=1:numel(obj)
        obj(i) = obj(i).remove(index);
    end
    return;
end

if any(index>obj.Num)
    error('Index is greater than the number of elements in the Union.');
end

% remove based on the index
obj.Set(index) = [];
   
end
