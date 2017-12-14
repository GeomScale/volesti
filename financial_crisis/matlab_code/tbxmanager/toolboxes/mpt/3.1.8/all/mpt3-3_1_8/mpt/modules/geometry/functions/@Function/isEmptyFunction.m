function ts = isEmptyFunction(obj)
%
%  ISEMPTY: Test if the object of the Function class contains a function handle. 
%  ==============================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      isempty(F)
%      F.isEmptyFunction
%    
%  
%  DESCRIPTION
%  -----------
%     The objects of Function class is considered as empty when it does not have a
%  function handle associated to it, i.e. the property Handle is empty.
%  
%  INPUT
%  -----
%     
%        
%          F Function object.                         
%            Class: Function                          
%              
%  
%  
%  SEE ALSO
%  --------
%     Function,  setHandle
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
 
 
ts = zeros(size(obj));
for i=1:numel(obj)
    ts(i) = isempty(obj(i).Handle);
end
