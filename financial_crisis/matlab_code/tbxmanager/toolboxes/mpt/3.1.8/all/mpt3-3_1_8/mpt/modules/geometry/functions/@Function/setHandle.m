function obj=setHandle(obj, h)
%
%  SETHANDLE: Assign function handle to existing Function object 
%  ==============================================================
%  
%  
%  SYNTAX
%  ------
%     
%      setHandle(F,@fun)
%      F.setHandle(@fun)
%    
%  
%  DESCRIPTION
%  -----------
%     Overwrites function handle of the Function object with a new one. This method
%  is suitable for specification of functions where the parameters are stored under
%  Data property. Look at examples for better explanation.
%  
%  INPUT
%  -----
%     
%        
%          F Existing Function object which function  
%            we want to overwrite.                    
%            Class: Function                          
%          F Representation of the new function given 
%            as function_handle class.                
%            Class: function_handle                   
%              
%  
%  
%  OUTPUT
%  ------
%     
%        
%          F Modified Function objects.               
%            Class: Function                          
%              
%  
%  
%  SEE ALSO
%  --------
%     AffFunction,  QuadFunction,  setHandle
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
 
 
narginchk(2, 2);

if numel(h)~=numel(obj)
    error('There must be %d handles present in a cell.',numel(obj));
end

if numel(obj)>1 && ~isa(h,'cell')
    error('Handles must be given in a cell.');
end

if numel(h)==1
    h = {h};
end

% make handles a column vector
h = h(:);

for i=1:numel(obj)
    
    if ~isa(h{i},'function_handle')
        error('Argument must a function handle.')
    end

    
    if strcmp('AffFunction',class(obj(i))) || strcmp('QuadFunction',class(obj(i)))
        error('setHandle: Can not change handle for affine or quadratic functions.');
    end
    
    obj(i).Handle = h{i};
    
end
end
