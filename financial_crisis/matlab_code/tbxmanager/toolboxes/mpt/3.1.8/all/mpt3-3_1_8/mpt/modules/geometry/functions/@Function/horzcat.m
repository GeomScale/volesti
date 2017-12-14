function y=horzcat(varargin)
%
%  HORZCAT: Horizontal concatenation for Function class. 
%  ======================================================
%  
%  
%  SYNTAX
%  ------
%     
%      F = horzcat(F1,F2)
%      F = [F1,F2]
%      F = [F1,F2,F3]
%    
%  
%  DESCRIPTION
%  -----------
%     The objects of Function class are concatenated into arrays only, no matrix
%  concatenation is allowed. Note that only the objects of the same class can be
%  concatenated into the same array.
%  
%  INPUT
%  -----
%     
%        
%          F1 Function object.                         
%             Class: Function                          
%          F2 Function object.                         
%             Class: Function                          
%               
%  
%  
%  OUTPUT
%  ------
%     
%        
%          F An array of Function objects.            
%            Class: Function                          
%              
%  
%  
%  SEE ALSO
%  --------
%     Function,  vertcat,  disp
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
 
 
e = cellfun(@(x)builtin('isempty',x),varargin,'UniformOutput',false);

% delete empty entries
varargin(cell2mat(e))=[];

% check if the sets are the same
s = cellfun(@class,varargin,'UniformOutput',false);
if length(unique(s))~=1
   error('Only the same objects can be concatenated.');
end

for i=1:length(varargin)
    varargin{i} = transpose(varargin{i}(:));
end

y = builtin('horzcat',varargin{:});
