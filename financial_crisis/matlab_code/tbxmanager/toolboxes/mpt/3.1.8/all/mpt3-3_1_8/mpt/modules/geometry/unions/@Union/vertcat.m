function y=vertcat(varargin)
%
%  VERTCAT: Vertical concatenation for union objecs. 
%  ==================================================
%  
%  
%  SYNTAX
%  ------
%     
%      U = [U1; U2]
%      S = vertcat(U1,U2)
%    
%  
%  DESCRIPTION
%  -----------
%     Overloaded method for vertical concatenation of unions. It is not possible to
%  concatenate objects of different type to the same array (e.g. PolyUnion and
%  Union). Similarly, it is not possible to concatenate into matrices, only vectors
%  are allowed.
%  
%  INPUT
%  -----
%     
%        
%          U1 Object of the Union class.               
%             Class: Union                             
%          U2 Object of the Union class.               
%             Class: Union                             
%               
%  
%  
%  OUTPUT
%  ------
%     
%        
%          U The array of the Union objects.          
%            Class: Union                             
%              
%  
%  
%  SEE ALSO
%  --------
%     horzcat
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
 
 
e = cellfun('isempty', varargin);
varargin(e)=[];

% check whether all arguments are of the same type
first_class = class(varargin{1});
for i = 2:numel(varargin)
	if ~isequal(class(varargin{i}), first_class)
		error('Only the objects of the same type can be concatenated.');
	end
end

% make sure that each argument is a column array
for i = 1:length(varargin)
	if size(varargin{i}, 2) > 1
		varargin{i} = varargin{i}(:);
	end
end

% concatenate arguments vertically
y = builtin('vertcat',varargin{:});

end
