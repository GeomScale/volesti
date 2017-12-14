function Q = getFacet(P, varargin)
%
%  GETFACET: Extract facet of the polyhedron specified by the inequality index. 
%  =============================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      Q = P.getFacet()
%      Q = P.getFacet(index)
%    
%  
%  DESCRIPTION
%  -----------
%     Extract the facet of the polyhedron P  specified by the inequality index. The
%  returned polyhedron Q  is formed by rewriting this ineqality as equality
%  constraint. The polyhedron P  is given as 
%              {        T                        T                       }
%          P = {x      a x<= b ,  i=1,..., m,   a   x= b   , j=1,..., m  }
%              {        i     i                  e,j    e,j            e }
%     Given the index, the polyhedron Q  is built as 
%      {        T                               T                       T          
%                                                }
%  Q = {x      a x<= b , forall i neq index,   a x = b , i in index,   a   x = b  
%                                  , j=1,..., m , }
%      {        i     i                         i     i                 e,j     e,j
%                                             e  }
%     The polyhedron P  must be given in its minimal representation (irredundant)
%  H-representation, otherwise an error is thrown.
%  
%  INPUT
%  -----
%     
%        
%          P     Polyhedron in H-representation           
%                Class: Polyhedron                        
%          index Index of a facet from polyhedron P       
%                which is less or equal than the number   
%                of hyperplanes defining P. If omited,    
%                all facets will be returned as a         
%                polyhedron array.                        
%                Class: double                            
%                  
%  
%  
%  OUTPUT
%  ------
%     
%        
%          Q Polyhedron Q  that represents            
%            lower-dimensional facet of the           
%            Polyhedron P.                            
%            Class: Polyhedron                        
%              
%  
%  
%  SEE ALSO
%  --------
%     isAdjacent
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
 
 
narginchk(1, 2);

% empty array
if numel(P)==0
    Q = Polyhedron;    
    return
elseif numel(P)>1
	error('This function does not support arrays of polyhedra. Use forEach().');
end

if ~P.irredundantHRep
    error('Polyhedron must be in its minimal representation. Use "minHRep()" to perform the redundancy elimination.');
end

n1 = size(P.H,1);
if nargin==2
	index = varargin{1}(:)';
	validate_indexset(index);
	if numel(index)>n1 || any(index>n1)
		error('Facet index set contains indices out of range.');
	end
else
	% no index provided => return all facets
	index = 1:n1;
end

Q = [];
for i = index
	% take complement constraints and form H-rep
	ic = setdiff(1:n1, i);

	% form new polyhedron
	Q = [Q; Polyhedron('H',P.H(ic,:),'He',[P.H(i,:); P.He])];
end

end
