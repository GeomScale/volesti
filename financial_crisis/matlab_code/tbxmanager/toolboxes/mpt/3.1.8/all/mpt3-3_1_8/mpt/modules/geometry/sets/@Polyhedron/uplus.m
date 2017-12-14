function Q = uplus(P)
%
%  UPLUS: Unitary plus for a polyhedron. 
%  ======================================
%  
%  
%  SYNTAX
%  ------
%     
%      Q = +P;
%    
%  
%  DESCRIPTION
%  -----------
%     Creates a new copy Q of the polyhedron P by copying all fields.
%  
%  INPUT
%  -----
%     
%        
%          P Polyhedron object.                       
%            Class: Polyhedron                        
%              
%  
%  
%  SEE ALSO
%  --------
%     uminus
%  

%  AUTHOR(s)
%  ---------
%     
%    
%   (c) 2010-2013  Colin Neil Jones: EPF Lausanne
%   mailto:colin.jones@epfl.ch 
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
 
 
Q(size(P)) = Polyhedron;

% copy polyhedron
for i=1:length(P)
    Q(i) = Polyhedron(P(i));
end



end
