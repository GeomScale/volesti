function T = triangulate(P)
%
%  TRIANGULATE: Triangulation of a polyhedron. 
%  ============================================
%  
%  
%  SYNTAX
%  ------
%     
%      T = P.triangulate;
%      T = triangulate(P)
%    
%  
%  DESCRIPTION
%  -----------
%     The polyhedron P  is split into an array of triangular T_i  polyhedra that
%  built a partition P =  T_i . The triangulate function works over
%  full-dimensional polyhedra, that are bounded and not empty.
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
%     minVRep,  computeVRep
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
 
 
error(P.rejectArray());

if ~P.isBounded || ~P.isFullDim || P.isEmptySet
    error('Only bounded, non-empty polyhedra in the full dimension can be triangulated.');
end
P.minVRep();
V = P.V;

% Triangulate V
%K = delaunayn(V);
K = delaunayn(V, {'Qt', 'Qbb', 'Qc', 'Qz'}); 
T(size(K,1))=Polyhedron;
for i=1:size(K,1)
    T(i) = Polyhedron(V(K(i,:),:));
end

% attach functions if defined
T.copyFunctionsFrom(P);

end
