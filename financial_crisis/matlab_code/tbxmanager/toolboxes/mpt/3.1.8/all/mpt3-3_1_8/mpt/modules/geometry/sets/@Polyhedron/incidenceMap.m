function iMap = incidenceMap(P)
%
%  INCIDENCEMAP: Compute the incidence map of this polyhedron. 
%  ============================================================
%  
%  
%  SYNTAX
%  ------
%     
%      iMap = P.incidenceMap
%      iMap = incidenceMap(P)
%    
%  
%  DESCRIPTION
%  -----------
%     Computes map describing containment of vertices in facets. Note: This
%  function computes both irredundant V and H-representations of the polyhedron and
%  can be time consuming.
%  
%  INPUT
%  -----
%     
%        
%          P Polyhedron in any format                 
%            Class: Polyhedron                        
%              
%  
%  
%  OUTPUT
%  ------
%     
%        
%          iMap       Incidence map or [] is P  is empty.      
%                     Class: struct                            
%          iMap.V     Vertices of P                            
%          iMap.R     Rays of P                                
%          iMap.H     Inequalities of P                        
%          iMap.He    Equalities of P                          
%          iMap.incVH Incidence of V  in H. incVH(i,j) = 1 if  
%                     [V(i,:);-1]*H(j,:)' = 0                  
%                     Class: sparse logical                    
%          iMap.incRH Incidence of R  in H. incRH(i,j) = 1 if  
%                     [R(i,:);0]*H(j,:)' = 0                   
%                     Class: sparse logical                    
%                       
%  
%  
%  SEE ALSO
%  --------
%     minVRep,  minHRep
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
 
 
global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

% deal with arrays
if numel(P)>1
	% return an array of structures
	iMap = P.forEach(@(elem) elem.incidenceMap);
    return
end

% pre-alocate output
iMap = struct('V', [], 'R', [], 'H', [], 'He', [], 'incRH', [], 'incVH', []);

% Need H-rep and V-rep
hRep = P.minHRep();
vRep = P.minVRep();

iMap.V  = vRep.V; 
iMap.R  = vRep.R;
iMap.H  = hRep.H; 
iMap.He = hRep.He;

nV = size(iMap.V,1);
nR = size(iMap.R,1);
nF = size(iMap.H,1);

% Vertex map is easy - just test for inclusion
iMap.incVH = sparse(abs(hRep.H*[vRep.V -ones(nV,1)]')' < MPTOPTIONS.abs_tol);

% A ray r is in a facet F if there exists a vertex v in F s.t. v+r in F 
iMap.incRH = spalloc(nR,nF,3*nR+3*nF);
for i=1:nF
  % Vertices in this facet
  v = vRep.V(iMap.incVH(:,i),:);
  
  % Test each ray
  for j=1:size(v,1)
    r = vRep.R + repmat(v(j,:),nR,1);
    iMap.incRH(:,i) = iMap.incRH(:,i) | abs(hRep.H(i,:)*[r -ones(nR,1)]')' < MPTOPTIONS.abs_tol;
  end
end

end
