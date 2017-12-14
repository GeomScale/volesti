function aff = affineHull(obj)
%
%  AFFINEHULL: Computes affine hull of a convex set. 
%  ==================================================
%  
%  
%  SYNTAX
%  ------
%     
%      He = affineHull(Set)
%      He = Set.affineHull
%    
%  
%  DESCRIPTION
%  -----------
%     Compute an implicitly-defined affine hull of the convex Set. The output is a
%  real matrix He that defines the affine set 
%                                       (    )             
%                           { x      He ( x  )= 0 }     (1)
%                                       ( 1  )             
%     If He is empty, then the affine hull is empty. The affine hull function for
%  general convex sets will only function for bounded sets. If you want the affine
%  hull of an unbounded set, then intersect your set with a large full-dimensional
%  box.
%  
%  INPUT
%  -----
%     
%        
%          Set Any object derived from the ConvexSet    
%              class, e.g. Polyhedron, YSet, ...        
%              Class: ConvexSet                         
%                
%  
%  
%  OUTPUT
%  ------
%     
%        
%          H The real matrix that defines the affine  
%            hull.                                    
%            Class: double                            
%              
%  
%  
%  SEE ALSO
%  --------
%     innerApprox,  outerApprox,  isEmptySet
%  

%  AUTHOR(s)
%  ---------
%     
%    
%   (c) 2010-2013  Colin Neil Jones: EPF Lausanne
%   mailto:colin.jones@epfl.ch 
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
 
 
global MPTOPTIONS

if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

% use forEach for arrays
error(obj.rejectArray());

% prepare output
aff = [];
if obj.isEmptySet
    return;
end

test  = eye(obj.Dim);
V = []; % Vertex discovered so far
while size(test,2) > 0
    % Choose a test direction in the null space of the basis
    t = test(:,1); test(:,1) = [];
    
    pos = obj.extreme(t);
    neg = obj.extreme(-t);
            
    if pos.exitflag == MPTOPTIONS.INFEASIBLE || pos.exitflag == MPTOPTIONS.ERROR || ...
            neg.exitflag == MPTOPTIONS.INFEASIBLE || neg.exitflag == MPTOPTIONS.ERROR
        error('Infeasible solution returned in non-empty set')
    end
    
    if pos.exitflag == MPTOPTIONS.UNBOUNDED || neg.exitflag == MPTOPTIONS.UNBOUNDED
        error('Can only compute the affine hull of bounded ConvexSets');
    end
    
    gap = pos.supp + neg.supp;
    if gap > MPTOPTIONS.abs_tol % This is a full-dimensional direction
        V = [V;pos.x';neg.x'];
        test = null(V(2:end,:)-ones(size(V,1)-1,1)*V(1,:));
    end
end

if ~isempty(V)
    n = null(V(2:end,:)-ones(size(V,1)-1,1)*V(1,:))';
    aff = [n n*V(1,:)'];
end

end
