function tf = isBounded(P)
%
%  ISBOUNDED: Test if a polyhedron is bounded. 
%  ============================================
%  
%  
%  SYNTAX
%  ------
%     
%      tf = P.isBounded
%      tf = isBounded(P)
%    
%  
%  DESCRIPTION
%  -----------
%     Return true if the polyhedron P  is bounded and false otherwise.
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
%          tf true if the polyhedron P  is bounded and 
%             false otherwise.                         
%             Class: logical                           
%             Allowed values:                          
%                                                      
%               true                                   
%               false                                  
%                                                      
%               
%  
%  
%  SEE ALSO
%  --------
%     isEmptySet,  isFullDim,  isBounded
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

% only vectors please
if ~isvector(P)
    error('Input argument must be vector.')
end

% deal with arrays of polyhedra
tf = false(size(P));
for i=1:length(P)
    % compute ChebyCenter only if P(i).Empty property has not been set
    if builtin('isempty',P(i).Internal.Bounded)        
        if P(i).isEmptySet
            % empty polyhedron is considered as bounded
            P(i).Internal.Bounded = true;
        elseif P(i).hasVRep
            P(i).Internal.Bounded = size(P(i).R,1) == 0;
        elseif P(i).hasHRep && ~P(i).isFullDim()
            He = P(i).affineHull();
            if ~isempty(He)
                % project the lower-dimensional polyhedron onto its affine hull
                % and check boundedness there (issue #112)
                P(i).Internal.Bounded = P(i).projectOnAffineHull(He).isBounded();
            else
                % note that we declare a set lower-dimensional if the radius of
                % its chebyshev's ball drops below a tolerance. that does not
                % imply there are equalities, that's why we also need to check
                % whether we in fact have a non-empty affine hull, otherwise we
                % would start to cycle (test_polyhedron_isbounded_17_pass)
                %
                % hence we check boundedness by computing the outer box
                % approximation
                P(i).outerApprox();
                P(i).Internal.Bounded = all(~isinf(P(i).Internal.lb)) & ...
                    all(~isinf(P(i).Internal.ub));
            end
            
        elseif P(i).hasHRep
            sol = P(i).chebyCenter;
            %P(i).Internal.Bounded = true;
            if sol.exitflag == MPTOPTIONS.UNBOUNDED || sol.r==Inf
                P(i).Internal.Bounded = false;
                continue;
            end
            % we cannot rely only on the Chebyshev center due to possible coplanar hyperplanes
            [n,m] = size(P(i).A);
            if m >= n,     % we have more variables than constraints: unbounded
                P(i).Internal.Bounded = false;
                continue;
            end
            nullA = null(P(i).A');
            colsNullA = size(nullA,2);
            rankA = n - colsNullA;
            
            if rankA < m,  % linearly dependant constraints
                P(i).Internal.Bounded = false;
                continue;
            end

            % we're using the Minkowski's theorem on polytopes:
            % Given a_1, ..., a_m unit vectors, and x_1, ... x_m > 0, there
            % exists a polytope having a_1,...,a_m as facets and x_1, ... x_m as
            % facet areas iff:
            %
            %  a_1 x_1 + ... + a_m x_m = 0
            %
            % Hence, checking boundedness of a polytope boils down to feasibility
            % LP.
            
            S.A = -nullA;
            S.b = -ones(n,1);
            S.f = zeros(1,colsNullA);
			S.Ae = []; S.be = []; S.lb = []; S.ub = []; S.quicklp = true;
            res = mpt_solve(S);
            
            if res.exitflag == MPTOPTIONS.OK % problem is feasible (unconstrained case)
                P(i).Internal.Bounded = true;
                % cannot actually happen (a_i are
                % not co-planar and sum(a_i x_i) =
                % 0 with x_i > 0)
                %  => polyhedron is bounded
            else
                % problem is infeasible => polyhedron is unbounded
                P(i).Internal.Bounded = false;
            end
            
        end
    end
    tf(i) = P(i).Internal.Bounded;
end
end
