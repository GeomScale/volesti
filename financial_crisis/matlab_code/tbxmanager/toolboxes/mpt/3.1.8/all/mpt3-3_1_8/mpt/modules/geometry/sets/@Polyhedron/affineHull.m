function [aff, H] = affineHull(P)
%
%  AFFINEHULL: Compute the affine hull of the Polyhedron. 
%  =======================================================
%  
%  
%  SYNTAX
%  ------
%     
%      hull = P.affineHull
%      hull = affineHull(P)
%    
%  
%  DESCRIPTION
%  -----------
%     Computes a minimum-rank matrix hull such that 
%                                       [     ]             
%                         P  { x | hull [ x   ]= 0 }     (1)
%                                       [ -1  ]             
%     Notes: 
%    
%     - If an H-rep is known, this function does not assume He is the affine hull,
%     but rather searches for implicit equalities. 
%     - The Polyhedron P is not modified after the call. 
%  
%  
%  INPUT
%  -----
%     
%        
%          P Polyhedron in any format.                
%            Class: Polyhedron                        
%              
%  
%  
%  OUTPUT
%  ------
%     
%        
%          hull Minimum-rank matrix representing the     
%               affine hull of P.                        
%               Class: double                            
%                 
%  
%  
%  SEE ALSO
%  --------
%     minVRep,  affineMap
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

% use P.forEach() for arrays
error(P.rejectArray());

if P.hasVRep
    if nargout < 2
        % Only the affine hull is requested:
        % Lift to a cone and compute affine set
        R = [P.V ones(size(P.V,1),1);P.R zeros(size(P.R,1),1)];
        aff = null(R)'; aff(:,end) = -aff(:,end);
        return
    elseif ~P.hasHRep
        % Hrep requested, compute Hrep first
        P.computeHRep();
    end
end

% Assume that the affine hull is more restrictive than He, and re-compute
%
% Relax all constraints and test one-by-one
H = P.H; He = P.He;
while 1
    done = true;
    for i = 1:size(H,1)
        Htmp = H;
        
        % Relax i+1:end constraints to make the polyhedron full-dim
        Htmp(i+1:end,end) = Htmp(i+1:end,end) + 1;
        
        % Test if the result is full-dimensional
        if ~fast_isFullDim(Htmp, He)
            % This constraint cannot be tightened => it's part of the affine hull
            He = [He;H(i,:)];
            H(i,:) = [];
            done = false;
            break
        end
    end
    if done, break; end
end

% Compute a min-rep of the affine hull of this polyhedron
aff = mpt_minAffineRep(He);

end
