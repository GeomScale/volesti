function Un = plus(U,Q)
%
%  PLUS: Minkowski addition for union of polyhedra 
%  ================================================
%  
%  
%  SYNTAX
%  ------
%     
%      U + W
%      U.plus(W)
%    
%  
%  DESCRIPTION
%  -----------
%     Computation of Minkowski addition for the union of polyhedra in the same
%  dimension. The algorithm proceeds in the following way: 
%    
%     1. Compute the Minkowski summation for each of the polyhedron contained in
%     the union to get. 
%     2. Compute the convex hull of the union. 
%     3. Compute the set difference between the convex hull and the union. 
%     4. Compute the set difference between the Minkowski sum for each polyhedron
%     and set obtained in the previous result. 
%    The result is a non-overlapping union of the polyhedra.
%  
%  INPUT
%  -----
%     
%        
%          U Union of polyhedra in the same           
%            dimension.                               
%            Class: PolyUnion                         
%          W Polyhedron to be summed with the union   
%            that is in the same dimension as the     
%            union.                                   
%            Class: Polyhedron                        
%              
%  
%  
%  SEE ALSO
%  --------
%     convexHull,  minus
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
 
 
global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

if ~isa(Q,'Polyhedron'),
    error('Second input argument must be a Polyhedron object!');
end
if numel(Q)~=1 
    error('Only single polyhedron "Q" is accepted.');
end
if numel(Q)>1
    if any(U.Dim~=[Q.Dim])
        error('The polyhedron array "Q" must be in the same dimension as the union.');
    end
end

% if Q is empty, create a new copy and quickly return
if Q.isEmptySet
    Un = PolyUnion(U.Set);
    Un.Internal = U.Internal;
    Un.Data = U.Data;
    return;
end


% MINKOWSKI SUMMATION ON UNIONS OF POLYHEDRA
%=============================================

% compute Minkowski sum for each polyhedron in the array
P = U.Set + Q;

% compute convex hull of the union
Phull=U.convexHull;

% compute minkowski summation for the convex hull
PM = plus(Phull,Q);

% compute the difference between the summation for the convex hull and P
T = PM \ P;

if numel(T)>0
    % discard low-dim polyhedra
    T(~T.isFullDim)=[];
        
    % R is final union
    R = PM \ T;
    
    % discard low-dim polyhedra
    R(~R.isFullDim)=[];
else
    R = PM;
end

% return new PolyUnion object
Un = PolyUnion(R);

% copy internal data including properties
Un.Internal = U.Internal;
Un.Data = U.Data;

% the new union does not have overlaps due to set difference operation
Un.Internal.Overlaps = false;

end
