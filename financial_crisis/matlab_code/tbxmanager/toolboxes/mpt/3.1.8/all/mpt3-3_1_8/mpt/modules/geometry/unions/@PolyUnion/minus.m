function Un = minus(U,Q)
%
%  MINUS: Pontryagin/Minkowski difference for union of polyhedra 
%  ==============================================================
%  
%  
%  SYNTAX
%  ------
%     
%      U - W
%      U.minus(W)
%    
%  
%  DESCRIPTION
%  -----------
%     Computation of Pontryagin or Minkowski difference for the union of polyhedra
%  in the same dimension. The algorithm for efficiently computing the Minkowski
%  difference between a union of polytopes and a polytope is based on a talk by S.
%  Rakovic and D. Mayne entitled Constrained Control ComputationsIt was the keynote
%  addressed at the GBT Meeting, London, November, 2002. The algorithm proceeds in
%  the following way: 
%    
%     1. Compute the convex hull of the union. 
%     2. Compute the Minkowski difference of each of the polyhedron from the convex
%     hull. 
%     3. Compute the set difference between the convex hull and the union. 
%     4. Compute the set difference between the Minkowski difference for each
%     polyhedron and the set obtained in the previous step. 
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
%     convexHull,  plus
%  

%  AUTHOR(s)
%  ---------
%     
%    
%   (c) 2010-2013  Martin Herceg: ETH Zurich
%   mailto:herceg@control.ee.ethz.ch 
%     
%    
%   (c) 2005  Mario Vasak: FER Zagreb
%   mailto:mario.vasak@fer.hr 
%     
%    
%   (c) 2003-2013  Michal Kvasnica: STU Bratislava
%   mailto:michal.kvasnica@stuba.sk 
%     
%    
%   (c) 2003-2005  Pascal Grieder: ETH Zurich
%   mailto:grieder@control.ee.ethz.ch 
%     
%    
%   (c) 2003  Mato Baotic: ETH Zurich
%   mailto:baotic@control.ee.ethz.ch 
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
    error('Second input argument must be a Polyhedron object');
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
if Q.isEmptySet || U.Num==0
    Un = PolyUnion(U.Set);
    Un.Internal = U.Internal;
    Un.Data = U.Data;
    return;
end

% MINKOWSKI DIFFERENCE ON UNIONS OF POLYHEDRA
%=============================================

% if we want to subtract full-dimensional polyhedron Q from any
% low-dimensional polyhedron- this will be empty, thus we remove any
% low-dim polyhedra first to simplify computations
if Q.isFullDim
    U.Set(~U.Set.isFullDim) = [];
end

% compute convex hull
Phull=U.convexHull;

% compute minkowski difference
PM = minus(Phull,Q);

% E is union of polytopes; sets inside Phull which are not covered by PA
E = Phull \ U.Set;

if numel(E)>0
    % discard low-dim polyhedra if all of U are full-dim
    E(~E.isFullDim)=[];

end
    
if numel(E)>0
    % flip Polytope
    QM = uminus(Q);
    
    % Minkowski addition on QM
    Emin = plus(E,QM);
    
    % Compute final polytopes
    R = PM \ Emin;
else
    R = PM;
end

% discard low-dim polyhedra if all of U are bounded
R(~R.isFullDim)=[];

% return new PolyUnion object
Un = PolyUnion(R);

% copy internal data including properties
Un.Internal = U.Internal;
Un.Data = U.Data;

% the new union does not have overlaps due to set-difference operation
if U.isFullDim
    Un.Internal.Overlaps = false;
end


end
