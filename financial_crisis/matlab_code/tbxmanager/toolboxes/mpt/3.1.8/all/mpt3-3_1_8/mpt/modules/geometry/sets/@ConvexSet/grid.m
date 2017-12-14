function X=grid(P, N)
%
%  GRID: Grid the convex set. 
%  ===========================
%  
%  
%  SYNTAX
%  ------
%     
%      x = grid(Set, N)
%      x = Set.grid(N)
%    
%  
%  DESCRIPTION
%  -----------
%     Gridding of the convex Set with respect to N linearly scaled gridding points.
%  The output x consist of points sorted vertically that belonging to the Set. The
%  principle of the algorithm is as follows: 
%    
%     1. Compute outer bounding hypercube. 
%     2. Grid the hypercube. 
%     3. Test each point for inclusion in the set, discarding those outside. 
%    Before running the algorithm, consider the number N of gridding points. If
%  this number is very large then it takes algorithm longer to grid the space
%  because an exhaustive search is done at the last step of the algorithm.
%  
%  INPUT
%  -----
%     
%        
%          Set Any object derived from the ConvexSet    
%              class, e.g. Polyhedron, YSet, ...        
%              Class: ConvexSet                         
%          N   The number of gridding point. It         
%              specifies with how many elements to      
%              scale the interval equidistantly.        
%                
%  
%  
%  OUTPUT
%  ------
%     
%        
%          x An array of points sorted vertically.    
%            The number of columns specifies the      
%            dimension.                               
%            Class: double                            
%              
%  
%  
%  SEE ALSO
%  --------
%     plot,  fplot,  contains
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
    MPTOPTIONS=mptopt;
end

narginchk(2, 2);
error(P.rejectArray());

if P.isEmptySet
    error('Empty set, there is nothing to be gridded here.');
end

if ~P.isBounded
    error('Can only grid bounded sets.'); 
end

% get bounding box
bbox = P.outerApprox;

lb = bbox.Internal.lb;
ub = bbox.Internal.ub;


% grid the state-space into equidistantly placed points
dimB = size(lb(:),1);
Xpoints = zeros(N, dimB);
for ii=1:dimB
    Xpoints(:,ii) = linspace(lb(ii),ub(ii),N)';
end

% generate all possible combinations of states
% one could use kron() here, but that one fails for high number of elements
n_states = dimB;
ZZ=[];
ZZ{n_states}=Xpoints(:,n_states);
for ii=n_states-1:-1:1,
    Zd=[];
    for jj=1:size(Xpoints,1),
        Zd=[Zd; repmat(Xpoints(jj,ii),length(ZZ{ii+1}),1)];
    end
    ZZ{ii}=Zd;
end
for ii=2:n_states,
    ZZ{ii}=repmat(ZZ{ii},length(ZZ{ii-1})/length(ZZ{ii}),1);
end
datapoints=[];
for ii=1:n_states,
    datapoints(:,ii) = ZZ{ii};
end
npoints = size(datapoints,1);
isin = false(npoints, 1);
for i = 1:npoints,
    isin(i) = P.contains(datapoints(i,:)');
end
X = datapoints(isin, :);



end
