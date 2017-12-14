function [index, details] = locatePoint(U,x)
%
%  LOCATEPOINT: Implementation of a graph search algorithm for a point location
%  ============================================================================
%  problem. 
%  =========
%  
%  
%  SYNTAX
%  ------
%     
%      [index, details] = locatePoint(U,x)
%    
%  
%  DESCRIPTION
%  -----------
%     Return an index of region in the polyunion U where the point x is contained.
%  If the point lies outside of the union of polyhedra, the output is empty. The
%  polyunion U must be returned from PLCP solver with an appropriate adjacency
%  list, otherwise the method is not applicable. The best performance is achieved
%  if the region exploration in PLCP solver has been done using breadth-first
%  search (BFS) where the regions have a leveled structure and the first region
%  lies in the middle of the partition. In this case, the graph traversal algorithm
%  proceeds through increasing levels outside from the first region until the
%  desired region is found. The set membership operation depends on the settings of
%  absolute tolerance that can be changed in Options settings.
%  
%  INPUT
%  -----
%     
%        
%          U Union of polyhedra returned from PLCP    
%            solver with an included adjacency list.  
%            Class: PolyUnion                         
%          x A point in the same dimension as         
%            PolyUnion given as real column vector.   
%            Class: double                            
%              
%  
%  
%  OUTPUT
%  ------
%     
%        
%          index   Index of a region where x  is contained, 
%                  otherwise an empty output [].            
%                  Class: double                            
%          details A structure with statistical information 
%                  numerical performance of the graph       
%                  traversal algorithm.                     
%                  Class: struct                            
%                  Allowed values:                          
%                                                           
%                    niter  Number of iterations.           
%                    noper  Total number of arithmetic      
%                     operations                            
%                     (multiplications+summations+compariso 
%                     ns).                                  
%                    multiplications  Multiplications       
%                     count.                                
%                    summations  Summation count.           
%                    comparisons  Comparisons count.        
%                                                           
%                    
%  
%  
%  SEE ALSO
%  --------
%     contains,  isInside
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
 
 
narginchk(2, 2);
% use U.forEach(@(u) u.contains(x)) to evaluate arrays
error(U.rejectArray());

% empty polyunion
if numel(U)==0 || U.Num<1
    index = [];
    details = [];
    return;
end

% check x
validate_realvector(x);
x = x(:);
if numel(x)~=U.Dim
    error('The point must have the same dimension as the union.');
end


% check if the adjacency list is present
if ~isfield(U.Internal,'adj_list')
    error('The union does not have an adjacency list. Please, use the polyunion output from PLCP solver that contains the adjacency list.');
end

% check the properties of the union
if isempty(U.Internal.Convex) || isempty(U.Internal.Overlaps) || ...
    isempty(U.Internal.Connected) || isempty(U.Internal.Bounded) || ...
    isempty(U.Internal.FullDim)        
    disp('Some of the required properties have not been determined for this union.')
    disp('The following properties are checked: convexity, overlaps, connectivity, boundedness, and full-dimensionality.');
    disp('This may take some time...');
end
if ~(U.isBounded && U.isFullDim && U.isConnected && U.isConvex && ~U.isOverlapping)
    error(['This method supports unions of polyhedra that are convex, non-overlapping, '... 
        'bounded, full-dimensional, connected, and come from PLCP solver with an adjacency list.']);
end

        
% call the graph traversal method
[index, details] = find_region(x,U.Set, U.Internal.adj_list);


end

function [index, details] = find_region(x,Parray,graph, index, details)
%
% For given point x, find a region index from the partition Parray of the
% explicit solution where the point lies.
%
% input:
%   x - value of parameter
%   Parray - an array of regions 
%   graph - represented by adjacency list
%

global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS=mptopt;
end

DEBUG=0;

if nargin<4
    index = 1;
    details.niter = 0;
    details.noper = 0;
    details.multiplications = 0;
    details.summations = 0;
    details.comparisons = 0;
end
if nargin<5
    details.niter = 0;
    details.noper = 0;
    details.multiplications = 0;
    details.summations = 0;
    details.comparisons = 0;
end

details.niter = details.niter + 1;

% extract polyhedron
P = Parray(index);

if DEBUG
    if index<=1
        plot(Parray); 
        hold on;
        text(x(1),x(2),'x');
    end
    xc = chebyCenter(P);
    text(xc.x(1),xc.x(2),num2str(index));
end  

% get direction
%P.normalize; % polytopes are normalized when returned from PLCP solver
d = P.A*x - P.b; % compute distance
direction = d > MPTOPTIONS.abs_tol;
pd = find(direction);

% count operations
n=P.Dim;
nv = size(P.H,1);
details.multiplications = details.multiplications + nv*n;
details.summations = details.summations + nv*n;
details.comparisons = details.comparisons + nv;
details.noper = details.multiplications + details.summations + details.comparisons;

if any(direction)
    % pick the direction with maximum distance over positive d
    [~,id]= max(d(pd));
    direction = false(size(direction));
    direction(pd(id))=true;
    
    % add comparisons
    details.comparisons = details.comparisons + numel(pd)-1;
    details.noper = details.multiplications + details.summations + details.comparisons;
    
    % next region to check
    index_new = graph{index}(direction);
    if length(index_new)>1
        % if there are more ways, pick the one which is not empty
        v = cellfun(@isempty,index_new);
        index_new = index_new(~v);
        if isempty(index_new)
            % point out of feasible area
            index = [];
            return;
        end
    end

    % if there are more choices, pick first    
    index_new=index_new{1};
    if numel(index_new)>1
        % if the facet neighbors to more regions, pick the first
        index_new=index_new(1);
    end
    
    if isempty(index_new)
        % point out of feasible area
        index = [];
    else
        % if x does not lie in this region, continue
        [index,details] = find_region(x,Parray,graph, index_new, details);
    end
    
end

end
