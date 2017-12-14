function kept = reduce(U)
%
%  REDUCE: Reduces the overlapping union to minimal number of sets. 
%  =================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      kept = U.reduce
%      kept = reduce(U)
%    
%  
%  DESCRIPTION
%  -----------
%     Simplifies the union of polyhedra by removing those regions who are
%  completely covered by others. Note that this algorithm is valid only for
%  overlapping unions.
%  
%  INPUT
%  -----
%     
%        
%          U    Union of polyhedra in the same           
%               dimension.                               
%               Class: PolyUnion                         
%          kept Vector of logical indices indicating if  
%               the region is non-redundant.             
%               Class: logical                           
%               Allowed values:                          
%                                                        
%                 true                                   
%                 false                                  
%                                                        
%                 
%  
%  
%  SEE ALSO
%  --------
%     merge,  isOverlapping
%  

%  AUTHOR(s)
%  ---------
%     
%    
%   (c) 2003-2013  Michal Kvasnica: STU Bratislava
%   mailto:michal.kvasnica@stuba.sk 
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

% deal with arrays
if numel(U)>1
    kept = cell(size(U));
    for i=1:numel(U)
        kept{i} = U(i).reduce;
    end
    return;
end


% if there is 0 or 1 set contained, return
if U.Num<=1
    kept = 1;
    return
end

% if regions do not overlap, we can exit quickly
if ~U.isOverlapping
	kept = 1:U.Num;
	return
end

% sort regions by size for better performance.
% reasoning: at each step we need to check whether polytope P(i) is covered
%            by all other remaining polytopes. it is more likely that small
%            regions will be kicked out at the beginning, reducing the size
%            of the polytope array we need to check P(i) against.

% get the radius of the polyhedra (if possible) to sort the polyhedra from
% the smallest 
r = zeros(1,U.Num);
for i=1:U.Num
    % compute chebycenter
    ip = U.Set(i).interiorPoint;
    if ~isempty(ip.r)
        r(i) = ip.r;
    else
        % no chebycenter available, estimate the radius based on the
        % average distance of the vertices to the center point
        d=U.Set(i).V-repmat(ip.x',size(U.Set(i).V,1),1);
        r(i) = mean(sqrt(sum(d.*d,2)));
    end
end
[r, index] = sort(r);


% prepare output
kept=true(1, U.Num);
k=1;

for selected = index,
    if MPTOPTIONS.verbose>=1
        fprintf('Region: %d/%d\n',k,U.Num);
    end
    k = k+1;
    regions_kept = find(kept);
    to_check = setdiff(regions_kept, selected);
    if isempty(to_check)
        % we kicked out all other regions, return
        break
    end
    % peform the containment test backwards because there is high
    % probability for overlap
    nc = numel(to_check);
    for j=nc:-1:1
        if U.Set(to_check(j)).contains(U.Set(selected))
            kept(selected) = false;
        end
    end
end


% write output
regions_kept = find(kept);
if isempty(regions_kept)
    U.Set = [];
else
    U.Set = U.Set(regions_kept);
end

end
