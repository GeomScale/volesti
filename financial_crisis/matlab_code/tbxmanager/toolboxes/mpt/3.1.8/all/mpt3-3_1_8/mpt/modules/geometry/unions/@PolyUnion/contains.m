function [isin, inwhich, closest] = contains(obj, x, fastbreak)
%
%  CONTAINS: Test if a point is contained inside the union of polyhedra in the same
%  ================================================================================
%  dimension. 
%  ===========
%  
%  
%  SYNTAX
%  ------
%     
%      [isin, inwhich, closest] = contains(U, x, fastbreak)
%      [isin, inwhich, closest] = U.contains(x)
%      [isin, inwhich, closest] = U.contains(x, fastbreak)
%    
%  
%  DESCRIPTION
%  -----------
%     Check if the point x is contained in any of the polyhedra stored in the union
%  U. The result it the logical statement if x in U  and false otherwise. If the
%  point is contained inside the union, indices of the corresponding polyhedra in
%  which the point lie are returned. If the point does not lie in the union, the
%  index of the region with the least distance to the point x  is returned.
%  
%  INPUT
%  -----
%     
%        
%          U         Single PolyUnion object that holds sets  
%                    polyhedra in the same dimension.         
%                    Class: PolyUnion                         
%          x         A point in the same dimension as the     
%                    union.                                   
%                    Class: double                            
%          fastbreak Do a quick stop in the consecutive       
%                    search when x is contained in the first  
%                    polyhedron it founds.                    
%                    Class: logical                           
%                    Allowed values:                          
%                                                             
%                      true                                   
%                      false                                  
%                                                             
%                    Default: false                           
%                      
%  
%  
%  OUTPUT
%  ------
%     
%        
%          isin    True if x in U  and false otherwise.     
%                  Class: logical                           
%                  Allowed values:                          
%                                                           
%                    true                                   
%                    false                                  
%                                                           
%          inwhich Indices of sets that contain x. If the   
%                  fastbreak option is turned on, a single  
%                  index is returned.                       
%                  Class: double                            
%          closest If the point is not contained inside the 
%                  union, this output indicates the index   
%                  of the set that is the closest to the    
%                  point x. Note: since this computation is 
%                  expensive, do not ask for the third      
%                  output unless you really need it.        
%                  Class: double                            
%                    
%  
%  
%  SEE ALSO
%  --------
%     feval
%  

%  AUTHOR(s)
%  ---------
%     
%    
%   (c) 2010-2013  Martin Herceg: ETH Zurich
%   mailto:herceg@control.ee.ethz.ch 
%     
%    
%   (c) 2003-2013  Michal Kvasnica: STU Bratislava
%   mailto:michal.kvasnica@stuba.sk 
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
 
 
narginchk(2, 3);
% use obj.forEach(@(u) u.contains(x)) to evaluate arrays
error(obj.rejectArray());

%% parse inputs
if nargin<3
	fastbreak = false;
end
check_hull_first = false;

isin = [];
inwhich = [];
closest = [];
if numel(obj)==0 || ( numel(obj)==1 && obj.Num==0 )
	return
end

%% validation
validate_realvector(x);
if size(x, 2)~=1 || numel(x)~=obj.Dim
	error('The point must be a %dx1 vector.', obj.Dim);
end

%% heuristics
if check_hull_first
	% check the convex hull first, maybe we could exit quickly.
	%
	% in theory, this is a great idea. in practice, it does not pay out.

	if numel(obj)==1 && isfield(obj.Internal, 'convexHull') && ...
			~isempty(obj.Internal.convexHull)
		H = obj.Internal.convexHull;
		if ~H.contains(x)
			% not in the convex hull, no point in evaluating further
			isin = false;
			if nargout==3
				% computing the closest region is slow anyhow, so there is
				% no harm in going via Union/contains
				[~, ~, closest] = obj.contains@Union(x, fastbreak);
			end
			return
		end
	end
end

%% search
% exploit Polyhedron/contains operating on arrays
c = obj.Set.contains(x, fastbreak);
isin = any(c);
% always return "inwhich" as a row vector
inwhich = find(c); 
inwhich = inwhich(:)';

%% find closest region if necessary
if nargout==3 && ~isin
	% computing the closest region is slow anyhow, so there is no harm in
	% going via Union/contains
	[~, ~, closest] = obj.contains@Union(x, true);
end

end
