function [isin, inwhich, closest] = contains(U, x, fastbreak)
%
%  CONTAINS: Test if a point is contained inside the union of convex sets. 
%  ========================================================================
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
%     Check if the point x is contained in any of the sets in the union U. The
%  result it a logical statement if x in U  and false otherwise. If the point is
%  contained inside the union, indices of sets in which the point lie are returned.
%  If the point does not lie in the union, the index of the region with the least
%  distance to the point x  is returned. All sets in the union must have the same
%  dimension, otherwise the evaluation cannot be done.
%  
%  INPUT
%  -----
%     
%        
%          U         Single union object that holds sets      
%                    derived from the ConvexSet class. If Uis 
%                    an array, use U.forEach().               
%                    Class: Union                             
%          x         A point in the same dimension as all the 
%                    sets in the union.                       
%                    Class: double                            
%          fastbreak Do a quick stop in the consecutive       
%                    search when x is contained in the first  
%                    set it founds.                           
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
error(U.rejectArray());

if numel(U)==0
	isin = [];
	inwhich = [];
	closest = [];
	return
end
if nargin<3
    fastbreak = false;
end

%% validation
validate_realvector(x);
if size(x, 2)~=1
	error('The point must be a column vector.');
end
nx = numel(x);
iscell_set = iscell(U.Set);
for i = 1:U.Num
	if ( iscell_set && U.Set{i}.Dim ~= nx ) || ...
			( ~iscell_set && U.Set(i).Dim ~= nx )
		error('All sets must be %d-dimensional.', nx);
	end
end

%% search
isin = false;
inwhich = [];
closest = [];
for i = 1:U.Num
	if ( iscell_set && U.Set{i}.contains(x) ) || ...
			(~iscell_set && U.Set(i).contains(x))
		isin = true;
		inwhich = [inwhich, i];
		if fastbreak
			return
		end
	end
end

%% find closest region if necessary
if ~isin && nargout>2
	d = Inf(1, U.Num);
	for i=1:U.Num
		if iscell_set
			s = distance(U.Set{i}, x);
		else
			s = distance(U.Set(i), x);
		end
		if ~isempty(s.dist)
			d(i) = s.dist;
		end
	end
	[~, closest] = min(d);
end

end
