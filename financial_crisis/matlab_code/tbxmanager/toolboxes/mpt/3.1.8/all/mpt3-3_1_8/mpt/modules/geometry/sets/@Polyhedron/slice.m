function S = slice(P, dims, values, varargin)
%
%  SLICE: Slice the polyhedron through given dimensions at specified values. 
%  ==========================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      S = P.slice(dims, values)
%      S = slice(P, dims, values)
%      S = P.slice(dims, values, 'keepDim', true/false)
%    
%  
%  DESCRIPTION
%  -----------
%     Compute the intersection of P  with a subspace spanning the dimensions dims
%  at given values. If the argument values are omitted, the value is assumed to be
%  zero, i.e. values = zeros(size(dims)). For a polyhedron given in
%  H-representation 
%                                                          
%                         P = x  |  Ax <= b,  A  x = b   , 
%                                              eq     eq   
%     the slice operation over dims at given values returns a polyhedron in a
%  reduced dimension P.Dim-length(dims) 
%                                                                  
%                  S =  x  |  A        x <= b - A         values . 
%                              (:,keep)          (:, dims)         
%     This corresponds to the default choice keepdim=false.
%    Alternatively, by invoking keepdim=true, the polyhedron S will be returned in
%  the same dimension as P 
%                                                                    
%               S =  x  |  Ax <= b,  A  x = b  , x(dims) == values . 
%                                     eq     eq                      
%  
%  
%  INPUT
%  -----
%     
%        
%          P      Polyhedron in any format                 
%                 Class: Polyhedron                        
%          dims   Dimensions to cut through                
%                 Class: double                            
%          values Set of values at which to compute the    
%                 slices.                                  
%                 Class: double                            
%                 Default: 0                               
%                   
%  
%  
%  OUTPUT
%  ------
%     
%        
%          S Polyhedron that represents the cut of    
%            the polyhedron P over the specified      
%            dimensions.                              
%            Class: Polyhedron                        
%              
%  
%  
%  SEE ALSO
%  --------
%     projection,  project
%  

%  AUTHOR(s)
%  ---------
%     
%    
%   (c) 2010-2013  Colin Neil Jones: EPF Lausanne
%   mailto:colin.jones@epfl.ch 
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
 
 
global MPTOPTIONS

narginchk(2, Inf);
if nargin<3
	values = zeros(size(dims));
end
if nargin<4
	options.keepDim = false;
end
if nargin>3
	% parse options
	ip = inputParser;
	ip.addParamValue('keepDim', false, @islogical);
	ip.parse(varargin{:});
	options = ip.Results;
end


%% deal with arrays
if numel(P)>1
	S = P.forEach(@(e) e.slice(dims, values, varargin{:}));
	return
end
        
%% validation
if P.isEmptySet
    error('Cannot slice empty polyhedra.');
end
% check dimensions
for i=1:numel(dims)
    validate_dimension(dims(i));
end
if any(dims>P.Dim)
    error('The second input cannot exceed dimension of the polyhedron.');
end
if numel(values)~=numel(dims)
	error('"values" must be a vector with %d element(s).', numel(dims));
end

%% computation

% require the H-representation (the getters computes it automatically if
% it does not exist)
if options.keepDim
	% keep dimension of the slice equal to the dimension of the polyhedron
	Z = zeros(numel(dims), P.Dim);
	for i = 1:numel(dims)
		Z(i, dims(i)) = 1;
	end
	S = Polyhedron('H', P.H, 'He', [P.He; [Z, values(:)]]);

else
	% return polyhedron in lower dimension
	keep_dims = setdiff(1:P.Dim, dims);
	remove_dims = dims;
	if isempty(P.He)
		% faster constructor
		S = Polyhedron(P.A(:, keep_dims), P.b-P.A(:, remove_dims)*values(:));
		
	else
		% check that the slice values satisfy equality constraints:
		%
		% find z s.t. P.Ae(:, keep)*z + P.Ae(:, remove)*values == P.be
		nz = length(keep_dims);
		lp.f = zeros(1, nz);
		lp.A = P.A(:, keep_dims);
		lp.b = P.b - P.A(:, remove_dims)*values(:);
		lp.Ae = P.Ae(:, keep_dims);
		lp.be = P.be - P.Ae(:, remove_dims)*values(:);
		lp.lb = []; lp.ub = [];
		lp.quicklp = true;
		sol = mpt_solve(lp);
		
		if sol.exitflag==MPTOPTIONS.OK
			% equality constraints are consistent
			S = Polyhedron('A', P.A(:, keep_dims), ...
				'b', P.b-P.A(:, remove_dims)*values(:), ...
				'Ae', P.Ae(:, keep_dims), ...
				'be', P.be-P.Ae(:, remove_dims)*values(:));
		else
			% equality constraints inconsistent => empty polyhedron of
			% correct dimension
			S = Polyhedron.emptySet(length(keep_dims));
		end
	end
	
	% slice functions
	for n = P.Functions.keys
		name = n{1};
		new = P.Functions(name).slice(dims, values);
		S.addFunction(new, name);
	end
end

end
