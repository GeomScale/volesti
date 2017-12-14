classdef QuadFunction < Function
%
%  QUADFUNCTION: Representation of quadratic functions in the form x'*H*x + F*x + g
%  ================================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      Q = QuadFunction(H,F,g)
%      Q = QuadFunction(H,F)
%      Q = QuadFunction(H)
%      Q = QuadFunction(H,F,g,Data)
%    
%  
%  DESCRIPTION
%  -----------
%     The QuadFunction class represents quadratic functions of the form f(x) =
%  x'*H*x + F*x + g  where H  is a real matrix, F  is a real matrix and gis a real
%  column vector. Dimensions of H, F  and g  must coincide such that the output is
%  a scalar.
%  
%  INPUT
%  -----
%     
%        
%          H    Real matrix representing the             
%               coefficients in the quadratic term H  in 
%               f(x) = x'*H*x + F*x + g .                
%               Class: double                            
%          F    Real matrix representing the             
%               coefficients in the linear term F  in    
%               f(x) = F*x + g .                         
%               Class: double                            
%          g    Real vector representing the affine      
%               terms g  in f(x) = F*x + g .             
%               Class: double                            
%          Data Any data related to the function.        
%                 
%  
%  
%  OUTPUT
%  ------
%     
%        
%          Q QuadFunction object.                     
%              
%  
%  
%  SEE ALSO
%  --------
%     Function,  AffFunction
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
 
 
	%
	% class for representing quadratic functions x'*H*x + F*x + g
	%
	% syntax: Q = QuadFunction(H,F,g)
	%         Q = QuadFunction(H,F)
	%         Q = QuadFunction(H)
	%         Q = QuadFunction(H,F,g,Data)
	
	properties (SetAccess=private)
		H % quadratic term
		F % linear term
		g % affine term
		D=0; % dimension of the domain
		R=0; % dimension of the range
	end
	properties(Dependent=true, SetAccess=private, Transient=true)
		weight % used to access the leading term in penalties
	end
	
	methods
		
		function weight = get.weight(obj)
			% QuadFunction.weight getter
			%
			% Returns the leading term of x'*H*x+F*x+g, i.e., H, when the
			% function is employed to represent a penalty
			weight = obj.H;
		end
		
		% Constructor
		function obj = QuadFunction(varargin)
			%
			% sets data for QuadFunction object
			% syntax: Q = QuadFunction(H,F,g)
			%         Q = QuadFunction(H,F)
			%         Q = QuadFunction(H)
			%         Q = QuadFunction(H,F,g,Data)
			%
			% for more details, type "help QuadFunction"
			
			global MPTOPTIONS
			if isempty(MPTOPTIONS)
				MPTOPTIONS = mptopt;
			end
			
			narginchk(1, 4);
			
			% check H
			Hm = varargin{1};
			if ~isa(Hm, 'sdpvar')
				validate_realmatrix(Hm);
			end
			if size(Hm,1)~=size(Hm,2)
				error('The matrix "H" must be square.');
			end
			
			% assign H
			if issparse(Hm)
				Hm = full(Hm);
			end
			obj.H = Hm;
			
			% get the dimension of the domain
			obj.D = size(Hm,1);
			
			% get the dimension of the range
			obj.R = 1;
			
			% only H provided
			if nargin==1
				obj.F = zeros(obj.R,obj.D);
				obj.g = zeros(obj.R,1);
			end
			
			% H, F provided
			if nargin>1
				% F is provided, check
				Fm = varargin{2};
				if ~isa(Fm, 'sdpvar')
					validate_realmatrix(Fm);
				end
				if size(Fm,1)~=obj.R
					error('The number of rows for matrix "F" must be %d.',obj.R);
				end
				if size(Fm,2)~=obj.D
					error('The number of columns for matrix "F" must be %d.',obj.D);
				end
				if issparse(Fm)
					Fm = full(Fm);
				end
				obj.F = Fm;
				obj.g = zeros(obj.R,1);
			end
			
			% H, F, g provided
			if nargin>2
				% g is provided, check
				gm = varargin{3};
				if ~isa(gm, 'sdpvar')
					validate_realvector(gm);
				end
				% make column vector
				if issparse(gm)
					gm = full(gm);
				end
				gm = gm(:);
				if length(gm)~=obj.R
					error('The vector "g" must be of the size %d.',obj.R);
				end
				obj.g = gm;
			end
			
			% Data provided
			if nargin>3
				obj.Data = varargin{4};
			end
			
			% full syntax
			obj.Handle = @obj.qf;
			
		end
		
		function status = eq(f, g)
			% Returns true if the functions are identical
			
			global MPTOPTIONS
			if isempty(MPTOPTIONS)
				MPTOPTIONS = mptopt;
			end
			
			if numel(f)~=numel(g)
				error('Matrix dimensions must agree.');
			elseif numel(f)>1
				% for arrays
				status = false(size(f));
				for i = 1:numel(f)
					status(i) = (f(i)==g(i));
				end
			else
				% for scalars
				
				% TODO: maybe we could alllow comparing AffF with QuadF that
				% has zero quadratic term
				status = isa(f, 'QuadFunction') && isa(g, 'QuadFunction') && ...
					(f.R==g.R) && (f.D==g.D) && ...
					norm(f.H-g.H)<MPTOPTIONS.zero_tol && ...
					norm(f.F-g.F)<MPTOPTIONS.zero_tol && ...
					norm(f.g-g.g)<MPTOPTIONS.zero_tol && ...
					isequal(f.Data, g.Data);
			end
		end
		
		function status = ne(f, g)
			% Returns true if two functions are not identical
			
			status = ~eq(f, g);
		end
		function new = slice(obj, dims, values)
			% Slice a quadratic function through given coordinates
			%
			% When f=x'*H*x + F*x+g, then f.slice(dims, values) produces a
			% new function fn=z'*Hn*z + Fn*z + gn where the "dims"
			% coordinates of "x" are replaced by fixed "values".
			narginchk(3, 3);
			
			% validation
			for i=1:numel(dims)
				validate_dimension(dims(i));
			end
			if any(dims>obj.D)
				error('Dimension must not exceed %d.', obj.D);
			end
			if numel(values)~=numel(dims)
				error('"values" must be a vector with %d element(s).', numel(dims));
			end
			values = values(:);
			% first re-order matrices H and F such that the dimensions to
			% be retained are the first, followed by those which should be
			% eliminated, i.e.:
			%         [H1 H2] [z]           [z]
			% [z' v']*[H3 H4]*[v] + [F1 F2]*[v] + g
			%
			% Note that we account for non-symmetric H-matrix
			keep = setdiff(1:obj.D, dims);
			neworder = [sort(keep(:)); dims(:)];
			Hr = obj.H(neworder, neworder);
			Fr = obj.F(:, neworder);
			nd = numel(keep);
			H1 = Hr(1:nd, 1:nd);
			H2 = Hr(1:nd, nd+1:end);
			H3 = Hr(nd+1:end, 1:nd);
			H4 = Hr(nd+1:end, nd+1:end);
			F1 = Fr(:, 1:nd);
			F2 = Fr(:, nd+1:end);
			
			% Fixing "u=values" gives the new function
			%   z'*H1*z + (F1 + v'*H2' + v'*H3)*z + ( v'*H4*v + F2*v + g)
			Hn = H1;
			Fn = F1 + values'*(H2' + H3);
			gn = values'*H4*values + F2*values + obj.g;
			new = QuadFunction(Hn, Fn, gn);
		end
		
	end
	methods (Hidden)
		function y=qf(obj,x)
			
			if ~isequal(size(x), [obj.D, 1])
				error('The input must be a %dx1 vector.', obj.D);
			end
			% output
			y = x'*obj.H*x + obj.F*x + obj.g;
		end
	end
end
