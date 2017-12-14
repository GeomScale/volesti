classdef NormFunction < Function
%
%  NORMFUNCTION: Class representing 1- or infinity-norm function. 
%  ===============================================================
%  
%  
%  SYNTAX
%  ------
%     
%      f = NormFunction(flag)
%      f = NormFunction(flag, Q)
%    
%  
%  DESCRIPTION
%  -----------
%     The common object for representing 1- and infinity- norm functions given as
%  f=||x||_p  where pin{1,oo}. The one norm is given as a sum of absolute values 
%                                         n       
%                                        --       
%                                        \        
%                                   f  = /   |x | 
%                                    1   --    i  
%                                        i=1      
%     and the infinity norm is given as the maximum over the absolute values 
%                                                      
%                             f   = max(|x |,...,|x |) 
%                              oo         1        n   
%     where n  is the dimension of the vector x. If the weighing matrix Q is
%  provided, then the product f=||Qx||_p  is considered. The weight Qdoes not need
%  to be square. Function value is always scalar.
%    2-norms are not supported because they are neither quadratic, nor piecewise
%  linear.
%    Do not use these objects in the user interface. Use OneNormFunction and
%  InfNormFunction objects instead.
%  
%  INPUT
%  -----
%     
%        
%          1 Flag indicating the type of the norm. It 
%            can be either 1 or Inf.                  
%            Class: double                            
%            Allowed values:                          
%                                                     
%              1                                      
%              Inf                                    
%                                                     
%          Q Weighing matrix where the number of      
%            columns determines the dimension of the  
%            vector argument.                         
%            Class: double                            
%              
%  
%  
%  OUTPUT
%  ------
%     
%        
%          f The NormFunction object.                 
%            Class: NormFunction                      
%              
%  
%  
%  SEE ALSO
%  --------
%     OneNormFunction,  InfNormFunction,  AffFunction,  QuadFunction
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
 
 
	%
	% class for representing norm-based functions
	%
	% syntax:
	%   f = NormFunction(1)      : f = norm(x, 1)
	%   f = NormFunction(Inf)    : f = norm(x, Inf)
	%   f = NormFunction(1, Q)   : f = norm(Q*x, 1)
	%   f = NormFunction(Inf, Q) : f = norm(Q*x, Inf)
	%
	% "Q" need not to be square. Function value is always scalar.
	%
	% 2-norms are not supported because they are neither quadratic, nor
	% piecewise linear.
	
	properties(SetAccess=private)
		weight=1; % weight (1 by default)
		type=1; % either 1 or Inf
		D=0; % dimension of the domain
		R=1; % dimension of the range, norms are always scalar-valued
	end
	
	methods
		
		% Constructor
		function obj = NormFunction(type, Q)
			%
			% syntax:
			%   f = NormFunction(1)      : f = norm(x, 1)
			%   f = NormFunction(Inf)    : f = norm(x, Inf)
			%   f = NormFunction(1, Q)   : f = norm(Q*x, 1)
			%   f = NormFunction(Inf, Q) : f = norm(Q*x, Inf)
			%
			% "Q" need not to be square. Function value is always scalar.
			%
			% 2-norms are not supported because they are neither quadratic,
			% nor piecewise linear.
			%
			% Do not use these objects in the user interface. Use
			% OneNormFunction and InfNormFunction objects instead.
			
			if nargin==0
				% when called from derived classes
				return
			end
			
			narginchk(1, 2);
			
			% validation of arguments is done in setters
			obj.type = type;
			if nargin==2
				obj.weight = Q;
			end
			
			obj.Handle = @(x) norm(obj.weight*x, obj.type);
		end
		
		function obj = set.weight(obj, Q)
			% obj.Q setter
			
			if isempty(Q)
				% restore norm(x, type)
				obj.weight = 1;
				obj.D = 0;
			elseif isscalar(Q)
				% restore unrestricted domain
				validate_realmatrix(Q);
				obj.weight = Q;
				obj.D = 0;
			else
				% domain is equal to number of columns of Q
				validate_realmatrix(Q);
				if issparse(Q)
					Q = full(Q);
				end
				obj.weight = Q;
				obj.D = size(Q, 2);
			end
		end
		
		function obj = set.type(obj, type)
			% obj.type setter
			
			if ~isnumeric(type) || ~ismember(type, [1, Inf])
				error('Norm type can only be either 1 or Inf.');
			end
			obj.type = type;
		end
		
		function display(obj)
			
			% TODO: inherit from Function, see how it's done in
			% AbstractController
            if numel(obj)>1
				fprintf('Array of %d norm functions\n',numel(obj));
				return
            end
			
            if numel(obj)==0
                disp('Empty function');
                return
            end

			if obj.D==0
				fprintf('%s-norm function\n', num2str(obj.type));
			else
				fprintf('%s-norm function in R^%i\n', num2str(obj.type), obj.D);
			end
		end
		
		function status = eq(f, g)
			% Returns true if the two functions are identical
			
			% TODO: move common code to Function/eq
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
				
				status = isa(f, 'NormFunction') && ...
					isa(g, 'NormFunction') && ...
					f.D==g.D && ...
					f.type==g.type && ...
					isequal(f.weight, g.weight);
			end
		end
		
		function status = ne(f, g)
			% Returns true if two functions are not identical
			
			% TODO: inherit from Function
			status = ~eq(f, g);
		end
		
	end
	
end
