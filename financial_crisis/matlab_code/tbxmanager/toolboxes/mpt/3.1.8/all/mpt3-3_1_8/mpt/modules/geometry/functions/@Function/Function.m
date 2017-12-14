classdef Function < handle & IterableBehavior
%
%  FUNCTION: Function representation for MPT 
%  ==========================================
%  
%  
%  SYNTAX
%  ------
%     
%      F = Function(@fun, data)
%      F = Function(@fun)
%      F = Function([],struct('p',2))
%    
%  
%  DESCRIPTION
%  -----------
%     The Function class represents functions handles that are associated to sets
%  in MPT. The class combines the function_handle with changeable parameters which
%  allows specifications of very general functions. When constructing the Function
%  object, the correctness of the given function is not tested. It is up to the
%  user to provide correct mapping, or test it via feval method for ConvexSet
%  class. Domain of the function is specified by the set it is associated to. If
%  the point lies outside of the domain, an error is thrown. The range of the
%  function is not known at the time of the construction of the Function object. It
%  can be determined after evaluation of the related function for given point. User
%  can associate any data with the function, including parameters of the function,
%  under the Datafield. These data can be modified anytime after the construction
%  of the object.
%  
%  INPUT
%  -----
%     
%        
%          Handle A function_handle that represents the    
%                 function. It can be an anonymous         
%                 function, e.g. f(x) = x^3  that          
%                 corresponds to syntax @(x)x^3  or it can 
%                 be a link to another file that evaluates 
%                 the expression, i.e. @file. Both of the  
%                 expressions are fine as long as the      
%                 given handle can be evaluated at the     
%                 time of the construction. If more        
%                 arguments (or parameters) are present in 
%                 the function, the arguments are          
%                 separated with a comma, e.g. f(x,y,z) =  
%                 x(y-z)^2  corresponds to @(x,y,z)        
%                 x*(y-z)^2 . For more info about          
%                 constructing handles see help for        
%                 function_handle class.                   
%                 Class: function_handle                   
%          Data   Any data related to the function. It can 
%                 be function parameters, variable bounds, 
%                 domain, range, measured data, symbolic   
%                 representation, etc. The data can be     
%                 changed after construction of the        
%                 Function object.                         
%                   
%  
%  
%  OUTPUT
%  ------
%     
%        
%          F Function object.                         
%            Class: Function                          
%              
%  
%  
%  SEE ALSO
%  --------
%     AffFunction,  QuadFunction,  setHandle
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
	% class for representing functions
	%
	% syntax: F = Function('Handle',@fun,'Data',any_data)
	%         F = Function(@fun)
	%
	%         F = Function('Data',struct('p',2));
	%         F.setHandle(@(x)*F.Data.p)
	
	
	%%
	properties(SetAccess = protected)
		Handle % function handle
		Internal % internal data
	end
	
	properties(SetAccess = public)
		Data; % additional user data
	end
	
	%%
	methods(Access = public)
		
		% Constructor
		function F = Function(Handle, Data)
			%
			% sets data for Function object
			%
			% syntax: F = Function(@fun, any_data)
			%         F = Function(@fun)
			%
			%         F = Function(struct('p',2));
			%         F.setHandle(@(x)*F.Data.p)
			%
			% for more details, type "help Function"
			
			if nargin==0
				% nothing to do, an empty object will be automatically
				% constructed
			elseif nargin==1
				if isa(Handle, 'function_handle')
					F.Handle = Handle;
				else
					F.Data = Handle;
				end
			elseif nargin==2
				F.Handle = Handle;
				F.Data = Data;
			end
		end
		
		function out = feval(obj, x)
			% evaluates the function at a point
			
			% For performance reasons we really don't want to declare "feval"
			% as a method. Instead, remove this method, and rename the
			% "Handle" property to "feval". Then "obj.feval(x)" works
			% correctly without any overhead due to method call.
			out = obj.Handle(x);
		end
		
		function new = slice(obj, dims, values)
			% Slice a function through given coordinates
			narginchk(3, 3);
			
			% A very general implementation is to restore the vector "x" by
			% x(dims)=values, x(keep)=z, where "keep" are the indices which
			% are _not_ sliced.
			%
			% This implementation supports arbitrary nonlinear functions
			% specified as function handles, but also provides slicing of
			% OneNormFunction and InfNormFunction objects.
			new = Function(@(z) obj.Handle(obj.restore_sliced_x(z, dims, values)));
		end
	end
	
	methods (Hidden)
		
		function obj = setInternal(obj,name,value)
			%
			% overwrites internal property for the Function object
			% (internal function)
			%
			% If we want to add the internal property from outside of
			% this class (e.g. inside the PLCP solver) e.g.
			%
			% obj.Internal.name = value,
			%
			% use the syntax:
			%         obj.setInternal('name',value)
			%
			% DO NOT USE THIS METHOD UNLESS YOU PERFECTLY KNOW WHAT
			% YOU ARE DOING
			%
			
			narginchk(3, 3);
			
			if ~ischar(name)
				error('Name must be a string.');
			end
			
			obj.Internal.(name) = value;
			
		end
	end
	
	methods (Static, Hidden)
		function x = restore_sliced_x(z, dims, values)
			% Restores "x" by fixing x(dims)=values and x(other)=z with
			% "other = setdiff(1:nx, dims)
			
			nx = numel(z)+numel(dims);
			keep = setdiff(1:nx, dims);
			x = zeros(nx, 1);
			x(keep) = z(:);
			x(dims) = values;
		end
	end
	
end
