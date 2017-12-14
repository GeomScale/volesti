classdef OneNormFunction < NormFunction
%
%  ONENORMFUNCTION: Class representing 1-norm function. 
%  =====================================================
%  
%  
%  SYNTAX
%  ------
%     
%      f = OneNormFunction(Q)
%    
%  
%  DESCRIPTION
%  -----------
%     The object for representing 1-norm function given as f=||Qx||_1. The function
%  is given as a sum of absolute values of the product y=Qx. 
%                                        n       
%                                       --       
%                                       \        
%                                   f = /   |y | 
%                                       --    i  
%                                       i=1      
%     where n  is the dimension of the vector y. The weight Qdoes not need to be
%  square. Function value is always scalar.
%  
%  
%  INPUT
%  -----
%     
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
%          f The OneNormFunction object.              
%            Class: OneNormFunction                   
%              
%  
%  
%  SEE ALSO
%  --------
%     InfNormFunction,  AffFunction,  QuadFunction
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
	% represents weighted 1-norm functions
	%
	% syntax:
	%   f = OneNormFunction(Q) : f = norm(Q*x, 1)
	%
	% "Q" need not to be square. Function value is always scalar.
	
	methods
		
		% Constructor
		function obj = OneNormFunction(Q)
			% Constructs a weighted 1-norm function object
			%
			% syntax:
			%   f = OneNormFunction(Q) : f = norm(Q*x, 1)
			%
			% "Q" need not to be square. Function value is always scalar.
			
			narginchk(1, 1);
			obj = obj@NormFunction(1, Q);
		end
		
	end
	
end
