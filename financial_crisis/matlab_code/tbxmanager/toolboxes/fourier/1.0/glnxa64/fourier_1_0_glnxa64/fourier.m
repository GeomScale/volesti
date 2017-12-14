function fourier
%
%  FOURIER: Fourier-Motzkin elimination algorithm. 
%  ================================================
%  
%  
%  SYNTAX
%  ------
%     
%      h = fourier(H,dim)
%      h = fourier(H,dim,tol,qtol)
%    
%  
%  DESCRIPTION
%  -----------
%     Implementation of the Fourier-Motzkin elimination algorithm in C. The
%  algorithm processes the system of linear inequalities over given dimensions to a
%  reduced system with the same solution set. This operation can be viewed as
%  projection on the specified dimensions.
%  
%  INPUT
%  -----
%     
%    
%      H    H-representation of a polyhedron, i.e.   
%           H=(A,b)  where Ax<=b  describe the       
%           inequalities.                            
%           Class: double                            
%      dim  Dimensions over which to project given   
%           as integers.                             
%           Class: double                            
%      tol  Numerical tolerance before two numbers   
%           are considered equal.                    
%           Class: double                            
%           Default: 1e-6                            
%      qtol Angle in degrees between two vectors     
%           before they are considered equal.        
%           Class: double                            
%           Default: 1e-2                            
%             
%  
%  
%  OUTPUT
%  ------
%     
%    
%      h New system of inequalities defined as    
%        h=(F,g)  where Fx<=g.                    
%        Class: double                            
%          
%  
%  
%  EXAMPLE(s)
%  ----------
%  
%  
%  Example 1
%  =========
%    % We have a polyhedron described in 3D with the following inequalities:
%       H = [-1 -1 -1 0; 3 -1 -1 1; -1 3 -1 2; -1 -1 3 3]; 
%       P = Polyhedron('H',H); 
%    % Plot the polyhedron
%       P.plot 
%    % Reduce the system of inequalies H by eliminating the first variable.
%    % This means that we want to perform the elimination over dimensions [2, 3].
%       h = fourier(H,[2 3]) 
%    % The reduced system has one variable left and is given in 2D.
%    % We can plot both sets to see the projection.
%       Q = Polyhedron('H',h); 
%       plot([P, Q]) 
%    
%  
%  SEE ALSO
%  --------
%     Polyhedron
%  

%  AUTHOR(s)
%  ---------
%     
%    
%   (c) 2010-2012  Colin Neil Jones: EPFL Lausanne
%   mailto:colin.jones@epfl.ch 
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
 
 
end