function tf = ne(P, S)
%
%  NEQ: Test if a polyhedron is not equal to another polyhedron. 
%  ==============================================================
%  
%  
%  SYNTAX
%  ------
%     
%      S ~= P
%    
%  
%  DESCRIPTION
%  -----------
%     Check if the polyhedron S is not equal to polyhedron P. The result it the
%  logical statement if P neq S  and false otherwise.
%  
%  INPUT
%  -----
%     
%        
%          P Polyhedron in any format                 
%            Class: Polyhedron                        
%          S Polyhedron in any format with the same   
%            dimension as P.                          
%            Class: Polyhedron                        
%              
%  
%  
%  OUTPUT
%  ------
%     
%        
%          tf True if P neq S  and false otherwise.    
%             Class: logical                           
%             Allowed values:                          
%                                                      
%               true                                   
%               false                                  
%                                                      
%               
%  
%  
%  SEE ALSO
%  --------
%     eq,  contains,  le,  lt,  gt,  ge
%  

%  AUTHOR(s)
%  ---------
%     
%    
%   (c) 2010-2013  Colin Neil Jones: EPF Lausanne
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
%  -------------------------------------------------------------------------------
%    
%      This document was translated from LaTeX by HeVeA (1).
%  ---------------------------------------
%    
%    
%   (1) http://hevea.inria.fr/index.html
 
 
tf = ~eq(P,S);

end
