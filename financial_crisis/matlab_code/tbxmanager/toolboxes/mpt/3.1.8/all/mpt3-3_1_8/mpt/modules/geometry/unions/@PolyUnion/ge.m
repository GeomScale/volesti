function status = ge(U1, U2)
%
%  GE: Test if a union of polyhedra contains another union. 
%  =========================================================
%  
%  
%  SYNTAX
%  ------
%     
%      U1 >= U2
%    
%  
%  DESCRIPTION
%  -----------
%     Check if the union of polyhedra U1 is a non-strict superset of another union
%  U2. The result it the logical statement if U1 >= U2 and false otherwise.
%  
%  INPUT
%  -----
%     
%        
%          U1 Union of polyhedra in the same           
%             dimension.                               
%             Class: PolyUnion                         
%          U2 Union of polyhedra in the same           
%             dimension.                               
%             Class: PolyUnion                         
%               
%  
%  
%  OUTPUT
%  ------
%     
%        
%          tf True if U1 >= U2 and false otherwise.    
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
%     le
%  

%  AUTHOR(s)
%  ---------
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
 
 
status = (U2 <= U1);

end
