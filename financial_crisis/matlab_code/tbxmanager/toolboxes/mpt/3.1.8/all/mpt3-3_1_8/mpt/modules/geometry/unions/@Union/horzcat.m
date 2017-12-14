function y=horzcat(varargin)
%
%  HORZCAT: Horizontal concatenation for union objecs. 
%  ====================================================
%  
%  
%  SYNTAX
%  ------
%     
%      U = [U1, U2]
%      S = horzcat(U1,U2)
%    
%  
%  DESCRIPTION
%  -----------
%     Overloaded method for horizontal concatenation of unions. It is not possible
%  to concatenate objects of different type to the same array (e.g. PolyUnion and
%  Union). Similarly, it is not possible to concatenate into matrices, only vectors
%  are allowed.
%  
%  INPUT
%  -----
%     
%        
%          U1 Object of the Union class.               
%             Class: Union                             
%          U2 Object of the Union class.               
%             Class: Union                             
%               
%  
%  
%  OUTPUT
%  ------
%     
%        
%          U The array of the Union objects.          
%            Class: Union                             
%              
%  
%  
%  SEE ALSO
%  --------
%     vertcat
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
 
 
y = vertcat(varargin{:});

end
