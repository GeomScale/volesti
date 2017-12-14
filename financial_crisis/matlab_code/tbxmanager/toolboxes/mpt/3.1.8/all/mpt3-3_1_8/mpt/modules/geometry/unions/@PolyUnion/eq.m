function ts = eq(U1,U2)
%
%  EQ: Returns true if the set covered by unions of polyhedra U_1  is the same as
%  ==============================================================================
%  the set covered by union U_2  and false otherwise. 
%  ===================================================
%  
%  
%  SYNTAX
%  ------
%     
%      ts = U1.eq(U2)
%      ts = U1 == U2
%    
%  
%  DESCRIPTION
%  -----------
%     Returns true if the union of polyhedra in the same dimension U1 covers the
%  same space as the union U2 which is in the same dimension. The result is a
%  logical statement if U_1 = U_2  holds true or not. The polyunion U1 can be an
%  array, whereas U2 is restricted to be a single polyunion object.
%  
%  INPUT
%  -----
%     
%        
%          U1 Union of polyhedra in the same           
%             dimension.                               
%             Class: PolyUnion                         
%          U2 Union of polyhedra in the same dimension 
%             as U1.                                   
%             Class: PolyUnion                         
%               
%  
%  
%  OUTPUT
%  ------
%     
%        
%          ts True if U1 == U2 and false otherwise.    
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
%     copy,  join
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
 
 
if ~isa(U2,'PolyUnion')
    error('The input argument must be "PolyUnion" object.');
end

% both polyhedra are empty arrays
if numel(U1)==0 && numel(U2)==0
    ts = true;
    return;
end
% of the polyhedra is empty array
if numel(U1)==0 || numel(U2)==0
    ts = false;
    return;
end

if numel(U2)~=1
    error('Only single "PolyUnion" object can be tested for equality.');
end

% deal with arrays
nU = numel(U1);
if numel(U1)>1
   ts = false(size(U1));
   for i=1:nU
        ts(i) = U1(i).eq(U2);
   end
   return;
end


if U1.Num==0 && U2.Num==0
    ts = true;
    return;
end
if U1.Num==0 || U2.Num==0
    ts = false;
    return
end

% check dimensions
dims = cell(nU, 1);
dims{:} = U1.Dim;

if any(diff([dims{:}]))
    error('The array of polyunions must be in the same dimension.');
end

if U1(1).Dim~=U2.Dim
    error('Unions must have the same dimension.');
end

% use the comparison based on set difference for arrays of polyhedra
ts = (U1.Set == U2.Set);


end
