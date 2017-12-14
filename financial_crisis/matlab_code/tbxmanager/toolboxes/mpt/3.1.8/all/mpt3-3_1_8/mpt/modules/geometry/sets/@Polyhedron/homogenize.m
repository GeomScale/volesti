function Pnew = homogenize(P, type)
%
%  HOMOGENIZE: Compute the homogenization of the given Polyhedron. 
%  ================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      H = P.homogenize(type)
%      H = homogenize(P, type)
%    
%  
%  DESCRIPTION
%  -----------
%     Compute the homogenization of the given Polyhedron. Parametrize the right
%  hand side of the inequalities/equalities that describe the polyhedron to get
%  homogenized system of equations Ax<= 0  and A_ex = 0 . Given Polyhedron 
%                                    n                      
%                         P = {x in R   | Ax <= b, A x = b }
%                                                   e     e 
%     the homogenization is 
%                          n        1                                
%                H = {xin R , t in R   | Ax - bt <= 0, A x - b t = 0}
%                                                       e     e      
%     where the t  is the lifting parameter. The dimension of the polyhedron H  is
%  by one greater than the dimension of P. If type = 'Hrep' or type = 'Vrep' is
%  specified, then the homogenization is returned in this form, otherwise the
%  returned type is equal to the type of P.
%  
%  INPUT
%  -----
%     
%        
%          P    Polyhedron in any format                 
%               Class: Polyhedron                        
%          type Desired type of the returned polyhedron  
%               Class: char                              
%               Allowed values:                          
%                                                        
%                 Hrep  Hyperplane representation.       
%                 Vrep  Vertex representation.           
%                                                        
%                 
%  
%  
%  OUTPUT
%  ------
%     
%        
%          H Homogenization of P polyhedron           
%            Class: Polyhedron                        
%              
%  
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
 
 
if nargin < 2 
    type = '';
end

%% deal with arrays
if length(P)>1
    Pnew(size(P)) = Polyhedron;
    for i=1:length(P)
        Pnew(i) = P(i).homogenize(type);
    end
    return;
end

if isempty(type)
    if P.hasHRep,
        type = 'hrep';
    else
        type = 'vrep';
    end
end


% empty polyhedron
if ~P.hasHRep && ~P.hasVRep
    Pnew = Polyhedron;
    return;
end


switch lower(type)
    
    case 'hrep'
        if ~P.hasHRep,
            P.minHRep();
        end
        Pnew = Polyhedron('H', [P.A -P.b zeros(size(P.H,1),1)], 'He', [P.Ae -P.be zeros(size(P.He,1),1)]);
    case 'vrep'
        if ~P.hasVRep,
            P.minVRep();
        end
        Pnew = Polyhedron('V', zeros(1,P.Dim+1), 'R', [P.V ones(size(P.V,1),1);P.R zeros(size(P.R,1),1)]);
    otherwise
        error('Unknown type %s', type);
end

end
