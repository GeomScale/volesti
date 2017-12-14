function display(P)
%
%  DISPLAY: Overload display for Polyhedron class. 
%  ================================================
%  
%  
%  SYNTAX
%  ------
%     
%      display(P)
%      P.display
%    
%  
%  DESCRIPTION
%  -----------
%     Default display for Polyhedron class.
%  
%  INPUT
%  -----
%     
%        
%          P Polyhedron object.                       
%            Class: Polyhedron                        
%              
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
 
 
if numel(P)==0
    fprintf('Empty polyhedron array.\n');
    return;
elseif numel(P) > 1
    fprintf('Array of %i polyhedra.\n', numel(P));
    return
end

if isempty(P.H_int) && isempty(P.V_int) && isempty(P.R_int)
    if size(P.He_int,1) == 0
        fprintf('Empty polyhedron in R^%i\n', P.Dim);
    else
        fprintf('Affine set with %i equations in R^%i\n', size(P.He,1), P.Dim);
    end
    return
end
fprintf('Polyhedron in R^%i with representations:\n', P.Dim);
fprintf('    H-rep ');
if P.hasHRep
    if P.irredundantHRep, fprintf('%-13s', '(irredundant)');
    else                  fprintf('%-13s', '(redundant)');
    end
    fprintf(' : Inequalities %3i | Equalities %3i\n', size(P.H,1), size(P.He,1));
else
    fprintf('%-13s : Unknown (call computeHRep() to compute)\n', '');
end
fprintf('    V-rep ');
if P.hasVRep
    if P.irredundantVRep, fprintf('%-13s', '(irredundant)');
    else                  fprintf('%-13s', '(redundant)');
    end
    fprintf(' : Vertices %3i | Rays %3i\n', size(P.V,1), size(P.R,1));
else
    fprintf('%-13s : Unknown (call computeVRep() to compute)\n', '');
end

% display attached functions (implemented in ConvexSet/displayFunctions)
P.displayFunctions;

end
