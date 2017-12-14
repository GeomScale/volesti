function approx = outerApprox(obj)
%
%  OUTERAPPROX: Computes outer bounding box of the convex set. 
%  ============================================================
%  
%  
%  SYNTAX
%  ------
%     
%      B = outerApprox(S)
%      B = S.outerApprox
%    
%  
%  DESCRIPTION
%  -----------
%     Compute the smallest axis-aligned hypercube that contains this set. The lower
%  and upper bounds of the hypercube are stored under Internal property, i.e.
%  Internal.lb for lower bound and Internal.ub for upper bound.
%  
%  INPUT
%  -----
%     
%        
%          S Any set derived from ConvexSet class,    
%            e.g. YSet or Polyhedron.                 
%            Class: ConvexSet                         
%              
%  
%  
%  OUTPUT
%  ------
%     
%        
%          B Bounding box B  described as Polyhedron  
%            in H-representation.                     
%            Class: Polyhedron                        
%              
%  
%  
%  SEE ALSO
%  --------
%     support,  separate
%  

%  AUTHOR(s)
%  ---------
%     
%    
%   (c) 2010-2013  Colin Neil Jones: EPF Lausanne
%   mailto:colin.jones@epfl.ch 
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
 
 
global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

if obj.isEmptySet
    approx = Polyhedron.emptySet(obj.Dim);
    approx.Internal.lb =  inf*ones(obj.Dim,1);
    approx.Internal.ub = -inf*ones(obj.Dim,1);
    return
end

% deal with arrays
no = numel(obj);
if no>1
    approx(size(obj)) = Polyhedron;
    for i=1:no
        approx(i) = obj(i).outerApprox;
    end
    return
end


if nargin < 2
    % Compute the support in all +- elementary directions.
    I   = eye(obj.Dim);
    pos = zeros(obj.Dim,1);
    neg = zeros(obj.Dim,1);
    for i=1:obj.Dim
        pos(i) =  obj.support( I(i,:)');
        neg(i) = -obj.support(-I(i,:)');
    end
    
    H = [eye(obj.Dim), pos; -eye(obj.Dim), -neg];
    
    approx = Polyhedron(H(:, 1:end-1), H(:, end));
    approx.Internal.lb = neg;
    approx.Internal.ub = pos;
else
    error('Polyhedral outer approximation not yet implemented');
end

% Input case 1: maxFacets, maxVertices, maxErr
% Input case 2: none


end
