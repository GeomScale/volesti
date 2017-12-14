function Pnew = mtimes(P, S)
%
%  MTIMES: Multiply two polyhedra, or a polyhedron with a matrix or scalar. 
%  =========================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      Q = P*S
%      Q = M*P
%      Q = a*P
%      Q = P*a
%    
%  
%  DESCRIPTION
%  -----------
%     
%    
%     1. Q = P*S where P and S are Polyhedra. Computes the product of the two
%     polyhedra: 
%                           Q = { (x,y)  |  x in P, y in S }
%   
%   
%     2. Q = M*P where P is a Polyhedron and M a matrix Computes the affine mapping
%     P.affineMap(M)  
%     3. Q = a*P or P*a where P is a Polyhedron and a is a scalar Computes the
%     scaling 
%                                  Q = {ax | x in P }
%    
%    Note that only single polyhedron S can be provided to compute a product.
%  
%  INPUT
%  -----
%     
%        
%          P Polyhedron in any format                 
%            Class: Polyhedron                        
%          S Polyhedron in any format                 
%            Class: Polyhedron                        
%          M Matrix of size n X P.Dim                 
%            Class: double matrix                     
%          a Scaling factor                           
%            Class: double                            
%              
%  
%  
%  OUTPUT
%  ------
%     
%        
%          Q Polyhedron: Either P*S, M*P  or a*P      
%            depending on the input                   
%            Class: Polyhedron                        
%              
%  
%  
%  SEE ALSO
%  --------
%     affineMap
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
 
 
global MPTOPTIONS

% Determine the type of the arguments
if isa(S, 'Polyhedron')
    if numel(S)>1
        error('Only one polyhedron "S" can be provided.');
    end
    S_class = 'Polyhedron';
elseif isnumeric(S) && isscalar(S)
    validate_realvector(S);
    S_class = 'scalar';
elseif isnumeric(S)
    validate_realmatrix(S);
    S_class = 'matrix';
else
    error('This type of the object "%s" is not supported as argument for "S".', class(S));
end

if isa(P, 'Polyhedron')
    P_class = 'Polyhedron';
elseif isnumeric(P) && isscalar(P)
    validate_realvector(P);
    P_class = 'scalar';
elseif isnumeric(P)
    validate_realmatrix(P);
    P_class = 'matrix';
else
    error('This type of the object "%s" is not supported as argument for "P".', class(P));
end

% deal with arrays if P is polyhedron
if strcmpi(P_class,'Polyhedron') && numel(P)>1
    Pnew(size(P)) = Polyhedron;
    for i=1:numel(P)
        Pnew(i) = mtimes(P(i),S);
    end
    return;
end

switch [P_class S_class]
    case ['Polyhedron' 'Polyhedron']
        if S.isEmptySet,
            Pnew = Polyhedron(P);
            return;
        end
        if P.isEmptySet,
            Pnew = Polyhedron(S);
            return;
        end
        
        % Compute the cartesian product of P x S
        % For now - just compute the convex hull...
        P.minHRep(); S.minHRep();
		if isempty(P.He_int) && isempty(S.He_int)
			Pnew = Polyhedron(blkdiag(P.A, S.A), [P.b;S.b]);
		else
			Pnew = Polyhedron('H', [blkdiag(P.A, S.A) [P.b;S.b]], 'He', [blkdiag(P.Ae, S.Ae) [P.be;S.be]]);
		end
        %%% TODO : Compute polyhedral product without taking convex hull
        
    case {['Polyhedron' 'scalar'], ['scalar' 'Polyhedron']}
        if isnumeric(P),
            alpha = P;
            poly = S;
        else
            alpha = S;
            poly = P;
		end
		
		if abs(alpha) < MPTOPTIONS.abs_tol
			% scaling with zero produces a singleton
			%
			% we deal with this explicitly as to correctly support R^n
			% (depends on resolution of issue #93)

			% obey representation of the input
			if poly.hasHRep
				Pnew = Polyhedron([eye(poly.Dim); -eye(poly.Dim)], ...
					zeros(2*poly.Dim, 1));
			else
				Pnew = Polyhedron(zeros(1, poly.Dim));
			end
			
		elseif poly.hasHRep
			if isempty(poly.He_int)
				Pnew = Polyhedron(poly.A, alpha*poly.b);
			else
				Pnew = Polyhedron('H',[poly.A alpha*poly.b],'He', poly.He);
			end
        else
            Pnew = Polyhedron('V',alpha*poly.V, 'R', poly.R);
        end
    case ['Polyhedron' 'matrix']
        Pnew = P.invAffineMap(S);
    case ['matrix' 'Polyhedron']
        Pnew = S.affineMap(P);
end

end
