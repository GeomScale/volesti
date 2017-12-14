function Pdiff = minus(P, S)
%
%  MINUS: Subtract a Polyhedron or a vector from a Polyhedron. 
%  ============================================================
%  
%  
%  SYNTAX
%  ------
%     
%      Q = P - S
%      Q = minus(P,S)
%      Q = P - x
%      Q = minus(P, x)
%    
%  
%  DESCRIPTION
%  -----------
%     Compute the Pontryagin difference of P  and S, or P  and x. 
%                    Q=P-S={ x in P  |  x+w in P forall w in S }
%     or 
%                              P-x = { y-x  |  y in P}
%  
%  
%  INPUT
%  -----
%     
%        
%          P Polyhedron in any format                 
%            Class: Polyhedron                        
%          S Polyhedron in any format                 
%            Class: Polyhedron                        
%          x Column vector of length P.Dim            
%            Class: Polyhedron                        
%              
%  
%  
%  OUTPUT
%  ------
%     
%        
%          Q Polyhedron Q = P-S or Q = P-x.           
%            Class: Polyhedron                        
%              
%  
%  
%  SEE ALSO
%  --------
%     plus
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
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end


% deal with arrays
if numel(P)>1
    Pdiff(size(P)) = Polyhedron;
    for i=1:numel(P)
        Pdiff(i) = minus(P(i),S);
    end
    return;
end

type = class(S);

switch type
    case 'Polyhedron'
        
        if numel(S)>1
            error('Only one polyhedron S is allowed.');
        end
                
        % if P is empty array
        if builtin('isempty',P)
            Pdiff = Polyhedron;
            return
        end
        
        if isEmptySet(S)
            Pdiff = P;
            return;
        end
        if P.Dim ~= S.Dim,
            error('P and S must have the same dimension'); 
        end
        
        % P-S = {x in P | x+w in P forall w in S}
        %%% TODO : Figure out how to compute Pontyargin diffs for V-rep
        P.minHRep();
        S.minHRep();
        
        % The affine hull of S must be a subset of the affine hull of P, or the
        % difference is empty.
        if rank([P.He;S.He], MPTOPTIONS.abs_tol) > rank(S.He, MPTOPTIONS.abs_tol)
            % empty polyhedron in the same dimension
			Pdiff = Polyhedron.emptySet(P.Dim);
            return
        end
        
        % special case P==S
        if P==S
            % empty polyhedron in the same dimension
            Pdiff = Polyhedron.emptySet(P.Dim);
            return;
        end
        
        % Subtract the support of S for each inequality of P
        Hn = [P.A P.b-S.support(P.A')];
        if any(Hn(:,end)==-Inf)
            % empty polyhedron in the same dimension
            Pdiff = Polyhedron.emptySet(P.Dim);
            return;
        end
        % remove remaining Inf rows
        Hn((Hn(:,end))==Inf,:) = [];
        Pdiff = Polyhedron('H', [P.H; Hn], 'He', P.He);
                       
    otherwise
        
        x = S(:);
        validate_realvector(x);
        if length(x) ~= P.Dim,
            error('Length of the vector must be P.Dim = %i.', P.Dim);
        end
        
        % Shift the polyhedron by x
        Pdiff = P.plus(-x);
end



end
