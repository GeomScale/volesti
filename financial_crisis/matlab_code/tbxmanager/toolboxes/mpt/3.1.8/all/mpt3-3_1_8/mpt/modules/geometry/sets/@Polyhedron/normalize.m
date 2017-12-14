function obj = normalize(obj)
%
%  NORMALIZE: Normalizes polyhedron in H-representation. 
%  ======================================================
%  
%  
%  SYNTAX
%  ------
%     
%      P.normalize
%      normalize(P)
%    
%  
%  DESCRIPTION
%  -----------
%     Normalize polyhedron by scaling each inequality a_ix<= b  and equality a_e,i
%  = b_e,i  such that each row has unitary norm, i.e. ||a_i||_2 = 1  and
%  ||a_e,i||_2 = 1 . The normalization routine is very useful to prevent from
%  numerical problems involved in badly scaled polyhedra.
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
    MPTOPTIONS=mptopt;
end

% deal with arrays
if numel(obj)>1
    for i=1:numel(obj)
        obj(i).normalize;
    end
    return
end


if obj.hasHRep
    
    A = obj.A;
    b = obj.b;

    % 2-norm of each facet
    n = matNorm(A);

	% normalize 0'*x<=+/-b to 0'*x<= +/- sign(b)
	% (correct sign is important as not to mess with trivial infeasibility)
	ZeroRows = (n<MPTOPTIONS.zero_tol);
	n(ZeroRows) = 1;
	b(ZeroRows) = sign(b(ZeroRows));

	% normalize each half-space (0*x<=b will be left intact)
	nA = A ./ repmat(n,1,size(A,2));
	nb = b ./ n;
	obj.H_int=[nA, nb];

    % normalize also equalities if present
    if ~isempty(obj.He_int)
        Ae = obj.Ae;
        be = obj.be;
        
        % 2-norm of each equality
        ne = matNorm(Ae);
        
        % if polyhedron was properly constructed, n should not contain zero
        % values
        if any(abs(ne)<MPTOPTIONS.zero_tol)
            nAe=Ae;
            nbe=be;
        else
            nAe = Ae ./ repmat(ne,1,size(Ae,2));
            nbe = be ./ repmat(ne,1,size(be,2));
        end
        obj.He_int=[nAe, nbe];
    end

end

end
