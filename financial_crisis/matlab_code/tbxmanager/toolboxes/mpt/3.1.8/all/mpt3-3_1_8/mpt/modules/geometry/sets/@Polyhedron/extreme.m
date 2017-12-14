function sol = extreme(obj, y)
%
%  EXTREME: Compute extremal point of a polyhedron in a given direction. 
%  ======================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      sol = P.extreme(y)
%    
%  
%  DESCRIPTION
%  -----------
%     P.extreme(y) solves the following problem: 
%                                                        
%                            J(y) = max   {y'x |  x in P}
%                                    x                   
%     and returns the optimizer / extreme point x.
%  
%  INPUT
%  -----
%     
%        
%          P Polyhedron in any format.                
%            Class: Polyhedron                        
%          y Direction to compute the extreme point.  
%            Class: double                            
%              
%  
%  
%  OUTPUT
%  ------
%     
%        
%          sol          Support of y  in P.                      
%                       Class: struct                            
%          sol.exitflag Integer value informing about the        
%                       termination status of the optimization.  
%                       Class: double                            
%          sol.x        An optimizer of max_x   {y'x |  x in P}  
%                       Class: double                            
%          sol.supp     Optimal value of max_x   {y'x |  x in P} 
%                                                                
%                       Class: double                            
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
 
 
global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

narginchk(2, 2);

if numel(obj)>1
	% deal with arrays
	sol = obj.forEach(@(x) x.extreme(y));

else
	% compute the extreme point
    validate_realvector(y);
    y=y(:);
    if length(y) ~= obj.Dim,
        error('The vector x must be of length %i.', obj.Dim);
    end
    
    lp   = obj.optMat;
    lp.f = obj.buildCost(-y(:)).f;
	lp.quicklp = true;
    opt = mpt_solve(lp);
    sol.exitflag    = opt.exitflag;
    sol.x       = [];
    sol.supp    = [];
    switch sol.exitflag
        case MPTOPTIONS.UNBOUNDED
            sol.supp = inf;
        case MPTOPTIONS.OK
            sol.x       = opt.xopt(1:obj.Dim);
            sol.supp    = -opt.obj;
    end
end

end
