function sol = shoot(obj, r, x0)
%
%  SHOOT: Maximize along a given ray within the polyhedron. 
%  =========================================================
%  
%  
%  SYNTAX
%  ------
%     
%      s = P.shoot(r, x0)
%      s = shoot(P, r, x0)
%    
%  
%  DESCRIPTION
%  -----------
%     Computes the solution to the optimization problem: 
%                                                            
%                        max    { alpha | x  + alpha r in P }
%                       alpha              0                 
%     The problem can be explained as to find the point that lies along the line
%  given by the ray rwhich starts from the point x0 such that it still lies inside
%  the set P. Typically, the solution of this set is the point that lies on the
%  boundary of the set, or can be Inf if the set is unbounded.
%  
%  INPUT
%  -----
%     
%        
%          P   Polyhedron in any format                 
%              Class: Polyhedron                        
%          dir Vector of size P.Dim                     
%              Class: double                            
%          x0  Vector of size P.Dim                     
%              Class: double                            
%              Default: Vector of all zeros             
%                
%  
%  
%  OUTPUT
%  ------
%     
%        
%          ret          Optimal solution                         
%                       Class: struct                            
%          ret.alpha    Maximum length of ray, [] if P  is       
%                       empty, Inf if unbounded.                 
%          ret.exitflag Integer value informing about the        
%                       termination status of the above          
%                       optimization problem.                    
%                       Class: double                            
%                       Allowed values:                          
%                                                                
%                         mptopt.OK                              
%                         mptopt.INFEASIBLE                      
%                         mptopt.UNBOUNDED                       
%                         mptopt.ERR                             
%                                                                
%          ret.x        Extreme point of ray : alpha r + x_0     
%                         
%  
%  
%  SEE ALSO
%  --------
%     extreme
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

validate_realvector(r);
r=r(:);

% check dimension
D = [obj.Dim];
if any(D(1)~=D)
    error('The polyhedron array must be in the same dimension.');
end

if nargin<3
    x0 = zeros(obj(1).Dim,1);
else
    validate_realvector(x0);
end

% deal with arrays
if numel(obj)>1
	% return an array of structures
	sol = obj.forEach(@(elem) elem.shoot(r,x0));
    return
end

% prealocate output
sol = struct('exitflag', [], 'x', [], 'alpha', []);

% check sizes of vectors
if length(r) ~= obj.Dim,
    error('The ray vector "r"  must have the length of %i.', obj.Dim);
end
if length(x0) ~= obj.Dim,
    error('The vector "x0" must have the length of %i.', obj.Dim);
end


if obj.Dim>0
    % Re-arrange pre-computed set description matrices to:
    %  A*blkdiag(r,I) [alpha;lam] <= b - A*[x0;0]
    S.A  = [obj.optMat.A(:,1:obj.Dim)*r obj.optMat.A(:,obj.Dim+1:end)];
    S.b  = obj.optMat.b-obj.optMat.A(:,1:obj.Dim)*x0;
    S.Ae = [obj.optMat.Ae(:,1:obj.Dim)*r obj.optMat.Ae(:,obj.Dim+1:end)];
    S.be = obj.optMat.be-obj.optMat.Ae(:,1:obj.Dim)*x0;
    
    S.lb = [-inf;obj.optMat.lb(obj.Dim+1:end)];
    S.ub = [ inf;obj.optMat.ub(obj.Dim+1:end)];
    S.f  = [-1;zeros(size(S.Ae,2)-1,1)];
    
    res = mpt_solve(S);
        
else
    res.exitflag = MPTOPTIONS.ERROR;
end

sol.exitflag = res.exitflag;
sol.x = [];
sol.alpha = Inf;

switch sol.exitflag
    case MPTOPTIONS.OK,
        sol.alpha = -res.obj;
        sol.x     = sol.alpha*r + x0;
%         if sol.alpha>=MPTOPTIONS.infbound
%             sol.alpha = Inf;
%         end
%     case MPTOPTIONS.UNBOUNDED
%         sol.alpha = Inf;
%     otherwise
%         sol.alpha = [];
%         sol.x = [];
end

end
