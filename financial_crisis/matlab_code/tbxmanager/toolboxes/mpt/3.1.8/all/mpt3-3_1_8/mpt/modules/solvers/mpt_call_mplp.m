function ret = mpt_call_mplp(S)
%
%  MPT_CALL_MPLP: A gateway function to MPLP solver (without errorchecks) 
%  =======================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      R = mpt_call_mplp(S)
%    
%  
%  DESCRIPTION
%  -----------
%     The function call to MPLP solver from Opt class. Note that this solver is not
%  capable of solving MPLP with the parameterized cost function, i.e. if there is
%  non-zero pF term. Using option settings for MPLP solver taken from MPT2.6.
%  
%  INPUT
%  -----
%     
%        
%          S              Object of the Opt class                  
%                         Class: Opt                               
%          S.H            Quadratic part of the objective          
%                         function.                                
%                         Class: double                            
%                         Default: 0                               
%          S.f            Linear part of the objective function.   
%                         Class: double                            
%          S.pF           Linear part of the objective function    
%                         for parameters.                          
%                         Class: double                            
%                         Default: 0                               
%          S.A            Linear part of the inequality            
%                         constraints Ax <= b + Btheta.            
%                         Class: double                            
%          S.b            Right hand side of the inequality        
%                         constraints Ax <= b + Btheta.            
%                         Class: double                            
%          S.pB           Right hand side of the inequality        
%                         constraints for parameters Ax <= b +     
%                         Btheta.                                  
%                         Class: double                            
%          S.Ae           Linear part of the equality constraints  
%                         A_ex=b_e + Etheta .                      
%                         Class: double                            
%                         Default: []                              
%          S.be           Right hand side of the equality          
%                         constraints A_ex=b_e + Etheta .          
%                         Class: double                            
%                         Default: []                              
%          S.pE           Right hand side of the equality          
%                         constraints for parameters A_ex=b_e +    
%                         Etheta .                                 
%                         Class: double                            
%                         Default: []                              
%          S.lb           Lower bound for the decision variables x 
%                         >= lb.                                   
%                         Class: double                            
%                         Default: []                              
%          S.ub           Upper bound for the decision variables x 
%                         <= ub.                                   
%                         Class: double                            
%                         Default: []                              
%          S.Ath          Linear part of the inequality            
%                         constraints A_thetatheta <= b_theta.     
%                         Class: double                            
%                         Default: []                              
%          S.bth          Right hand side of the inequality        
%                         constraints A_thetatheta <= b_theta.     
%                         Class: double                            
%                         Default: []                              
%          S.M            Linear matrix involved in LCP.           
%                         Class: double                            
%                         Default: []                              
%          S.q            Right hand side vector involved in LCP.  
%                         Class: double                            
%                         Default: []                              
%          S.Q            Linear matrix involved in parametric     
%                         formulation of LCP.                      
%                         Class: double                            
%                         Default: []                              
%          S.n            Number of decision variables.            
%                         Class: double                            
%          S.d            Number of parameters.                    
%                         Class: double                            
%          S.m            Number of inequalities in Ax <= b +      
%                         Btheta.                                  
%                         Class: double                            
%          S.me           Number of equalities in A_ex=b_e +       
%                         Etheta.                                  
%                         Class: double                            
%          S.problem_type A string specifying the problem to be    
%                         solved                                   
%                         Class: char                              
%                         Default: []                              
%          S.varOrder     Order of variables if the problem was    
%                         processed by YALMIP first.               
%                         Class: double                            
%                         Default: []                              
%          S.Internal     Internal property of Opt class.          
%                         Class: struct                            
%                         Default: []                              
%          S.recover      Affine map for MPLP problems if there    
%                         were any equalities present and have     
%                         been removed by eliminateEquations       
%                         method.                                  
%                         Class: struct                            
%          S.recover.Y    Matrix of the affine map x = Yy + th(    
%                          theta                                   
%                            1    ) . The map is from the          
%                         optimization variables involed in the    
%                         reduced MPLP  y(theta) |->> x  to the    
%                         original MPLP.                           
%                         Class: double                            
%                         Default: []                              
%          S.recover.th   Matrix of the affine map x = Yy + th(    
%                          theta                                   
%                            1    ) . The map is from the          
%                         optimization variables involed in the    
%                         reduced MPLP  y(theta) |->> x  to the    
%                         original MPLP.                           
%                         Class: double                            
%                         Default: []                              
%                         Default: []                              
%                           
%  
%  
%  OUTPUT
%  ------
%     
%        
%          R                           result structure                         
%                                      Class: struct                            
%          R.xopt                      Optimal solution                         
%                                      Class: PolyUnion                         
%          R.mplpsol                   Structure with the solution as returned  
%                                      by MPLP solver.                          
%                                      Class: struct                            
%          R.mplpsol.Pn                Array of polytopes in MPT2 format.       
%                                      Class: polytope                          
%          R.mplpsol.Fi                Cell array of matrices of the control    
%                                      law given as x=F_itheta + G_i.           
%                                      Class: cell                              
%          R.mplpsol.Gi                Cell array of matrices of the control    
%                                      law given as x=F_itheta + G_i.           
%                                      Class: cell                              
%          R.mplpsol.activeConstraints Index set of active constraints          
%                                      Class: cell                              
%          R.mplpsol.Phard             Feasible domain.                         
%                                      Class: polytope                          
%          R.mplpsol.details           More details about the solution and the  
%                                      computation.                             
%                                      Class: struct                            
%                                        
%  
%  
%  SEE ALSO
%  --------
%     Opt,  mpt_solvemp
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
    MPTOPTIONS = mptopt;
end

if ~isa(S,'Opt')
    error('mpt_call_mplp: Input argument must be an "Opt" class.');
end

if ~strcmpi(S.problem_type,'LP')
    error('mpt_call_mplp: MPLP solver does not solve %s problems!',S.problem_type);
end

% store the original problem
opt = S.copy;
% make a copy of the original problem before we modify it via
% eliminateEquations()
S = S.copy();

if S.me>0   
    % eliminate equations first
    S.eliminateEquations;
end

% In validation of Opt class there are prepreprocessing functions
% for tightening the bounds on the parametric solution. One of the task performs
% extraction of lower and upper bounds on the variables from equation 
% G*U <= W + E*th and puts them into separate fields. We need to put these 
% bounds back.

Matrices.G = S.A;
Matrices.W = S.b;
Matrices.E = S.pB;
ilb = (S.lb==-Inf) | (S.lb<-MPTOPTIONS.infbound);
iub = (S.ub==Inf)  | (S.ub>MPTOPTIONS.infbound);
if any(~ilb)
    % put ones at the positions where there is ub
    L = -eye(S.n);
    Matrices.G = [Matrices.G; L(~ilb,:)];
    Matrices.W = [Matrices.W; -S.lb(~ilb)];
    Matrices.E = [Matrices.E; zeros(nnz(~ilb),S.d)];
end
if any(~iub)
    % put ones at the positions where there is ub
    U = eye(S.n);
    Matrices.G = [Matrices.G; U(~iub,:)];
    Matrices.W = [Matrices.W; S.ub(~iub)];
    Matrices.E = [Matrices.E; zeros(nnz(~iub),S.d)];
end

Matrices.H = S.f;
Matrices.F = S.C;
Matrices.bndA = S.Ath;
Matrices.bndb = S.bth;
if any(S.pF(:))
    % add parameterized cost
    Matrices.D = S.pF;
end

if MPTOPTIONS.verbose >= 1
    disp('Calling mpt_mplp_26 with default options...')
end
start_time = clock;
[r.Pn,r.Fi,r.Gi,r.activeConstraints,r.Phard,r.details]=mpt_mplp_26(Matrices);

% Re-order variables if this came from YALMIP
P = speye(opt.n);
if ~isempty(opt.varOrder)
    P = P(opt.varOrder.requested_variables,:);
end

% convert @polytope to @Polyhedron
reg = toPolyhedron(r.Pn);

for i=1:length(reg)
    
    % add only not empty regions (with region_tol)
    if ~reg(i).isEmptySet
        %regions.add(reg);
        
        % y = Fi*th + Gi = [Fi Gi]*[th;1]
        % x = Y*y + T*[th;1];
        % x = (Y*[Fi Gi]+T)*[th;1]
        
        if opt.me
            % Compute affine mapping from parameter to primal        
            T = S.recover.Y*[r.Fi{i} r.Gi{i}] + S.recover.th;
        else
            T = [r.Fi{i} r.Gi{i}];
        end
        
        % compute primal variables
        Lprimal = P*T;
        reg(i).addFunction(AffFunction(Lprimal(:,1:end-1),Lprimal(:,end)),'primal');
        
        % compute the objective value
        Y = T(:,1:end-1);
        R = T(:,end);

        % add parameterized cost theta'*pF*x (issue #122)
        lt = R'*opt.pF + opt.f'*Y + opt.C;
        at = opt.f'*R + opt.c;
        reg(i).addFunction(AffFunction(lt,at),'obj');

    end
end

ret.xopt = PolyUnion('Set',reg,'Domain', toPolyhedron(r.Phard),...
	'Convex',true,'Overlaps',false,'Bounded',true,'Fulldim',true,'Connected',true);
ret.xopt.setInternal('convexHull', toPolyhedron(r.Phard));
ret.mplpsol = r;
if ret.xopt.Num>0
	ret.exitflag = MPTOPTIONS.OK;
	ret.how = 'ok';
else
	ret.exitflag = MPTOPTIONS.INFEASIBLE;
	ret.how = 'infeasible';
end
ret.stats.solveTime = etime(clock, start_time);

end
