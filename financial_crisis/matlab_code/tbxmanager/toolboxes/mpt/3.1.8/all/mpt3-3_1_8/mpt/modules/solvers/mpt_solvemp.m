function R = mpt_solvemp(S)
%
%  MPT_SOLVEMP: A gateway function to solve parametric optimization problems
%  =========================================================================
%  (without errorchecks) 
%  ======================
%  
%  
%  SYNTAX
%  ------
%     
%      R = mpt_solvemp(S)
%    
%  
%  DESCRIPTION
%  -----------
%     The main routine for fast calls to parametric optimization solvers. In fact,
%  it is a subroutine of Opt as a part of solve method. The Opt class serves as
%  general wrapper for preprocessing the data involved in optimization, including
%  necessary error checks. Once the data are valid, then are passed to mpt_solvemp
%  function that calls the appropriate parametric solver. It is assumed that
%  MPLP/MPQP entering this function are of the form 
%                                  1  T           T                 
%                            min   - x Hx+Ftheta+f x            (1) 
%                                  2                                
%                           s.t.   Ax <= b + Btheta             (2) 
%                                                                   
%                                  A x = b  + E                 (3) 
%                                   e     e                         
%                                  lb <= x <= ub                (4) 
%                                                                   
%                                  A     theta <= b             (5) 
%                                   theta          theta            
%     where the matrices H, F, A, A_e  A_theta, B, E, and vectors f, b, b_e,
%  b_theta, lb, ub are the problem data. Vector x  represents decision variables
%  and theta  are parameters. The PLCP must be given as: 
%                                w - Mz = q + Qtheta          (6) 
%                                             w >= 0          (7) 
%                                             z >= 0          (8) 
%                                      T                          
%                              w(theta) z(theta) = 0          (9) 
%                                                                 
%                             A     theta <= b               (10) 
%                              theta          theta               
%                                                            (11) 
%     where the matrices M, Q, A_theta, and vectors q, b_theta  are the problem
%  data, then z, w  are the decision variables and theta  are the parameters. The
%  routine mpt_solve processes data from any of the above optimization problems by
%  passing it to the appropriated solver. Following solvers are supported: 
%    
%     - PLCP is the default solver, called automatically when invoking parametric
%     optimization. It can solve MPQP/MPLP and PLCP problems. 
%     - MPLP is the solver from MPT2.6 which is still supported for development
%     purposes. 
%     - MPQP is the solver from MPT2.6 which is still supported for development
%     purposes. 
%    Note that this function must contain all properties of Opt class that have
%  been properly validated in order to perform a correct call to given parametric
%  solver. It is recommended to use Opt.solve method instead.
%  
%  INPUT
%  -----
%     
%        
%          S                  Object of the Opt class                  
%                             Class: Opt                               
%          S.H                Quadratic part of the objective          
%                             function.                                
%                             Class: double                            
%                             Default: []                              
%          S.f                Linear part of the objective function.   
%                             Class: double                            
%          S.pF               Linear part of the objective function    
%                             for parameters.                          
%                             Class: double                            
%                             Default: []                              
%          S.A                Linear part of the inequality            
%                             constraints Ax <= b + Btheta.            
%                             Class: double                            
%          S.b                Right hand side of the inequality        
%                             constraints Ax <= b + Btheta.            
%                             Class: double                            
%          S.pB               Right hand side of the inequality        
%                             constraints for parameters Ax <= b +     
%                             Btheta.                                  
%                             Class: double                            
%          S.Ae               Linear part of the equality constraints  
%                             A_ex=b_e + Etheta .                      
%                             Class: double                            
%                             Default: []                              
%          S.be               Right hand side of the equality          
%                             constraints A_ex=b_e + Etheta .          
%                             Class: double                            
%                             Default: []                              
%          S.pE               Right hand side of the equality          
%                             constraints for parameters A_ex=b_e +    
%                             Etheta .                                 
%                             Class: double                            
%                             Default: []                              
%          S.lb               Lower bound for the decision variables x 
%                             >= lb.                                   
%                             Class: double                            
%                             Default: []                              
%          S.ub               Upper bound for the decision variables x 
%                             <= ub.                                   
%                             Class: double                            
%                             Default: []                              
%          S.Ath              Linear part of the inequality            
%                             constraints A_thetatheta <= b_theta.     
%                             Class: double                            
%                             Default: []                              
%          S.bth              Right hand side of the inequality        
%                             constraints A_thetatheta <= b_theta.     
%                             Class: double                            
%                             Default: []                              
%          S.M                Linear matrix involved in LCP.           
%                             Class: double                            
%                             Default: []                              
%          S.q                Right hand side vector involved in LCP.  
%                             Class: double                            
%                             Default: []                              
%          S.Q                Linear matrix involved in parametric     
%                             formulation of LCP.                      
%                             Class: double                            
%                             Default: []                              
%          S.n                Number of decision variables.            
%                             Class: double                            
%          S.d                Number of parameters.                    
%                             Class: double                            
%          S.m                Number of inequalities in Ax <= b +      
%                             Btheta.                                  
%                             Class: double                            
%          S.me               Number of equalities in A_ex=b_e +       
%                             Etheta.                                  
%                             Class: double                            
%          S.problem_type     A string specifying the problem to be    
%                             solved                                   
%                             Class: char                              
%                             Default: []                              
%          S.vartype          A string array reserved for              
%                             MPMILP/MPMIQP.                           
%                             Class: char                              
%                             Default: "                               
%          S.solver           S string specifying which solver should  
%                             be called.                               
%                             Class: char                              
%                             Default: []                              
%          S.isParametric     Logical scalar indicating that the       
%                             problem is parametric.                   
%                             Class: double or logical                 
%                             Default: 1                               
%          S.varOrder         Order of variables if the problem was    
%                             processed by YALMIP first.               
%                             Class: double                            
%                             Default: []                              
%          S.Internal         Internal property of Opt class.          
%                             Class: struct                            
%                             Default: []                              
%          S.recover          Affine map for MPLP/MPQP problems after  
%                             transformation to LCP.                   
%                             Class: struct                            
%          S.recover.uX       Matrix of the affine map x = uX(         
%                              w                                       
%                              z  ) + uTh(                             
%                              theta                                   
%                                1    ) . The map is from the          
%                             optimization variables involed in LCP    
%                             w(theta),z(theta) |->> x  and in the     
%                             original LP/QP.                          
%                             Class: double                            
%                             Default: []                              
%          S.recover.uTh      Matrix of the affine map x = uX(         
%                              w                                       
%                              z  ) + uTh(                             
%                              theta                                   
%                                1    ) . The map is from the          
%                             optimization variables involed in LCP    
%                             w(theta),z(theta) |->> x  and in the     
%                             original LP/QP.                          
%                             Class: double                            
%                             Default: []                              
%          S.recover.lambdaX  Matrix of the affine map x = lambdaX(    
%                              w                                       
%                              z  ) + lambdaTh(                        
%                              theta                                   
%                                1    ) . The map is from the          
%                             optimization variables involed in LCP    
%                             w(theta),z(theta) |->> lambda  and the   
%                             Lagrangian multipliers in the original   
%                             LP/QP.                                   
%                             Class: double                            
%                             Default: []                              
%          S.recover.lambdaTh Matrix of the affine map x = lambdaX(    
%                              w                                       
%                              z  ) + lambdaTh(                        
%                              theta                                   
%                                1    ) . The map is from the          
%                             optimization variables involed in LCP    
%                             w(theta),z(theta) |->> lambda  and the   
%                             Lagrangian multipliers in the original   
%                             LP/QP.                                   
%                             Class: double                            
%                             Default: []                              
%                             Default: []                              
%                               
%  
%  
%  OUTPUT
%  ------
%     
%        
%          R           result structure                         
%                      Class: struct                            
%          R.xopt      Optimal solution with the associated     
%                      functions for optimizer, multipliers and 
%                      the objective value.                     
%                      Class: PolyUnion                         
%          R.exitflag  An integer value that informs if the     
%                      result was feasible (1), or otherwise    
%                      (different from 1)                       
%                      Class: double                            
%          R.how       A string that informs if the result was  
%                      feasible ('ok'), or if any problem       
%                      appeared through optimization            
%                      Class: char                              
%          R.solveTime Information about the time that elapsed  
%                      during the computation in seconds.       
%                      Class: double                            
%          R.stats     Further details from the parametric      
%                      solver.                                  
%                      Class: struct                            
%                        
%  
%  
%  SEE ALSO
%  --------
%     mpt_solve
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
 
 
narginchk(1, 1);

if ~isa(S,'Opt')
    error('mpt_solvemp: Input argument must be an "Opt" class.');
end

% check parameters
if S.d<1
    error('mpt_solvemp: Problem does not contain parametric inputs.');
end

% call specific solver
switch upper(S.solver)
           
    case {'PLCP'}
        
        % call PLCP solver
        R = mpt_call_plcp(S);

    case {'ENUMPLCP'}
        
        % call ENUMPLCP solver
        R = mpt_call_enum_plcp(S);
        
    case {'MPQP'}
        
        % call MPQP solver
        R = mpt_call_mpqp(S);
        
    case {'MPLP'}
        
        % call MPLP solver
        R = mpt_call_mplp(S);
        
    case {'ENUMPQP'}
        % enumeration-based pQP solver
        R = mpt_enum_pqp(S, struct('regions', true));

    case {'RLENUMPQP'}
        % enumeration-based pQP solver
        R = mpt_enum_pqp(S, struct('regions', false));

    otherwise
        
        % unknown solver
        error('mpt_solvemp: Solver not recognized.');
end
        
    
    
end
