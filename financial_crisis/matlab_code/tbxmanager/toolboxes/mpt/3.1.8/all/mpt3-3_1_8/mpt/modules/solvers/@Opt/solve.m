function [ret,opt] = solve(opt, varargin)
%
%  SOLVE: The main routine for solving optimization problems 
%  ==========================================================
%  
%  
%  SYNTAX
%  ------
%     
%      result = problem.solve
%      result = problem.solve(th)
%      result = solve(problem)
%    
%  
%  DESCRIPTION
%  -----------
%     The main routine for solving non-parametric (LP, QP, MILP, MIQP) and
%  parametric problems (MPLP, MPQP, PLCP). The result is returned in the
%  appropriate format depending on the problem. The Opt class serves as general
%  wrapper for preprocessing the data involved in optimization, including necessary
%  error checks. Once the data are valid, then are passed to mpt_solve or
%  mpt_solvemp function that calls the appropriate solver without any errorchecks.
%    For parametric problems it is possible to solve the problem for a particular
%  value of the parameters theta  if provided as an argument th.
%  
%  INPUT
%  -----
%     
%        
%          problem Object of the Opt class that defines the 
%                  optimization problem to be solved.       
%                  Class: Opt                               
%                    
%  
%  
%  OUTPUT
%  ------
%     
%        
%          result          Structure with the date that represents  
%                          the solution to given problem. For       
%                          non-parametric problems the solution is  
%                          returned with the following fields.      
%                          Class: struct                            
%          result.xopt     Optimal solution for primal variables.   
%                          Class: double                            
%          result.obj      Objective value.                         
%                          Class: double                            
%          result.lambda   Lagrangian multipliers                   
%                          Class: double                            
%          result.exitflag An integer value that informs if the     
%                          result was feasible (1), or otherwise    
%                          (different from 1)                       
%                          Class: double                            
%          result.how      A string that informs if the result was  
%                          feasible ('ok'), or if any problem       
%                          appeared through optimization            
%                          Class: char                              
%          result          Structure with the date that represents  
%                          the solution to given problem. For       
%                          parametric problems the solution is      
%                          returned with the following fields.      
%                          Class: struct                            
%          result.xopt     Optimal solution for variables w, z      
%                          from PLCP reformulation if the problem   
%                          was solved using PLCP solver. If the     
%                          original problem was given as MPLP/MPQP, 
%                          then this field returns also primal,     
%                          dual variables, and the objective value  
%                          for given problem. The solution is given 
%                          as a collection of polyhedra in the same 
%                          dimension with certain properties and is 
%                          given as PolyUnion class. The function   
%                          data associated to each variable is      
%                          stored under Function, in particular in  
%                          res.Set(i).Funcj where the index i       
%                          corresponds to i-th region and index j   
%                          to j-th function.                        
%                          Class: PolyUnion                         
%          result.exitflag An integer value that informs if the     
%                          result was feasible (1), or otherwise    
%                          (different from 1)                       
%                          Class: double                            
%          result.how      A string that informs if the result was  
%                          feasible ('ok'), or if any problem       
%                          appeared through optimization            
%                          Class: char                              
%          result.stats    Other details from the computation that  
%                          might be of interest, such as the number 
%                          of pivots, elapsed time, etc.            
%                          Class: struct                            
%                            
%  
%  
%  SEE ALSO
%  --------
%     Opt,  mpt_solve,  mpt_solvemp
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

% deal with arrays
ret=cell(size(opt));
if numel(opt)>1
    for i=1:numel(opt)
        ret{i} = opt(i).solve;
    end
    return
end

% Determine problem type
if opt.isParametric && nargin==1
    % for parametric solvers use "mpt_solvemp" which depends
    % on OPT and POLYHEDRON classes
    ret = mpt_solvemp(opt);
else
    % non-parametric solvers can be called directly.
	%
	% alternatively, we can solve a parametric problem for a particular
	% value of the parameter
    ret = mpt_solve(opt, varargin{:});
end

end
