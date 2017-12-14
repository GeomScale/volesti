function mpt_demo_lti5
%
%  MPT_DEMO_LTI5: Demostration of problem formulation using additional properties. 
%  ================================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      mpt_demo_lti5
%    
%  
%  DESCRIPTION
%  -----------
%     Formulation of MPC problem using terminal cost, terminal set constraints, and
%  move blocking constraint.
%  
%  SEE ALSO
%  --------
%     mpt_demo_lti1,  mpt_demo_lti2,  mpt_demo_lti3,  mpt_demo_lti4
%  

%  AUTHOR(s)
%  ---------
%     
%    
%   (c) 2003-2013  Michal Kvasnica: STU Bratislava
%   mailto:michal.kvasnica@stuba.sk 
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
 
 
close all

% define an LTI prediction model
lti = LTISystem('A', [1 1; 0 1], 'B', [1; 0.5]);

% define the MPC controller
horizon = 5;
ctrl = MPCController(lti, horizon);

% define quadratic penalties
ctrl.model.x.penalty = QuadFunction(eye(2));
ctrl.model.u.penalty = QuadFunction(1);

% add a terminal set constraint (see help SystemSignal/filter_terminalSet)
ctrl.model.x.with('terminalSet');
ctrl.model.x.terminalSet = Polyhedron('lb', [-1; -1], 'ub', [1; 1]);

% add an LQR terminal penalty (see help SystemSignal/filter_terminalPenalty)
lqr_penalty = ctrl.model.LQRPenalty();
ctrl.model.x.with('terminalPenalty');
ctrl.model.x.terminalPenalty = lqr_penalty;

% add a move-blocking constraint (the last 3 moves are to be constant)
ctrl.model.u.with('block');
ctrl.model.u.block.from = ctrl.N-2;
ctrl.model.u.block.to = ctrl.N;

% obtain the optimal control input
x0 = [-4; 0];
[Uonl, feasible] = ctrl.evaluate(x0)

% we can also ask for full open-loop predictions:
[~, ~, openloop] = ctrl.evaluate(x0)

% plot the open-loop predictions
ctrl.model.plot()

end
