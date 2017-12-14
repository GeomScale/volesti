function mpt_demo_lti2
%
%  MPT_DEMO_LTI2: Demonstrates online MPC for LTI system 
%  ======================================================
%  
%  
%  SYNTAX
%  ------
%     
%      mpt_demo_lti2
%    
%  
%  DESCRIPTION
%  -----------
%     Demonstrates online MPC for LTI system.
%  
%  SEE ALSO
%  --------
%     mpt_demo_lti1,  mpt_demo_lti3,  mpt_demo_lti4,  mpt_demo_lti5
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

% Define an LTI dynamics x^+ = A x + B u
% (note that we can ommit the output equation)
A = [1 1; 0 1];
B = [1; 0.5];
lti = LTISystem('A', A, 'B', B);

% Define an MPC controller using "lti" as the prediction model
ctrl = MPCController(lti);

% Set the prediction horizon
ctrl.N = 10;

% Add constraints on predicted states
ctrl.model.x.min = [-5; -5];
ctrl.model.x.max = [5; 5];

% Add constraints on predicted control inputs
ctrl.model.u.min = -1;
ctrl.model.u.max = 1;

% Use quadratic state penalty with identity weighting matrix
W = eye(2);
ctrl.model.x.penalty = QuadFunction(W);

% Set quadratic input penalty with identity weighting matrix
W = 1;
ctrl.model.u.penalty = QuadFunction(W);

% Obtain the optimal control input for a given initial condition
x0 = [-4; 0];
u = ctrl.evaluate(x0);

% We can also ask for the open-loop predictions:
[u, feasible, openloop] = ctrl.evaluate(x0);
openloop

% Plot the open-loop trajectories
ctrl.model.plot()

end
