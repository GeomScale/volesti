function mpt_demo_lti3
%
%  MPT_DEMO_LTI3: Demonstrates simulation of the closed-loop system 
%  =================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      mpt_demo_lti3
%    
%  
%  DESCRIPTION
%  -----------
%     Demonstrates simulation of the closed-loop system.
%  
%  SEE ALSO
%  --------
%     mpt_demo_lti1,  mpt_demo_lti2,  mpt_demo_lti4,  mpt_demo_lti5
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
 
 
A = [1 1; 0 1];
B = [1; 0.5];
C = [1 0];
D = 0;
model = LTISystem('A', A, 'B', B, 'C', C, 'D', D);

% Next we define an MPC controller:
horizon = 5;
ctrl = MPCController(model, horizon);

% Specify the MPC problem (constraints and penalties):
ctrl.model.x.min = [-5; -5];
ctrl.model.x.max = [5; 5];
ctrl.model.u.min = -1;
ctrl.model.u.max = 1;
ctrl.model.x.penalty = OneNormFunction(eye(model.nx)); % 1-norm type penalty
ctrl.model.u.penalty = InfNormFunction(eye(model.nu)); % Inf-norm type penalty

% Create the closed-loop system:
loop = ClosedLoop(ctrl, model);

% Simulate the closed loop from a given initial condition
x0 = [-4; 0];
N_sim = 20;
data = loop.simulate(x0, N_sim);

% Plot the simulated state trajectories
plot(0:N_sim, data.X);

end
