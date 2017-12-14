function mpt_demo_lti4
%
%  MPT_DEMO_LTI4: Construction of explicit controller for LTI system 
%  ==================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      mpt_demo_lti4
%    
%  
%  DESCRIPTION
%  -----------
%     Construction of explicit controller for LTI system
%  
%  SEE ALSO
%  --------
%     mpt_demo_lti1,  mpt_demo_lti2,  mpt_demo_lti3,  mpt_demo_lti5
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
% Use quadratic cost function
ctrl.model.x.penalty = QuadFunction(eye(model.nx));
ctrl.model.u.penalty = QuadFunction(eye(model.nu));

% Finally, we can convert the controller to an explicit form:
disp('Generating the explicit solution:');
expctrl = ctrl.toExplicit()

% Compare optimal solutions
x0 = [-4; 0];
disp('Optimal control input obtained by evaluating the explicit solution:');
Uexp = expctrl.evaluate(x0)

disp('Optimal control input obtained by solving the optimization problem on-line:');
Uonl = ctrl.evaluate(x0)

close all
% plot the explicit optimizer
expctrl.feedback.fplot();
title('PWA representation of the optimal control input')

% plot the value function
figure
expctrl.cost.fplot();
title('Explicit cost function');

% plot the regions
figure
expctrl.partition.plot()
title('Regions of the polyhedral partition');


end
