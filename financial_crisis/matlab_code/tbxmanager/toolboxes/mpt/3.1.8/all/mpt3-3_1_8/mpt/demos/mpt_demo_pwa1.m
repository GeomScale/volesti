function mpt_demo_pwa1
%
%  MPT_DEMO_PWA1: Demonstration for modeling PWA systems 
%  ======================================================
%  
%  
%  SYNTAX
%  ------
%     
%      mpt_demo_pwa1
%    
%  
%  DESCRIPTION
%  -----------
%     Demonstration for modeling PWA system comprising of two linear time-invariant
%  models: 
%                                {                        
%                                { A x + Bu if  x  <= 0   
%                            +   {  1            1        
%                           x  = {                        
%                                { A x + Bu if  x  >= 0   
%                                {  2            1        
%                            y = Cx + Du                  
%  
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
 
 
B = [0; 1]; C = [1 0]; D = 0;

% First dynamics:
alpha = -pi/3;
A1 = 0.8*[cos(alpha) -sin(alpha); sin(alpha) cos(alpha)];

dyn1 = LTISystem('A', A1, 'B', B, 'C', C, 'D', D);
% We need to tell that dynamics #1 should be active if x_1<=0:
dyn1.setDomain('x', Polyhedron([1 0], 0));

% Second dynamics:
alpha = pi/3;
A2 = 0.8*[cos(alpha) -sin(alpha); sin(alpha) cos(alpha)];

dyn2 = LTISystem('A', A2, 'B', B, 'C', C, 'D', D);
% Region of validity of the dynamics (x_1>=0):
dyn2.setDomain('x', Polyhedron([-1 0], 0));

% Create the PWA description using an array of LTI systems:
pwa = PWASystem([dyn1 dyn2]);

% Optionally we can set constraints:
pwa.x.min = [-10; -10];
pwa.x.max = [10; 10];
pwa.y.min = -10;
pwa.y.max = 10;
pwa.u.min = -1;
pwa.u.max = 1;

% Define an on-line MPC controller for such a system
horizon = 2;
onl_ctrl = MPCController(pwa, horizon);
% Set panalties used in the cost function:
onl_ctrl.model.x.penalty = OneNormFunction(10*eye(2));
onl_ctrl.model.u.penalty = OneNormFunction(1);

% Construct the explicit solution
exp_ctrl = onl_ctrl.toExplicit();

% Obtain the closed-loop optimizers for a particular initial condition
x0 = [-4; 0];
Uonl = onl_ctrl.evaluate(x0)
Uexp = exp_ctrl.evaluate(x0)

end
