function mpt_demo_deployment_onlineMPC
%
%  MPT_DEMO_DEPLOYMENT_ONLINEMPC: Application of online MPC controller with the
%  ============================================================================
%  help of Simulink interface 
%  ===========================
%  
%  
%  SYNTAX
%  ------
%     
%      mpt_demo_deployment_onlineMPC
%    
%  
%  DESCRIPTION
%  -----------
%     Demonstration of real-time control using online MPC controller. The demo
%  relies on standalone LCP solver that solves the optimization problem given as
%  linear-complementarity problem (LCP). The LCP solver must be present in the
%  installation path including the extended version lcprtw that contains
%  pre-compiled libraries for linking with OpenWatcom compiler. The demo has been
%  tested under Real-Time Windows target and XPC target on Windows 32-bit platform.
%  
%     Note that this demo can be compiled and run on Windows 32-bit platform for
%  Matlab r2012a!
%    Deployment steps: 
%    
%     1. Generate the online controller that has the desired properties. 
%     2. Export the online controller to YALMIP using toYalmip method of
%     MPCController class. 
%     3. Formulate a parametric optimization problem and specify the feedback
%     variables that are used as parameters. 
%     4. Create an instance of Opt class and transform the optimization problem to
%     LCP using qp2lcp method. 
%     5. Create a Simulink scheme with the S-Function block that links to LCP
%     solver. 
%     6. In the Simulink scheme, choose code generation options and pick rtwin.tlc
%     as the system target file  that corresponds to Real-Time Windows target or
%     xpctarget.tlc for XPC target. 
%     7. Press "CTRL+B" that executes the code generation and compiles the code. 
%     8. In the Simulink scheme choose "Simulation->External" option and press
%     "Connect To Target". 
%     9. Start the simulation to verify the controller in real-time. 
%  
%  
%  SEE ALSO
%  --------
%     mpt_demo_deployment_explicitMPC
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
 
 
A=[   0.5403   -0.8415; 0.8415    0.5403];
B=[ -0.4597; 0.8415];
C=[1 0];
D=0;

% linear discrete-time model with sample time 1
sys = ss(A,B,C,D,1);

model = LTISystem(sys);

% set constraints on output
model.y.min = -10;
model.y.max = 10;

% set constraints on input
model.u.min = -1;
model.u.max = 1;

% weights on states/inputs
model.x.penalty = QuadFunction(eye(2));
model.u.penalty = QuadFunction(1);

% terminal set
Tset = model.LQRSet;

% terminal weight
PN = model.LQRPenalty;

% add terminal set and terminal penalty
model.x.with('terminalSet');
model.x.terminalSet = Tset;
model.x.with('terminalPenalty');
model.x.terminalPenalty = PN;


% MPC controller
ctrl = MPCController(model,5);

% export variables to YALMIP
v = ctrl.toYALMIP;

% initial condition is the parameter theta
theta = v.variables.x(:,1); 

% extract remaining variables
states = v.variables.x(:,2:end);
inputs = v.variables.u;
outputs = v.variables.y;

% transform the problem to parametric LCP
problem = Opt(v.constraints,v.objective,theta,[inputs(:); states(:); outputs(:)]);
problem.qp2lcp;

% open the simulink scheme
mpt_demo_rtw_onlinempc

% export the problem data to workspace for simulation to run
assignin('base','problem',problem);

end
