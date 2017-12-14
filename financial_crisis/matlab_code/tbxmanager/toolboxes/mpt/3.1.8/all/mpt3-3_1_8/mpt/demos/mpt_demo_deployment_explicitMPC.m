function mpt_demo_deployment_explicitMPC
%
%  MPT_DEMO_DEPLOYMENT_EXPLICITMPC: Application of explicit MPC controller with the
%  ================================================================================
%  help of Simulink interface 
%  ===========================
%  
%  
%  SYNTAX
%  ------
%     
%      mpt_demo_deployment_explicitMPC
%    
%  
%  DESCRIPTION
%  -----------
%     Demonstration of real-time control using explicit MPC controller. This demo
%  shows how to design an explicit controller for usage in real-time control in
%  Real-Time Workshop. The code generation and compilation for real-time has been
%  tested under Real-Time Windows target on Windows 32-bit and 64-bit platforms.
%    Deployment steps: 
%    
%     1. Generate the explicit controller that has the desired properties. 
%     2. Export the explicit controller to C-code using mpt_exportToC function. 
%     3. Create a Simulink scheme with the S-Function block that represents
%     generated controller. 
%     4. In the Simulink scheme, choose code generation options and pick rtwin.tlc
%     as the system target file (e.g. for Real-Time Windows target). 
%     5. In the code generation options, go to "Custom Code" tab and put in "Source
%     files" the absolute path to generated files. 
%     6. Press "CTRL+B" that executes the code generation and compiles the code. 
%     7. In the Simulink scheme choose "Simulation->External" option and press
%     "Connect To Target". 
%     8. Start the simulation to verify the controller. 
%  
%  
%  SEE ALSO
%  --------
%     mpt_demo_deployment_onlineMPC
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

% explicit controller
ectrl = ctrl.toExplicit;

% export explicit controller to C
dir_name = 'rtw_explicitMPC';
file_name = 'EMPCcontroller';
ectrl.exportToC(file_name,dir_name);

% compile the S-function
p = [pwd,filesep,dir_name,filesep];
mex(['-LargeArrayDims -I',dir_name],[p,file_name,'_sfunc.c'],[p,file_name,'.c']);

% open the simulink scheme
mpt_demo_rtw_explicitmpc

% set the source files in the Simulink scheme
str = sprintf('"%s%s_sfunc.c" "%s%s.c"',p,file_name,p,file_name);
set_param('mpt_demo_rtw_explicitmpc','CustomSource',str);
set_param('mpt_demo_rtw_explicitmpc/Controller','FunctionName',[file_name,'_sfunc']);

end
