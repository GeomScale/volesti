function mpt_demo_lti1
%
%  MPT_DEMO_LTI1: Simulation of LTISystem 
%  =======================================
%  
%  
%  SYNTAX
%  ------
%     
%      mpt_demo_lti1
%    
%  
%  DESCRIPTION
%  -----------
%     Simulation of LTI systems.
%  
%  SEE ALSO
%  --------
%     mpt_demo_lti2,  mpt_demo_lti3,  mpt_demo_lti4,  mpt_demo_lti5
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
lti = LTISystem('A', A, 'B', B, 'C', C, 'D', D);

% Set the initial state of the system
lti.initialize([1; 1.5]);

% Ask for the states
disp('Initial state:');
x = lti.getStates()

% Update the system's state using some control action
disp('State update using u=0.5:');
u = 0.5;
lti.update(u); % this updates the internal state
x = lti.getStates() % ask for the updated states


% Value of the state update can also be directly obtained from update().
disp('Next state update using u=-0.6:');
u = -0.6;
next_x = lti.update(u)

disp('System output using the last known state:');
y = lti.output()

end
