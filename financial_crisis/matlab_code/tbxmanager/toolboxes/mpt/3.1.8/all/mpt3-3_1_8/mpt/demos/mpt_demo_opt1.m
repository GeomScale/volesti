function mpt_demo_opt1
%
%  MPT_DEMO_OPT1: Demonstration for using Opt interface 
%  =====================================================
%  
%  
%  SYNTAX
%  ------
%     
%      mpt_demo_opt1
%    
%  
%  DESCRIPTION
%  -----------
%     Demonstration for using Opt interface for solving optimization problems
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
 
 
close all

%% formulate LP problem min f'*x s.t. A*x<=b
f = randn(8,1);
A = randn(25,8);
b = 5*rand(25,1);

disp('Formulate LP')
problem1 = Opt('f',f,'A',A,'b',b)

disp('Solve LP')
res1 = problem1.solve

pause

%% formulated MPLP problem min f'*x s.t. A*x<=b+B*theta
% put bounds on the parameters -1<theta<1

disp('Formulate MPLP')
problem2 = Opt('f',f,'A',A,'b',b,'pB',ones(25,1),'Ath',[-1;1],'bth',[1;1])

disp('Solve MPLP')
res2 = problem2.solve

disp('Solution is stored as "primal", we can plot it');
res2.xopt.fplot('primal','show_set',true,'LineWidth',3)
title('Primal solution')

pause

%% formulate MPLP problem using MPT2-solver
disp('Formulate MPLP with MPT2-MPLP solver')
problem3 = Opt('f',f,'A',A,'b',b,'pB',ones(25,1),'Ath',[-1;1],'bth',[1;1],'solver','MPLP')

disp('Solve MPLP using MPT2')
res3 = problem3.solve

disp('Plot the objective value')
res3.xopt.fplot('obj','show_set',true,'LineWidth',3)
title('Objective value')

pause

%% formulate problem using YALMIP
disp('Formulate problem using YALMIP');
% Model data
A = [0.5 -1;1 0];
B = [1;0];
nu = 1; % Number of inputs

% MPC data
Q = eye(2);
R = 1;
N = 4;

% Initial state
x0 = sdpvar(2,1);

% setup the problem
u = sdpvar(nu,N);

constraints = [];
objective = 0;
x = x0;
for k = 1:N
 x = A*x + B*u(k);
 objective = objective + norm(Q*x,1) + norm(R*u(k),1);
 constraints = [constraints, -1 <= u(k)<= 1, -5<=x<=5];
end

% formulate the problem
problem4 = Opt(constraints, objective, x0, u)


disp('Solve the problem');
res4 = problem4.solve

disp('Plot the partition')
res4.xopt.plot

pause


%% formulate problem using MPT2
disp('Formulate problem by as in MPT2.');

sysStruct.A= [1 1; 0 1];
sysStruct.B= [1; 0.5];
sysStruct.C= [1 0; 0 1];
sysStruct.D= [0;0];

%set constraints on output
sysStruct.ymin = [-5; -5];
sysStruct.ymax = [5; 5];

%set constraints on input
sysStruct.umin = -1;
sysStruct.umax = 1;

% problem formulation
probStruct.norm=2;
probStruct.Q=eye(2);
probStruct.R=1;
probStruct.N=5;
probStruct.subopt_lev=0;

Matrices = mpt_constructMatrices(sysStruct,probStruct);

problem5 = Opt(Matrices)


disp('Solve the problem');
res5 = problem5.solve

disp('Plot he objective function')
res5.xopt.fplot('obj')

end
