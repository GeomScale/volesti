function mpt_demo_functions1
%
%  MPT_DEMO_FUNCTIONS1: Demonstration of functions associated to sets 
%  ===================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      mpt_demo_functions1
%    
%  
%  DESCRIPTION
%  -----------
%     Demonstration of functions associated to sets.
%  
%  SEE ALSO
%  --------
%     mpt_demo_functions2
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


%% constructing general functions
disp('Create Function object f(x) = x');
F1 = Function(@(x)x)

pause

disp(' ');disp(' ');
disp('Create Function object f(x) = x(1)-x(2).^3 ');
F2 = Function(@(x) x(1)-x(2).^3 )


pause

disp(' ');disp(' ');
disp('Create Function with parameter K, f(x) = -log(K*x)');
disp('Since the parameter value may change, we first create the  object with the parameter K.')
F3 = Function([],struct('K',eye(2)))

pause

disp(' ');disp(' ');
disp('Once the object has been created, we can assign the handle and refer to already stored parameter K.')
F3.setHandle(@(x) -log(F3.Data.K*x))

disp(' ');disp(' ');
disp('We can change the value of the parameter any time later.');
F3.Data.K = 2*eye(2)

pause
%% constructing linear and affine functions
disp(' ');disp(' ');
disp('Construct object with an affine map f(x) = 6*x + 1')
L1 = AffFunction(6,1)

pause

disp(' ');disp(' ');
disp('Construct object with an affine map f(x) = -x(1)+x(2) + 1')
L2 = AffFunction([-1,1],1)

disp(' ');disp(' ');
disp('Vector function f(x) = eye(2)*[x1;x2] + [1;2]')
L3 = AffFunction(eye(2),[1;2])

disp(' ');disp(' ');
disp('Linear function f(x) = eye(5)*x');
L4 = AffFunction(eye(5))

pause
%% constructing quadratic functions
disp(' ');disp(' ');
disp('Quadratic map f(x) = 0.5x^2 + 1')
Q1 = QuadFunction(0.5,1)

pause

disp(' ');disp(' ');
disp('Quadratic map f(x) = x^2 - 4')
Q2 = QuadFunction(2,-4)

pause

disp(' ');disp(' ');
disp('Quadratic map f(x) = x(1)^2+x(2)^2 + 1')
Q3 = QuadFunction(eye(2),[0,0],1)

pause

%% assign function to a set
disp(' ');disp(' ');
disp('Function can be attached to a set after construction.');
disp('Construct the polyhedron first');
P1 = Polyhedron('lb',-1,'ub',1)
pause
disp('Add a quadratic function under the name "a".');
P1.addFunction(QuadFunction(4,-1),'a');

pause
disp(' ');disp(' ');
disp('Add an affine function under the name "b".');
P2 = Polyhedron('V',[-1 1;1 1; -1 -1])
P2.addFunction(AffFunction(-eye(2),[-1;2]),'b')

disp(' ');disp(' ');
disp('You can assign multiple functions to a set under different names.')
P3 = Polyhedron('lb',[-1;-1],'ub',[1;1]);
P3.addFunction(Function(@(x)x.^2-x.^3+1),'gain');
P3.addFunction(AffFunction(randn(2)),'power')

pause

disp(' ');disp(' ');
disp('Multiple functions can be assigned only at separate calls.');
P4 = Polyhedron('V',randn(6,2));
L(1) = AffFunction(randn(2),randn(2,1));
L(2) = AffFunction(randn(2),randn(2,1));
P4.addFunction(L(1),'a')
P4.addFunction(L(2),'b')

pause

%% plotting
disp(' ');disp(' ');
disp('Plot the function over the set is based on the function name');
P3.fplot('gain')
pause
P3.fplot('power')

pause
disp(' ');disp(' ');
disp('Plot the first element of the vector valued function "a"');
P4.fplot('a','position',1,'color','y')

pause
disp(' ');disp(' ');
disp('Plot the second element of the vector valued function "b" based on the index.');
P4.fplot('b','position',2,'color','m')

pause

%% evaluation
disp(' ');disp(' ');
disp('Evaluate the function based on the name');
P3.feval([1;0],'power')

disp(' ');disp(' ');
disp('Evaluate the function "gain"');
P3.feval([-1;0],'gain')

end
