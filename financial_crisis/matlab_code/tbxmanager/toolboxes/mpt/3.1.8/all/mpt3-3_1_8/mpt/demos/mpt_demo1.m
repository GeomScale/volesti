function mpt_demo1
%
%  MPT_DEMO1: Demonstration of basic usage of the geometric library 
%  =================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      mpt_demo1
%    
%  
%  DESCRIPTION
%  -----------
%     Basic usage of the new interface to geometric library
%  
%  SEE ALSO
%  --------
%     mpt_demo_sets1,  mpt_demo_functions1,  mpt_demo_unions1
%  

%  AUTHOR(s)
%  ---------
%     
%    
%   (c) 2010-2013  Colin Neil Jones: EPF Lausanne
%   mailto:colin.jones@epfl.ch 
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
 
 
clear;

%% Create polytope:
P = Polyhedron('V', randn(10,2));

% Plotting
clf;
P.plot('color','b','alpha',0.3,'linewidth',1,'linestyle','--');
axis off

% set(gcf, 'renderer', 'zbuffer');
%print -dpng '2dPolytope.png'

% Double description created and stored automatically
P
%% Create Polyhedron
disp('Polyhedron: Press any key to continue...')
pause

P = Polyhedron('V', randn(10,2), 'R', randn(1,2));
clf;
P.plot;

disp('Polyhedron: Press any key to continue...')
pause


P = Polyhedron('V', randn(50,3), 'R', randn(1,3));
clf;
P.plot;
axis vis3d
figure(1);


%% Lower-dimensional polyhedra
disp('Low-dimensional polyhedron: Press any key to continue...')
pause

P = Polyhedron('H', [randn(30,3) ones(30,1)], 'He', [randn(1,3) 0]);
clf;
P.plot('alpha',0.3);


%% Boxes
disp('Box: Press any key to continue...')
pause

P = Polyhedron('lb', -rand(3,1), 'ub', rand(3,1), 'H', [randn(10,3) 1.5*ones(10,1)]);
clf;
P.plot('alpha',0.3);


%% Polyhedron queries
disp('Project onto polyhedron: Press any key to continue...')
pause

% Project onto polyhedron
P = Polyhedron(randn(10,2), ones(10,1));
clf;
P.plot('color','g');
hold on;
axis square;


x = 5*randn(2,1);
sol = P.project(x);
pplot([x sol.x]','bo');
fprintf('Distance: %.2f\n', sol.dist);


% Separate point
disp('Separate point: Press any key to continue...')
pause

sep = P.separate(x);
v = axis;
s = Polyhedron('He', sep, 'lb', [v(1);v(3)],'ub',[v(2);v(4)]);
s.plot;

% Interior point
disp('Interior point: Press any key to continue...')
pause

sol = P.interiorPoint;
pplot(sol.x, 'ro');


%% Addition
disp('Addition: Press any key to continue...')
pause

P = Polyhedron(3*randn(10,2));
Q = Polyhedron('H', [randn(10,2) ones(10,1)]);

clf;
P.plot; hold on;
Q.plot;

disp('Press any key to continue...')
pause

PQ = P+Q;
PQ.plot('alpha',0.1,'color','r');

%% Lower-dimensional addition 
disp('lower-dimensional addition: Press any key to continue...')
pause

for i=1:5
  P(i) = Polyhedron('V', [0 0 0;randn(1,3)]);
end
clf;
plot(P,'linewidth',2); hold on;

disp('Press any key to continue...')
pause

Q = Polyhedron;
for i=1:5
  Q = Q + P(i);
end
Q.plot('alpha',0.1,'color','b');
axis vis3d

%% Standard

% Can compute distances, test inclusion, convex hulls, test equality,
% incidence maps, intersections, boundedness, full-dim, 

%% Convex sets
disp('Convex sets: Press any key to continue...')
pause

x = sdpvar(2,1);
T = eye(2);
F = [x <= 0.1 ; x'*T'*T*x <= 1];
opt = sdpsettings('solver','sedumi','verbose',0);
Y = YSet(x, F, opt);

clf;
Y.plot('alpha',0.3,'color','b','grid',50);
hold on;

fprintf('Press any key to test point separation...\n');
pause

x = [-2;-2];
sol = Y.project(x);
pplot([x sol.x]','bo');

% Separate point
disp('Separate point: Press any key to continue...')
pause

sol = Y.separate(x);
v = axis;
s = Polyhedron('He', sol, 'lb', [v(1);v(3)],'ub',[v(2);v(4)]);
s.plot('alpha',0.2);

% Compute outer approximation
fprintf('Press any key to test outer approximation...\n');
pause

O = Y.outerApprox;
O.plot('alpha',0.1,'color','b');
axis(axis*1.1);

%% Create set of Polytopes
disp('Create set of polytopes: Press any key to continue...')
pause

Ps = PolyUnion;
for i=1:5
  Ps.add(Polyhedron(randn(10,2)) + 5*randn(2,1));
end
for i=1:5
  Ps.add(Polyhedron('H', [randn(10,2) ones(10,1)], 'He', [randn(1,2) 0]) + 5*randn(2,1));
end
for i=1:5
  Z = Polyhedron;
for j=1:3, Z = Z + Polyhedron([0 0;randn(1,2)]); end
  Ps.add(Z+5*randn(2,1));
end

clf;
Ps.plot;

%% Create complex
disp('Triangulation: Press any key to continue...')
pause

T = Polyhedron('H', [randn(30,2) ones(30,1)]);
T = triangulate(T);

clf;
T.plot;

%% Polyhedral functions
disp('Polyhedral functions: Press any key to continue...')
pause

F = Polyhedron('V', 2*randn(10,2));
F.addFunction(Function(@(x) sin(x(1))*cos(x(2))),'func1');

clf;
F.fplot;

%% PolyUnion
disp('Functions over sets: Press any key to continue...')
pause

tic
for i=1:5
  F(i) = Polyhedron('V', 2*randn(10,2)) + 5*randn(2,1);
end
F.addFunction(Function(@(x) sin(x(1))*cos(x(2))),'wave');
Ps = PolyUnion(F);

clf;
Ps.fplot;
toc


%% Parametric solutions
disp('Parametric solutions');
disp('Solving MPLP');
pause

n = 2;
m = 1;
N = 5;

A = [1 1; 0 1];
B = [1; 0.5];
xlb = -10*ones(n,1);
xub = 10*ones(n,1);
ulb = -5*ones(m,1);
uub = 5*ones(m,1);
R = 2*eye(m);
Q = 0.2*eye(n);

% MPQP
% formulate problem using Yalmip
x = sdpvar(n,N,'full');
u = sdpvar(m,N-1,'full');   
cost = 0;
F = [];
for i=1:N-1
  F = F + (x(:,i+1) == A*x(:,i) + B*u(:,i));
  F = F + (xlb <= x(:,i) <= xub);
  F = F + (ulb <= u(:,i) <= uub);
  
  if i > 1
    cost = cost + x(:,i)'*Q*x(:,i);
  end
  cost = cost + u(:,i)'*R*u(:,i);
end
F = F + [xlb <= x(:,end) <= xub];
cost = cost + x(:,end)'*Q*x(:,end);

% solve using MPT2.6
fprintf('===> SOLVING USING MPQP <===\n');
t0=cputime;
% change globally the solver
mptopt('pqpsolver','MPQP');

% construct problem using MPQP solver
problem1 = Opt(F, cost, x(:,1), u(:));

% solve
res1 = problem1.solve;
t1 = cputime-t0;
fprintf('MPQP solution took %.2f seconds.\n',t1);

% solve using PLCP
fprintf('\n\n===> SOLVING USING PLCP <===\n');
t0=cputime;
% change globally the solver
mptopt('pqpsolver','PLCP');

% call problem constructor
problem2 = Opt(F, cost, x(:,1), u(:));

% solve
res2 = problem2.solve;
t2 = cputime-t0;
fprintf('PLCP solution took %.2f seconds.\n\n',t2);

fprintf('Regions: PLCP vs MPT_MPQP : %i vs %i\n', res2.xopt.Num, res1.xopt.Num);
fprintf('Time   : PLCP vs MPT_MPQP : %.2fs vs %.2fs\n', t2, t1);

if n <= 3
  figure(1);
  plot(res1.xopt)
  title('PLCP');
  axis tight;
  figure(2);
  plot(res2.xopt);
  title('MPQP');
  axis tight;
end

disp('Solving MPLP');
pause


% MPLP
cost = 0;
F = [];
for i=1:N-1
  F = F + (x(:,i+1) == A*x(:,i) + B*u(:,i));
  F = F + (xlb <= x(:,i) <= xub);
  F = F + (ulb <= u(:,i) <= uub);
  
  if i > 1
    cost = cost + norm(Q*x(:,i),1);
  end
  cost = cost + norm(R*u(:,i),1);
end
F = F + [xlb <= x(:,end) <= xub];
cost = cost + norm(Q*x(:,end),1);

% solve using MPT2.6
fprintf('===> SOLVING USING MPLP <===\n');
t0 = cputime;
mptopt('plpsolver','MPLP');
problem3 = Opt(F, cost, x(:,1), u(:));
res3 = problem3.solve;
t1 = cputime - t0;
fprintf('MPLP solution took %.2f seconds.\n',t1);


% solve using PLCP
fprintf('\n\n===> SOLVING USING PLCP <===\n');
t0 = cputime;
mptopt('plpsolver','PLCP');
problem4 = Opt(F, cost, x(:,1), u(:));
res4 = problem4.solve;
t2 = cputime - t0;
fprintf('PLCP solution took %.2f seconds.\n\n',t2);

fprintf('Regions: PLCP vs MPT_MPLP : %i vs %i\n', res4.xopt.Num, res3.xopt.Num);
fprintf('Time   : PLCP vs MPT_MPLP : %.2fs vs %.2fs\n', t2, t1);

if n <= 3
  figure(1);
  plot(res4.xopt)
  title('PLCP');
  axis tight;
  figure(2);
  plot(res3.xopt);
  title('MPLP');
  axis tight;
end

disp('The end.')


end
