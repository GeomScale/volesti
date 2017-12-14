function mpt_demo_functions2
%
%  MPT_DEMO_FUNCTIONS2: Demonstration of functions over unions of polyhedra 
%  =========================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      mpt_demo_functions2
%    
%  
%  DESCRIPTION
%  -----------
%     Demonstration of functions over unions of polyhedra.
%  
%  SEE ALSO
%  --------
%     mpt_demo_functions1
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


%% constructing union of triangular polyhedra
disp('Create random polyhedron')
P = 10*ExamplePoly.randVrep

disp(' '); disp(' ');
disp('Triangulate the polyhedron to get a complex.');
T = P.triangulate

disp(' '); disp(' ');
disp('For each of the polyhedron, assign affine function with the same name')
for i=1:numel(T)
    T(i).addFunction(AffFunction(eye(2),[-1;1]),'phi');
end

pause
disp(' '); disp(' ');
disp('Construct the polyunion object U.')
U = PolyUnion('Set',T,'FullDim',true,'Bounded',true,'Overlaps',false,'Convex',true)

pause
disp(' '); disp(' ');
disp('Plot the function over the polyhedra')
U.fplot

pause,
close all

%% construct overlapping union
disp('Create 3 random polyhedra.');
for i=1:3
    Q(i) = ExamplePoly.randVrep+5*rand(2,1);
end

pause
disp('Assign two quadratic functions to each of the polyhedra.')
for i=1:3
    Q(i).addFunction(QuadFunction(eye(2),randn(1,2),randn(1)),'alpha')
    Q(i).addFunction(QuadFunction(eye(2),randn(1,2),randn(1)),'beta')
end
pause

disp(' '); disp(' ');
disp('Create union without specifying any properties.');
PU = PolyUnion(Q)

pause
disp(' '); disp(' ');
disp('Plot the functions over polyhedra.');
PU.fplot('beta')

pause
disp('Plot the functions over polyhedra with some properties.');
PU.fplot('beta','show_set',true)


end
