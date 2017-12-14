function slide = mpt_demo_sets3
%
%  MPT_DEMO_SETS3: Demo that illustrates operations on polyhedra 
%  ==============================================================
%  
%  
%  SYNTAX
%  ------
%     
%      mpt_demo_sets3
%    
%  
%  DESCRIPTION
%  -----------
%     Demo that illustrates geometric operations with polyhedra. In particular, it
%  is shown how to 
%    
%     - convert between H- and V-representation 
%     - intersection of polyhedra 
%     - comparisons between polyhedra 
%     - set algebra: affine maps, set summation, set difference 
%  
%  
%  SEE ALSO
%  --------
%     mpt_demo_sets1,  mpt_demo_sets2
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
 
 
if nargout<1,
  playshow mpt_demo_sets3
else
    %%intro

    slide(1).code={
        'clc; cla;'
        'axis([-4 4 -2 2]); axis off; text(-4,0.5,''Geometric manipulation with Polyhedra'',''FontSize'',16);'};
    slide(1).text={
        'This tour presents the geometric manipulations with polyhedra.'
        ''
        '  - Conversion between H- and V-representation'
        '  - Intersections of polyhedra'
        '  - Comparisons of polyhedra'
        '  - Set algebra: affine maps, set summation, set difference'
    };

    slide(2).code = {
        'P = Polyhedron([-1 -2; -0.3 0.4; 0.5 -0.4; 0.5 0.5],[1.7;1.4; -1.1; -0.6]);'
        'Q = Polyhedron([-1 1; 0.5 1; 1 -0.4]);'
        'cla; axis on; hold on; grid on;'
        'plot(P,''color'',''r''); plot(Q,''color'',''b'');'
        'plot(Polyhedron(''He'',[-1 -2 1.7]));'
        'plot(Polyhedron(''He'',[-0.3 0.4 1.4]));'
        'plot(Polyhedron(''He'',[0.5 -0.4 -1.1]));'
        'plot(Polyhedron(''He'',[0.5 0.5 -0.6]));'
        'plot(-1,1,''k.'',0.5, 1,''k.'',1,-0.4,''k.'',''MarkerSize'',20);'
        'axis([-4 2 -1 2]); hold off'
        'P.H'
        'Q.V'
    };
    slide(2).text = {
      'Construct two polyhedra in 2D.'
      'The red polyhedron P is given in H-representation:'
      '>> P = Polyhedron([-1 -2; -0.3 0.4; 0.5 -0.4; 0.5 0.5],[1.7;1.4; -1.1; -0.6]);'
      ''
      'Each hyperplane defining the polyhedron P is shown as a line and can be accessed using'
      '>> P.H '
      ''
      'The blue polyhedron Q is given in V-representation:'
      '>> Q = Polyhedron([-1 1; 0.5 1; 1 -0.4]);'
      ''
      'Each vertex defining the polyhedron Q is shown as a black dot and can be accessed using'
      '>> Q.V'
    };

    %extreme and hull
    slide(3).code={
        'P.V;'
        'Q.H;'
        'cla; axis on; hold on;'
        'plot(P,''color'',''r''); plot(Q,''color'',''b'');'
        'plot(Polyhedron(''He'',[-1 -2 1.7]));'
        'plot(Polyhedron(''He'',[-0.3 0.4 1.4]));'
        'plot(Polyhedron(''He'',[0.5 -0.4 -1.1]));'
        'plot(Polyhedron(''He'',[0.5 0.5 -0.6]));'
        'plot(-1,1,''k.'',0.5, 1,''k.'',1,-0.4,''k.'',''MarkerSize'',20);'
        'plot(Polyhedron(''He'',Q.H(1,:)));'
        'plot(Polyhedron(''He'',Q.H(2,:)));'
        'plot(Polyhedron(''He'',Q.H(3,:)));'
        'plot(P.V(1,1),P.V(1,2),''k.'',P.V(2,1),P.V(2,2),''k.'',P.V(3,1),P.V(3,2),''k.'',P.V(4,1),P.V(4,2),''k.'',''MarkerSize'',20);'
    };
    slide(3).text={    
        'It is always possible to go from V to H representation and vice versa.'
        'Please not that this operation is computationally very demanding and may'
        'consume time for larger input data.'
        ''
        'To compute V-representation of polyhedron P, just refer for V, R getters'
        '>> P.V'
        ''
        'To compute H-representation of polyhedron Q, refer to H, He getters or A, b, Ae, be'
        '>> Q.H'
        ''
        'Both polyhedra have now defined H- and V-representation.'
        'This can be verified with the method e.g. on the polyhedron P'
        '>> P.hasHRep'
        '>> P.hasVRep'
    };

    % minimal representation
    slide(4).code={
        'cla; axis off;'
        'G = Polyhedron(randn(20,2),rand(20,1))'
        'G.minHRep'
        'G.H'
        'H = Polyhedron(randn(20,2))'
        'H.minVRep'
        'H.V'        
    };
    slide(4).text={    
        'Redundant H-representation means, that the number of hyperplanes'
        'defining the polyhedron is not minimal. To compute the minimal'
        'H-representation of a polyheron, use the "minHRep" method.'        
        ''
        'Define a polyhedron with 20 random hyperplanes.'
        '>> G = Polyhedron(randn(20,2),rand(20,1))'
        ''
        'Compute its  minimal (nonredundant) H-representation'
        '>> G.minHRep'
        ''
        'You can see that the number of hyperplanes have been reduced by calling the H-property'
        '>> G.H'
        ''
        'Redundant V-representation means, that the number of extreme points'
        'and rays is not minimal. To compute the minimal V-representation use'
        '"minVRep" method.'
        ''
        'Define a polyhedron with 20 random vertices.'
        '>> H = Polyhedron(randn(20,2))'
        ''
        'Compute its minimal (nonredundant) V-representation'
        '>> H.minVRep'
        ''
        'The resulting number of vertices can be accessed via V-property'
        '>> H.V'
    };


%% intersection
    slide(5).code = {
        'cla; axis on; hold on;'
        'plot(P,''color'',''r''); plot(Q,''color'',''b'');'
        'T = Polyhedron([-3, -1; -3.5, 0; -2 2])'
        'plot(T,''color'',''g'');'
        'axis([-4 2 -1 2]); hold off'
        'L = P.intersect(T)'
        'L.isEmptySet'
        'K = Q.intersect(T)'
        'K.isEmptySet'
    };
    slide(5).text = {
        'We define a new polyhedron T in V-representation:'
        '>> T = Polyhedron([-3, -1; 1, 4; 0 -2]);'
        ''
        'The polyhedron T is shown in green and it intersects the polyhedron P.'
        'The intersection can be computed as'
        '>> L = P.intersect(T)'
        ''
        'The intersected set L is not empty and we can verify that using'
        '>> L.isEmptySet'
        ''
        'On the contrary, the intersection with the polyhedron Q is empty'
        '>> K = Q.intersect(T)'
        '>> K.isEmptySet'
    };

   % intersection of two unbounded polyhedra
   slide(6).code = {
       'P1 = Polyhedron(''lb'',[0;0])'
       'P2 = Polyhedron([8 -1;-3 -1],[0; 0])'
       'S = P1.intersect(P2)'
       'cla; axis on; hold on;'
       'plot(P1,''color'',''r'');'
       'plot(P2,''color'',''b'');'
       'plot(S,''color'',''pink'');'
       'axis([-1 2 -1 2]);'
       'S.isFullDim'
       'S.isEmptySet'
       'S.isBounded'
       };
   slide(6).text = {
     'Consider the following unbounded polyhedra.'
     ''
     'Unbounded polyhedron P1'
     '>> P1 = Polyhedron(''lb'',[0;0]);'
     ''
     'Unbounded polyhedron P2'
     '>> P2 = Polyhedron([8 -1;-3 -1],[0; 0])'
     ''
     'The intersection is a lower-dimensional polyhedron and is plotted in pink color'
     '>> S = P1.intersect(P2)'
     ''
     'We can ask for other properties of the set S'
     '>> S.isFullDim'
     '>> S.isEmptySet'
     '>> S.isBounded'
   };

    %% affine maps
    % scaling
    slide(7).code={
    'cla; axis on;'
    'R = Polyhedron(''lb'',[1, 1],''ub'',[3, 2])'
    'R1 = R*0.6'
    'hold on; plot(R,''color'',''red'');'
    'plot(R1,''color'',''green'');'
    'axis([0 4 0 4]);'
    };
    slide(7).text={
        'Polyhedra can be scaled, rotated and shifted. This is achieved using affine maps.'
        'Consider a rectangle'
        '>> R = Polyhedron(''lb'',[1, 1],''ub'',[3, 2]);'
        ''
        'Scaling is achieved by multiplying a polyhedron with a constant:'
        '>> R1 = R*0.6'
        ''
        'It does not matter if the multiplication is from the left or right.'
    };

    % rotate and scale
    slide(8).code = {
        'A = [1 2; 0.5 -1]'
        'R2 = A*R'
        'cla; hold on; plot(R,''color'',''red'');'
        'plot(R2,''color'',''green'');'
        'axis([0 8 -2 3]);'
    };
    slide(8).text = {
      'Rotation and scaling is achieved by multiplication with a matrix from the left.'
      'For instance, the following matrix'
      '>> A = [1 2;0.5 -1]'
      '>> R2 = A*R'
      'changes the polyhedron to a parallelogram R2.'
    };

    % shifting
    slide(9).code = {
        'R3 = R+[4;3];'
        'R4 = R-[4;3];'
        'cla; hold on;'
        'plot(R,''color'',''red'');'
        'plot(R3,''color'',''green'');'
        'plot(R4,''color'',''blue'');'
        'axis([-5 9 -3 6]);'
        };
    slide(9).text = {
      'Addition or subtraction of a polyhedron with a point is called shifting.';
      'Add the point [4;3] to the rectangle R'
      '>> R3 = R + [4;3]'
      ''
      'Subtract the point point [4;3]'
      '>> R4 = R - [4;3]'      
      ''
      'The rectangle R3 is shown in green and R4 in blue color.'
    };

    slide(10).code = {
        'cla; hold on;'
        'S1 = R.affineMap([1 0.5])'
        'S1.Dim'
        'S2 = R.affineMap([1 0.5; -0.5 2])'
        'S3 = R.affineMap([1 0.5; -0.5 2; 0.4 -1.3])'
        'plot(R,''color'',''red'');'
        'plot(S1,''color'',''green'');'
        'axis([1 4 -1 4]);'
    };
    slide(10).text = {
        'Depending on the transformation matrix A (n x d) three cases'
        'of affine transformation are distinguished for a polyhedron in dimension d:'
        ' If n < d then this is projection'
        ' If n  > d then this is a lifting'
        ' If n == d then this is rotation/skew.'
        ''
        'Example of projection, the set S is shown in green.'
        '>> S1 = R.affineMap([1 0.5])'
        ''
        'Note that in the projection operation the dimension has reduced to 1.'
        '>> S1.Dim'
        };

    % rotation/skew
    slide(11).code = {
        'cla; hold on;'
        'plot(R,''color'',''red'');'
        'plot(S2,''color'',''blue'');'
        'S2.Dim'
    };
    slide(11).text = {
       'Example of rotation/skew, the set S is shown in blue.'
       '>> S2 = R.affineMap([1 0.5; -0.5 2]);'
       ''
       'The dimension of the set is preserved'
       '>> S2.Dim'
    };

    % lifting
    slide(12).code = {
        'cla; hold on;'
        'plot(R,''color'',''red'');'
        'plot(S3,''color'',''magenta'');'
        'S3.Dim'
    };
    slide(12).text = {
        'Example of lifting.'
        'The set S is shown in magenta and the dimension has increased to 3.'        
        '>> S3 = R.affineMap([1 0.5; -0.5 2; 0.4 -1.3])'
        '>> S3.Dim'
    };


    %relational operators
    slide(13).code={
        'cla; view(2); hold on;'
        'plot(0.8*R,''color'',''red'');'
        'plot(1.2*R,''color'',''blue'');'
        'plot(0.8*R+[0.5;0.8],''color'',''green'');'
        'axis([0 4 0 3]);'
        '0.8*R <= 1.2*R'
        '0.8*R+[0.5;0.8] <= 1.2*R'
    };
    slide(13).text={
        'One can use a set of overloaded operators to perform comparisons on polyhedra.'
        ''
        'P2 == P1  Tests if two polyhedra are equal'
        'P2 <= P1  Tests if P2 is a subset of P1'
        'P2 >= P1  Tests if P1 is a subset of P2'
        ''
        'Is polyhedron 0.8*R (in red) a subset of 1.2*R (in blue)?'
        '>> 0.8*R <= 1.2*R'
        ''
        'Is polyhedron 0.8*R+[0.5;0.8] (in green) a subset of 1.2*R (in blue)?'
        '>> 0.8*R+[0.5;0.8] <= 1.2*R'
    };

    % convex hull
    slide(14).code = {
        'P(1) = Polyhedron([-1 -1;1 0;0 1; -1 0])'
        'P(2) = Polyhedron([-0.5 0; 2 0.5; 2 -0.5])'
        'H = PolyUnion(P).convexHull'
        'cla; hold on; plot(P);'
        'plot(H,''wire'',true,''linewidth'',3);'
        'axis([-2 3 -2 2]);'
    };
    slide(14).text = {
        'The convex hull of an union of polyhedra in the same dimension can be computed using the "convexHull" method.'
        'Put two polyhedra in an array'
        '>> P(1) = Polyhedron([-1 -1;1 0;0 1; -1 0])'
        '>> P(2) = Polyhedron([-0.5 0; 2 0.5; 2 -0.5])'
        ''
        'Compute the convex hull'
        '>> H = PolyUnion(P).convexHull '
    };

    %minkowski addition
    slide(15).code={
        'cla; hold on;'
        'P1 = Polyhedron(''H'',[0.5 0 1.9; 1.5 -2 -0.4;-1.7 0.8 0.4; -0.5 -0.4 -1.1]);'
        'P2 = Polyhedron(''lb'',[0; 0],''ub'',[0.5 1]);'
        'M = P1 + P2'        
        'plot(P1,''color'',''red'');'
        'plot(P2,''color'',''yellow'');'
        'plot(M,''wire'',true,''linewidth'',3);'
        'axis([-1 5 -1 10]);'        
    };
    slide(15).text={
        'Let us assume the two polytopes as plotted on your screen. P1 is plotted in red.'
        '>> P1 = Polyhedron(''H'',[0.5 0 1.9; 1.5 -2 -0.4;-1.7 0.8 0.4; -0.5 -0.4 -1.1]);'
        ''
        'Polyhedron P2 is plotted in yellow.'
        '>> P2 = Polyhedron(''lb'',[0; 0],''ub'',[0.5 1]);'
        'Then, the Minkowski sum of these polytopes can be calculated by typing:'
        ''
        '>> M = P1 + P2'
        ''
        'The Minkowski sum M is depicted in wired frame.'
    };

    %pontryagin difference
    slide(16).code={
        'cla;clc;'
        'P1 = Polyhedron(''lb'',[-1;-1],''ub'',[1;1]);'
        'P2 = Polyhedron(''lb'',[-0.2;-0.2],''ub'',[0.2;0.2]);'
        'N = P1-P2'
        'hold on;'
        'plot(P1,''color'',''red'');'
        'plot(P2,''color'',''yellow'');'
        'plot(N,''wire'',true,''linewidth'',3);'
        'axis([-2 2 -2 2]);'
    };
    slide(16).text={
        'Similarly, one can compute Pontryagin difference by using the overloaded minus (-) operator:'
        'Define polyhedron P1 (depicted in red)'
        '>> P1 = Polyhedron(''lb'',[1;1],''ub'',[1;1])'
        ''
        'and the polyhedron P2 (in yellow)'
        '>> P2 = Polyhedron(''lb'',-[0.2;0.2],''ub'',[0.2;0.2]);'
        ''
        'Compute the Pontryagin difference in wired frame.'
        '>> N = P1 - P2',
    };

    %pontryagin difference with sets
    slide(17).code={
        'cla;clc;'
        'P1 = Polyhedron([-4 2; -2 4; 0 2; 0 -2; -2 -4; -4 -2]);'
        'P2 = -P1;'
        'P3 = Polyhedron([eye(2); -eye(2)], [1;1;1;1]*0.5);'
        'S = [P1 P2]-P3'        
        'hold on;'
        'plot(P1,''color'',''red'');'
        'plot(P2,''color'',''yellow'');'
        'plot(P3,''color'',''blue'');'
        'plot(S,''wire'',true,''linewidth'',3);'
        'axis([-5 5 -5 5])'
        'hold on;'
    };
    slide(17).text={
        'Pontryagin difference also accepts Polyhedron arrays in the same dimension, as shown in this example:'
        ''
        'Polyhedron P1 is shown in red color.'
        '>> P1 = Polyhedron([-4 2; -2 4; 0 2; 0 -2; -2 -4; -4 -2]);'
        '>> P2 = -P1'
        'Mirrors polytope P1 around the origin and is shown in yellow.'
        ''
        'Polyhedron P3 is depicted with blue color.'
        '>> P3 = Polyhedron(''lb'',[-0.5; 0.5],''ub'',[0.5;0.5]);'
        '>> S = [P1 P2] - P3'
        ''
        'In this case, the solution S consists of a finite number of polytopes and it is shown in wired frame.'
    };

    %set difference
    slide(18).code={
        'clc;cla;'
        'P1=Polyhedron([-3 3; 0 3; -3 -3; 0 -3]);'
        'P2=Polyhedron([0 0; 5 0; 5 -3; 0 -3]);'
        'P3=Polyhedron([eye(2); -eye(2)], [1;1;1;1]);'
        'P4=Polyhedron([-1 -2; 2 -3; -1 -3]);'
        'hold on'
        'plot([P1,P2]);'
        'plot([P3,P4],''wire'',true);'
    };
    slide(18).text={
        'To show the computation of set differences, let us define following polytopes:'
        ''
        '>> P1=Polyhedron([-3 3; 0 3; -3 -3; 0 -3]);'
        '>> P2=Polyhedron([0 0; 5 0; 5 -3; 0 -3]);'
        '>> P3=Polyhedron(''lb'',[-1;-1],''ub'',[1;1]);'
        '>> P4=Polyhedron([-1 -2; 2 -3; -1 -3]);'
        'P1 and P2 are depicted in full color and P3 and P4 are drawn in wireframe.'
        ''
        'The set difference between set S1 = [P1 P2] and set S2 = [P3 P4] can be then computed as follows:'
        ''
        '>> S1 = [P1 P2]'
        '>> S2 = [P3 P4]'
        '>> Sdiff = S1 \ S2'
        ''
        'The set difference is displayed on the next slide...'
    };

    slide(19).code={
        'clc;cla; '
        'hold on'
        'S1 = [P1, P2];'
        'S2 = [P3, P4];'
        'Sdiff = S1\ S2;'
        'plot(Sdiff);'
    };
    slide(19).text={
        '...continuation of the previous slide.'
        ''
        '>> S1 = [P1 P2]'
        '>> S2 = [P3 P4]'
        '>> Sdiff = S1 \ S2'
    };

   % slicing the zonotope
   slide(20).code = {
       'clc; cla;'
       'B = Polyhedron(''lb'',[-1;-1;-1],''ub'',[1;1;1]);'
       'A = [1, -2, 0.3; 0.2 0 -5; -1, -2.4 0];'
       'v = [0;0;5];'
       'S = A*B + v;'
       'Sn = Polyhedron(''Ae'',[0 0 1],''be'',6);'
       'Sm = S.intersect(Sn);'
       'plot(S,''color'',[0.9 0.9 0.9],Sm,''color'',''darkblue'',''linewidth'',3,''wire'',true,''linestyle'',''--'');'
       'C = S.slice(3,6);'
       ''
       'hold on;'
       'plot(C,''color'',''lightgreen'');'
       'axis([-5 5 -10 10 0 10]);'
   };
   slide(20).text = {
    'Cut trough a three-dimensional polytope.'
    'Create a box B in 3D'
    '>> B = Polyhedron(''lb'',[-1;-1;-1],''ub'',[1;1;1]);'
    ''
    'Compute an affine map of the box B with a matrix A and shift it to a point [0;0;5].'
    '>> A = [1, -2, 0.3; 0.2 0 -5; -1, -2.4 0];'
    '>> v = [0;0;5];'
    '>> S = A*B + v;'
    ''
    'Create a cut through the polytope S by fixing the third dimension x3 = 6 that is shown in wired frame.'
    '>> C = S.slice(3,6);'
    ''
    'The resulting cut C is plotted in light green color.'
   };
 
%     % set difference in lower-dimensional polyhedron
%    slide(20).code = {
%        'B = Polyhedron(''lb'',[-1;-1],''ub'',[1;1]);'
%        'A = Polyhedron(''Ae'',[1 -1],''be'',0);'
%        'S = B\A;'
%        'cla; hold on;'
%        'plot([B,A]);'       
%        'plot(S,''wire'',true,''linestyle'',''--'',''linewidth'',3);'       
%        'axis([-2 2 -2 2]);'
%    };
%    slide(20).text = {
%     'Compute the set difference between a polytope and lower-dimensional polyhedron.'
%     'Create a box B'
%     '>> B = Polyhedron(''lb'',[-1;-1],''ub'',[1;1]);'
%     ''
%     'Affine set that intersects B'
%     '>> A = Polyhedron(''Ae'',[1 -1],''be'',0);'
%     ''
%     'The set difference results in two polytopes plotted in wired frame.'
%     '>> S = B \ A'
%         
%    };
  
    %chebyball
    slide(21).code={
        'cla; view(2); hold on;' 
        'T = Polyhedron([-1 -3; 1 -0.5; -0.5 3]);'
        'T.computeHRep'
        'data = T.chebyCenter;'
        'x = sdpvar(2,1);'        
        'S = YSet(x,norm(x-data.x,2)<data.r);'
        'plot(T,''color'',''yellow'');'
        'plot(S,''wire'',true);'
        'axis([-2 2 -4 4]);'
    };
    slide(21).text={
        'For polyhedra in H-representation it is possible to inscribe a Chebyshev ball.'
        'Consider the following polyhedron'
        '>> T = Polyhedron([-1 -3; 1 -0.5; -0.5 3]);'
        ''
        'Polyhedron must be in H-representation.'
        '>> T.computeHRep'
        ''
        'The Chebyshev ball is characterized by the center x and the radius r.'
        'The data of the ball can be obtained via'
        '>> data = T.chebyCenter'
    };

    %isinside
    slide(22).code={
        'cla; hold on'
        'T.contains([0.1; 0.2])'
        'plot(T,''color'',''yellow'');'
        'plot(0.1,0.2,''k.'',''MarkerSize'',20);'
        'axis([-2 2 -4 4]);'
        };
    slide(22).text={
        'To find out if a given point lies in the interior of some polytope, one can use the method "contains"',
        ''
        '>> T.contains([0.1; 0.2])'
        ''
        'Returns 1 (true) if a points (x1=0.1, x2=0.2) lies in interior of Polyhedron T.'
    };

    %end
    slide(23).code={
        'cla; axis on; hold on'
        'P1 = ExamplePoly.randHrep;'
        'P2 = ExamplePoly.randVrep;'
        'plot(P1,''color'',''blue'');'
        'plot(P2,''color'',''red'');'
        'methods(Polyhedron)'
        'axis([-5 5 -5 5]);'
    };
    slide(23).text={
        'Functions that generate random polyhedra:',
        ''
        'Random H-representation (shown in blue):'
        '>> P1 = ExamplePoly.randHrep'
        ''
        'Random V-representation (shown in red):'
        '>> P2 = ExamplePoly.randVrep'
        ''
        'To see all methods on Polyhedra type'
        '>> methods(Polyhedron)'
        ''
        'For further information consult the individual help files, e.g.'
        '>> help Polyhedron/contains'
        '>> doc Polyhedron/contains'
    };



end

end

