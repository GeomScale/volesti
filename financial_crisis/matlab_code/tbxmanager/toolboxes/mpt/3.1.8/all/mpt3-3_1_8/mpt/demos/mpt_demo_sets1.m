function slide=mpt_demo_sets1
%
%  MPT_DEMO_SETS1: Demonstrates construction and basic properties of the Polyhedron
%  ================================================================================
%  object 
%  =======
%  
%  
%  SYNTAX
%  ------
%     
%      mpt_demo_sets1
%    
%  
%  DESCRIPTION
%  -----------
%     This demo presents the basic usage for the Polyhedron object. In particular
%  it is shown how to: 
%    
%     - create Polyhedron objects 
%     - access data stored in the Polyhedron object 
%     - query basic properties of the Polyhedron set 
%  
%  
%  SEE ALSO
%  --------
%     mpt_demo_sets2,  mpt_demo_sets3
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
  playshow mpt_demo_sets1
else
    %intro

    slide(1).code={
        'clc; cla; axis([-4 4 -2 2]);'
        'axis([-4 4 -2 2]); grid off; axis off; text(-4,0.5,''Basic manipulation with Polyhedra'',''FontSize'',16);'};
    slide(1).text={
        'This tour presents the basics for working with polyhedra. In particular:'
        ''
        '  - Creation of Polyhedron objects'
        '  - Accessing internal information stored in the Polyhedron object'
        '  - Querying basic properties of the Polyhedron set'
    };

    %% Creating polyhedra

    slide(2).code={
        '' };
    slide(2).text={
        'A Polyhedron is defined as an intersection of a finite number of '
        'hyperplanes (H-representation) or as an convex union of finite'
        'number of points (V-representation)'
        ''
        'In general, a Polyhedron object can be initialized in two ways:'
        '- providing the intersecting hyperplanes given as inequalities A*x<=b'
        '  and equalities Ae*x=b '
        '- providing the set of extreme points (vertices) V and rays R'
    };

    %H representation

    slide(3).code={
        'disp(''Polyhedron P1:'');'
        'P1=Polyhedron([eye(2); -eye(2)],[1;1;1;1])'
        'axis on; grid on; hold on; plot(P1,''color'',''r''); axis([-2 2 -2 2]);'
        'xlabel(''x_1''); ylabel(''x_2''); hold off; '
    };
    slide(3).text={
        'Let us now create a Polyhedron object using the H-representation:',
        ''
        '>> A = [eye(2); -eye(2)];'
        '>> b = [1; 1; 1; 1];'
        '>> P1 = Polyhedron(A, b)'
        ''
        'The A and b matrices essential state that:'
        'x1<=1, x2<=1, x1>=-1, x2>=-1'
        ''
        'The created Polyhedron object is plotted on your screen. Since the'
        'polyhedron is bounded, it is often referred to as polytope.'
    };

    %V representation

    slide(4).code={
        'fprintf(''\nPolytope P2:\n'');'
        'P2=Polyhedron([0.2 -0.2; 0.2 0.4; -0.3 0.4; -0.3 -0.5])'
        'cla; hold on; plot(P2,''color'',''b''); axis([-2 2 -2 2]);' 
        'xlabel(''x_1''); ylabel(''x_2''); hold off; '
     };
    slide(4).text={
        'Now we define a polytope in V-representation having the following vertices:'
        'V1 = [0.2 -0.2], V2=[0.2 0.4], V3=[-0.3 0.4], V4=[-0.3 -0.5]'
        ''
        'This can be done by typing:'
        '>> V = [0.2 -0.2; 0.2 0.4; -0.3 0.4; -0.3 -0.5]'
        '>> P2 = Polyhedron(V)'
        ''
        'The resulting polytope is plotted on your screen.'
    };

    % lb, ub
    slide(5).code={
        'P3 = Polyhedron(''lb'',-1,''ub'',2);'
        'cla; hold on; plot(P3,''color'',''g'',''linewidth'',3); axis([-2 2 -2 2]);'
        'P4 = Polyhedron(''lb'',[-2,-1.5],''ub'',[1,-0.5])'
        'plot(P4,''color'',''g''); hold off'
    };
    slide(5).text={
        'Polyhedra can be created using lower and upper bounds on the variables.'
        'For instance, the interval -1 <= x <= 2 determines a polyhedron in 1D.'
        'To create such polyhedron we specify "lb" and "ub" arguments:'
        '>> P3 = Polyhedron(''lb'',-1,''ub'',2)'
        ''
        'Polyhedron P4 in 2D:'
        '>> P4 = Polyhedron(''lb'',[-2,-1.5],''ub'',[1,-0.5])'
    };
    % equalities
    slide(6).code = {
        'P5 = Polyhedron(''Ae'',[1 0.5],''be'',0)'
        'P6 = Polyhedron(''A'',[-1 0.5 2; -3 -1.5 0.2; 0.4 -1.8 -1],''b'',[2 3 1],''ub'',[1 1 2],''Ae'',[1 2 0],''be'',0)'
        'cla; hold on; plot([P5, P6]); hold off'
    };
    slide(6).text = {
        'The polyhedron can contain equalities in the form Ae*x=be.'
        'The polyhedron containing only equalities is referred to as affine set'
        'and the polyhedron containing both equalities and inequalities is '
        'referred to as lower-dimensional polyhedron.'
        ''
        'Create an affine set x1=-0.5x2 in 2D'
        '>> P5 = Polyhedron(''Ae'',[1 0.5],''be'',0)'
        ''
        'Create a lower dimensional polyhedron in 3D with '
        'equality constrains x1 = -2*x2'
        '>> P6 = Polyhedron(''A'',[-1 0.5 2; -3 -1.5 0.2; 0.4 -1.8 -1],...'
        '''b'',[2,3,1],''ub'',[1,1,2],''Ae'',[1 2 0],''be'',0)'    
    };
    % unbounded in H
    slide(7).code = {
        'P7 = Polyhedron([1 -1;-2.5 1],[0.5;0.5])'
        'cla; hold on; view(2); plot(P7); axis([-2 2 -2 2]); hold off'
    };
    slide(7).text = {
        'Polyhedra can be unbounded. In the sequel we construct unbounded'
        'polyhedron using H-representation, x1-x2<0.5, -2.5x1+x2<0.5'
        ''
        '>> P7 = Polyhedron([1 -1;-2.5 1],[0.5;0.5])'       
    };
    % unbounded in V
    slide(8).code = {
        'v1 = [-1, 0]; v2 = [1, 0.5], v3 = [1, -0.5];'
        'V = [v1; v2; v3];'
        'R = [1, 0];'
        'P8 = Polyhedron(''V'',V,''R'',R);'
        'cla; hold on; plot(P8); axis([-2 7 -2 2]); hold off'
    };
    slide(8).text = {
        'Construct unbounded polyhedron in V-representation using vertices and rays.'
        'We use three points to specify the polyhedron'
        ''
        '>> v1 = [-1,0]; v2 = [1, 0.5], v3 = [1, -0.5]; '
        '>> V = [v1; v2; v3];'
        ''
        'and one ray that shows the direction in x1 axis in which the polyhedron is unbounded'
        '>> R = [1, 0];'
        ''
        'The polyhedron is constructed as'
        '>> P8 = Polyhedron(''V'',V,''R'',R)'
    };            
    % general constructor
    slide(9).code = {
        'A = randn(5,2); b = rand(5,1);'
        'P9 = Polyhedron(A, b);'
        'V = randn(5,2);'
        'P10 = Polyhedron(V);'
        'cla; hold on; plot([P9, P10]); hold off; axis([-5 5 -5 5]);'
    };
    slide(9).text = {
        'In general it is possible to construct a polyhedron object '
        'using parameter-value pairs. The parameters can be the matrices'
        'of the inequality constraints "A", "b", or equality constraints'
        '"Ae", "be", or lower/upper bounds "lb", "ub", or vertices "V",'
        'rays "R".'
        ''
        'For quick construction of polyhedron objects is always'
        'recommended to use the short syntax:'
        ''
        '>> A = randn(5,2); b = rand(5,1);'
        '>> P9 = Polyhedron(A, b);'
        ''
        'or'
        '>> V = randn(5,2);'
        '>> P10 = Polyhedron(V);'
    };
    
    % creating arrays
    slide(10).code={
        'P11 = Polyhedron([eye(2); -eye(2)],[1.5; 1.5; 0; 0])'
        'S1 = [P1 P2 P11];'
        'S2 = [P3, P4];'
        'cla; hold on; plot([S1,S2]); axis([-2 2 -2 2]);hold off;'
    };
    slide(10).text={
        'Polyhedra can be concatenated into arrays:'
        ''
        '>> P11 = Polyhedron([eye(2); -eye(2)], [1.5; 1.5; 0; 0])'
        '>> S1 = [P1 P2 P11];'
        '>> S2 = [P3; P4];'
    };

    %indexing of arrays
    slide(11).code={
        'axis off; cla; fprintf(''\nIndexing a Polyhedron array:\n'');'
        'S = [S1; S2];'
        'R1=S(1), fprintf(''\n'');'
        'R2=S(1:3), fprintf(''\n'');'
        'R3=S([3 2]), fprintf(''\n'');'
        'R4=S(end), fprintf(''\n'');'
        };
    slide(11).text={
        'Subindexing can be used to access individual polyhedra'
        'stored in an array S=[S1, S2]'
        ''
        '>> R1=S(1)'
        'Returns the first element of polyhedron array S.'
        ''
        '>> R2=S(1:3)'
        'Returns polyhedra stored at positions 1 to 3 in S.'
        ''
        '>> R3=S([3 2])'
        'Returns polyhedra stored at positions 3 and 2. It is '
        'the same as typing R3=[S(3) S(2)].'
        ''
        '>> R4=S(end)'
        'Returns the final element of S. It is identical to R4=S(length(S)).'
    };

    %% Accessing internal elements
    % H, He
    slide(12).code={
        'P = Polyhedron([eye(2); -eye(2)], [1;1;1;1]);'
        'A = P.A'
        'b = P.b'
        'Ae = P.Ae'
        'be = P.be'
        };
    slide(12).text={
        'The H-representation of a polytope can be accessed by '
        'referring to H, He-properties. The matrices H, He are composed as'
        'H = [A, b] and He = [Ae, be].'
        ''
        '>> P = Polyhedron(''A'',[eye(2); -eye(2)],''b'',[1;1;1;1],''Ae'',[0.2 -0.4],''be'', 0)'
        '>> H = P.H'
        '>> He = P.He'
        ''
        'It is also possible to extract individual matrices A, b, Ae, be:'
        '>> A = P.A'
        '>> b = P.b'
        '>> Ae = P.Ae'
        '>> be = P.be'
    };
    % V, R
    slide(13).code={
        'Q = Polyhedron([eye(2); -eye(2)], [1;1;1;1]);'
        'V = Q.V'
        'R = Q.R'
        };
    slide(13).text={
        'The V-representation of a polytope can be accessed by '
        'referring to V, R-properties. The vertices V are stored under'
        'V-property and rays under R-property'
        ''
        '>> Q = Polyhedron(''V'',[0.4 -3 5;-3 1.6 0.5; 2 -4 -0.7],''R'',[-0.4, 0.6, 0.5])'
        '>> V = Q.V'
        '>> R = Q.R'
    };
    % H-V rep
    slide(14).code={
        'P.hasHRep'
        'P.hasVRep'
        'Q.hasHRep'
        'Q.hasVRep'
        };
    slide(14).text={
        'We can ask if the polyhedron is stored under H- or V- representation'
        'using the following functions:'
        ''
        '>> P.hasHRep'
        '>> P.hasVRep'
        '>> Q.hasHRep'
        '>> Q.hasVRep'
    };
    % Dim
    slide(15).code = {
        'P.Dim'
        'Q.Dim'
    };
    slide(15).text = {
        'The dimension of the polyhedron can be accessed by referring to '
        'Dim property.'
        ''
        '>> P.Dim'
        '>> Q.Dim'
    };
    slide(16).code = {
        'P1 = Polyhedron(''A'',[1 -2],''b'',0.5,''Data'',struct(''alpha'',1,''beta'',-pi/3));'
        'P1.Data'
        'P1.Data.theta = {''name1'',''name2''};'
        'P1.Data.c = [-0.4, 0.3, 0.2, -0.3, -0.2, 0.5];'
        'P1.Data'
    };
    slide(16).text = {
        'It is possible to store any data with the Polyhedron object'
        'which is stored under "Data" property.'
        ''
        'One way of storing data is to put them as argument to "Data" parameter, i.e.'
        '>> P1 = Polyhedron(''A'',[1 -2],''b'',0.5,''Data'',struct(''alpha'',1,''beta'',-pi/3));'
        ''
        'Data can be accessed by referring to "Data" property'
        '>> P1.Data'
        ''
        'Second way of attaching any data to Polyhedron object'
        'is by direct assignment to "Data" property.'
        '>> P1.Data.theta = {''name1'',''name2''};'        
        '>> P1.Data.c = [-0.4, 0.3, 0.2, -0.3, -0.2, 0.5];'
    };

    %% Querying basic properties
    % 1D
    slide(17).code = {
        'R1 = Polyhedron([1 1],1);'
        'cla; axis on; hold on; plot(R1); hold off; axis([-2 2 -2 2]);'
        'R1.isEmptySet'
        'R1.isBounded'
        'R1.isFullDim'
        'R1.Dim'
    };
    slide(17).text = {
    'Create inequality set  x1+x2 <= 1'
    '>> R1 = Polyhedron([1 1],1)'
    ''
    'Investigate the properties of the set R1:'
    'Is the set empty?'
    '>> R1.isEmptySet'
    ''
    'Is the set bounded?'
    '>> R1.isBounded'
    ''
    'Is the set full-dimensional?'
    '>> R1.isFullDim'
    ''
    'What is the dimension of the set?'
    '>> R1.Dim'
    };

    % 2D
    slide(18).code = {
        'A = [1 1;-1 1];'
        'b = [1;1];'
        'H = [A,b];'
        'R2 = Polyhedron(''H'',H);'
        'cla; hold on; plot(R2); hold off; axis([-2 2 -2 2]);'
        'R2.isEmptySet'
        'R2.isBounded'
        'R2.isFullDim'
        'R2.Dim'        
    };
    slide(18).text = {
    'Create inequality set x1+x2<=1, -x1+x2<=1'
    '>> A = [1 1;-1 1];'
    '>> b = [1;1];'
    '>> H = [A,b];'
    '>> R2 = Polyhedron(''H'',H)'
    ''
    'Investigate the properties of the set R2:'
    'Is the set empty?'
    '>> R2.isEmptySet'
    ''
    'Is the set bounded?'
    '>> R2.isBounded'
    ''
    'Is the set full-dimensional?'
    '>> R2.isFullDim'
    ''
    'What is the dimension of the set?'
    '>> R2.Dim'        
    };

    % 2D lowdim
    slide(19).code = {
        'R3 = Polyhedron(''A'',[1 1;-1 1;2 -5],''b'',[1;1;3],''Ae'',[0.5 -3],''be'',0.5)'
        'cla; hold on; plot(R3); hold off; axis([-2 2 -2 2]);'
        'R3.isEmptySet'
        'R3.isBounded'
        'R3.isFullDim'
        'R3.Dim'        
    };
    slide(19).text = {
    'Create set based on inequalities x1+x2<=1, -x1+x2<=1, x1-x2<=2 and equalities 0.5*x1-3*x2=0.5'
    '>> R3 = Polyhedron(''A'',[1 1;-1 1;2 -5],''b'',[1;1;3],''Ae'',[0.5 -3],''be'',0.5)'
    ''
    'Investigate the properties of the set R3:'
    'Is the set empty?'
    '>> R3.isEmptySet'
    ''
    'Is the set bounded?'
    '>> R3.isBounded'
    ''
    'Is the set full-dimensional?'
    '>> R3.isFullDim'
    ''
    'What is the dimension of the set?'
    '>> R3.Dim'        
    };
    
    % 3D
    slide(20).code = {
        'R4 = Polyhedron([1 -1 1; -1 -1 -1; -1 1 -1])'
        'cla; hold on; plot(R4); hold off; axis([-2 2 -2 2]);'
        'R4.isEmptySet'
        'R4.isBounded'
        'R4.isFullDim'
        'R4.Dim'        
    };
    slide(20).text = {
    'Create the set out of the vertice V, v1=[1,-1,1], v2=[-1,-1,-1], v3=[-1, 1, -1].'
    '>> R4 = Polyhedron([1 -1 1; -1 -1 -1; -1 1 -1])'
    ''
    'Investigate the properties of the set R4:'
    'Is the set empty?'
    '>> R4.isEmptySet'
    ''
    'Is the set bounded?'
    '>> R4.isBounded'
    ''
    'Is the set full-dimensional?'
    '>> R4.isFullDim'
    ''
    'What is the dimension of the set?'
    '>> R4.Dim'        
    };
 
    % 3D unbounded
    slide(21).code = {
        'R5 = Polyhedron(''V'',[1 -1 1; -1 -1 -1; -1 1 -1],''R'',[0,0.5,0])'
        'cla; hold on; plot(R5); hold off; axis([-2 2 -2 7]);'
        'R5.isEmptySet'
        'R5.isBounded'
        'R5.isFullDim'
        'R5.Dim'        
    };
    slide(21).text = {
    'Create the set out of the vertice V, v1=[1,-1,1], v2=[-1,-1,-1], v3=[-1, 1, -1], and ray [0,0.5,0]'
    '>> R5 = Polyhedron(''V'',[1 -1 1; -1 -1 -1; -1 1 -1],''R'',[0,0.5,0])'
    ''
    'Investigate the properties of the set R5:'
    'Is the set empty?'
    '>> R5.isEmptySet'
    ''
    'Is the set bounded?'
    '>> R5.isBounded'
    ''
    'Is the set full-dimensional?'
    '>> R5.isFullDim'
    ''
    'What is the dimension of the set?'
    '>> R5.Dim'        
    };


end

end
