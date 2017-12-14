function slide = mpt_demo_unions1
%
%  MPT_DEMO_UNIONS1: Demo that illustrates working with unions of convex sets 
%  ===========================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      mpt_demo_unions1
%    
%  
%  DESCRIPTION
%  -----------
%     Demo that illustrates working with unions of convex sets. It will be shown
%  how to 
%    
%     - construct Union objects and 
%     - add and remove sets from the Union object. 
%  
%  
%  SEE ALSO
%  --------
%     mpt_demo_unions2
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
    playshow mpt_demo_unions1
else
    %intro
    
    slide(1).code={
        'clc; cla; axis([-4 4 -2 2]);'
        'axis([-4 4 -2 2]); grid off; axis off; text(-4,0.5,''Basic manipulation with unions of convex sets'',''FontSize'',16);'};
    slide(1).text={
        'This tour presents the basics for working with union of convex sets.'
        ''
        '  - Construction of Union objects'
        '  - Adding and removing sets from the Union object'
        };
    
    %% Creating union
    % two YSets in same dimension
    slide(2).code={
        'x = sdpvar(2,1);'
        'triangle = YSet(x,[1.7*x(1)-1*x(2)<=5.7; -x(1)+4*x(2)<=7; -x(1)-1.6*x(2)<=-4.1])'
        'box = YSet(x,[ 0<=x(1)<=2; 0<=x(2)<=3])'
        'U = Union([box, triangle])'
        'axis on; grid on; hold on;'
        'plot(U,''colororder'',''random'');'
        'axis([-1 6 -1 4]);'
     };
    slide(2).text={
        'Union is defined a collection of finite number of convex sets.'
        ''
        'We can create unions out of any objects derived from ConvexSet class.'
        'For instance, YSet objects and Polyhedron objects can be contained in the same union object.'
        ''
        'Create a triangle and box sets as YSet objects.'
        '>> x = sdpvar(2,1);'
        '>> triangle = YSet(x,[1.7*x(1)-1*x(2)<=5.7; -x(1)+4*x(2)<=7; -x(1)-1.6*x(2)<=-4.1])'
        '>> box = YSet(x,[ 0<=x(1)<=2; 0<=x(2)<=3])'
        ''
        'Create an union of box and triangle by constructing Union object'
        '>> U = Union([box, triangle])'
     };
     % two polyhedra in different dimension
     slide(3).code={
        'P1 = Polyhedron(''lb'',[0;0;0],''ub'',[1;1;1])'
        'P2 = P1.projection([1:2]) + [3;-1]'
        'U1 = Union([P1,P2])'
        'cla; grid on; hold on;'
        'U1.plot(''colororder'',''random'');'
        'axis([-1 5 -2 2 0 1]);'
     };
    slide(3).text={
        'The sets contained in the Union object does not need to have the same dimension.'
        ''
        'Create array of two polyhedra in different dimension'
        '>> P(1) = Polyhedron(''lb'',[0;0;0],''ub'',[1;1;1])'
        '>> P(2) = P1.projection([1:2]) + [3;-1]'
        ''
        'Put the polyhedra to an Union object'
        '>> U1 = Union(P)'
    };

    % polyhedron and YSet in union
    slide(4).code = {
        'P3 = ExamplePoly.randVrep'
        'x = sdpvar(2,1);'
        'F = norm(x)<=1;'
        'C = YSet(x,F);'
        'U2 = Union(P3)'
        'U2.add(C)'
        'cla; hold on;'
        'U2.plot(''colororder'',''random''); view(2);'
        'axis([-5 5 -5 5]);'
    };
    slide(4).text = {
        'Here we show how to create union from objects that are of different class.'
        ''
        'Create random polyhedron'
        '>> P3 = ExamplePoly.randVrep'
        ''
        'Create a circle as YSet object'
        '>> x = sdpvar(2,1);'
        '>> F = norm(x)<=1;'
        '>> C = YSet(x,F);'
        ''
        'To create Union object first we put together objects of the same class.'
        '>> U2 = Union(P3)'
        ''
        'Secondly, we use the following method to add the set C to the Union.'
        '>> U2.add(C)'
    };
    
    % constructing union of triangular polyhedra
    slide(5).code = {
        'P = 10*Polyhedron(''V'',randn(17,2)) + [10;0]'
        'T = P.triangulate'
        'U = Union(T(1:2))'
        'cla; hold on;'
        'U.plot(''colororder'',''random'');'
        'axis([-50 50 -50 50]);'
     };
    slide(5).text = {
        'In this slide we show how to add and remove sets from the Union object.'
        ''
        'Create random polyhedron'
        '>> P = 10*Polyhedron(''V'',randn(17,2)) + [10;0]'
        ''
        'Triangulate the polyhedron to get an array of polyhedra.'
        '>> T = P.triangulate'
        ''
        'Construct the Union object U out of the first two polyhedra.'
        '>> U = Union(T(1:2))'
        '' 
        'The union is plotted on this slide.'
     };
    slide(6).code = {
    'U.add(T(3:end))'
    'cla; hold on;'
    'U.plot(''colororder'',''random'');'
    };
    slide(6).text = {
    'We can add the remaining polyhedra to union U with the "add" method.'
    '>> U.add(T(3:end))'
    };

    slide(7).code = {
      'U.Set'
      'U.remove([2,3])'
      'cla; hold on;'
      'U.plot(''colororder'',''random'');'
    };
    slide(7).text = {
    'We can also remove sets from the union using the "remove" method.'
    'Removal of sets is based on the index of sets as they are stored in the array.'
    'To see the array of stored sets, we need to access the "Set" property:'
    '>> U.Set'
    'To remove the set from the union, call "remove" method with the indices of sets to remove.'
    '>> U.remove([2,3])'
    };

    slide(8).code = {
      'U.add(C)'
      'cla; hold on;'
      'U.plot(''colororder'',''random'');'
    };
    slide(8).text = {
    'Or add the sets of the different class, e.g. YSet.'
    ''
    'Add the circle C from previous slides.'
    '>> U.add(C)'
    };
    slide(9).code = {
        'cla; axis off;'
        };
    slide(9).text = {
      'The stored sets can be accessed via referring to "Set" property.'
      '>> U.Set'
      ''
      'The total number of sets that are contained in the union can be found by'
      '>> U.Num'
    };
 
    
end

end
