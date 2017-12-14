function slide = mpt_demo_unions2
%
%  MPT_DEMO_UNIONS2: Demo that illustrates working with unions of polyhedra 
%  =========================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      mpt_demo_unions2
%    
%  
%  DESCRIPTION
%  -----------
%     Demo that illustrates working with unions of polyhedra. It will be shown how
%  to: 
%    
%     - construct PolyUnion objects, 
%     - query basic properties from the union, and 
%     - perform geometric operations with unions of polyhedra. 
%  
%  
%  SEE ALSO
%  --------
%     mpt_demo_unions1
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
    playshow mpt_demo_unions2
else
    %intro
    
    slide(1).code={
        'clc; cla; axis([-4 4 -2 2]);'
        'axis([-4 4 -2 2]); grid off; axis off; text(-4,0.5,''Basic manipulation with unions of polyhedra'',''FontSize'',16);'};
    slide(1).text={
        'This tour presents the basics for working with union polyhedra.'
        ''
        '  - Construction of PolyUnion objects'
        '  - Querying basic properties of the union.'
        '  - Operations with unions of polyhedra'
        };
    
    
    %% union of two polyhedra
    slide(2).code = {
        'P1 = Polyhedron(''V'',[-0.25, 0.8; 0.6, -0.1; 0.2, 0.7; 0.3, -0.6; -0.6, 0.2])'
        'P2 = Polyhedron(''V'',[-0.8, -0.3; 2.2, -0.4; 0.6, -1.3; 0.9, -1.3; 1.12, 0.8;-1 -0.9])'
        'U1 = PolyUnion([P1,P2])'
        'cla; axis on; grid on; hold on;'
        'U1.plot'
        'axis([-3 3 -3 3]);'
    };
    slide(2).text = {
        'The PolyUnion objects represents union of polyhedra.'
        'Only the polyhedra in the same dimension can be put together to form PolyUnion object.'
        '>> P1 = Polyhedron(''V'',[-0.25, 0.8; 0.6, -0.1; 0.2, 0.7; 0.3, -0.6; -0.6, 0.2])'
        '>> P2 = Polyhedron(''V'',[-0.8, -0.3; 2.2, -0.4; 0.6, -1.3; 0.9, -1.3; 1.12, 0.8;-1 -0.9])'
        ''
        'Put the polyhedra to an PolyUnion object'
        '>> U1 = PolyUnion([P1,P2])'
    };
 
    %% adding polyhedra to PolyUnion
    slide(3).code = {
        'P3 = Polyhedron(''lb'',[-3;-3],''ub'',[-1;0])'
        'U1.add(P3)'
        'cla; hold on;'
        'U1.plot;'
        'axis([-5 5 -5 5]);'
    };
    slide(3).text = {
        'The PolyUnion object was created without specifying any properties, '
        'just by converting an array of polyhedra in the same dimension to PolyUnion object.'
        'This union of polyhedra is very general and we can easily add and remove sets from the union.'
        ''
        'Create a polyhedron P3 that is to be added to the union U1'
        '>> P3 = Polyhedron(''lb'',[-3;-3],''ub'',[-1;0])'
        ''
        'We can add the polyhedron P3 to the PolyUnion using the following approach'
        '>> U1.add(P3)'
     };
    slide(4).code = {
        'U1.add([ExamplePoly.randVrep, ExamplePoly.randVrep])'
        'cla; hold on;'
        'U1.plot;'
    };
    slide(4).text  = {
        'We can add multiple polyhedra to the union.'
        'For instance, in the following we add two randomly generated polyhedra.'
        '>> U1.add([ExamplePoly.randHrep, ExamplePoly.randVrep])'
    };
    % removing polyhedra
    slide(5).code = {
        'U1.Set'
        'U1.Num'
        'U1.remove([4,5])'
        'cla; hold on;'
        'U1.plot;'
    };
    slide(5).text = {
     'The sets can be removed as well from the PolyUnion.'
     'Removing is based on the index of the set as it is stored in array.'
     'The polyhedra are stored under "Set" property of the PolyUnion object.'
     '>> U1.Set'
     ''
     'To find out the number of sets in the union, type'
     '>> U1.Num'
     ''
     'Based on the number of sets we can remove the sets from the union.'
     'For instance, to remove the last two added polyhedra type'
     '>> U1.remove([4,5])'
    };
    
    %% properties of polyunion
    slide(6).code = {
        'cla; axis off;'
        'U1.isConvex'
        'U1.isOverlapping'
        'U1.isFullDim'
        'U1.isConnected'
        'U1.isBounded'        
    };
    slide(6).text = { 
        'We can investigate the properties of the polyunion similarly as for polyhedra.'
        'Is the union convex?'
        '>> U1.isConvex'
        ''
        'Has the union overlapping regions?'
        '>> U1.isOverlapping'
        ''
        'Is the union built only from full-dimensinal polyhedra?'
        '>> U1.isFullDim'
        ''
        'Is the union connected?'
        '>> U1.isConnected'
        ''
        'Is the union build only from bounded polyhedra?'
        '>> U1.isBounded'
        ''
        'Note that it is computationally expensive to determined these properties,'
        'so use these methods with care.'
    };

    %% create polyunion with certain properties
    slide(7).code = {
        'P4 = Polyhedron(''ub'',-1)'
        'P5 = Polyhedron(''lb'',-1,''ub'',1)'
        'P6 = Polyhedron(''lb'',1,''ub'',2)'
        'U2 = PolyUnion(''Set'',[P4,P5,P6],''Overlaps'',false,''Convex'',true)'
        'cla; axis on; grid on; hold on;'
        'U2.plot(''LineWidth'',3);'
        'axis([-3 3 -1 1]);'
        };
    slide(7).text = {
        'To avoid expensive computation for querying properties of the union of polyhedra,'
        'it is possible to assign these properties at the construction.'
        ''
        'For instance, create polyhedra that do not overlap and form convex union.'
        '>> P4 = Polyhedron(''ub'',-1)'
        '>> P5 = Polyhedron(''lb'',-1,''ub'',1)'
        '>> P6 = Polyhedron(''lb'',1,''ub'',2)'
         ''
        'Create PolyUnion object with these particular properties.'
        '>> U2 = PolyUnion(''Set'',[P4,P5,P6],''Overlaps'',false,''Convex'',true)'
    };
    
    slide(8).code = {
        'try, U2.add(Polyhedron(''lb'',0)); catch me; disp(me.message); end;'
        'U2.add(Polyhedron(''lb'',2))'
        'cla; hold on;'
        'U2.plot(''LineWidth'',3);'
    };
    slide(8).text = {
        'We cannot add polyhedra to the union that violate the assigned properties.'
        '>> try'
        '      U2.add(Polyhedron(''lb'',0))'
        '   catch me'
        '      disp(me.message);'
        '   end'
        ''
        'We can add only polyhedra that do not violate the properties.'
        '>> U2.add(Polyhedron(''lb'',2))'
    };
    
    %% operations on polyunions
    slide(9).code = {
        'G = Polyhedron(''V'',randn(26,2)) + [5;0]'
        'T = G.triangulate'
        'U3 = PolyUnion(''Set'',T,''Convex'',true,''Overlaps'',false,''FullDim'',true,''Bounded'',true,''Connected'',true)'
        'clc; cla; hold on;'
        'U3.plot;'
        'axis([0 10 -5 5]);'
        'U3.Internal'
        };
    slide(9).text = {
        'The following properties can be assigned at the construction'
        ' - Convex (is the union of polyhedra convex?)'
        ' - Overlaps (do the polyhedra overlap?)'
        ' - Bounded (are the polyhedra bounded?)'
        ' - FullDim (are the polyhedra full-dimensional?)'
        ' - Connected (are the polyhedra connected?)'
        ''
        'Now we show a case where all of these properties can be assigned at the construction.'
        'Generate random bounded polyhedron.'
        '>> G = Polyhedron(''V'',randn(26,2)) + [5;0]'
        ''
        'Triangulate the set.'
        '>> T = G.triangulate'
        ''
        'Create a polyunion object with certain properties that we know about T.'
        '>> U3 = PolyUnion(''Set'',T,''Convex'',true,''Overlaps'',false,''FullDim'',true,''Bounded'',true,''Connected'',true)'
        ''
        'The assigned properties can be accessed by referring to "Internal" property of the PolyUnion object,i.e.'
        '>> U3.Internal'
    };
    
    % merging
    slide(10).code = {
      'U3.merge;'
      'cla; hold on;'
      'U3.plot;'
    };
    slide(10).text = {
        'Sometimes the union of polyhedra can be simplified to another union that contains fewer number of regions.'
        'One approach is referred to as merging and can be invoked via following method'
        '>> U3.merge'
        ''
        'Note that this operation is computationally demanding and not all unions can be simplified.'
    };
     
    slide(11).code = {
        'P(1) = Polyhedron(''V'',[ -1.5,-0.5; 0.4, -1.3; 1.3, 0.3; -2, 0; -0.2,0.8; 2.5, -1.9]);'
        'P(2) = 2*P(1);'
        'P(3) = P(2)+[1;2];'
        'U4 = PolyUnion(P)'
        'U4.reduce'
        'cla; hold on;'
        'U4.plot;'
        'axis([-5 7 -5 5]);'
    };
    slide(11).text = {
        'Another way of simplifying the union of polyheda is using the "reduce" method.'
        'This method eliminates the regions that are completely covered by other regions.'
        'To demonstrate this function, create a polyunion object without any properties.'
        '>> P(1) = Polyhedron(''V'',[ -1.5,-0.5; 0.4, -1.3; 1.3, 0.3; -2, 0; -0.2,0.8; 2.5, -1.9]);'
        '>> P(2) = 2*P(1);'
        '>> P(3) = P(2)+[1;2];'
        '>> U4 = PolyUnion(P)'
        ''
        'Reduce the union by'    
        '>> U4.reduce'
     };
    % bounding box
    slide(12).code = {
        'B = U4.outerApprox'
        'B.plot(''wire'',true,''LineWidth'',2,''LineStyle'',''--'');'
     };
    slide(12).text = {
        'Very useful method for working with unions of polyhedra is to find outer bounding box approximation.'
        'This can be achieved using "outerApprox" method:'
        '>> B = U4.outerApprox'
        ''
        'The outer approximation using the bounding box is shown in dashed wired frame.'
    };
    % convexHull
    slide(13).code = {
        'H = U4.convexHull'
        'H.plot(''wire'',true,''LineWidth'',2,''LineStyle'',''-.'');'
    };
    slide(13).text = {
        'Tighter approximation of the union can be achieved using "convexHull" method.'
        'The convex hull is defined as the smallest convex set that approximates the union from outsied.'
        'To compute the convex hull of the union use'
        '>> H = U4.convexHull'
        ''
        'The outer approximation using the convex hull is shown in dash-dotted wired frame.'
     };

end
end
