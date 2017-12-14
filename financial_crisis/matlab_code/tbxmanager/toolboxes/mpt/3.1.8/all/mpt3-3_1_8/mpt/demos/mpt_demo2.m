function slide=mpt_demo2
%
%  MPT_DEMO2: Tour through visualization capabilities of the toolbox 
%  ==================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      mpt_demo2
%    
%  
%  DESCRIPTION
%  -----------
%     Demonstration how to plot basic sets and explanation of various plotting
%  options.
%  
%  SEE ALSO
%  --------
%     mpt_demo1,  mpt_demo_sets1,  mpt_demo_sets2
%  

%  AUTHOR(s)
%  ---------
%     
%    
%   (c) 2003-2013  Michal Kvasnica: STU Bratislava
%   mailto:michal.kvasnica@stuba.sk 
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
  playshow mpt_demo2
else
    blue='b';
    red='r';
    slide(1).code={
        'clc; cla; axis([-4 4 -2 2]);',
        'axis([-4 4 -2 2]); grid off; axis off; text(-2.5,0.5,''Visualization tools'',''FontSize'',16);'
        };
    slide(1).text={
        'The Topic of this tutorial is to show some basic and advanced features of the MPT toolbox in terms of visualization.',
    };

    slide(2).code={
        'cla; axis on; grid on; hold on;'
        'axis([-2 2 -2 2]);'        
        'P = Polyhedron(''lb'',[-1;-1],''ub'',[1;1])'
        'plot(P);'
        'xlabel(''x_1''); ylabel(''x_2''); hold off; '
        };
    slide(2).text={
        'Let us first create a single Polyhedron:'
        ''
        '>> P = Polyhedron(''lb'',[-1;-1],''ub'',[1;1])'
        ''
        'The basic command for plotting is:'
        '>> plot(P)'
        ''
        'This function plots Polyhedron P in default color scheme "mpt".'
    };

    slide(3).code={
        'cla; hold on;'
        'plot(P,''color'',''b'');'
        'hold off;'
        };
    slide(3).text={
        'To plot a Polyhedron in a different color, use the following call:'
        ''
        '>> plot(P, ''color'',''b'')'
        ''
        'This will plot the Polyhedron in blue color.'
    };

    slide(4).code={
        'cla; hold on'
        'plot(P, ''color'',''lightgreen'');'
        'axis([-2 2 -2 2]);'
        'c = charToColor(''lightgreen'')'
     };
    slide(4).text={
        'The color can be specified using a string with the color name, i.e.'
        ''
        '>> plot(P, ''color'',''lightgreen'')'
        ''
        'Or using the numerical values of matlab RGB color scheme that varies between 0 and 1.'
        '>> plot(P, ''color'',[ 0.56 0.93 0.56]);'
        ''
        'The numerical value of a specific color can be obtained with the help of'
        '>> c = charToColor(''lightgreen'')'
    };


    slide(5).code={
        'cla; hold on;'
        'plot(P,''wire'',true);'
    };
    slide(5).text={
        'To draw a Polyhedron in wireframe (i.e. no color fill), pass the "wire" option as true'
        ''
        '>> plot(P,''wire'',true)'
    };

    slide(6).code={
        'cla; hold on;'
        'plot(P,''wire'',true,''linewidth'',3,''linestyle'',''--'');'
    };
    slide(6).text={
        'The line thickness as well as the line style can be changed.'
        'The additional options are passed as parameter-value pairs to the "plot" function.'
        ''
        'Plot the polyhedron in wired frame, in thick dashed line.'
        '>> plot(P,''wire'',true,''linewidth'',3,''linestyle'',''--'')'
    };


    slide(7).code={
        'cla; hold on;'
        'P1 = Polyhedron([-4 2; -2 4; 0 2; 0 -2; -2 -4; -4 -2]);'
        'P2 = Polyhedron([eye(2); -eye(2)], [1.5; 3; 0; 0.5]);'
        'P3 = Polyhedron([2.1 3.2;3.2 1; 3.2 -2.2;1.5 -1.1;1.5 2]);'
        'plot(P,P1,P2,P3);'
        'axis([-5 5 -5 5]);'
        };
    slide(7).text={
        'The plot function also accepts multiple Polyhedron arguments as shown by this example:'
        ''
        '>> P1 = Polyhedron([-4 2; -2 4; 0 2; 0 -2; -2 -4; -4 -2]);'
        '>> P2 = Polyhedron([eye(2); -eye(2)], [1.5; 3; 0; 0.5]);'
        '>> P3 = Polyhedron([2.1 3.2;3.2 1; 3.2 -2.2;1.5 -1.1;1.5 2]);'
        '>> plot(P, P1, P2, P3)'
        'or'
        '>> S = [P P1];'
        '>> plot(S, P2, P3);'
    };

    slide(8).code = {
        'cla; hold on;'
        'plot(P1,''color'',''green'',P2,''color'',''green'',''alpha'',0.5,P3,''color'',''green'',''alpha'',0.2);'
    };
    slide(8).text = {
        'Any visualizable property can be changed seperately for each polyhedron using parameter-value approach.'
        ''
        'Here we change the transparency by modifying the option "alpha" which can take values between 0 and 1.'
        '0 means completely transparent and 1 is opaque.'
        '>> plot(P1,''color'',''green'',P2,''color'',''green'',...'
                 '        ''alpha'',0.5,P3,''color'',''green'',''alpha'',0.2);'
    };

 
    slide(9).code={
        'cla; hold on;'
        'Q1 = Polyhedron([-1 8; -6 -2; -2.5 -0.6]);'
        'Q2 = Polyhedron(''V'',[-1 0;1 4; 2.3 1],''R'',[0.5 0.5]);'
        'plot(Q1, ''Color'',''olive'',''Marker'',''x'',''MarkerSize'',12,Q2,''Marker'',''o'',''MarkerSize'',16);'
        'axis([-8 12 -3 12]);'
        };
    slide(9).text={
        'Plot two sets in different markings and in different size.'
        ''
        'The set Q1 is a polytope.'
        '>> Q1 = Polyhedron([-1 8; -6 -2; -2.5 -0.6]);'
        ''
        'The set Q2 is unbounded polyhedron.'
        '>> Q2 = Polyhedron(''V'',[-1 0;1 4; 2.3 1],''R'',[0.5 0.5]);'
        ''
        '>> plot(Q1, ''Color'',''olive'',''Marker'',''x'',''MarkerSize'',12,Q2,''Marker'',''o'',''MarkerSize'',16)'
    };

    slide(10).code={
        'cla; hold on;'
        'x = sdpvar(2,1);'
        'S = YSet(x, [2*x(1)-0.4*x(2) <= 0.5; norm(2*x-[1;2],2)<=1]);'
        'plot(S,''color'',''blue'',''alpha'',0.5,''linewidth'',3,''linestyle'',''-.'');'
        'axis([-0.5 1 0 2]);'
        };
    slide(10).text={
        'The same options for plotting are valid for convex sets reprented using "YSet" object.'
        ''
        '>> x = sdpvar(2,1);'
        '>> S = YSet(x, [2*x(1)-0.4*x(2) <= 0.5; norm(2*x-[1;2],2)<=1]);'
        ''
        'This plots set S in blue, transparent and with edges in dash dotted style.'
        '>> plot(S,''color'',''blue'',''alpha'',0.5,''linewidth'',3,''linestyle'',''-.'');'
    };
 
    slide(11).code={
        'cla; hold on;'
        'P3D = Polyhedron(''lb'',[-1, -2, -3],''ub'',[1, 2, 3]);'
        'plot(P3D,''color'',''chocolate'');'
        'axis([-4 4 -4 4 -4 4]);'
        'xlabel(''x_1''); ylabel(''x_2''); zlabel(''x_3''); '
    };
    slide(11).text={
        'Plotting of Polyhedrons in higher dimensions:'
        ''
        'The MPT toolbox provides a possibility to visualize polyhedra in 3D. To show it, we create a rectangle:'
        ''
        '>> P3D = Polyhedron(''lb'',[-1, -2, -3],''ub'',[1, 2, 3]);'
        '>> plot(P3D,''color'',''chocolate'');'
    };

    slide(12).code={
        'cla; hold on; grid on;'
        'Q1 = Polyhedron(rand(8,3));'
        'Q2 = Polyhedron(rand(10,3));'
        'Q3 = Polyhedron(rand(9,3));'
        'Q4 = Polyhedron(rand(11,3));'
        'Q5 = Polyhedron(rand(13,3));'
        'S1 = [Q1 Q2 Q5]; S2 = [Q4 Q3];'
        'plot(S1,S2);'
        'axis([0 1 0 1 0 1]);'
    };
    slide(12).text={
        'Again, the plot command can handle multiple input arguments also in 3D, as shown by the following example:',
        ''
        '>> Q1 = Polyhedron(rand(8,3));'
        '>> Q2 = Polyhedron(rand(10,3));'
        '>> Q3 = Polyhedron(rand(9,3));'
        '>> Q4 = Polyhedron(rand(11,3));'
        '>> Q5 = Polyhedron(rand(13,3));'
        ''
        '>> S1 = [Q1 Q2 Q5]; S2 = [Q4 Q3];'
        '>> plot(S1, S2);'
        ''
        'Individual coloring can be specified by the user.'
    };

    slide(13).code={
        'cla; hold on;'
        'plot(Q4.slice(2,0.5));'
        'axis([0 1 0 1]);'
        'xlabel(''x_1''); ylabel(''x_2'');'
    };
    slide(13).text={
        'It is sometimes useful (and also necessary for dimensions > 3) to plot a section through a Polyhedron.'
        ''
        'This can be done by slicing the polyhedron.'
        'Plot a cut through Polyhedron along the second dimension at value x2=0.5'
        ''
        '>> plot(Q4.slice(2,0.5));'
        ''
        'If the values at which to cut are omitted, zero is assuming as default value.'
    };

    slide(14).code={
        'cla; hold on; view(2);'
        'P4D = Q1*Q2;'
        'plot(P4D.projection(1:2),''color'',''plum'',''alpha'',0.6);'
        'axis([0 1 0 1]);'
    };
    slide(14).text={
        'Other option with more dimensional polyhedra is plotting projections on given axis.'
        ''
        '>> P4D = P1*P2';
        ''
        'Plot the projection of the polyhedron P4D in dimensions 1 and 2 in "plum" color with 0.6 transparency.'
        '>> plot(P4D.projection(1:2),''color'',''plum'',''alpha'',0.6);'
        ''
        'Note that projection is a computationally demanding operation.'
    };
 
    % partition - colormap1
    slide(15).code={
        'cla; hold on; grid on;'        
        'Z = Polyhedron(rand(20,2));'
        'T = Z.triangulate; '
        'U = PolyUnion(T); '
        'U.plot(''colormap'',''colorcube'');'
    };
    slide(15).text={
        'The most interesting part in visualisation is plotting of partitions.'       
        'Here we can apply different colormaps to distinguish the plots.'
        ''
        'Create a partition in 2D by triangulating a randomly generated polytope.'
        '>> Z = Polyhedron(rand(20,2));'
        '>> T = Z.triangulate; '
        ''
        'Put the polyhedra to a PolyUnion object that is devoted to represent partitions.'
        '>> U = PolyUnion(T); '
        ''
        'Plot the partition in "colorcube" colormap"'
        '>> U.plot(''colormap'',''colorcube'');'
    };

    % partition - colormap2
    slide(16).code={
        'cla; hold on; grid on;'
        'U.plot(''colormap'',''prism'')'
    };
    slide(16).text={
        'Change the colormap and plot the partition in different scheme.'
        'For instance, now we select "prism" colormap.'
        ''
        '>> U.plot(''colormap'',''prism'');'
        ''
        'The default colormap is "mpt".'
    };

    % partition - colororder
    slide(17).code={
        'cla; hold on; grid on;'
        'U.plot(''colororder'',''random'')'
    };
    slide(17).text={
        'By default, the colors of each polyhedron are plotted in a fixed order.'
        'The fixed order of colors can be changed to random. This is achieved by:'
        ''
        '>> U.plot(''colororder'',''random'');'
        ''
        'For other help on plotting of sets, see the corresponding documentation.'
    };

    % the end slide
    slide(18).code={
        'cla; axis([-2 2 -2 2]); grid off; axis off; text(-0.5,0.5,''The end'',''FontSize'',16); title(''''); ',
    };
    slide(18).text={
        ''
    };

end
