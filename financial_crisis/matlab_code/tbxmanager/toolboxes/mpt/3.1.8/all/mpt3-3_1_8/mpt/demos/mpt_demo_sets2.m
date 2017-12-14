function slide = mpt_demo_sets2
%
%  MPT_DEMO_SETS2: Construction and basic properties of sets created in Yalmip 
%  ============================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      mpt_demo_sets2
%    
%  
%  DESCRIPTION
%  -----------
%     Slideshow on how to construct and work with YSet objects. In particular, it
%  is demonstrated how to: 
%    
%     - Create Yalmip YSet objects. 
%     - Access internal information stored in the YSet object. 
%     - Query basic properties of the YSet object. 
%  
%  
%  SEE ALSO
%  --------
%     mpt_demo_sets1,  mpt_demo_sets3
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
    playshow mpt_demo_sets2
else
    %intro
    
    slide(1).code={
        'clc; cla; axis([-4 4 -2 2]);'
        'axis([-4 4 -2 2]);axis off; text(-4,0.5,''Basic manipulation with convex sets'',''FontSize'',16);'
        };
    slide(1).text={
        'This tour presents basic usage with general convex sets, in particular'
        ''
        '  - Creation of Yalmip "YSet" objects'
        '  - Accessing internal information stored in the YSet object'
        '  - Querying basic properties of the YSet object'
        };
    
    % 2D set
    slide(2).code = {
        'x = sdpvar(2,1);'
        'F1 = [-1 <= x <= 1];'
        'S1 = YSet(x,F1);'
        'axis on; grid on; hold on; plot(S1); axis([-2 2 -2 2]); hold off;'
        'xlabel(''x_1''); ylabel(''x_2'');'
    };
    slide(2).text = {
        'Construction general convex sets follows the YALMIP syntax. '
        'First, one has to define variables that define the set.'
        'For instance, continuous variable with 2 elements in a column'
        'can be declared as'
        ''
        '>> x = sdpvar(2,1);'                
        'Once the variable has been declared, the set can be defined by'
        'concatenating various constraints, e.g. specify lower and upper'
        'bounds'
        ''
        '>> F1 = [-1 <= x <= 1];'
        ''
        'This type of constraints are typical for YALMIP and MPT3 interfaces'
        'them using YSet objects. To construct YSet object use the following'
        'syntax'       
        ''
        '>> S1 = YSet(x,F1)'
    };

    % 2D set
    slide(3).code = {
        'x = sdpvar(2,1);'
        'F2 = [ 1*x(1)-0.9*x(2) <= 2.5; 0.6*x(1)+1.3*x(2) >= -2.3; 0.2*x(1)+0.9*x(2) <= 4.5];'
        'S2 = YSet(x,F2);'
        'cla; axis on; hold on; plot(S2); axis([-30 10 -4 12]); hold off;'
        'xlabel(''x_1''); ylabel(''x_2'');'
    };
    slide(3).text = {
        'Here we construct a polyhedron as YSet object by concatenating few inequalities.'
        ''
        'Define the variable first'
        '>> x = sdpvar(2,1);'                
        ''
        'Now we create linear constraint set'
        '>> F2 = [ 1*x(1)-0.9*x(2) <= 2.5; 0.6*x(1)+1.3*x(2) >= -2.3; 0.2*x(1)+0.9*x(2) <= 4.5]'
        ''
        'and construct YSet object as'
        ''
        '>> S2 = YSet(x,F2)'
    };

    % 2D quadratic function
    slide(4).code = {
        'x = sdpvar(2,1);'
        'F3 = [ x(2)>= 1.5*(x(1)+5).^2-5*x(1)-20; -3.5*x(1)+2.4*x(2) <= 28];'
        'S3 = YSet(x,F3);'
        'cla; axis on; hold on; plot(S3,''grid'',100); axis([-6 0 -1 12]); hold off;'
        'xlabel(''x_1''); ylabel(''x_2'');'
    };
    slide(4).text = {
        'Using YSet object we can plot more general sets. For instance, consider'
        'function f(x) = 1.5*(x+5)^2-5*x-20 which is a quadratic function. '
        'We can create a convex set by intersecting this function with a hyperplane.'
        ''        
        'First we create the nonlinear constraint'
        '>> x = sdpvar(2,1);'  
        '>> F3 = [ x(2)>= 1.5*(x(1)+5).^2-5*x(1)-20;'
        ''
        'and add one hyperplane'
        '>> F3 = [F3; -3.5*x(1)+2.4*x(2) <= 28]'
        '>> S3 = YSet(x,F3)'
        ''
        'Displaying of the general sets takes more time than with '
        'polyhedra because the sets are sampled so be patient.'
    };

    % array of two sets
    slide(5).code = {
        'y = sdpvar(1);'
        'S4 = YSet(y, [-5 <= y <= 1]);'
        'x = sdpvar(2,1);'
        'S5 = YSet(x, [x(2).^2+x(1)+2.1*x(2)<=1; -x(1)+0.5*x(2).^2-2.1*x(2)<=1]);'
        'cla; axis on; hold on; plot([S4,S5],''grid'',100); axis([-6 6 -2 2]); hold off;'        
        };
    slide(5).text = {
      'YSet objects can be grouped to arrays. For instance, we can group 1D object'
      'with 2D object.'
      ''
      'Create 1D set:'
      '>> y = sdpvar(1);'
      '>> S4 = YSet(y, [-5 <= y <= 1])'
      ''
      'Create 2D set:'
      'x = sdpvar(2,1);'
      '>> S5 = YSet(x, [x(2).^2+x(1)+2.1*x(2)<=1; -x(1)+0.5*x(2).^2-2.1*x(2)<=1])'
      ''
      'Put the sets S1, S2 in an array'
      '>> S = [S4, S5];'
    };
    
    slide(6).code = {
        'cla; axis off;'
        'X = sdpvar(2);'
        'A = [1, 2; -0.5 1.6];'
        'F4 = [A''*X + X*A <= -eye(2)];'
        'S6 = YSet(X(:), F4);'
    };
    slide(6).text = {
        'It is possible to construct a convex set using linear matrix inequalities.'
        'For instance, the inequality A''*X + X*A <= -I'
        'can be written as'
        ''
        '>> X = sdpvar(2)'
        '>> A = [1, 2; -0.5 1.6];'
        ''
        'Create the constraint set'
        '>> F4 = [A''*X + X*A <= -eye(2)]'
        ''
        'To construct YSet object we need to pass the variables as vector, i.e.'
        '>> S6 = YSet(X(:), F4) '
        ''
        'As the object is in the dimension 4, it cannot be plotted.'
    };

    slide(7).code = {
    'x = sdpvar(2,1);';
    'F5 = [ 3*norm(x-[1;2])<=5 ];'
    'S7 = YSet(x,F5);'
    'cla; axis on; hold on; S7.plot; axis([-1 3 0 4]); hold off; '
    };

    slide(7).text = {
    'Create a circle centered at the point [1;2].'
    ''
    '>> x = sdpvar(2,1);'
    '>> F5 = [ 3*norm(x-[1;2])<=5 ]'
    '>> S7 = YSet(x,F5)'            
    };

    % intersect the circle with the box
    slide(8).code = {
        'x = sdpvar(2,1);'
        'box = ( [0;0.5] <= x <= [1; 2] );'
        'circle = ( 2*norm(x-1) <= 1 );'
        'S8 = YSet(x,box + circle);'
        'cla; axis on; hold on; S8.plot; axis([0.2 1.2 0 2]); hold off;'
    };

    slide(8).text = {
        'We can concatenate multiple constraint sets to one set.'
        ''
        'For instance, define a box in 2D:'
        '>> x = sdpvar(2,1);'
        '>> box = ( [0;0.5] <= x <= [1; 2] )'
        ''
        'Define circle as '
        '>> circle = ( 2*norm(x-1) <= 1 )'
        ''
        'Create a set that contains both box and circle constraints.'
        '>> S8 = YSet(x,box + circle)'
    };

    % MPC constraints
    slide(9).code = {
        'A = [2 -1;1 0];'
        'B = [1;0];'
        'nu = 1;'
        ''
        'Q = eye(2);'
        'R = 1;'
        'N = 2;'
        'x0 = [2;1];'
        ''
        'u = sdpvar(nu,N);'
        'constraints = [];'
        'objective = 0;'
        'x = x0;'
        ''
        'for k = 1:N,'
        '       x = A*x + B*u(k);'
        '       objective = objective + norm(Q*x,1) + norm(R*u(k),1);'
        '       constraints = [constraints, -1 <= u(k)<= 1, -5<=x<=5];'
        'end;'
        ''
        'S9 = YSet(u,constraints);'
        'cla; axis on; hold on; S9.plot(''color'',''olive''); axis([-1.2 1.2 -1.2 1.2]); hold off;'
    };

    slide(9).text = {
        'This example shows how to model sets arising in MPC.'
        'Model data:'
        '>> A = [2 -1;1 0];'
        '>> B = [1;0];'
        ''
        'Number of inputs'
        '>> nu = 1;'
        ''
        'MPC tuning parameters'
        '>> Q = eye(2);'
        '>> R = 1;'
        '>> N = 2;'
        ''
        'Initial state'
        '>> x0 = [2;1];'
        ''
        'Setup of optimization problem'        
        '>> u = sdpvar(nu,N);'
        '>> constraints = [];'
        '>> objective = 0;'
        '>> x = x0;'
        '>> for k = 1:N'
        '       x = A*x + B*u(k);'
        '       objective = objective + norm(Q*x,1) + norm(R*u(k),1);'
        '       constraints = [constraints, -1 <= u(k)<= 1, -5<=x<=5];'
        '   end'
        ''
        'Create the convex arising from MPC constraints'
        '>> S9 = YSet(u,constraints)'            
    };


    %% Access internal properties
    slide(10).code = {
        'cla; axis off;'
        'x = sdpvar(6,1);'
        'S = YSet(x, x >= 0);'
        'S.vars'
        'S.constraints'
    };
    slide(10).text = {
      'We can access internal data stored under YSet objects using "dot" syntax.'
      'The data are stored under "vars" and "constraints" properties.'
      ''
      'Construct simple YSet object:'
      '>> x = sdpvar(6,1)'
      '>> S = YSet(x, x >= 0)'
      ''
      'To access its elements, refer to the following properties:'
      '>> S.vars'
      '>> S.constraints'           
    };

    slide(11).code = {
        'cla; axis off;'
        'x = sdpvar(3,1);'
        'S = YSet(x, x.^2 <= 1);'
        'S.Data.a = [1, 2];'
        'S.Data.b = {[0.1 -0.5],[0.12 -0.51],[0.11 -0.53]};'
        'S.Data.name = ''friction'';'
    };
    slide(11).text = {
      'It is possible to store with the YSet as well as with Polyhedron object arbitrary data.'
      'For this purpose there is "Data" property.'
      ''
      'To construct an object follow the default procedure'
      '>> x = sdpvar(3,1)'
      '>> S = YSet(x, x.^2 <= 1)'
      ''
      'Any data can be assigned to "Data" property after the object has been created'
      '>> S.Data.a = [1, 2]'
      '>> S.Data.b = {[0.1 -0.5],[0.12 -0.51],[0.11 -0.53]}'
      '>> S.Data.name = ''friction'''
    };

    
    %% properties of the Set
    slide(12).code = {
        'x = sdpvar(10,1);'
        'Y1 = YSet(x, x>=0);'
        'Y1.isBounded'
    };
    slide(12).text  = {
        'We can investigate basic properties of the convex set such as boundedness and emptyness.'
        'For this purpose there are two functions "isBounded" and "isEmptySet".'
        'For instance the set'
        ''
        '>> x = sdpvar(10,1);'
        '>> Y1 = YSet(x, x>=0);'
        ''
        'is obviously unbounded and we can query this information using'
        '>> Y1.isBounded'
    };
    slide(13).code = {
        'x = sdpvar(2,1);'
        'F = [ randn(12,2)*x <= rand(12,1) ];'
        'Y2 = YSet(x,F);'
        'cla; axis on; hold on; Y2.plot; hold off;'
        'Y2.isEmptySet'
    };
    slide(13).text = {
        'For sets that contain many constraints it is useful to know if this set is empty or not.'
        'Assume that we have a set formed by 12 inequalities:'
        ''
        '>> x = sdpvar(2,1);'
        '>> F = [ randn(12,2)*x <= rand(12,1) ]'
        '>> Y2 = YSet(x,F);'
        ''
        'To find out if the set Y2 is empty or not, use'        
        '>> Y2.isEmptySet'
        ''
        'Is the set bounded?'
        '>> Y2.isBounded'
        ''
        'What is the dimension of the set?'
        'Y2.Dim'
    };

    % more inequalities
    slide(14).code = {
       'x = sdpvar(2,1);'
       'constraints = [ x(1)-3*x(2)<=4; 4*x(2)-0.5*x(1)<=1; 0 <= x <=3 ];'
       'Y3 = YSet(x,constraints);'
       'cla; axis on; hold on; Y3.plot;axis([-0.2 3.2 -0.2 1]); hold off;'
       'Y3.isBounded'
       'Y3.Dim'        
    };
    slide(14).text = {
       'Create a set comprising of three inequalities:'
       '>> x = sdpvar(2,1);'
       '>> constraints = [ x(1)-3*x(2)<=4; 4*x(2)-0.5*x(1)<=1; 0 <= x <=3 ];'
       '>> Y3 = YSet(x,constraints)'
       ''
       'Is the set bounded?'
       '>> Y3.isBounded'
       ''
       'What is the dimension of the set?'
       '>> Y3.Dim'
    };    
    slide(15).code = {
        'x = sdpvar(2,1);'
        'constraints = [-5 <= x <= 5; x(1) -2*x(2) == 1];'
        'Y4 = YSet(x,constraints);'
        'cla; axis on; hold on; Y4.plot; axis([-0.5 5.5 -1 0]); hold off;'
        'Y4.isEmptySet'
        'Y4.isBounded'
    };
    slide(15).text = {
        'We can add equality constraints to a convex set.'        
        '>> constraints = [-5 <= x <= 5; x(1) -2*x(2) == 1];'
        '>> Y4 = YSet(x,constraints)'
        ''
        'Is this set empty?'
        '>> Y4.isEmptySet'
        ''
        'Is the set bounded?'
        '>> Y4.isBounded'
    };
      
end

end
