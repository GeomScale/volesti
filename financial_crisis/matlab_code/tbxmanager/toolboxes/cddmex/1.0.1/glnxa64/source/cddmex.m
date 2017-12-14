%==========================================================================
% Help file for CDDMEX
%
% CDDMEX is a Matlab MEX file interface (authors: Mato Baotic, Fabio Torrisi),
% obtained by linking CDD library CDDLIB-093 (author: Komei Fukuda), which
% allows calls to some of CDD functions from within MATLAB.
%
% CDD is a software package written by K. Fukuda which boasts a wide array
% of efficient algorithms for many polytope manipulations as well as a fast
% LP solver. The full CDD library is free and availible from:
% http://www.cs.mcgill.ca/~fukuda/soft/cdd_home/cdd.html
%
% Version: cddmex v1.0
%          cddlib v0.94g
%
% See the bottom of this file for more details.
%==========================================================================
%
% Syntax: outstruct=cddmex(action,instruct)
%
% action:  'hull'             Convex hull of a V polyhedron
%          'extreme'          Vertex/ray enumeration of an H polyhedron
%          'reduce_h'         Minimal H representation (of an H polyhedron)
%          'reduce_v'         Minimal V representation (of a V polyhedron)
%          'solve_lp'         Solve a Linear Program using Criss-Cross method
%          'solve_lp_DS'      Solve a Linear Program using Dual-Simplex method
%          'version'          Version of the mex function
%
% Further info: see CDDMEX.C
%
% Not documented functions:
%          'copy_v'           Simple copy of an H polyhedron
%          'copy_h'           Simple copy of a V polyhedron
%          'v_hull_extreme'   ???
%          'adj_extreme'      Extreme points and adjecancy list of an H polyhedron
%          'find_interior'    Interior point of an H polyhedron
%
% Not compiled functions:
%          'adjacency'        Adjacency list of a V polyhedron
%                             Note that V has to be in a minimal representation!
%
%
% Example: Find the extreme points/rays of P={x: A1*x=B1, A2*x<=B1}
%
%    A1=[1 1 1];B1=[1];
%    A2=[eye(3);-eye(3)];B2=[1;1;1;2;2;2];
%    H=struct('A',[A1;A2],'B',[B1;B2],'lin',(1:size(B1,1))');
%    % H.lin=indices of equality constraints
%    V=cddmex('extreme',H);
%    % The rows of V.V are the extreme points of P
%
% Example: Find the convex hull of v1,v2,v3
%
%    V=struct('V',[v1';v2';v3']);
%    H=cddmex('hull',V);
%
% Example: Find minimal H representation
%
%    A1=[-1 0; 0 -1; 1 0; 0 1; 1 1];B1=[0;0;1;1;1];
%    H1=struct('A',A1,'B',B1);
%    [Hred]=cddmex('reduce_h',H1);
%    %  Hred.A * x <= Hred.B  is a minimal H representation
%    [Hred,ind_redrows]=cddmex('reduce_h',H1);
%    %  ind_redrows are the indices of redundant rows
%    % Note: if H1 is empty polyhedron, then Hred
%    % containes irreducibly incosistent set (IIS)
%
% Example: Find minimal V representation
%
%    V1=[0 0; 0 1; 1 0; 1 1; 0.2 0.8; 0.7 0.3];
%    V1=struct('V',V1);
%    % each row of V1.V is a vertex
%    [Vred]=cddmex('reduce_v',V1);
%    %  Vred.V  is a minimal V representation
%    [Vred,ind_redvert]=cddmex('reduce_v',V1);
%    %  ind_redvert are the indices of redundant vertices
%
%
% Example: Solve a Linear Program
% 		Syntax
%          OUT=cddmex('solve_lp',IN)
% 		
% 		solves problem
%          min	IN.obj*x
%           x
%          s.t.	IN.A*x <= IN.B
%          	IN.A(IN.lin,:) x == IN.B(IN.lin)
% 		
% 		where
%           IN  -	input structure
%          	IN.A	- constraints matrix
%          	IN.B	- constraints matrix
%          	IN.obj	- objective function
%          	IN.lin	- indices of equality constraints
%          	
%           OUT -	structure with solution (similiar to lpsolve)
%          	OUT.xopt	- primal solution
%          	OUT.lambda	- dual solution
%          	OUT.how		- Fukuda's code for LP solution status
%          			  0 = dd_LPSundecided,
%          			  1 = dd_Optimal 
%          			  2 = dd_Incosistent
%          			  3 = dd_DualIncosistent
%          			  4 = dd_StrucIncosistent
%          			  5 = dd_StrucDualIncosistent
%          			  6 = dd_Unbounded
%          			  7 = dd_DualUnbounded
%          	OUT.objlp	- optimal value
%
%       Note: CrissCross algorithm is used when solving LP.
%          	
% 	Example
% 		objective=[-1 1];
% 		A1=[1 1; -1 0; 0 -1];
% 		B1=[1; 0; 0];
%        	IN=struct('obj',objective,'A',A1,'B',B1);
% 			
% 		OUT = cddmex('solve_lp',IN);
% 		
% 		OUT =
% 		    xopt: [1 0]
% 		    lambda: [-1 0 -2]
% 		    how: 1
% 		    objlp: -1
%
%
% Note: the projection function is only available in CDDMEX.M. A projection
% can be obtained here by (1) vertex enumeration, (2) projection of vertices,
% (3) hull of projected vertices. Example:
%
%    H=struct('A',A1,'B',B1);
%    V=cddmex('extreme',H);
%    Vproj.V=V.V(:,1:2); % Project over the first two coordinates
%    Hproj=cddmex('hull',Vproj); % Get the projection
%
%
%==============================================================================
%
% /* The cdd library cddlib-093a was written by
%    Komei Fukuda, fukuda@ifor.math.ethz.ch
%    Version 0.93, August 11, 2003
%    Standard ftp site: ftp.ifor.math.ethz.ch, Directory: pub/fukuda/cdd
% */
%
% /* cddlib.c : C-Implementation of the double description method for
%    computing all vertices and extreme rays of the polyhedron 
%    P= {x :  b - A x >= 0}.  
%    Please read COPYING (GNU General Public Licence) and
%    the manual cddlibman.tex for detail.
% */
%
% /*  This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation; either version 2 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program; if not, write to the Free Software
%     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
% */
%===============================================================================
%
% History:
% cddmex.c v0.1, An initial Matlab MEX wrapper for cdd library written by
% 		 Fabio Torrisi (torrisi@control.ee.ethz.ch) and
% 		 Mato Baotic (baotic@control.ee.ethz.ch),
%		 Zurich, December 2002.
%
% cddmex.c v1.0, Revision update by Mato Baotic, Zurich, October 2003.


