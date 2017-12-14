% ECOS - Embedded COnic Solver.
%
% Self-dual homogeneous embedding interior point implementation for optimization
% over linear or second-order cones. ECOS does not support semi-definite
% cones (feel free to contact the developers if you wish to add SDP support 
% to ECOS.
%
%   [x,y,info,s,z] = ECOS(c,G,h,dims) Solves a pair of primal and dual
%   cone programs
% 
%        minimize    c'*x
%        subject to  G*x + s = h
%                    s >= 0
%
%        maximize    -h'*z
%        subject to  G'*z + c = 0
%                    z >= 0.
%
%   The inequalities are with respect to a cone K defined as the Cartesian
%   product of N+1 cones:
% 
%        K = K_0 x K_1 x .... x K_N.
% 
%     The first cone, K_0, is the nonnegative orthant of dimension dims.l.
%     The next N cones are second order cones of dimension dims.q(1), ...,
%     dims.q(N), where a second order cone of dimension m is defined as
% 
%         { (u0, u1) in R x R^{m-1} | u0 >= ||u1||_2 }.
% 
%     INPUT arguments:
% 
%         c is a dense column vector of size n
% 
%         dims is a struct with the dimensions of the components of cone K.
%         It has two fields.
%         - dims.l, the dimension of the nonnegative orthant C_0, with l>=0.
%         - dims.q, a row vector of N integers with the dimensions of the 
%           second order cones K_1, ..., K_N. (N >= 0 and q(i) >= 3.)
% 
%         G is a sparse matrix of size (m,n), where
% 
%             m = dims.l + dims.q(1) + ... + dims.q(N).
%               = dims.l + sum(dims.q)
% 
%         Each column of G describes a vector
% 
%             v = ( v_0, v_1, ..., v_N )
% 
%         in V = R^dims.l x R^dims.q(1) x ... x R^dims.q(N)
%         stored as a column vector
% 
%             [ v_0; v_1; ...; v_N ].
% 
%         h is a dense column vector of size m, representing a vector in V,
%         in the same format as the columns of G.
%         
%
%
%   [x,y,info,s,z] = ECOS(c,G,h,dims,A,b) Solves a pair of primal and 
%   dual cone programs
% 
%        minimize    c'*x
%        subject to  G*x + s = h
%                    A*x = b
%                    s >= 0
%
%        maximize    -h'*z - b'*y
%        subject to  G'*z + A'*y + c = 0
%                    z >= 0.
%
%      where c,G,h,dims are defined as above, and A is a sparse matrix of 
%      size (p,n), and b is a dense matrix of size (p,1).
%       
%      It is assumed that rank(A) = p and rank([A; G]) = n.
%
% 
%   [x,y,info,s,z] = ECOS(c,G,h,dims,opts) and 
%   [x,y,info,s,z] = ECOS(c,G,h,dims,A,b,otps) are as above, with the struct 
%   otps used to control settings of the solver. Use the function ECOSOPTIMSET
%   to obtain a default initialization, and see help ECOSOPTIMSET for more 
%   details.
% 
%
% Details on ECOS can be found at http://embotech.com/ECOS and in the paper 
%    Alexander Domahidi, Eric Chu, Stephen Boyd. "ECOS: An Embedded
%    Conic Solver." In proceedings of European Control Conference (ECC), 
%    pp. 3071-3076, Zurich, Switzerland, July 2013."
%
% More details are given in A. Domahidi's PhD Thesis (Chapter 9): 
%    http://e-collection.library.ethz.ch/view/eth:7611
%
%
% (c) A. Domahidi, ETH Zurich & embotech GmbH, Zurich, Switzerland, 2012-15.
%
% The branch and bound module is (c) Han Wang, Stanford University.
%
% LICENSE: ECOS is distributed under GPLv3. For commercial licenses and 
%          professional support contact embotech at ecos@embotech.com.
%
% See also ECOSQP ECOSOPTIMSET ECOS_LICENSE
