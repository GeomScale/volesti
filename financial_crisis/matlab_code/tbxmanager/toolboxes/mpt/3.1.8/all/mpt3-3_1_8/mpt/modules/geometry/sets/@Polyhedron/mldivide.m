function res = mldivide(P,Pn,noconstruction)
%
%  MLDIVIDE: Set difference between polyhedra 
%  ===========================================
%  
%  
%  SYNTAX
%  ------
%     
%      R = P \ Q
%      R = mldivide(P,Q)
%    
%  
%  DESCRIPTION
%  -----------
%     Function computes the set difference between polyhedron P  and Qwhich can be
%  both single polyhedra or arrays of polyhedra in the same dimension The set
%  difference operation is defined as 
%                       R = P \ Q ={ x  |  xin P  , xnotin Q }
%     where the output R  can comprise of multiple polyhedra.
%  
%  INPUT
%  -----
%     
%        
%          P Polyhedron in any format                 
%            Class: Polyhedron                        
%          Q Polyhedron in any format                 
%            Class: Polyhedron                        
%              
%  
%  
%  OUTPUT
%  ------
%     
%        
%          R Polyhedron (or array) R = P? .           
%            Class: Polyhedron                        
%              
%  
%  
%  SEE ALSO
%  --------
%     plus,  minus
%  

%  AUTHOR(s)
%  ---------
%     
%    
%   (c) 2003  Mato Baotic: ETH Zurich
%   mailto:baotic@control.ee.ethz.ch 
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
 
 
global MPTOPTIONS
if isempty(MPTOPTIONS)
   MPTOPTIONS = mptopt;
end
if nargin<3
	% boolean flag indicating that we can exit quickly if P\Pn is non-empty
	noconstruction=false;
end

validate_polyhedron(P);
validate_polyhedron(Pn);

% improve performance by discarding any empty sets
Pn = Pn(find(~Pn.isEmptySet()));

if numel(Pn)==0 || numel(P)==0
	% 1) P=emptyset => P\Pn = emptyset
	% 2) Pn=emptyset => P\Pn = P
	res = P;

elseif numel(P)==1 && numel(Pn)==1
	% set difference between two polyhedra (full- or lower-dimensional)
	res = single_setdiff(P, Pn, noconstruction);

elseif numel(Pn)==1
	% array \ polyhedron can be computed by applying single_setdiff() to
	% each element of the array
	%
	% supports both full- and lower-dimensional polyhedra
	res = [];
	for i = 1:numel(P)
		R = single_setdiff(P(i), Pn, noconstruction);
		res = [res R];
	end

elseif numel(P)==1 && P.isFullDim() && all(Pn.isFullDim())
	% polyhedron \ array, all fully dimensional
	%
	% special fast implementation via regiondiff
	res = vanilla_regiondiff(P, Pn, noconstruction);
	
else
	% array \ array or polyhedron \ array with lowdim elements
	res = P;
	for i = 1:numel(Pn)
		res = mldivide(res, Pn(i), noconstruction);
	end
	
end

% kick out empty polyhedra
if numel(res)>0
	res = res(find(~res.isEmptySet()));
end
if numel(res)==0 && numel(P)>0
	res = Polyhedron.emptySet(P(1).Dim);
end

% TODO: either return a copy of P without functions, or make sure that
% every polyhedron of the set difference has the functions of the input
% propagated to it.

end

%%
function R = vanilla_regiondiff(P, Q, noconstruction)
%
% fast regiondiff for full-dimensional polyhedra
%
% vanilla implementation of Alg. 2 from
% Baotic: Polytopic Computation in Constrained Optimal Control,
% Automatika, 50(3-4), 119-134(2009)

R = [];
k = 1;
intersects = zeros(1, numel(Q));
for i = 1:length(Q)
	if Q(i).Dim ~= P.Dim
		error('MPT:Polyhedron', 'All polyhedra must be in the same dimension.');
	end
	
	% all Q's must be in minimal H-representation
	Q(i).minHRep();
	
	% we only need to consider those Q's which intersect with P. moreover,
	% this function assumes that all polyhedra are fully dimensional,
	% that's why we are only interested in full-dimensional intersections
	intersects(i) = fast_isFullDim([P.H; Q(i).H]);
end
% kick out elements that do not intersect with P
Q = Q(find(intersects));
if numel(Q)==0
	% P\emptyset = P
	R = P;
	return
end

PH = P.H;
QkH = Q(k).H;
for j = 1:size(QkH, 1)
	% S = P.intersect(negation of the j-th face of Q)
	An = [PH(:, 1:end-1); -QkH(j, 1:end-1)];
	bn = [PH(:, end); -QkH(j, end)];
	if fast_isFullDim([An, bn])
		% fully dimensional intersection exists
		S = Polyhedron(An, bn);
		if k < length(Q)
			R = [R vanilla_regiondiff(S, Q(k+1:end), noconstruction)];
		else
			if noconstruction
				R = S;
				return
			else
				R = [R S];
			end
		end
	end
	% P = P.intersect(j-th face of Q)
	PH = [PH; QkH(j, :)]; 
end

end

%%
%---------------------------------------------------------
function res = single_setdiff(Y, R, noconstruction)
% set difference between two polyhedra
%
% both input arguments can either be lower- or fully-dimensional
%
% the output is in minimal H-representation unless noconstruction=true

global MPTOPTIONS

% quickly check several trivial special cases:
if Y.isEmptySet()
	% special case: if Y is empty set, then Y\R = empty set
	res = Polyhedron.emptySet(Y.Dim);
	return

elseif Y.Dim ~= R.Dim
	error('MPT:Polyhedron', 'All polyhedra must be in the same dimension.');

elseif R.isEmptySet() || ~Y.doesIntersect(R)
	% special cases where Y\R = Y
	% 1) R is empty set
	% 2) Y and R do not intersect
	res = Polyhedron(Y);
	return
	
end

% identify implicit equalities
Y.minAffineRep();
R.minAffineRep();

Y.minHRep();
R.minHRep();
res = [];
R_A = R.A;
R_b = R.b;
Y_A = Y.A;
Y_b = Y.b;

if Y.isFullDim() && R.isFullDim()
	% both are fully dimensional, just flip half-spaces and check
	% full dimensionality of the results
	for i = 1:length(R_b)
		A = [-R_A(i, :); R_A(1:i-1, :); Y_A];
		b = [-R_b(i); R_b(1:i-1); Y_b];
		if fast_isFullDim([A, b])
			Ri = Polyhedron(A, b);
			if noconstruction
				% exit quickly if we are only checking coverage
				res = Ri;
				return
			else
				res = [res Ri];
			end
		end
	end

elseif Y.isFullDim() || (~isempty(R.Ae) && ~isempty(Y.Ae) && ...
		rank(null(R.Ae)) < rank(null(Y.Ae)))
	% 1) Y is full-dim, R is low-dim
	% 2) both are low-dim, but R has a smaller dimension of
	%    the affine hull
	%
	% unless we support open half-space, the set
	% difference is equal to Y
	res = Polyhedron(Y); % copy constructor
	
elseif R.isFullDim() || (~isempty(R.Ae) && ~isempty(Y.Ae) && ...
		rank(null(R.Ae)) > rank(null(Y.Ae)) )
	% 1) Y is lower-dimensional, but R is full dimensional
	% 2) both are lower dimensional, but R has higher affine dimension
	%
	% open half-spaces of R need to be shifted and emptienies has to be
	% checked 
	Y_He = Y.He;
	shift_tol = MPTOPTIONS.rel_tol*2;
	for i = 1:length(R_b)
		A = [-R_A(i, :); R_A(1:i-1, :); Y_A];
		b = [-R_b(i)-shift_tol; R_b(1:i-1); Y_b];
		if ~fast_isEmptySet([A b], Y_He)
			% shift the half-space back
			b(1) = b(1)+shift_tol;
			Ri = Polyhedron('H', [A b], 'He', Y_He);
			if noconstruction
				% exit quickly if we are only checking coverage
				res = Ri;
				return
			else
				res = [res Ri];
			end
		end
	end
		
elseif rank([R.He; Y.He], MPTOPTIONS.zero_tol) > max(size(R.He, 1), size(Y.He, 1))
	% both are lower dimensional, but their affine hulls do not intersect,
	% hence Y\R=Y
	res = Polyhedron(Y);
	
elseif isempty(R.He) && isempty(Y.He) && Y<=R
	% both are lower-dimensional, but with no affine subspace, most
	% probably a vertex
	res = [];

else
	% Y and R are both lower-dimensional and affdim(Y)<=affdim(R)
	
	% project R and Y on the fully-dimensional null space of R
	F = null(R.Ae);
    x0 = R.Ae\R.be;
	Y_Z = Polyhedron('H', [Y_A*F, Y_b-Y_A*x0], 'He', [Y.Ae*F, Y.be-Y.Ae*x0]);
	R_Z = Polyhedron(R_A*F, R_b-R_A*x0);

	% comptute set difference where R_Z is now fully dimensional
	res = single_setdiff(Y_Z, R_Z, noconstruction);
	
	% project back on the original space
	res_new = [];
	for i = 1:numel(res)
		if ~res(i).isEmptySet()
			A = res(i).A;
			b = res(i).b;
			Ri = Polyhedron('H', [A*pinv(F), b+A*pinv(F)*x0], ...
				'He', Y.He);
			res_new = [res_new Ri];
		end
	end
	res = res_new;
end

if numel(res)==0
	% make sure the output is always a polyhedron in the dimension of Y
	res = Polyhedron.emptySet(Y.Dim);
end

end
