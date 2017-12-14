function [V, cells] = mpt_voronoi(S, varargin)
% Computes the Voronoi diagram of a set of points
%
% The Voronoi diagram of a set of points S = [s_1, ..., s_n] is a
% polyhedral partition consisting of polyhedra P_i, i = 1,...,n, where
%   P_i = { x | dist(x, s_i) <= dist(x, s_j), \forall i \ne j }
% where dist(a, b) is the Euclidian distance of two points
%
% Syntax:
%   V = mpt_voronoi(S)
%   V = mpt_voronoi(S, 'bound', B)
%   [V, P] = mpt_voronoi(S, 'bound', B)
%
% Inputs:
%   S: seed points stored column-wise
%   B: artificial bounding of the voronoi cells as a Polyhedron object
%
% Outputs:
%   V: all Voronoi cells as a PolyUnion (only non-empty cells)
%   P: all cells as an array of Polyhedron objects (includes potentially
%      empty cells)

% at least one input please
narginchk(1, Inf);
if ~isa(S, 'double') || isempty(S)
	error('First input must be a non-empty set of points.');
end
[nx, n_points] = size(S);
if n_points<2
	error('More than one point please.');
end

% parse and validate options
ip = inputParser;
ip.addParamValue('bound', [], @validate_polyhedron);
ip.parse(varargin{:});
options = ip.Results;
if ~isempty(options.bound)
	if options.bound.Dim~=nx
		error('The bound must be a polyhedron in R^%d.', nx);
	elseif options.bound.isEmptySet
		error('The bound must not be an empty set.');
	end
end

% create cells
cells(1:n_points) = Polyhedron;
for i = 1:n_points
	A = zeros(n_points-1, nx);
	b = zeros(n_points-1, 1);
	idx = 1;
	for j = setdiff(1:n_points, i)
		A(idx, :) = 2*(S(:, j)-S(:, i))';
		b(idx) = S(:, j)'*S(:, j) - S(:, i)'*S(:, i);
		idx = idx + 1;
	end
	C = Polyhedron(A, b);
	if ~isempty(options.bound)
		C = C.intersect(options.bound);
	end
	C.Data.voronoi.seed = S(:, i);
	cells(i) = C;
end

V = PolyUnion(cells);
% by definition the cells do not overlap
V.setInternal('Overlaps', false);
% by definition the union of cells is convex
V.setInternal('Convex', true);
% by definition the cells are connected
V.setInternal('Connected', true);

end
