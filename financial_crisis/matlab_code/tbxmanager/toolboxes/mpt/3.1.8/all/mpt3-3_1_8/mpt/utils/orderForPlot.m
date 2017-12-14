function I = orderForPlot(v, n)
% Compute an ordering of the points v, which lie on a hyperplane with
% normal n s.t. they form a closed surface as a patch object.
%

narginchk(2, 2);

validate_realmatrix(v);
validate_realvector(n);
if numel(n)<3
    error('The second argument should be a real vector with at least 3 elements.');
end

% Projective basis
n = n(:) / norm(n);
B = null(n')';

% Make sure that B satisfied the right-hand rule
% so that all faces have an outward facing orientation
B(2,:) = cross(n,B(1,:));

% An interior point
v0 = mean(v);

% Shift points to contain 0
v = v - ones(size(v,1),1)*v0;

% Sort by angle
th = atan2(B(1,:)*v', B(2,:)*v');
[th,I] = sort(th);
