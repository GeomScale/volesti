function msg = validate_vector(x, nx, label)
% Validates that "x" is an (nx x 1) vector
%
% Returns a non-empty string if "x0" does not match dimension
% of the state vector. returns an empty vector if everything is
% ok.

if nargin<3
	label = 'point';
end
if ~isequal(size(x), [nx 1])
	msg = sprintf('The %s must be a %dx1 vector.', label, nx);
else
	msg = '';
end

end
