function y=validate_realinfvector(v)
%
% check if the argument is a real vector, otherwise throw an error
%
% empty argument (v=[]) is considered as valid
if nargin~=1
    error('validate_realinfvector: One argument is required.');
end

if isnumeric(v) && length(size(v))<=2 && min(size(v))<=1 && ~any(isnan(v(:))) && isreal(v)
    y=true;
else
    error('Input argument must be a real vector.');
end
end
