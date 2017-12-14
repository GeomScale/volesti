function y=validate_realvector(v)
%
% check if the argument is a real vector, otherwise throw an error
%
% empty argument (v=[]) is considered as valid

if isnumeric(v) && (isvector(v) || isempty(v)) && all(isfinite(v)) && isreal(v)
    y=true;
else
    error('Input argument must be a real vector.');
end
end
