function y=validate_dimension(v)
%
% Check if the argument is a proper dimension, otherwise throw an error
%

if isreal(v) && isscalar(v) && isfinite(v) && v>0 && mod(v, 1)==0
    y = true;
else
    error('Input argument is a not valid dimension.');
end
end
