function y=validate_logicalscalar(v)
%
% check if the argument is a logical scalar, otherwise throw an error
%
if nargin~=1
    error('validate_logicalscalar: One argument is required.');
end
    if islogical(v) && isscalar(v)
        y=true;
    else
        error('Input argument must be a logical scalar.');
    end

end
