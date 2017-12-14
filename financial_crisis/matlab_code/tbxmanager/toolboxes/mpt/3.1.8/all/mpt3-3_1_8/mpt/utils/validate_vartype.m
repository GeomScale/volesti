function y=validate_vartype(varargin)
%
% check if the argument belongs to C-continuous, I-integer, B-binary,
% S-semicontinuous, N-semiinteger variable
%
narginchk(1, 1);

% find first character different from C, I, B, S, N
if size(varargin{:},1)>size(varargin{:},2)
    str = transpose(varargin{:}); % must be row
else
    str = varargin{:}; % must be row
end
c = regexpi(str,'[^CIBSN]', 'once');
if ~isempty(c)
    error('validate_vartype: Only C, I, B, S, N characters are allowed.');
else
    y = true;
end
