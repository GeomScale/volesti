function ts = validate_marker(s)
%
% checks the validity of marker
%

narginchk(1, 1);

if ~ischar(s)
    error('Argument must be a char.');
end

% prepare output
ts = true;

list = {'.','o','x','+','*','s','d','v','^','<','>','p','h','none'};

index = find(strcmp(s,list),1);

if isempty(index)
    error('Marker style must be one of these ".", "o", "x", "+", "*", "s", "d", "v", "^", "<", ">", "p", "h" or "none".');
end

end
