function ts = validate_linestyle(s)
%
% checks the validity of linestyle
%

narginchk(1, 1);

if ~ischar(s)
    error('Argument must be a char.');
end

% prepare output
ts = true;

list = {'-',':','-.','--','none'};

index = find(strcmp(s,list),1);

if isempty(index)
    error('Style of the line must be one of these "-", ":", "-.", "--", or "none".');
end

end
