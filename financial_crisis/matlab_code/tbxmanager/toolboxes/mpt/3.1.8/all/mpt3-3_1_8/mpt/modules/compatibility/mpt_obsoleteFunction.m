function mpt_obsoleteFunction(name)
% Warn a user that function 'name' is obsolete

if nargin==0
	% automatically determine name of the caller
	s = dbstack;
	name = s(2).name;
end

fprintf('Function %s is obsolete and will be removed in a future MPT version.\n', name);
