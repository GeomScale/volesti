function v = struct2vars(s, v)
% extracts sdpvar objects from all fields of a structure
% "s" and concatenates them into a single column vector
% "v". dives recursively into fields of "s".

if isa(s, 'sdpvar')
	v = [v; s(:)];
elseif isstruct(s)
	n = fieldnames(s);
	for j = 1:length(n)
		v = struct2vars(s.(n{j}), v);
	end
end

end
