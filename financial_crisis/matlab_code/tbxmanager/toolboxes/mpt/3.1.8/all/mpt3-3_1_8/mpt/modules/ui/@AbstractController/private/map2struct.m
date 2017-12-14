function [s, inits] = map2struct(m, inits, master_component, master_key)
% converts a containers.Map to a structure (operates recursively)

if isa(m, 'containers.Map')
	k = m.keys;
	s = [];
	for i = 1:length(k)
		v = m(k{i});
		if nargin<3
			master_component = k{i};
		end
		if nargin<4
			local_master_key = k{i};
		else
			local_master_key = [master_key '.' k{i}];
		end
		if ~isempty(v)
			[new_s, inits] = map2struct(v, inits, master_component, local_master_key);
			s.(k{i}) = new_s;
		end
	end
else
	s = m.var;
	if m.parametric
		% include this variable into list of initial conditions
		v_init.component = master_component;
		v_init.name = master_key;
		v_init.var = s;
		v_init.dims = size(s);
		if isempty(inits)
			inits = v_init;
		else
			inits(end+1) = v_init;
		end
	end
end

end
