function h = plot_internal(obj, idx, varargin)
% Plots a single union
%
% This is an internal helper called from Union/plot. Must take color index
% as the second input.

h = [];
if iscell(obj.Set)
	for k = 1:obj.Num
		hk = plot(obj.Set{k}, varargin{:}, 'array_index', idx);
		idx = idx + 1;
		h = [h; hk];
	end
else
	% plot the whole set at once
	h = plot(obj.Set, varargin{:}, 'array_index', idx);
end

end
