function [objects, options] = parsePlotOptions(main_class, varargin)
% Parses options for MPT3's plotting functions
%
% [objects, options] = parsePlotOptions(main_class, ...
%                   P1, 'opt1', val1, 'opt2', val2, P2, P3, 'opt3', 'val3')
%
% objects = {P1, P2, P3}
% options{1} = {'opt1', val1, 'opt2', val2}
% options{2} = {}
% options{3} = {'opt3', 'val3}
%
% Objects in "obj" are recognized if isa(P, main_class) is true.

% This is an internal helper.

obj_pos = find(cellfun(@(x) isa(x, main_class), varargin));
nobj = length(obj_pos);
objects = cell(1, nobj); % stack of objects to plot
options = cell(1, nobj); % stack of options for each object
obj_pos = [obj_pos length(varargin)+1];
for i = 1:nobj
	objects{i} = varargin{obj_pos(i)};
	options{i} = varargin(obj_pos(i)+1:obj_pos(i+1)-1);
end

end
