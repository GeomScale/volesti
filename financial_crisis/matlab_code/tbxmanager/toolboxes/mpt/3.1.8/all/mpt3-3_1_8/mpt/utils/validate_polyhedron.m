function y=validate_polyhedron(varargin)
%
% check if the argument is of a polyhedron class, otherwise throw an
% error
%
if nargin<1
    error('validate_polyhedron: At least One argument is required.');
end
if ~isvector(varargin)
    error('validate_polyhedron: Only arrays of polyhedra are supported.');
end
for i=1:length(varargin)
    if isa(varargin{i},'Polyhedron')
        y=true;
    else
        error('Input argument must be a "Polyhedron" class.');        
    end
end
end
