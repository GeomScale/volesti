function h_all=mpt_plotPWQ(Pn,Ai,Bi,Ci,meshgridpoints,Options)
%MPT_PLOTPWQ Plots a PWQ function defined over polyhedral partition
%
% THIS FUNCTION IS OBSOLETE!
%
% Use PolyUnion/fplot() instead.

mpt_obsoleteFunction;

narginchk(4, 6);

if ~isa(Pn,'polytope')
    error('mpt_plotPWQ: First input argument MUST be a polytope');
end

Q = toPolyhedron(Pn);
Q.removeAllFunctions();
for i = 1:numel(Q)
	Q(i).addFunction(QuadFunction(Ai{i}, Bi{i}, Ci{i}), 'pwq');
end
h_all=Q.fplot;
