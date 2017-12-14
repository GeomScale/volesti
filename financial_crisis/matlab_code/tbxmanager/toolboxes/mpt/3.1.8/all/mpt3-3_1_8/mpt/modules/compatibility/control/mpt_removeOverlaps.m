function out = mpt_removeOverlaps(in, varargin)
% Wrapper around PolyUnion/min()

% convert inputs into an array of PolyUnion objects
PUs = mpt_mpsol2pu(in);

if numel(in)>1
	% remove overlaps based on cost
	NoOverlaps = PUs.min('obj');
else
	NoOverlaps = PUs;
end

% convert back to MPT2 format
out = [];
out.sysStruct = in{1}.sysStruct;
out.probStruct = in{1}.probStruct;
out.details = in{1}.details;
out.Pfinal = polytope(NoOverlaps.Domain);
out.Pn = polytope(NoOverlaps.Set);
out.Fi = cell(1, NoOverlaps.Num);
out.Gi = cell(1, NoOverlaps.Num);
out.Ai = cell(1, NoOverlaps.Num);
out.Bi = cell(1, NoOverlaps.Num);
out.Ci = cell(1, NoOverlaps.Num);
for i = 1:NoOverlaps.Num
	f = NoOverlaps.Set(i).Functions('primal');
	c = NoOverlaps.Set(i).Functions('obj');
	out.Fi{i} = f.F;
	out.Gi{i} = f.g;
	out.Ai{i} = [];
	out.Bi{i} = c.F;
	out.Ci{i} = c.g;
end
out.dynamics = zeros(1, NoOverlaps.Num);
out.overlaps = 0;
out.convex = 0;

end
