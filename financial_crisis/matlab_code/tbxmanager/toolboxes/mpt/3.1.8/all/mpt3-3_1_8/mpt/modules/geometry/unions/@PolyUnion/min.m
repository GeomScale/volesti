function out = min(PUs, func, coefficient)
% Minimum of PolyUnions 'inPUs' using function 'func'
%
% This function is basically MPT3 implementation of mpt_removeOverlaps.%
%
% If coefficient=-1, this function picks the maximum instead of the minimum

% Note: This is a slow, but compact implementation.

global MPTOPTIONS
if isempty(MPTOPTIONS)
	MPTOPTIONS = mptopt;
end

if nargin < 2 || isempty(func)
	% if no function is specified, take the first one
	if length(PUs(1).listFunctions)>1
		error('Please specify which function to use for comparison.');
	else
		fnames = PUs(1).listFunctions;
		if isempty(fnames)
			error('The object has no functions.');
		end
		func = fnames{1};
	end
end

if ~PUs(1).hasFunction(func)
	error('Couldn''t find function "%s".', func);
end

% make sure all other polyunions define the function
for i = 2:length(PUs)
	if ~PUs(i).hasFunction(func)
		error('Polyunion #%d must define function "%s".', i, func);
	end
end

if nargin<3
	% coefficient=1 means we are looking for the maximum
	% coefficient=-1 will search for the maximum
	coefficient = 1;
end

% TODO: make sure all regions are in H-representation
%code goes here

% make sure that the comparing function is a scalar PWA function
% TODO: extend this check to all regions of all polyunions
fun = PUs(1).Set(1).getFunction(func);
if ~isa(fun, 'AffFunction')
	error('Only PWA functions can be used for comparison.');
elseif fun.R~=1
	error('Only scalar-valued functions can be handled.');
end

% domains of input unions will be preserved
domains = cat(1, PUs.Domain);

% deal with overlaps within of individual polyunions
newPUs = [];
for i = 1:numel(PUs)
	if PUs(i).isOverlapping
		% treat regions of overlapping partitions as one-region polyunions
		for j = 1:PUs(i).Num
			newPUs = [newPUs, PolyUnion(PUs(i).Set(j))];
		end
	else
		% add the full polyunion if it does not contain overlaps
		newPUs = [newPUs, PUs(i)];
	end
end
PUs = newPUs;

nR = 0;
Pfinal = []; % list of regions

% list of intersecting partitions
PartitionsIntersect = true(numel(PUs));
for i = 1:numel(PUs)
	for j = setdiff(1:numel(PUs), i)
		if numel(PUs(i).Domain)==1 && ...
				numel(PUs(j).Domain)==1 && ...
				~PUs(i).Domain.intersect(PUs(j).Domain).isFullDim()
			% partitions do not intersect, no reason to check regions
			%
			% TODO: support domains composed of multiple polyhedra
			PartitionsIntersect(i, j) = false;
		end
	end
end

% map of intersecting regions
t = clock;
fprintf('Computing the intersection map...');
RegionsIntersect = cell(numel(PUs), numel(PUs));
for i = 1:numel(PUs)-1
    for j = i+1:numel(PUs)
        if ~PartitionsIntersect(i, j), continue, end
        map = false(PUs(i).Num, PUs(j).Num);
        for ir = 1:PUs(i).Num
            for jr = 1:PUs(j).Num
                status = PUs(i).Set(ir).doesIntersect(PUs(j).Set(jr), 'fully');
                map(ir, jr) = status;
            end
        end
        RegionsIntersect{i, j} = map;
    end
end
fprintf(' done in %.2f seconds\n', etime(clock, t));

for ipart = 1:numel(PUs)

	fprintf('Union %d (out of %d)...\n', ipart, numel(PUs));

	for ireg = 1:numel(PUs(ipart).Set)

        % list of regions in which the function is better than in the
        % 'ipart/ireg' region
        better_regions = [];

		for jpart = setdiff(1:numel(PUs), ipart)
		
            if ~PartitionsIntersect(ipart, jpart)
                % partitions do not intersect, no reason to check regions
                continue
            elseif ipart < jpart
                dointersect = RegionsIntersect{ipart, jpart};
            else
                dointersect = RegionsIntersect{jpart, ipart}';
            end
			
			for jreg = 1:numel(PUs(jpart).Set)
				
				% do the regions intersect?
                if dointersect(ireg, jreg)
					% compare function in the intersection
                    Q = PUs(ipart).Set(ireg).intersect(PUs(jpart).Set(jreg));
					
					ifun = PUs(ipart).Set(ireg).getFunction(func);
					jfun = PUs(jpart).Set(jreg).getFunction(func);
					Fj = jfun.F;
					gj = jfun.g;
					Fi = ifun.F;
					gi = ifun.g;
					
					% allow to compute maximum by setting coeffcient=-1
					Fdiff = coefficient*(Fj-Fi);
					gdiff = coefficient*(gj-gi);

					include_Q = false;
					if norm(Fdiff)+norm(gdiff) <= MPTOPTIONS.abs_tol
						% function is identical in both regions
						if ipart < jpart
							% this makes sure that we remove only one such
							% region 
							include_Q = true;
						end
						
					else
						% functions are different, find which part of
						% region 'ireg' has the smaller function values
						
						% TODO: correcly compare piecewise constant
						% functions (probably requires lifting)
						if norm(Fdiff) <= MPTOPTIONS.zero_tol
							% lift by one dimension to deal with constant
							% cost
							infbox = Polyhedron('lb', 0, ...
								'ub', MPTOPTIONS.infbound);
							Ireg = PUs(ipart).Set(ireg)*infbox;
							Jreg = PUs(jpart).Set(jreg)*infbox;
							Qreg = Ireg.intersect(Jreg);
							if ~Qreg.isFullDim
								error('Regions should intersect.');
							end
							Q = Polyhedron('H', [Qreg.H; Fdiff gdiff 0], 'He', Qreg.He);
							if Q.isFullDim
								Q = Q.projection(1:Q.Dim-1).minHRep();
								include_Q = true;
							end
						else
							Q = Polyhedron('H', [Q.H; Fdiff -gdiff], 'He', Q.He);
							if Q.isFullDim
								% the dominating subregion is fully
								% dimensional, include it
								include_Q = true;
							end
						end
						
					end
					
					if include_Q
						better_regions = [better_regions Q];
					end

				end % non-empty intersection
				
			end % jreg
			
		end % jpart
		
		if numel(better_regions) > 0
			% get all subregions of 'ireg' where the function is minimal
			
			Ri = PUs(ipart).Set(ireg) \ better_regions;
			
			if ~Ri.isFullDim
				% If the set difference is empty, the other regions with lover
				% function values are covering 'ireg'. 
				continue
				
			else
				% Otherwise we add all subregions of Ri

				% TODO: Polyhedron/mldivide should keep the
				% functions

				% make a copy before removing functions, otherwise
				% PUs(ipart).Set(ireg) could be affected
				Ri = Polyhedron(Ri);
				Ri.removeAllFunctions();
				for mm = 1:numel(Ri)
					if Ri(mm).isFullDim
						% we keep this region
						nR = nR+1;
						% copy all function from PUs(ipart).Set(ireg)
						Ri(mm).copyFunctionsFrom(PUs(ipart).Set(ireg));
					end
				end
				Pfinal = [Pfinal Ri];
			end

		else
			nR = nR+1;
			Pfinal = [Pfinal PUs(ipart).Set(ireg)];
		end
		
	end % ireg

end % ipart

if isa(Pfinal, 'double') && isempty(Pfinal)
	out = PolyUnion;
else
	out = PolyUnion('Set', Pfinal', 'overlaps', false, 'Domain', domains);
end

end
