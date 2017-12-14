function U = merge(obj, function_name, optimal, optimal_merging)
%
%  MERGE: Greedy merging of polyhedra 
%  ===================================
%  
%  
%  SYNTAX
%  ------
%     
%      U.merge
%      merge(U)
%    
%  
%  DESCRIPTION
%  -----------
%     Simplifies the union of polyhedra by merging the neighboring polyhedra if
%  their union is convex. The algorithm cycles through the regions and checks if
%  any two regions form a convex union. If so, the algorithm combines them in one
%  region, and continues checking the remaining regions. To improve the solution,
%  multiple merging loops can be enabled in options.
%  
%  INPUT
%  -----
%     
%        
%          U Union of polyhedra in the same           
%            dimension.                               
%            Class: PolyUnion                         
%              
%  
%  
%  SEE ALSO
%  --------
%     reduce
%  

%  AUTHOR(s)
%  ---------
%     
%    
%   (c) 2010-2013  Martin Herceg: ETH Zurich
%   mailto:herceg@control.ee.ethz.ch 
%     
%    
%   (c) 2003-2013  Michal Kvasnica: STU Bratislava
%   mailto:michal.kvasnica@stuba.sk 
%     
%    
%   (c) 2005  Frank J. Christophersen: ETH Zurich
%   mailto:fjc@control.ee.ethz.ch 
%     
%    
%   (c) 2005  Tobias Geyer: ETH Zurich
%   mailto:geyer@control.ee.ethz.ch 
%  
%  

%  LICENSE
%  -------
%    
%    This program is free software; you can redistribute it and/or modify it under
%  the terms of the GNU General Public License as published by the Free Software
%  Foundation; either version 2.1 of the License, or (at your option) any later
%  version.
%    This program is distributed in the hope that it will be useful, but WITHOUT
%  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
%  FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%    You should have received a copy of the GNU General Public License along with
%  this library; if not, write to the  Free Software Foundation, Inc.,  59 Temple
%  Place, Suite 330,  Boston, MA 02111-1307 USA
%  -------------------------------------------------------------------------------
%    
%      This document was translated from LaTeX by HeVeA (1).
%  ---------------------------------------
%    
%    
%   (1) http://hevea.inria.fr/index.html
 
 
global MPTOPTIONS

%% parsing of inputs
narginchk(1, Inf);
if nargin==1
	function_name = '';
	optimal_merging = false;
elseif nargin==2
	optimal_merging = false;
elseif nargin==3
	function_name = '';
	optimal_merging = optimal;
	optimal = 'optimal';
elseif nargin==4
	if ~isequal(lower(optimal), 'optimal')
		error('Unrecognized option "%s".', num2str(optimal));
	end
else
	error('Wrong number of input arguments.');
end

%% deal with arrays
if numel(obj)>1
	for i = 1:numel(obj)
		U(i) = U(i).merge(function_name, 'optimal', optimal_merging);
	end
	return
end

%% validation
if ~ischar(function_name)
	error('The function name must be a string.');
elseif ~isempty(function_name) && ~obj.hasFunction(function_name)
	error('No such function "%s" in the object.', function_name);
elseif isempty(function_name) && ~isempty(obj.listFunctions)
	error('Function name must be specified.');
end


%% merging
if nargout==0
	% no outputs = in-place merging
	U = obj;
else
	% output requested = do not modify the source object
	U = obj.copy();
end
% if there is 0 or 1 set contained, return
if U.Num<=1
    return
end

if isempty(function_name)
	% merging of regions (not considering functions)
	if optimal_merging
		U.Set = sub_optimalMerging(U);
	else
		U.Set = sub_greedyMerging(U);
	end
else
	% We are going to merge against a specified function.
	Rmerged = [];
	
	% Get list of regions in which the function has the same expression
	[uF, uMap] = U.Set.uniqueFunctions(function_name);
	if MPTOPTIONS.verbose >= 0
		fprintf('%d regions with %d unique functions.\n', U.Num, length(uF));
	end
	for i = 1:length(uF)
		% For each subset of "identical" regions, create a new union and
		% merge it
		regions = find(uMap==i);
		R = U.Set(regions);
		if MPTOPTIONS.verbose >= 0
			nbefore = length(R);
			fprintf('Function #%d: %d -> ', i, nbefore);
		end
		if length(R)>1
			% Make a copy of the regions, since we are going to remove
			% functions
			Ui = PolyUnion(Polyhedron(R));

			% To prevent PolyUnion/reduce from determining overlapping
			% status, we copy the information already provided to us.
			Ui.Internal.Overlaps = U.Internal.Overlaps;
			
			% Remove all functions, since the pure merge() code introduces
			% new regions, but does not automatically attach functions to
			% them.
			Ui.removeAllFunctions();
			Ui.merge('optimal', optimal_merging);
			R = Ui.Set;
			% Reattach functions
			for j = 1:numel(R)
				% Note that all functions but "FuncName" are meaningless
				% after merging.
				R(j).copyFunctionsFrom(U.Set(regions(1)));
				% Replace the main function
				R(j).addFunction(uF(i), function_name);
			end
		end
		Rmerged = [Rmerged; R];
		if MPTOPTIONS.verbose >= 0
			nafter = length(R);
			if nafter < nbefore
				fprintf('%d (reduction by %d)\n', nafter, nbefore-nafter);
			else
				fprintf('%d\n', nafter);
			end
		end
	end
	if MPTOPTIONS.verbose >= 0
		fprintf('Reduction by %.0f%% from %d to %d regions.\n', ...
			(U.Num-length(Rmerged))/U.Num*100, U.Num, length(Rmerged));
	end
	U.Set = Rmerged;
end

if numel(U.Set) < numel(U.Domain)
	% if the merged set is simpler than the original domain, use the
	% former as the new domain
	U.Domain = U.Set;
end

end

%% optimal merging subroutine
%------------------------------------------------------------------
function N = sub_optimalMerging(U)
% merges regions optimally

Hull = U.convexHull();
Complement = Hull \ U.Set;
% 0 = overlapping regions
% Inf = non-overlapping regions
Algorithm = 0;
N = mpt_optMerge(U.Set, struct('PAdom', Hull, 'PAcompl', Complement, ...
	'algo', Algorithm));

end

%% greedy merging subroutine
%------------------------------------------------------------------
function N = sub_greedyMerging(U)
% greedy merging

global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end
Options.trials = MPTOPTIONS.modules.geometry.unions.PolyUnion.merge.trials;

% try to reduce the union first
U.reduce;

% compute list of neighbours
Ri.BC = sub_buildBClist(U.Set);
% store polyhedra inside Ri
Ri.Pn = U.Set;
% put the total number of polyhedra inside Ri
Ri.nR = U.Num;

% initialize the current best solution
Ri_best.nR = inf;
Ri_min = inf; Ri_max = 0;

% loop of trials
trial = 1;
while trial <= Options.trials
    if MPTOPTIONS.verbose>=2,
        fprintf('Trial %i/%i:\n', trial, Options.trials); 
    elseif MPTOPTIONS.verbose>=1,
        fprintf('Trial %i/%i: %i\n', trial, Options.trials, Ri.nR);
    end
    Ri_mer = sub_merge(Ri);
    if MPTOPTIONS.verbose>=1,
        fprintf(' --> %i', Ri_mer.nR); 
    end
    
    % update minimum and maximum
    Ri_min = min(Ri_min, Ri_mer.nR);
    Ri_max = max(Ri_max, Ri_mer.nR);
    if Ri_mer.nR < Ri_best.nR
        Ri_best = Ri_mer;
        trial = trial + 1;
        if MPTOPTIONS.verbose>=1,
            fprintf(' (new minimum)\n'); 
        end
    else
        trial = trial + 1;
        if MPTOPTIONS.verbose>=1,
            fprintf('\n'); 
        end
    end;
    if Ri_min<=1,
        break
    end
end;

N = Ri_best.Pn;

if MPTOPTIONS.verbose>=1,
    fprintf('  ==> min: %i  max: %i\n', Ri_min, Ri_max');
end

end

%-----------------------------------------------------------------------
function Ri = sub_merge(Ri)

global MPTOPTIONS

% start the trial with the permutation K
iterateAgain = 1;
iter = 0;
while iterateAgain
        
    iterateAgain = 0;
    iter = iter+1;
    if MPTOPTIONS.verbose>=1, 
        fprintf('  iteration %i: merging %i --> ', iter, Ri.nR); 
    end;

    % indicators whether region k has been (1) included in some union or not (0)
    issorted=zeros(Ri.nR,1);
    newRi.nR=0;
    newRi.Pn = [];
    
    % find a random permutation of the indices (of polyhedra)
    K = randperm(Ri.nR);

    old2new = K;
    for k = K
        
        Pc = Ri.Pn(k);
        BCc=setdiff(Ri.BC{k},0);
        changed=0;
        if issorted(k)
            continue;
        end
        firstloop=1;
        while firstloop || changed
            firstloop=0;
            changed=0;
            BCc_old=BCc;
            BCc_old=BCc_old(randperm(length(BCc_old))); % permute BBc_old
            for ind_l=1:length(BCc_old), %1:Ri.nR
                ind2=BCc_old(ind_l);

                if(ind2==k || issorted(ind2))
                    continue;
                end
                
                PU = PolyUnion([Pc, Ri.Pn(ind2)]);
                how = PU.isConvex;
                %[Pu,how] = union([Pc Ri.Pn(ind2)],Options);
                if how
                    Pc = PU.Internal.convexHull;
                    BCc=union(BCc, setdiff(Ri.BC{ind2},0));
                    if MPTOPTIONS.verbose>=2,
                        disp(['regions ' num2str([k ind2]) ' are joined']);
                    end;
                    issorted(ind2)=1;
                    old2new(ind2)=newRi.nR+1;
                    changed=1;
                    iterateAgain = 1;
                end
            end
        end
        newRi.nR=newRi.nR+1;
        old2new(k)=newRi.nR;
        newRi.Pn = [newRi.Pn Pc];
        newRi.BC{newRi.nR}=BCc;
        issorted(k)=1;
    end
    for i=1:newRi.nR
        newRi.BC{i}=unique(old2new(newRi.BC{i}));
    end

    if MPTOPTIONS.verbose>=1,
        percent = (Ri.nR-newRi.nR)/Ri.nR * 100;
        fprintf('%i (%2.1f percent)\n', newRi.nR, percent);
    end;

    Ri = newRi;
    clear newRi;

end;


end

%-----------------------------------------------------------------------
function BC = sub_buildBClist(Pn)
% Inputs:  array of polyhedra in the same dimension
%
% Outputs: BC: list of neighbours

global MPTOPTIONS
nR = length(Pn);

if MPTOPTIONS.verbose>=1
    fprintf('Computing list of neighbors...\n');
end

% M(i, j)=1 means that Pn(i) and Pn(j) are neighbors
M = zeros(nR);
for i = 1:nR-1
	for j = i+1:nR
		if ~M(i,j) && ~M(j,i)
			if Pn(i).isNeighbor(Pn(j));
				M(i,j)=1;
				M(j,i)=1;
			end
		end
	end
end

BC = cell(1,nR);
for i = 1:nR
	BC{i} = find(M(i, :));
end

end
