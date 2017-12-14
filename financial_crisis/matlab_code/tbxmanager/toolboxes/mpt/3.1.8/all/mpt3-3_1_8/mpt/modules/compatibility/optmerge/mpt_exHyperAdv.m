%===============================================================================
%
% Title:        mpt_exHyper
%                                                             
% Project:      Optimal merging of polyhedra with the same PWA dynamics
%
% Input:        Hi, Ki: definition of the polyhedral partition
%               color.Reg: color of each region
%               color.Table: regions with the same color 
%               Opt.lpsolver:
%               Opt.verbose: 0: silent
%                            1: verbose only important information
%                            2: verbose everything
%               Opt.cdd: 0: don't use cdd
%                        1: use cdd
%
% Output:       M: set of markings. The entries have the following meaning: 
%                  -1:   polyhedron is on the - side of the hyperplane; 
%                        the hyperplane is a facet of the polyhedron
%                  -0.5: polyhedron is on the - side of the hyperplane;
%                        the hyperplane is redundant
%                   0:   polyhedron is cut by the hyperplane
%                  +0.5: polyhedron is on the + side of the hyperplane;
%                        the hyperplane is redundant
%                  +1:   polyhedron is on the + side of the hyperplane;
%                        the hyperplane is a facet of the polyhedron
%               color: as for the input
%               HypArr: hyperplane arrangement Hi*[xr; ur] = Ki
%               dom: domain (redundant hyperplanes have been removed)
%
% Authors:      Tobias Geyer <geyer@control.ee.ethz.ch>, Fabio D. Torrisi

% (C) 2004 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2003-2004 Tobias Geyer, Automatic Control Laboratory, ETH Zurich,
%          geyer@control.ee.ethz.ch

% ---------------------------------------------------------------------------
% Legal note:
%          This program is free software; you can redistribute it and/or
%          modify it under the terms of the GNU General Public
%          License as published by the Free Software Foundation; either
%          version 2.1 of the License, or (at your option) any later version.
%
%          This program is distributed in the hope that it will be useful,
%          but WITHOUT ANY WARRANTY; without even the implied warranty of
%          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%          General Public License for more details.
% 
%          You should have received a copy of the GNU General Public
%          License along with this library; if not, write to the 
%          Free Software Foundation, Inc., 
%          59 Temple Place, Suite 330, 
%          Boston, MA  02111-1307  USA
%
% ---------------------------------------------------------------------------

% History:      date        ver.    subject 
%               2003.03.11  1.0.0   initial version 
%               2003.05.15  1.1     outsource part of the code into HyperList
%                                   this is more elegant but not necessarily faster
%               2003.06.22  1.2     works now also for non-convex polyhedral partitions
%               2003.07.01  1.2.1   works now also if initial hyperplane arrangement is empty
%               2003.07.02  1.2.2   mark redundant hyperplanes as -0.5 / 0.5
%               2003.08.14  1.3     increased computational efficiency
%                                   (less LPs for extending hyperplanes)
%               2003.08.18  1.4     avoid usage of hyperList (--> faster)
%               2003.08.23  1.5     only based on colors (not dynamics anymore);
%                                   ensure that redundant hyperplanes have a +/- 0.5
%               2004.04.02  2.0     rewritten for MPT
%               2004.04.23  2.1     numerically more stable,
%                                   accepts now also overlapping polyhedra,
%                                   uses hyparr to compute markings
%                                   allows to cluster and simplify hyperplanes
%               2004.06.16  2.2     options to set maximal size of HA
%
% Requires:     tg_polyreduce,
%               polyinnerball, 
%               lpsolve
%
% Contact:      Tobias Geyer
%               Automatic Control Laboratory
%               ETH Zentrum, CH-8092 Zurich, Switzerland
%
%               geyer@control.ee.ethz.ch
%
%               Comments and bug reports are highly appreciated
%
%===============================================================================

function [M, color, HypArr, dom] = mpt_exHyperAdv(Hi, Ki, PA, color, PAhull, Opt)

% According to Ziegler, we define the markings in the following way:
% Example: delta = [-1 -1 +1;
%                   -1 +1 +1] 
%          implies the existence of 2 hyperplanes generating 3 polyhedra: 
%               1. A x <= B;
%               2. A(1,:) <= B(1),  A(2,:) >= B(2);
%               3. A(1,:) >= B(1),  A(2,:) >= B(2);
% This is according to the definition of the paper, but contrary to hys2pwa!


% tolerances: 
% -------------------------------------------------------------------------

% tol.facet.coeff: when dealing with coefficients of hyperplanes
if ~isfield(Opt, 'tol_facet_coeff')
    Opt.tol_facet_coeff = 1e-6;     
end;

% tol.facet.round: when rounding facet coefficients to this precision
if ~isfield(Opt, 'tol_facet_round')
    Opt.tol_facet_round = 1e-9;    
end;

% tol.polyreduce: when using polyreduce 
if ~isfield(Opt, 'tol_polyreduce')
    Opt.tol_polyreduce = 1e-6;      
end;

% tol.minR: minimal required radius of Chebycheff ball
if ~isfield(Opt, 'tol_minR')
    Opt.tol_minR = 1e-5;            
end;

% tol.simplifyHA: 1-norm of cluster of hyperplanes from center 
% clustering is done in order to simplify the hyperplane arrangement
% and to thus reduce the complexity of the solution. Clustering is
% beneficial, if we have many hyperplanes that are almost the same. 
% Then we can take a weighted average of such a cluster.
% When the tolerance is 0, no clustering is performed
if ~isfield(Opt, 'tol_simplifyHA')
    Opt.tol_simplifyHA = 0;
end;

% maximal size of hyperplane arrangement
if ~isfield(Opt, 'maxHA')
    Opt.maxHA = inf;
end;

if ~isfield(color, 'Reg')
    col = color;
    color = [];
    color.Reg = col;
    clear col
end;

% remarks:
% if we have double hyperplanes: reduce tol_round
% if we are missing some hyperplanes: increase tol_round


% Overview of Algorithm:
% ----------------------
% missing...


if Opt.verbose>1,
    disp(' ');
    disp(' *** extend hyperplanes ***');
    disp(' ');
end


% -1.) restrict polyhedra to domain and remove flat polyhedra
% -------------------------------------------------------------------------

if Opt.verbose>1, 
    disp('  restrict polyhedra to domain and remove flat polyhedra');
end;

keep = [];
for i = 1:length(PA)
%     if dointersect(PA(i), PAhull)
%         keep(end+1) = i;
%     end
    PP = PA(i).intersect(PAhull);
	rc = PP.chebyCenter.r;
    if rc >= Opt.tol_minR, 
        PA(i) = PP;
        keep(end+1) = i; 
    end;
end;
if Opt.verbose>1,
    fprintf('  ... removed %i flat polyhedra\n', length(PA)-length(keep));
end;
PA = PA(keep);
Hi = {};
Ki = {};
for i = 1:length(PA)
	PA(i).minHRep();
	Hi{i} = PA(i).A;
	Ki{i} = PA(i).b;
end
color.Reg = color.Reg(keep);



% 0.) norm all facets (of the polyhedra and of the domain) using the 2-norm
% -------------------------------------------------------------------------

if 1 == 0
    
    if Opt.verbose>1, 
        disp('  remove redundant facets of polyhedra');
    end;

    keep = [];
    for i=1:length(Ki)
        nr_old = length(Ki{i});
        [Hi{i}, Ki{i}, polyEmpty, keptrows] = tg_polyreduce(Hi{i}, Ki{i}, Opt.tol_polyreduce, Opt.lpsolver, 1);
        if polyEmpty, error('polyhedron is empty'); end;
        if Opt.verbose, 
            if nr_old-length(keptrows) > 0 & Opt.verbose>1, 
                fprintf('    #%i: removed %i redundant facets\n', i, nr_old-length(keptrows)); 
            end;
        end;
    end;
    
    if Opt.verbose>1
        disp('  norm facets of polyhedra and domain'); 
    end;
    
    dim = size(Hi{1}, 2);   % dimension (number of variables per hyperplane)
    
    % (i) polyhedra:
    for i=1:length(Ki)
        n = sqrt( sum(Hi{i}.*Hi{i}, 2) );       % vector of 2-norms
        Hi{i} = Hi{i} .* repmat(1./n, 1, dim);
        Ki{i} = Ki{i} ./ n;
    end;

end;


if length(Ki) < 2
    if Opt.verbose>1, 
        fprintf('  terminating mpt_exHyperAdv because number of polyhedra smaller than 2\n');
    end;
    
    M = [];
    color = [];
    HypArr.Hi = [];
    HypArr.Ki = [];
    dom.A = [];
    dom.B = [];
    
    return
end;


% 1.) collect all facets
% -------------------------------------------------------------------------

% Collect all (normed) facets of the polyhedra (whose first non-zero 
% element is positive). Sort them. This yields Hyp.

if Opt.verbose>0
    fprintf('  collect all facets\n');
end;

% dimension of state-input space
nxu = size(Hi{1}, 2);

% collect all facets
% determine length of this collection (to allocate memory for it)
Nr = 0; % number of rows
for i = 1:length(Ki)
    Nr = Nr + length(Ki{i});
end;

% allocate memory for collection, 
% where A*x = B is in the collection given by [A B]
FacColl = NaN*ones(Nr, nxu+1);

% collect the facets    
pos = 1;    % position
for i = 1:length(Ki)
    % facets
    A = Hi{i};
    B = Ki{i};
    
    % determine the signs of the first non-zero entries of A
    sign_a = NaN*ones(length(B), 1);
    for j = 1:length(B)
        n = 1;
        while abs(A(j,n)) <= Opt.tol_facet_coeff
            n = n+1;
            if n>nxu, error('a is all zero'); end;
        end;
        % the n-th entry differs now by more than tol from 0
    
        % determine its sign
        sign_a(j) = sign(A(j,n));  
        if abs(sign_a(j)) <= Opt.tol_facet_coeff, error('sign is zero'); end;
    end;
    
    % write 
    Ind = pos:pos+length(Ki{i})-1;
    pos = pos+length(Ki{i});
    FacColl(Ind,:) = [A B].*repmat(sign_a(:), 1, nxu+1);
end;

% round FacColl to a certain tolerance
% remark: This avoids, that sortrow runs into numerical problems, like
% sorting like that:
%   -0.00000000000000   1.00000000000000   3.75000000000000
%                   0   1.00000000000000 -10.00000000000000
%                   0   1.00000000000000 -10.00000000000000
%                   0   1.00000000000000   3.75000000000000
FacColl = round(FacColl/Opt.tol_facet_round) * Opt.tol_facet_round;

% sort the collection
FacCollSorted = sortrows(FacColl);

% get indices of facets without duplicates
KeepInd = 1:length(FacCollSorted);  % keep facets with these indices
for i = 1:length(FacCollSorted)-1
    if abs(FacCollSorted(i,:)-FacCollSorted(i+1,:)) <= Opt.tol_facet_coeff, 
        KeepInd(find(KeepInd==i+1)) = [];
    end;
end;

% this gives us the preliminary hyperplane arrangement
Hyp.A = FacCollSorted(KeepInd,1:nxu);
Hyp.B = FacCollSorted(KeepInd,nxu+1);

if Opt.verbose>0, fprintf('  ... found %i different hyperplanes\n', length(Hyp.B)); end;




% 2.) check size of HA
% -------------------------------------------------------------------------

% We check here the number of the hyperplanes in the HA. If the HA is 
% absolutely (=hopelessly) too large, we stop the code here to avoid a never 
% terminating computation.
% We expect, that the master program (like mpt_optMergeDivCon) then divides
% the intractable problem into subproblems.

if length(Hyp.B) > Opt.maxHA
    fprintf('  terminating mpt_exHyperAdv because HA is too large\n');
    
    M = [];
    color = [];
    HypArr.Hi = Hyp.A;
    HypArr.Ki = Hyp.B;
    dom.A = [];
    dom.B = [];
    
    return
end;



% 2.) collect information about relative position of facets
% -------------------------------------------------------------------------

if Opt.verbose>0
    fprintf('  collect information about relative position of facets\n');
end;

% Go through all hyperplanes of the hyperplane arrangment Hyp and determine 
% the relative position of each polyhedron with respect to that hyperplane 
% (cuts, on - or on + side, a facets / no facets). This matrix M has as many 
% rows as hyperplanes and as many columns as polyhedra. 
% The entries have the following meaning: 
%	-1:   polyhedron is on the - side of the hyperplane; 
%         the hyperplane is a facet of the polyhedron
%	+1:   polyhedron is on the + side of the hyperplane;
%         the hyperplane is a facet of the polyhedron
%   NaN:  else (no facet or cuts through)
% example: M(4,:) = [-1 NaN NaN 1 1 NaN]: 
%      We have 6 polyhedra. 
%      The polyhedron #1 is on the - side of the 4-th hyperplane, the 
%      polyhedra #4 and #5 are on the + side. The 4-th hyperplane is a
%      facet of the polyhedra #1, #4 and #5.

% go through all hyperplanes
M = NaN*ones(length(Hyp.B), length(Ki));
for i=1:length(Hyp.B)
    
    %if Opt.verbose, fprintf(' #%i', i); end;
    
    % consider the i-th hyperplane a*x = b
    a = Hyp.A(i,:);
    b = Hyp.B(i);
    
    % go through all polyhedra
    for p = 1:length(Ki)
        % is the i-th hyperplane a facet of the p-th polyhedron?
        % all facets and hyperplanes are normed. 
        % thus, we can simply compare them ...
        H = repmat([Hyp.A(i,:) Hyp.B(i)], length(Ki{p}), 1);
        sameInd = find(sum(abs([Hi{p} Ki{p}] - H), 2) < Opt.tol_facet_coeff);
        oppoInd = find(sum(abs([Hi{p} Ki{p}] + H), 2) < Opt.tol_facet_coeff);
        if length(sameInd) > 1, warning('polyhedron has redundant facets'), end;
        if length(oppoInd) > 1, warning('polyhedron has redundant facets'), end;
       
        if sameInd
            % the facet with index 'sameInd' is the same as the i-th
            % hyperplane
            same = 1;
        elseif oppoInd
            % the facet with index 'oppoInd' is the opposite of the i-th
            % hyperplane
            same = -1;
        else
            same = 0;
        end
        if length(sameInd) > 1, sameInd = sameInd(1); end;
        if length(oppoInd) > 1, oppoInd = oppoInd(1); end;
        if sameInd & oppoInd, error('facet can not point in two directions'); end;
        
        if same
            % yes, it is
            % --> add the position of the polyhedron to Hpos
            M(i,p) = (-1)*same;      % -1, +1 (non-redundant)
        end;
        
    end;
end;




% 3.) remove hyperplanes that are only facets of polyhedra with same color
% -------------------------------------------------------------------------

% which hyperplanes can we remove?
% (hyperplanes that touch only polyhedra with the same color) 
keepInd = [];
for i = 1:size(M,1)
    h = M(i,:);     % markings of i-th hyperplane

    % indices of polyhedra that share this hyperplane as facet
    indFac = find(abs(h)>0.75); 
    
    % check if all these polyhedra have the same color
    % if yes, remove hyperplane (if no, keep it)
    if length(unique(color.Reg(indFac)))>1 % number of colors is > 1
        keepInd(end+1) = i;
    end;
end;

% make sure that at least one hyperplanes is left
% (in order to avoid an error later)
if isempty(keepInd), 
    keepInd = [1:min(8,length(Hyp.B))];
end;

% remove 'unnecessary' hyperplanes from the arrangement
Hyp.A = Hyp.A(keepInd,:);
Hyp.B = Hyp.B(keepInd);

% build reduced hyperplane arrangement and remove duplicates
Mred = M(keepInd,:);
[Mdummy, I, J] = unique(round(Mred)','rows');

% adapt color accordingly
% map: Mred(:,i) ---I(i)---> M(:,I(i))
% map: M(:,i) ---J(i)---> Mred(:,J(i))
colorRed.Reg = color.Reg(J);

if Opt.safeMode
    % check that only markings with the same color are mapped to a new marking
    for i=1:size(Mred,2)
        ind = find(J==i);
        % old markings with indices ind are mapped on new marking with index i
        % do they all have the same color?
        if length(unique(color.Reg(ind))) > 1,
            error('mpt_exHyperAdv: markings with different color are mapped onto each other');
        end;
    end;
end;

if Opt.verbose>1,
    fprintf('  ... removed %i ''unnecessary'' hyperplane(s) \n', size(M,1)-length(keepInd));
    fprintf('  ... removed %i marking(s)\n', size(M,2)-size(Mred,2));
end;




% 3.5.) check size of HA
% -------------------------------------------------------------------------

% If the Options 'maxHA' has been set, we check here the number of the 
% hyperplanes in the HA. If the HA is too large, we stop the code here to
% avoid a very time consuming and probably never terminating computation.
% We expect, that the master program (like mpt_optMergeDivCon) then divides
% the intractable problem into subproblems.

if length(Hyp.B) > Opt.maxHA
    if Opt.verbose>0, 
        fprintf('  terminating mpt_exHyperAdv because HA is too large\n');
    end;
    
    M = [];
    color = [];
    HypArr.Hi = Hyp.A;
    HypArr.Ki = Hyp.B;
    dom.A = [];
    dom.B = [];
    
    return
end;



% 4.) analyze and simplify hyperplane arrangement
% -------------------------------------------------------------------------

if Opt.tol_simplifyHA > 0

    if Opt.verbose>0,
        fprintf('  find clusters (with 1-norm %0.5f) of similar hyperplanes\n', Opt.tol_simplifyHA);
    end
    
    % cluster hyperplanes (clusters can have size of 1)
    HA = [Hyp.A Hyp.B];
    HypClust = [];
    unclustInd = 1:length(Hyp.B);
    while length(unclustInd) > 0
        % start a new cluster by moving the first unclustered hyperplane
        % to the new cluster
        HypClust{end+1} = unclustInd(1);
        unclustInd(1) = [];
        
        % check all unclustered hyperplanes and move them to the cluster
        % if the distance (1-norm) from the center (linear weighted) of the cluster
        % to the new hyperplane is smaller than a threshold
        foundNew = 1;
        while foundNew & length(unclustInd) > 0
            foundNew = 0;
            
            clustCenter = sum(HA(HypClust{end},:),1) / length(HypClust{end});
            for i = 1:length(unclustInd)
                dis = norm(clustCenter-HA(unclustInd(i),:));
                if dis < Opt.tol_simplifyHA
                    HypClust{end}(end+1,:) = unclustInd(i);
                    unclustInd(i) = [];
                    foundNew = 1;
                    break;
                end;            
            end;
        end;
    end;
    
    if Opt.verbose>1
        fprintf('  ... clusters in hyperplane arrangement\n');
    end;
    numClust = 0;
    for c=1:length(HypClust)
        if Opt.verbose>1
            fprintf('     %2i-th cluster:', c);
            for i=HypClust{c}(:)'
                fprintf(' [ ');
                fprintf('%0.6f ', Hyp.A(i,:));
                fprintf('] = %0.6f\n', Hyp.B(i));
                if i~=HypClust{c}(end), fprintf('                   '); end;
            end;
        end;
        if length(HypClust{c}) >1, numClust = numClust + 1; end;
    end;
    if Opt.verbose>1,
        fprintf('  ... found %i non-trivial clusters\n', numClust);
    end
    
    % replace old hyperplane arrangement by clustered hyperplane arrang.
    HAnew = [];
    for c=1:length(HypClust)
        HAnew(end+1,:) = sum(HA(HypClust{c},:),1) / length(HypClust{c});
    end;
    if Opt.verbose>1,
        fprintf('  ... reduced hyperplane arrangement by %i hyperplanes\n', size(HA,1)-size(HAnew,1));
    end
    
    % have we found a new hyperplane arrangement?
    clustPerf = size(HA,1) ~= size(HAnew,1);
    if clustPerf
        Hyp.A = HAnew(:,1:end-1);
        Hyp.B = HAnew(:,end); 
    end;
    
else
    
    if Opt.verbose>1,
        fprintf('  clustering disabled\n');
    end
    clustPerf = 0;
    
end;



% 5.) get markings of hyperplane arrangement and their colors
%     use hyparr
%     restrict markings to domain given by PAhull
% -------------------------------------------------------------------------

if Opt.verbose>0,
    fprintf('  derive markings of hyperplane arrangement with %i hyperplanes\n', length(Hyp.B));
end

% get hyperplane arrangement 
% note: * M has {-1, 0, 1} entries
%       * radii of Chebycheff balls are required to be larger than Opt.minR (speed-up)
%       * restrict search for markings to domain (speed-up)
%       * for the domain, we use a first guess of domain
%         (we know that the domain is a subset of the convex hull)
dom.H = PAhull.A;
dom.K = PAhull.b;
dom.A = PAhull.A;
dom.B = PAhull.b;
dom.constr = 1;

Opt.verbose = Opt.verbose + 1; 
M = (-1)*mpt_hyparr2(Hyp, dom, dom, Opt);
Opt.verbose = Opt.verbose - 1;

if Opt.verbose>1
    fprintf('  ... hyperplane arrangement has %i markings\n', size(M,2));
    fprintf('  derive colors of cells in hyperplane arrangement\n');
end;



% remove polyhedra and their markings with radii < Opt.minR
% and determine colors and errors
keepInd = [];
colorNew.Reg = [];
errorStat = [];
for i = 1:size(M,2)
    % build polyhedron  
    facetInd = find( abs(M(:,i)) > 0.75 );
    H = [(-1)*diag(M(facetInd,i)) * Hyp.A(facetInd,:); dom.A];
    K = [(-1)*diag(M(facetInd,i)) * Hyp.B(facetInd)  ; dom.B];
       
    % Chebycheff ball in polyhedron
    [x, r] = polyinnerball(H, K);
    if Opt.verbose==2, fprintf('  polyhedron #%i: R=%1.4f\n', i, r); end;
    
    if r >= Opt.tol_minR, 
                
        % get color
        col = NaN;
        if ~clustPerf
            % we can use simple approach:
            % in which original polyhedra does the 'center' of the marking lie?
            [isin, inwhich, closest] = PA.isInside(x);
            if isin
                colStat = unique(color.Reg(inwhich));
                if length(colStat) > 1
                    if Opt.verbose>1,
                        disp('mpt_exHyperAdv: center of ''marking'' lies in polyhedra with different colors');
                    end
                    errorStat(end+1,:) = r;
                else
                    col = colStat(1);
                end;
            else
                col = color.Reg(closest);
            end;
        else
            % we have to use a more complicated approach:
            % find intersections of polyhedra with the one corresponding to the marking
            % one row of statistics = [radius, color];
            
            % fast approach
            stat = [];
            for j = 1:length(Ki)
                [x, r] = polyinnerball([H; Hi{j}], [K; Ki{j}]);
                if r >= Opt.tol_minR
                    stat(end+1,:) = [r, color.Reg(j)];
                end;
            end;
            
            % take the color of the biggest intersection
            % *** doesn't work when intersections of different colors have
            % some radii ***
            %stat = sortrows(stat, 1);
            %col = stat(end,2);
            
            if ~isempty(stat)
                % find for each occuring color the biggest intersection
                C = unique(stat(:,2))';
                statMax = [];
                for c = C
                    rMax = max(stat(find(stat(:,2)==c),1));
                    statMax(end+1,:) = [rMax c];
                end;
                statMax = sortrows(statMax, 1);
                
                % choose:
                % 1.) if we have only one color: take this color
                % 2.) if we have more than one: check if the difference in the
                % radii is significant. If not, we cannot choose. Issue a
                % warning and remove the marking (--> don't care)
                if size(statMax,1) < 1 & Opt.verbose>1
                    disp('mpt_exHyperAdv: cell does not interesect with any polyhedron')
                elseif size(statMax,1) == 1
                    col = statMax(end,2);
                else
                    rDiff = abs( (statMax(end,1)-statMax(end-1,1))/statMax(end-1,1) );
                    if rDiff < 0.05
                        if Opt.verbose>1,
                            fprintf('mpt_exHyperAdv: cannot determine color of cell (normed rDiff = %4.1e) --> remove cell\n', rDiff)
                        end
                        errorStat(end+1,:) = statMax(end,:);
                    else
                        col = statMax(end,2);
                        errorStat(end+1,:) = statMax(end-1,:);
                    end;
                end;
            else
                % intersections are empty or too smal
                % --> no color (don't care)
            end;
            
        end;
        
        if ~isnan(col)
            % keep marking and assign color
            keepInd(end+1) = i;
            colorNew.Reg(end+1) = col;
        end;
        
    end;
end;
if Opt.verbose>0,
    fprintf('  ... removed %i polyhedra from hyperplane arrangement with radii below %2.2e or undeterminable color\n', size(M,2)-length(keepInd), Opt.tol_minR);
    %fprintf('  ... keep %i polyhedra in hyperplane arrangement\n', length(keepInd));
end;
if clustPerf & Opt.verbose>1,
    fprintf('  ... max error due to reduction of hyper. arrangement (radius of polyhedron with wrong color): %1.5f\n', max(errorStat(:,1)));
end;
if Opt.plot & ~isempty(errorStat)
    figure; clf;
    hist(errorStat(:,1));
    xlabel('radii of polyhedra with wrong color');
    ylabel('number of polyhedra');
    title('histogramm of error due to reduction of hyperplane arrangement');
end;
M = M(:,keepInd);




% 6.) remove domain
% -------------------------------------------------------------------------

% Remove hyperplanes that are bounding the domain. These are hyperplanes
% with lines in M that have only -1 and -0.5 or 1 and 0.5 entries. Adapt 
% Hyp accordingly.

if Opt.verbose>0,
    fprintf('  remove hyperplanes bounding the domain\n');
end;

% get indices of hyperplanes that do not belong to the domain
% (hyperplanes with different signs and possibly also 0-elements)
KeepInd = 1:size(M,1);      % keep hyperplanes with these indices
RemInd = [];                % remove hyperplanes with these indices
for i = 1:size(M,1)
    if all(M(i,:)<-0.25) | all(M(i,:)>0.25) 
        % hyperplane bounds only domain
        KeepInd(find(KeepInd==i)) = [];
        RemInd(end+1) = i;
    end;
end;

% add the hyperplanes with indices in RemInd to the domain
if ~isempty(RemInd)
    sign_a = sign(M(RemInd,1));
    dom.A = [dom.A; Hyp.A(RemInd,:)*(-1).*repmat(sign_a, 1, nxu)];
    dom.B = [dom.B; Hyp.B(RemInd)  *(-1).*sign_a];
end;
if Opt.verbose>1,
    fprintf('  ... moved %i hyperplanes from the hyperplane arrangement to the domain\n', length(RemInd));
end;

% and remove redundant hyperplanes from the domain
oldL = length(dom.B);
[dom.A, dom.B, isemptypoly] = tg_polyreduce(dom.A, dom.B, Opt.tol_polyreduce, Opt.lpsolver, 1);
if isemptypoly, error('domain is empty'); end;
if Opt.verbose>1,
    fprintf('  ... removed %i hyperplanes from the domain\n', oldL-length(dom.B));
end;

% remove the hyperplanes from Hyp and M
Hyp.A = Hyp.A(KeepInd,:);
Hyp.B = Hyp.B(KeepInd);
M     = M(KeepInd,:);



% 7.) display and format outputs
% -------------------------------------------------------------------------

if Opt.verbose>0,
    fprintf('  ==> hyperplane arrangement contains %i hyperplanes and %i markings\n', length(Hyp.B), size(M,2)); 
end
    
if (Opt.verbose == 2) & (length(Hyp.B) < 100)
    % show all hyperplanes
    fprintf('\n')
    fprintf('  hyperplane arrangement\n');
    for i=1:length(Hyp.B)
        fprintf('      %i-th hyperplane: [ ', i);
        fprintf('%0.6f ', Hyp.A(i,:));
        fprintf('] = %0.6f\n', Hyp.B(i));
    end;
end;

if Opt.verbose == 2
    disp(' ')
    disp('  domain:')
    disp(dom.A)
    disp(dom.B)
end;

HypArr.Hi = Hyp.A;
HypArr.Ki = Hyp.B;

color = colorNew;

% complement color information
% color.Reg: color of each region
% color.Table: regions with the same color 
colorSorted = sort(unique(color.Reg));
for c = colorSorted; 
    color.Table{c} = find(color.Reg==c);
end;

return
