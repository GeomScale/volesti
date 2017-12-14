%===============================================================================
%
% Title:        Merge5
%                                                             
% Project:      Optimal merging of polyhedra with the same PWA dynamics
%
% Description:  Merges the white polyhedra recursively by dividing the 
%               polyhedron given by m_poly into two polyhedra given by 
%               m_poly1 and m_poly2 and branching on both of them.
%               The procedure is drastically sped up by using bound 
%               techniques and by only branching on hyperplanes that are 
%               supporting hyperplanes of black polyhedra.
%
% Inputs:       m_poly: marking of the polyhedron we are dealing with. 
%                       The components of the marking vector are in 
%                       {-1, 0, 1}. A 0 means, that the corresponding 
%                       hyperplane is not fixed.
%               M_white: markings of the white polyhedra
%               M_black: markings of the black polyhedra
%               level: level in the search tree
%               n_curr: current number of merged white polyhedra
%               n_up: upper bound on the number of merged white polyhedra
%
% Outputs:      M_mer: markings of the merged white polyhedra
%               nodes: number of nodes visited 
%
% Authors:      Tobias Geyer <geyer@control.ee.ethz.ch>, Fabio D. Torrisi

% History:      date        ver.    subject 
%               2002.11.xx  1.0.0   initial version 
%               2003.03.03  1.0.1   simplified updating of M_mer and n_up
%               2003.07.02  1.1     branch first over hyperplanes that separate 
%                                   as many white from black polyhedra as possible
%                                   speed up distribution of polyhedra
%               2003.08.25  1.2     improved branching rule (exploiting +/- 0.5 markings)
%               2003.09.07  1.3     remove black polyhedra, that are not contained 
%                                   in the envelope of white polyhedra
%
% Requires:     -
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

function [M_mer, nodes] = mpt_merge5(m_poly, M_white, M_black, level, n_curr, n_up, maxNodes)  

if nargin < 7
    maxNodes = inf;
end;

show = 0;

nodes = 1;

% get the envelope of the white polyhedra
% the i-th element of the envelope is 
%   -1 if all white polyhedra have a -0.5 or -1 marking as i-th element
%   +1 if all white polyhedra have a 0.5 or 1 marking as i-th element
%   0  else
% to speed up things, we derive the envelope by refining m_poly, i.e. we
% check if 0-elements of m_poly can be set to +/-1.
m_env = m_poly; 
poly_zeroInd = find(m_poly==0);
for i=poly_zeroInd'
    if all(round(M_white(i,:))==-1)
        m_env(i) = -1;
    elseif all(round(M_white(i,:))==+1)
        m_env(i) = +1;
    end;
end;

% remove black polyhedra that are not in the envelope
env_fixedInd = find(m_env~=0);
if env_fixedInd
    rem_ind = find(any( round(M_black(env_fixedInd,:)) ~= repmat(m_env(env_fixedInd), 1, size(M_black,2)) ));
    M_black(:,rem_ind) = [];
end;


if show
    for i=1:level, fprintf('             '); end;
    fprintf('  %i: %s (n_curr=%i, n_up=%i)', level, int2str(m_env)', n_curr, n_up);
end;

% if there are no white and no black polyedra, ... 
if isempty(M_white) & isempty(M_black)
    % ... we have found a non-existing one
    M_mer = [];
    if show, fprintf(' ...x\n'); end;
    
    % if there are no white polyhedra but only black ones, ...
elseif isempty(M_white) & ~isempty(M_black)
    % ... we have found one or more black polyhedra
    M_mer = [];
    if show, fprintf(' ...b\n'); end;
    
    % if there are white polyedra and no black ones, ...
elseif ~isempty(M_white) & isempty(M_black)
    % ... we have found one or more white polyhedra
    M_mer = m_env;
    if show, fprintf(' ...w\n'); end;
    
else
    % we have white and black polyhedra in our current polyhedron
    if show, fprintf('\n'); end;
    
    % index vector of open hyperplanes (0=fixed, 1=open)
    % (i.e. not yet fixed hyperplanes have marking 0)
    open_hyp = ( m_env == 0 );
    
    % index vector of non-redundant hyperplanes of black polyhedra 
    % (0=redundant, 1=non-redundant)
    % a hyperplane is redundant, if it has only -0.5 or 0.5 entries in all markings
    non_red = ( sum(abs(M_black), 2) > 0.5*size(M_black, 2) );
    
    % index set of hyperplanes on which we want to branch, i.e. hyperplanes, that
    % * haven't been fixed yet and
    % * are in terms of the black polyhedra non-redundant (this we avoid 
    % hyperplanes that cut only through white polyhedra and no facet of a 
    % black polyhedron is a subset of such a hyperplane)
    branch_hyp = find( open_hyp+non_red >= 2 );
    
    if length(branch_hyp) == 0
        warning('we should not be here.')
        % as all hyperplanes have been fixed and there is more than one bad guy
    end;
    
    
    % we branch in the following order:
    % 1.) on hyperplanes, that separate two non-connected groups of white 
    %     polyhedra. We are only interested in the first hyperplane that
    %     has these properties. Thus we can cut our problem into two.
    % if we don't have such a hyperplane, we branch 
    % 2.) on hyperplanes, such that in one halfspace are only white or 
    %     black polyhedra. If there are in one halfspace only white (black)
    %     polyhedra, we want to have a hyperplane with as many white
    %     (black) polyhedra as possible.
    % if we don't have such hyperplanes, we branch 
    % 3.) on hyperplanes, that are separating as many white from black
    %     polyhedra as possible (thus number of white facets plus number of
    %     black facets is the criterion). For this, we need redundancy
    %     information (+/- 0.5 markings).
    
    if length(branch_hyp) > 2
        [T, found_sep] = branchHeuristic_Sep(branch_hyp, M_white, M_black);
        if ~found_sep
            [T, found_only] = branchHeuristic_Only(branch_hyp, M_white, M_black);
%             if ~found_only
%                 T = branchHeuristic_asMany(branch_hyp, M_white, M_black);
%             end;
        end;
        branch_hyp = T(:,1);
    end;
    
    
    M_mer = []; % matrix keeping the markings of the merged polyhedra 
    % which are currently best (minimal)
    
    % Branch over all open hyperplanes
    % and put the resulting merged polyhedra in the structure M_mer{}.
    % Bound if n_curr >= n_up, where n_curr is the current number of merged
    % polyhedra and n_up is an upper bound.
    while ~isempty(branch_hyp) & (n_curr < n_up) & (nodes < maxNodes) %& (level < 5) & (nodes < 1000)
         
        % get current hyperplane
        hyp = branch_hyp(end);
        branch_hyp(end) = [];
        
        % polyhedron1: fix the marking to -1
        m_poly1 = m_poly; m_poly1(hyp) = -1;
        
        % polyhedron2: fix the marking to 1
        m_poly2 = m_poly; m_poly2(hyp) = +1;
        
        % distribute the white and black polyhedra to m_poly1 and m_poly2
        M_white1 = M_white( :, find(round(M_white(hyp,:)) == -1));
        M_black1 = M_black( :, find(round(M_black(hyp,:)) == -1));
        M_white2 = M_white( :, find(round(M_white(hyp,:)) == +1));
        M_black2 = M_black( :, find(round(M_black(hyp,:)) == +1));
        
        %plotInterm(M_white1, M_black1, M_white2, M_black2, hyp)
        
        % branch
        [M_mer1, nodes1] = mpt_merge5(m_poly1, M_white1, M_black1, level+1, n_curr, n_up, maxNodes);
        
        if n_curr+size(M_mer1,2) < n_up
            [M_mer2, nodes2] = mpt_merge5(m_poly2, M_white2, M_black2, level+1, n_curr+size(M_mer1,2), n_up, maxNodes);
        else 
            M_mer2 = []; nodes2 = 0;
        end
        % comment: It is not 100% clear, if the restriciton above would interfere
        % with the optimality of the algorithm. It saves hardly time, so better
        % keep it switched off.
        
        % update the nodes counter
        nodes = nodes + nodes1 + nodes2;
        
        % display
        if show==2
            disp('results from branching:')
            disp([M_mer1 M_mer2])
        end;
        
        % check if we have found a better branching result, 
        % i.e. if the new branching result [M_mer1 M_mer2] has less resulting 
        % merged white polyhedra then the current one M_mer.
        % If M_mer is empty, i.e. trial==1, we update in any case
        if isempty(M_mer) | (size([M_mer1 M_mer2], 2) < size(M_mer, 2))
            % update merged white polyhedra
            M_mer = [M_mer1 M_mer2];
            
            % try to lower upper bound 
            n_up = min( n_up, n_curr + size(M_mer, 2) ); 
        end;
        
    end;
    
    if show
        for i=1:level, fprintf('             '); end;
        fprintf('  --> %s (l=%i, n=%i)\n', int2str(M_mer')', level, nodes); 
    end;
    
end;







function [T, found] = branchHeuristic_Only(branch_hyp, M_white, M_black);

% 1.) on hyperplanes, such that in one halfspace are only white or 
%     black polyhedra. If it has in one halfspace only white (black)
%     polyhedra, we want to have a hyperplane with as many white
%     (black) polyhedra as possible.

% build table T, where
%       each row holds an entry
%       columns: 1=hyp, 2=SepHyp * #otherColor, 3=max(#white&black facets), 
T = NaN*ones(length(branch_hyp), 2);

% for 1.)
for i=1:length(branch_hyp)
    hyp = branch_hyp(i);
    
    % we know, that all white and black polyhedra are within m_poly
    % thus, it is enough, to find the number of white and black
    % polyhedra that have at position 'hyp' a + or - entry. This gives
    % us the number of white and black polyhedra on both sides of the
    % hyp-th hyperplane.
    n_white1 = length(find(round(M_white(hyp,:)) == -1));
    n_black1 = length(find(round(M_black(hyp,:)) == -1));
    
    n_white2 = length(find(round(M_white(hyp,:)) == +1));
    n_black2 = length(find(round(M_black(hyp,:)) == +1));
    
    % is there any hyperplane, that has in one of its
    % halfspaces only white or black polyhedra?
    % If it has in one halfspace only white (black) polyhedra, we want to 
    % have a hyperplane with as many white (black) polyhedra as possible.
    SepHyp = (n_white1==0)*n_black1 + (n_black1==0)*n_white1 + (n_white2==0)*n_black2 + (n_black2==0)*n_white2;
    
    % update table
    T(i,:) = [hyp, SepHyp];
end;

% sort table according to branching rule 
T = sortrows(T, 2);

% if we have found hyperplanes, that are separating, remove all the
% non-separating hyperplanes
if T(end,2) >= 1
    % remove all non-separating hyperplanes
    %T(T(:,2)==0, :) = [];
   
    found = 1;
else
    found = 0;
end;

% randomize the order of the non-separating hyperplanes
ind = find(T(:,2)==0);
rand_ind = ind(randperm(length(ind)));
T(ind,:) = T(rand_ind,:);

return


function [T] = branchHeuristic_asMany(branch_hyp, M_white, M_black);

% 2.) on hyperplanes, that are separating as many white from black
%     polyhedra as possible (thus number of white facets plus number of
%     black facets is the criterion). For this, we need redundancy
%     information (+/- 0.5 markings)

% try 2.)
T = NaN*ones(length(branch_hyp), 2);

for i=1:length(branch_hyp)
    hyp = branch_hyp(i);
    % we want to find a hyperplane that separates as many white from
    % black polyhedra.
    
    % find all the white polyhedra that have the hyperplane as a facet
    whiteFacetInd = find( abs(M_white(hyp,:))>0.75 );
    
    % take these white polyhedra and invert the hyperplane. How many of
    % these new polyhedra are black?
    % go through the black polyhedra
    n_wbFacets = length(whiteFacetInd);
    for w=whiteFacetInd
        m_new = M_white(:,w);
        m_new(hyp) = m_new(hyp)*(-1);
        if find(all( round(M_black) == round(m_new)*ones(1,size(M_black,2)) ))
            n_wbFacets = n_wbFacets+1;
        end;
    end;
    
    % update table
    T(i,:) = [hyp, n_wbFacets];
end;

% sort table according to branching rule 
T = sortrows(T, 2);

return



function [T, found] = branchHeuristic_Sep(branch_hyp, M_white, M_black);

% search for hyperplane that cut the problem into two.
% Thus, if we have at least two groups of non-connected white
% polyhedra, we search for a hyperplane, that separates two
% non-connected groups. We are only interested in the first hyperplane that
% has these properties. Thus we can cut our problem into two.

M_white = round(M_white);
M_black = round(M_black);

T = [];

for hyp=branch_hyp(:)'
    % distribute the polyhedra
    M_white1 = M_white( :, find(M_white(hyp,:) == -1));
    M_white2 = M_white( :, find(M_white(hyp,:) == +1));
    
    % go through all white polyhedra on the '-' side. Each neighbour on the
    % '+' side must be a black polyhedra, thus must not be contained in
    % M_white2. If this is fulfilled for all polyhedra in M_white1, we have
    % indeed a hyperplane that cuts our problem into two.
    cutting = 1;
    for w=1:size(M_white1,2)
        m_white1 = M_white1(:,w);
        m_new = m_white1;
        m_new(hyp) = m_new(hyp)*(-1);
        
        % old and slow
%         if find(all( M_white2 == m_new*ones(1,size(M_white2,2)) ))
%             cutting = 0;
%             break;
%         end;
        
        % new and apparently faster
        ind = 1:size(M_white2,2);
        % we only need to check for the hyperplanes in branch_hyp
        for j=branch_hyp(:)'
            ind = ind(find( m_new(j) == M_white2(j,ind) ));
            if isempty(ind), break; end;
        end;
        if ind
            cutting = 0;
            break;
        end;
        
    end;
    
    if cutting
        T = hyp;
        break;
    end;
end;

if T
    found = 1;
else
    found = 0;
end;
return