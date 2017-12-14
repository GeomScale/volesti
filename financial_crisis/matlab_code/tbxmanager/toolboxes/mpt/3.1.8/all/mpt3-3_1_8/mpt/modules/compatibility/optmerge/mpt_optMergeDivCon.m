%===============================================================================
%
% Title:    mpt_optMergeDivCon
%                                                             
% Project:  Optimal merging of a set of polyhedra using divide and
%           conquer strategy
%
% Input:    PA: polyhedral array
%           Opt: Options with the following fields
%              .color: color of polyhedra
%               This is an integer 1,2,... The algorithm will merge
%               polyhedra with the same color. If color is not specified, it
%               is assumed to be 1 for all polyhedra.
%               .PAdom: 
%               the set of polyhedra (the so called domain) the problem is 
%               defined in. If not given, the algorithm assumes the domain is 
%               given by the hull (or envelope if the hull computation fails) 
%               of PA.
%               .PAcompl:  
%               polyhedral array that is within the domain and not 
%               in PA - it is the so called complement. If not specified, 
%               the algorithm assumes PAcompl is empty. Then the 
%               complement is filled up optimally by the merging algorithm.
%               Refer also to the second remark below.
%              .verbose: 
%               0: silent (=default)
%               1: verbose only important information
%               2: verbose everything
%              .algo: 
%               0: overlapping merging (boolean minimization)
%               [1, inf]: optimal merging, with the given upper bound on nodes
%               (default is inf)
%              .maxHA: [1, inf]: maximal number of hyperplanes for which we
%               can compute the markings for sub1 or 2 (default inf)  
%              .maxHA_sub1and2: [1, inf]: maximal number of hyperplanes for
%               which we can compute the markings for the combined problem
%               of sub1 and 2 (default maxHA+10) 
%              .divideSub: 
%               0: use an artificial hyperplane parallel to axis of the
%                  coordinate with the longest edge of the bounding box to
%                  divide a problem into 2 subproblems (fast) 
%               1: go through all hyperplanes, choose the one that minimizes
%                  the cost abs(# GL in Sub1 - # GL in Sub2) + # GL in Sub1 + # GL in Sub2 (slow) (=default)
%              .tol.simplifyHA: 1-norm of cluster of hyperplanes from center 
%               clustering is done in order to simplify the hyperplane
%               arrangement and to thus reduce the complexity of the
%               solution. Clustering is beneficial, if we have many
%               hyperplanes that are almost the same.  Then we can take a
%               weighted average of such a cluster. When the tolerance is
%               0, no clustering is performed.
%
% Output:       PAmer: merged polyhedra array
%               colorMer: colors with the following fields
%                 .Reg:   color of each region
%                 .Table: regions with the same color 
%
% Author:       Tobias Geyer <geyer@control.ee.ethz.ch>
%
% Remarks:      some important remarks and explanations:
%               * the input polyhedral array may be overlapping
%                   First, we identify all existing hyperlanes, remove unecessary ones,
%                   and compute then the markings in the hyperplane arrangement. Then,
%                   we identify the colors of the cells in the arrangement by chosing
%                   one cell, getting the center of it and then identifying the 
%                   polyhedron in PA that includes the center. Thus overlaps do not matter.
%               * the input polyhedral array may contain small holes
%                   Holes may always occur (also in the hyperplane arrangement based on
%                   which we merge). Holes are assigned to 'don't cares' and are filled
%                   with neighbouring polyhedra. Acutally, to reduce complexity, we remove
%                   very small cells (say with radii < 1e-5) from the arrangement. These
%                   are then 'don't cares'. 
%               * we keep in the hyperplane arrangment only hyperplanes that are facets
%                 of polyhedra with different colors
%                   These hyperplanes separate the colors. Shouldn't we keep also the other
%                   hyperplanes? No, as we are allowing for overlapping polyhedra (merging
%                   based on boolean minimization), we only need these 'dividing' hyperplanes.
%                   Adding additional hyperplanes does not reduce the number of resulting
%                   polyhedra. Proof?
%               * we plan to later allow to also simplify the hyperplanes, i.e. to 
%                   aggregate similar hyperplanes in only one.
%
% History:      date        ver.    subject 
%               2004.06.16  1.0     initial version based on mpt_iterMerge
%               2004.07.14  1.1     divide problem into subproblems using 'best' hyperplane
%               2005.07.20  1.2     improved help 
%
% Requires:     mpt_exHyperAdv,
%               espresso or merge5,
%               plotPaCol,
%               mpt2cell,
%               MPT toolbox
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

function [PAmer, colorMer] = mpt_optMergeDivCon(PA, Opt)

global mptOptions

% complement options
% --------------------------------------------------------------------

if nargin < 2,
    Opt = [];
end;
if nargin < 1, 
    disp('mpt_optMergeDivCon   version 1.2   Tobias Geyer 2002-2005');
    return;
end;

% build full (original) problem
if ~isfield(Opt, 'color'), 
    color.Reg = ones(1,length(PA));
else
    color.Reg = Opt.color(:)';
end;

if ~isfield(Opt, 'plot')
    Opt.plot = 0;
end;
if ~isfield(Opt, 'algo'), 
    Opt.algo = 0;
end;
if ~isfield(Opt, 'verbose'), 
    Opt.verbose = 0; 
end;
if isfield(Opt, 'PAdom');
    PAdom = Opt.PAdom;
else
    PAdom = [];
end;
if isfield(Opt, 'PAcompl');
    PAcompl = Opt.PAcompl;
else
    PAcompl = [];
end;
if ~isfield(Opt, 'maxHA');
    Opt.maxHA = 200;
end;
if ~isfield(Opt, 'maxHA_sub1and2');
    Opt.maxHA_sub1and2 = Opt.maxHA + 10;
end;
if ~isfield(Opt, 'divideSub');
    Opt.divideSub = 1;
end;
if ~isfield(Opt, 'letter');
    Opt.letter = '';
end;
fname = tempname;
Opt.logfile = '';
% if ~isfield(Opt, 'logfile');
%     Opt.logfile = [fname '_log.txt'];
% end;
if ~isfield(Opt, 'infile');
    Opt.infile = [fname '_in.txt'];
end;
if ~isfield(Opt, 'outfile');
    Opt.outfile = [fname '_out.txt'];
end;
if ~isfield(Opt, 'tol_simplifyHA')
    Opt.tol_simplifyHA = 0;       
end;


Opt.safeMode = 0;   
Opt.lpsolver = mptOptions.lpsolver;
if ismac
    Opt.espressoCmd = which('espresso.mac');
elseif ispc
    Opt.espressoCmd = which('espresso.exe');
else
    Opt.espressoCmd = which('espresso.linux');
end



% prepare
% --------------------------------------------------------------------
% call mpt_hyparr2 with no input arguments => it will clear it's persistent
% variable, otherwise we get errors
mpt_hyparr2;


% check colors:
% each polyhedron must have a color which is integer
if length(PA) ~= length(color.Reg),
    error('the number of polyhedra and colors must be the same')
end;
if any(round(color.Reg)~=color.Reg), 
    error('the colors must be integers')
end;
colorSorted = sort(unique(color.Reg));
if min(colorSorted) < 1 or max(colorSorted) > size(colorSorted)
    error('colors need to be 1, 2, ...')
end;    

% if the user has not specified the domain,
% we quickly compute here the hull or the envelope
if isempty(PAdom)
    try
        disp('computing convex hull')
        PAdom = hull(PA);
    catch
        disp('error when computing convex hull, trying to compute envelope instead');
        PAdom = envelope(PA);
    end;
    disp('Warning: you have not specified the set of polyhedra (the so called');
    disp('domain) the problem is defined in. We assume the domain is given by');
    disp('the hull (or envelope if the hull computation fails).');
end;

% parts, that are within the domain and not in PA, are the so called complement
% Often, the complement is empty. If not specified, we assume it is empty.
% Then the complement is filled up optimally by the merging algorithm.
if isempty(PAcompl)
    disp('Warning: you have not specified the set difference between the set of');
    disp('polyhedra to be merged and the domain. We thus assume it is empty.');
    disp('Note that a non-empty not specified complement is filled up by the');
    disp('merging algorithm optimally!');
else
    % add the complement to the polyhedra array PA and associate an
    % auxiliary color to the complement
    colorSorted = sort(unique(color.Reg));
    maxCol = colorSorted(end);
    PA = horzcat(PA, PAcompl);
    color.Reg = [color.Reg, ones(1,length(PAcompl))*(maxCol+1)];
end;

% is there more than one color? 
% If not, there is nothing to be merged...
if length(unique(color.Reg)) == 1
    if Opt.verbose>0,
        disp('only one color: nothing to merge - giving back the domain')
    end
    PAmer = PAdom;
    colorMer.Reg = color.Reg;
    c = color.Reg(1);
    colorMer.Table{c} = [1:length(PAdom)];
    return;
end;

if ~isempty(Opt.logfile)
    % open logfile
    Opt.fname = [Opt.logfile '.log'];
    Opt.fid = fopen(Opt.fname, 'wt');
    fprintf(Opt.fid, ['logfile of mpt_optMergeDivCon\n\n']);
end;

% start to measure the time
t_total = cputime;



% merge using divide and conquer technique
% --------------------------------------------------------------------

Prob.PA = PA;
%Prob.color = color.Reg;
Prob.color = color;     % thus we have Prob.color.Reg
Prob.PAdom = PAdom;
ProbMer = mergeDivCon(Prob, Opt, 0);


% format output
PAmer = ProbMer.PA;
colorMer = ProbMer.color;


% some final outputs / statistics
% --------------------------------------------------------------------

nr_org = length(PA);
nr_mer = length(PAmer);

% if Opt.verbose>0,
%     disp(' ')
%     disp('********************************************************')
%     fprintf('total computation time: %0.2f s\n', cputime-t_total);
%     fprintf('# regions of original PWA model:          %i\n', nr_org);
%     fprintf('# regions after merging:                  %i\n', nr_mer);
%     fprintf('--> reduction by %i regions or %0.2f percent\n', nr_org-nr_mer, (nr_org-nr_mer)/nr_org*100);
%     disp(' ')
% end
% 
% if Opt.verbose>0 & ~isempty(Opt.logfile)
%     fprintf(Opt.fid, '\ntotal computation time: %0.2f s\n', cputime-t_total);
%     fprintf(Opt.fid, '# regions of original PWA model:          %i\n', nr_org);
%     fprintf(Opt.fid, '# regions after merging:                  %i\n', nr_mer);
%     fprintf(Opt.fid, '--> reduction by %i regions or %0.2f percent\n', nr_org-nr_mer, (nr_org-nr_mer)/nr_org*100);
%     fclose(Opt.fid);
% end;

return








function ProbMer = mergeDivCon(Prob, Opt, level)
% Sub.PA: polyhedral partition
% Sub.color: color information
   
%     % map colors from Sub.color (e.g. [1 3 7 8]) to successive integer
%     % sequence (e.g. [1 2 3 4]), i.e. 'flatten' the color information
%     colorMap_Flat2Full = unique(Prob.color);
%     colorMap_Full2Flat = [];
%     for i = 1:length(colorMap_Flat2Full)
%         colorMap_Full2Flat(colorMap_Flat2Full(i)) = i;
%     end;
           
    % check colors:
    % each polyhedron must have a color which is integer
    if length(Prob.PA) ~= length(Prob.color.Reg),
        error('the number of polyhedra and colors must be the same')
    end;
    if any(round(Prob.color.Reg)~=Prob.color.Reg), 
        error('the colors must be integers')
    end;
        
    % is there more than one color? 
    % If not, there is nothing to be merged...
    if length(unique(Prob.color.Reg)) == 1
        if Opt.verbose>0,
            disp('only one color: nothing to merge - giving back the domain')
        end
        
        if ~isempty(Opt.logfile)
            for i = 1:level, fprintf(Opt.fid, '   '); end;
            fprintf(Opt.fid, '%i --> 1 (t=0.0s) only 1 color\n', length(Prob.PA));
        end;
        
        ProbMer.PA = Prob.PAdom;
        
        c = Prob.color.Reg(1);
        ProbMer.color.Reg = c * ones(1,length(Prob.PAdom));
        colorMer.Table{c} = [1:length(Prob.PAdom)];
        return;
    end;
    
    % turn polyhedral array in halfspace representation
%     Prob.P.Hi = []; Prob.P.Ki = [];
%     for i = 1:length(Prob.PA)
%         [H, K] = double(Prob.PA(i));
%         Prob.P.Hi{i} = H;
%         Prob.P.Ki{i} = K;
%     end;
    %[Prob.P.Hi, Prob.P.Ki] = mpt2cell(Prob.PA);
    [Prob.P.Hi, Prob.P.Ki] = double(Prob.PA);
   
    % start to measure the time
    t_total = cputime;
    
    
    
    % derive hyperplane arrangement with markings
    % if necessary, split problem into subproblems
    % --------------------------------------------------------------------
    
    % derive hyperplane arrangement with markings
    OptexHyper = Opt;
    [M, color, GL, Dom] = mpt_exHyperAdv(Prob.P.Hi, Prob.P.Ki, Prob.PA, Prob.color.Reg, Prob.PAdom, OptexHyper);

    if isempty(M)
    
        % too many hyperplanes - computation stopped.
        % divide problem into subproblems
        
        [Psub1, Psub2] = deriveSubproblems(Prob, GL, Opt);
        
        % divide problem into 2 subproblems
        Sub1.PA = []; Sub2.PA = [];
        ind1 = []; ind2 = [];
        for i = 1:length(Prob.PA)
            %if i/500 == round(i/500), fprintf('    poly %i/%i\n', i, length(Prob.PA)); end;
            P1 = Prob.PA(i) & Psub1;
            P2 = Prob.PA(i) & Psub2;
            if isfulldim(P1) 
                ind1(end+1) = i;
                Sub1.PA = horzcat(Sub1.PA, P1);
            end;
            if isfulldim(P2) 
                ind2(end+1) = i;
                Sub2.PA = horzcat(Sub2.PA, P2);
            end;
        end;
        
        % colors for subproblem
        Sub1.color.Reg = Prob.color.Reg(ind1);
        Sub2.color.Reg = Prob.color.Reg(ind2);
        
        % update domains for subproblems
        Sub1.PAdom = Prob.PAdom & Psub1;
        Sub2.PAdom = Prob.PAdom & Psub2;    
        
        if Opt.plot
            figure; clf; hold on;
            plotPaCol(Sub1.PA, Sub1.color.Reg);
            plotPaCol(Sub2.PA, Sub2.color.Reg);
        end;
        
        % branch
        if Opt.verbose>1,
            fprintf('*** dividing problem into 2 subproblems: %i --> %i and %i ***\n', length(Prob.PA), length(Sub1.PA), length(Sub2.PA));
        end
        if ~isempty(Opt.logfile)
            for i = 1:level, fprintf(Opt.fid, '   '); end;
            fprintf(Opt.fid, [num2str(level) ': ']);
            fprintf(Opt.fid, '%i = %i + %i\n', length(Prob.PA), length(Sub1.PA), length(Sub2.PA));
        end;

        if Opt.verbose>1,
            disp(' ')
            disp('************************* subproblem 1 *************************')
        end
        t_start = cputime;
        SubMer1 = mergeDivCon(Sub1, Opt, level+1);
        if ~isempty(Opt.logfile)
            for i = 1:level+1, fprintf(Opt.fid, '   '); end;
            fprintf(Opt.fid, 'sub1: ');
            fprintf(Opt.fid, '%i --> %i (t=%0.1fs)\n', length(Sub1.PA), length(SubMer1.PA), cputime-t_start);
        end;
        
        if Opt.verbose>1,
            disp(' ')
            disp('************************* subproblem 2 *************************')
        end
        t_start = cputime;
        SubMer2 = mergeDivCon(Sub2, Opt, level+1);
        if ~isempty(Opt.logfile)
            for i = 1:level+1, fprintf(Opt.fid, '   '); end;
            fprintf(Opt.fid, 'sub2: ');
            fprintf(Opt.fid, '%i --> %i (t=%0.1fs)\n', length(Sub2.PA), length(SubMer2.PA), cputime-t_start);
        end;
        
        if Opt.verbose>1,
            disp(' ')
            disp('************************* subproblem 1 and 2 *************************')
        end
        
        % combine merged subproblems 
        Com = [];
        
        Com.PA = horzcat(SubMer1.PA, SubMer2.PA); 

        Com.color.Reg = [SubMer1.color.Reg, SubMer2.color.Reg];
        Com.color.Table = [];
%         Col = 1:max(length(SubMer1.color.Table), length(SubMer2.color.Table));
%         for col = Col
%             ind1 = SubMer1.color.Table{col};
%             ind2 = SubMer2.color.Table{col} + length(SubMer1.PA);
%             Com.color.Table{col} = [ind1, ind2];
%         end;
        
        Com.PAdom = Prob.PAdom;
        
        % the two combined subproblems are now our new problem
        Prob = Com; clear Com
        
        % turn polyhedral array in halfspace representation
        Prob.P.Hi = []; Prob.P.Ki = [];
        for i = 1:length(Prob.PA)
            [H, K] = double(Prob.PA(i));
            Prob.P.Hi{i} = H;
            Prob.P.Ki{i} = K;
        end;        

        % is there more than one color? 
        % If not, there is nothing to be merged...
        if length(unique(Prob.color.Reg)) == 1
            if Opt.verbose>0,
                disp('only one color: nothing to merge - giving back the domain')
            end
            
            if ~isempty(Opt.logfile)
                for i = 1:level, fprintf(Opt.fid, '   '); end;
                fprintf(Opt.fid, '%i --> 1 (t=0.0s) only 1 color\n', length(Prob.PA));
            end;
            
            ProbMer.PA = Prob.PAdom;
            
            c = Prob.color.Reg(1);
            ProbMer.color.Reg = c * ones(1,length(Prob.PAdom));
            colorMer.Table{c} = [1:length(Prob.PAdom)];

            return;
        end;
        
         % derive hyperplane arrangement with markings
        OptexHyper = Opt;
        OptexHyper.maxHA = Opt.maxHA_sub1and2;
        [M, color, GL, Dom] = mpt_exHyperAdv(Prob.P.Hi, Prob.P.Ki, Prob.PA, Prob.color.Reg, Prob.PAdom, OptexHyper);
        
    end;
    
    if isempty(M)
        if ~isempty(Opt.logfile)
            for i = 1:level, fprintf(Opt.fid, '   '); end;
            fprintf(Opt.fid, 'HA too large - exiting\n');
        end;
        
        %ProbMer.PA = Prob.PA;
        %ProbMer.color = Prob.color;
        ProbMer = Prob;
        return
    end;
        
    % now we have:
    % M, color, GL, Dom
    
    % plot
    if Opt.plot
        plotPolyMark(M, color, GL, Dom)
    end;
    
    if Opt.verbose == 2
        disp(' ')
        disp('  table specifying the regions with the same PWA dynamics (or color):')
        for i=1:length(color.Table), fprintf('  color %i:',i), fprintf(' %i', color.Table{i}), fprintf('\n'), end;
    end;
    
    
    
    % merge
    % --------------------------------------------------------------------
    
    if Opt.verbose>1, 
        disp(' ');
        disp(' **** merge **** ');
        disp(' ')    
    end;
    
    Col = [];
    for i = 1:length(color.Table)
        if ~isempty(color.Table{i}), Col(end+1) = i; end;
    end;
    
    % initialize structure PAmer and colorMer
    PAmer = [];
    colorMer = [];
    colorMer.Reg = [];
    for c=Col, colorMer.Table{c} = []; end;
    
    % merge the polyhedra of each color individually
    for c = Col
        
        % merge polyhedra with indices color.Table{c} using the set of 
        % markings M referring to the hyperplanes GL and restrict plotting 
        % to the domain P
        M_mer = mergePol(color.Table{c}, M, GL, Dom, Opt);
        
        % write the merged polyhedra together with the PWA dynamic into the new
        % structure Pmer{f}
        for i=1:length(M_mer)
            
            % get the corresponding polyhedron 
            % (defined by the hyperplane arrangement)
            Hi = (-1)*diag(M_mer{i}) * GL.Hi;
            Ki = (-1)*diag(M_mer{i}) * GL.Ki;
            
            % remove the zero rows (rows where M_mer == 0)
            Hi( find(M_mer{i}==0), : ) = [];
            Ki( find(M_mer{i}==0), : ) = [];
            
            % get minimal polyhedral description
            Pmer = polytope([Hi; Dom.A], [Ki; Dom.B]) & Prob.PAdom;
            
            % store the corresponding polyhedron 
            if isfulldim(Pmer), 
                % add the polyhedron
                PAmer = horzcat(PAmer, Pmer);
                
                % add the color
                colorMer.Reg(end+1) = c;
                colorMer.Table{c}(end+1) = length(PAmer);
            else
                if Opt.verbose>1,
                    warning('white polyhedron is (almost) empty');
                end
            end;
            
        end;
        
    end
    
    if length(PAmer) < length(Prob.PA)
        % give back merged problem
        ProbMer.PA = PAmer;
        ProbMer.color = colorMer;
    else
        % give back original problem
        ProbMer = Prob;
        if Opt.verbose>0,
            fprintf('merged problem larger than original problem --> give original problem back')
        end
    end;
    
    
    % plot
    if Opt.plot
        figure; plotPaCol(ProbMer.PA, ProbMer.color.Reg)
    end;

    
    
    % some final outputs / statistics
    % --------------------------------------------------------------------
    
    nr_org = length(Prob.PA);
    nr_mer = length(ProbMer.PA);

    if Opt.verbose>0,
        disp(' ')
        fprintf('total computation time: %0.2f s\n', cputime-t_total);
        fprintf('# regions of original PWA model:          %i\n', nr_org);
        fprintf('# regions after extension of hyperplanes: %i\n', size(M,2));
        fprintf('# regions after merging:                  %i\n', nr_mer);
        fprintf('--> reduction by %i regions or %0.2f percent\n', nr_org-nr_mer, (nr_org-nr_mer)/nr_org*100);
        disp(' ')
    end
    
    delete(Opt.infile);
    delete(Opt.outfile);
    if ~isempty(Opt.logfile),
        delete(Opt.logfile);
    end
return    





function [p, r] = prefix( N )
	% computes the prefix and the remainder of an hyperplane arrangement N, 
	% where the prefix is the set of rows that have all the same elements
	% and the remainder are the other rows
	p = find( all( N - N(:,1)*ones(1,size(N,2)) == 0, 2));
	r = 1:size(N,1); r(p) = [];
	p=p(:); r=r(:);
return



  
function M_mer = mergePol(I, M, GL, Dom, Opt);
    % merge polyhedra 
    % with indices I using the set of markings M referring to the hyperplanes 
    % GL and restrict plotting to the domain Dom

    % The regions which we want to merge share possibly some hyperplanes 
    % (or markings) which we call prefix. Within this prefix, the regions
    % can be marked using only the remainder. 
    % The regions within the prefix which we want to merge, are called
    % white polyhedra, whereas all the other regions within the prefix are
    % called black polyhedra.
    
    % get the prefix index 
    [pre_ind, rem_ind] = prefix( round(M(:, I)) );
    
    % get the prefix marking
    pre_m = round(M(pre_ind, I(1)));
    
    % the white and black polyhedra have markings within the prefix; the 
    % prefix is removed to reduce complexity (=we keep only the remainder)
    m_white = M(rem_ind, I);
    
    % black candidate polyhedra are all the other regions
    black_cand = 1:size(M, 2);
    black_cand(I) = [];
    
    % the black polyhedra are the black candidate polyhedra having the 
    % same prefix 'pre_m'
    black = find( all( pre_m*ones(1,length(black_cand)) == round(M(pre_ind, black_cand)), 1) );
    m_black = M(rem_ind, black_cand(black));
    
    % --> now we have:
    %   * the prefix (both the index 'pre_ind' and the marking 'pre_m'),
    %   * the markings 'm_white' of the white polyhedra which we want to 
    %     merge (\in {-1, -0.5, 0, 0.5, 1}), and
    %   * the markings 'm_black' of the black polyedra within the prefix
    %     (\in {-1, -0.5, 0, 0.5, 1})

    if Opt.verbose>1,
        fprintf('  %i white, %i black --> ', size(m_white,2), size(m_black,2))
    end
    t = cputime;
    
    % if there are no black polyhedra ...
    if isempty(m_black)
        % ... the remainder has only don't cares or 0
        m_mer = rem_ind*0;
        nodes = 0;
        
    % if we have only one white region ...
    elseif size(m_white, 2) == 1
        % ... we do not need to merge anything
        m_mer = m_white;
        nodes = 0;
        
    else
        if Opt.safeMode
            % check data:
            % do any markings occur both in m_white and m_black? That would
            % be very bad...
            for i=1:size(m_black,2)
                if find(all( m_white == m_black(:,i)*ones(1,size(m_white,2)) ))
                    error('marking is both white and black');
                end;
            end;
        end;
        
        %Opt.algo = 0;
        if Opt.algo > 0
            % optimal merging based on branch and bound
            % resulting in non-overlapping polyhedra
            [m_mer, nodes] = mpt_merge5(rem_ind*0, m_white, m_black, 0, 0, inf, Opt.algo);
        else
            % optimal merging based on boolean minimization
            % resulting in overlapping polyhedra
            
            m_white = round(m_white);
            m_black = round(m_black);
            
            % check markings
            if find(m_white==0), error('unexpected 0 in M_on'); end;
            if find(m_black==0), error('unexpected 0 in M_off'); end;
            
            % for large problems: pre-merge 
            if size(m_black,2) > inf %1000
                % perform a quick distance 1 merge on the black markings
                writeEspressoOnOff(Opt.infile, m_black, m_white);
                tt=evalc(['!' Opt.espressoCmd ' -Dd1merge ' Opt.infile ' > ' Opt.outfile]);
                m_black = readEspressoOne(Opt.outfile);
            end;
            
            % for large problems: pre-merge 
            if size(m_white,2) > inf %1000
                % perform a quick distance 1 merge on the white markings
                writeEspressoOnOff(Opt.infile, m_white, m_black);
                tt=evalc(['!' Opt.espressoCmd ' -Dd1merge ' Opt.infile ' > ' Opt.outfile]);
                m_white = readEspressoOne(Opt.outfile);
            end;
            
            % merge enforcing minimality of polyhedra (by using the option -Dexact)
            writeEspressoOnOff(Opt.infile, m_white, m_black);
            %eval(['!' Opt.espressoCmd ' -Dexact ' Opt.infile ' > ' Opt.outfile]);
            tt=evalc(['!' Opt.espressoCmd ' ' Opt.infile ' > ' Opt.outfile]);
            m_mer = readEspressoOne(Opt.outfile);
            
            nodes = 0;
        end;
    end;
    
    if Opt.verbose>1,
        fprintf('%i white (%0.5fs, %i nodes)\n', size(m_mer,2), cputime-t, nodes)
    end

    % add the prefix to the merged polyhedra
    M_mer = {};
    for i=1:size(m_mer,2)
        M_mer{i}(rem_ind) = m_mer(:,i);
        M_mer{i}(pre_ind) = pre_m;
    end;
return



function [Psub1, Psub2] = deriveSubproblems(Prob, GL, Opt)

% get hyperplane H*[x;u] = K to divide the problem into 2 subproblems
if Opt.divideSub == 0

    % approach 0: 
    % use a new hyperplane parallel to axis of the
    % coordinate with the longest edge of the bounding box
    % this can be done very quickly
    
    % bounding box of full problem
    BB.box = bounding_box(Prob.PAdom); 
    BB.vert = extreme(BB.box);
    BB.min = min(BB.vert);
    BB.max = max(BB.vert);
    BB.len = BB.max - BB.min;
    
    % get hyperplane parallel to axis of the
    % coordinate with the longest edge of the bounding box
    [dummy, i] = max(BB.len);
    dim = dimension(Prob.PA(1));
    H = zeros(1,dim); H(i) = 1;
    K = BB.min(i) + 0.5*BB.len(i);
    
else
    
    % approach 1:
    % go through all hyperplanes, choose the one that minimizes the cost
    % cost = abs(# GL in Sub1 - # GL in Sub2) + # GL in Sub1 + # GL in Sub2
    % this is computationally expensive (in particular for large problems)
    
    stat = [];
    
    % set of hyperplanes we can use to divide the problem into 2 subproblems
    hSet = 1:length(GL.Ki);
    hSet = hSet(randperm(length(hSet)));
    
    % restrict number of test hyperplanes according to size of HA
    if length(hSet) > 500, hSet = hSet(1:10); 
    elseif length(hSet) > 20, hSet = hSet(1:20); 
        %elseif length(hSet) > 40, hSet = hSet(1:round(length(hSet)/2));
    end;
    % set bound on maximal cost maxCost    
    maxCost = 1.2*length(GL.Ki); 

    % build statistics 'stat'
    bestCost = inf;
    while (length(hSet) > 0) & (bestCost > maxCost)
        % choose one hyperplane randomly and remove it from the set
        hSet = hSet(randperm(length(hSet)));
        h = hSet(1);
        hSet(find(hSet==h)) = [];
        
        fprintf('  testing hyperplane #%i: ', h);
        
        % get hyperplane
        H = GL.Hi(h,:); K = GL.Ki(h);
        
        % get halfspaces
        Psub1 = polytope(H, K);
        Psub2 = polytope(-H, -K);
        
        % divide problem into 2 subproblems
        Sub1.PA = []; Sub2.PA = [];
        ind1 = []; ind2 = [];
        for i = 1:length(Prob.PA)
            [PolyA, PolyB] = double(Prob.PA(i));
            m = intersectHP1(PolyA, PolyB, H, K, Opt.lpsolver, 1e-8);
            if m == -1
                % polyhedron is in subproblem 1
                ind1(end+1) = i;
                Sub1.PA = horzcat(Sub1.PA, Prob.PA(i));
            elseif m == +1
                % polyhedron is in subproblem 2
                ind2(end+1) = i;
                Sub2.PA = horzcat(Sub2.PA, Prob.PA(i));
            else
                % polyhedron lies in both subproblems
                P1 = Prob.PA(i) & Psub1;
                P2 = Prob.PA(i) & Psub2;
                if isfulldim(P1) 
                    ind1(end+1) = i;
                    Sub1.PA = horzcat(Sub1.PA, P1);
                end;
                if isfulldim(P2) 
                    ind2(end+1) = i;
                    Sub2.PA = horzcat(Sub2.PA, P2);
                end;
            end;
        end;
        
        if isempty(Sub1.PA) | isempty(Sub2.PA)
            stat(end+1,:) = [h, 0, 0, inf];
            
            fprintf('hyperplane is not suitable (as it is bounding)\n');
        else
            % halfspace description of subproblems
            %[Sub1.P.Hi, Sub1.P.Ki] = mpt2cell(Sub1.PA);
            [Sub1.P.Hi, Sub1.P.Ki] = double(Sub1.PA);
            [Sub2.P.Hi, Sub2.P.Ki] = double(Sub2.PA);
            
            % colors for subproblem
            Sub1.color.Reg = Prob.color.Reg(ind1);
            Sub2.color.Reg = Prob.color.Reg(ind2);
            
            % update domains for subproblems
            Sub1.PAdom = Prob.PAdom & Psub1;
            Sub2.PAdom = Prob.PAdom & Psub2;  
            
            % get hyperplane arrangement (without the markings)
            OptGetGL = Opt;
            OptGetGL.maxHA = 0;
            OptGetGL.verbose = 0;
            
            [Sub1.M, Sub1.color.Reg, Sub1.GL] = mpt_exHyperAdv(Sub1.P.Hi, Sub1.P.Ki, Sub1.PA, Sub1.color.Reg, Sub1.PAdom, OptGetGL);

            % compute only size of HA in subproblem 2 if the subproblem 1
            % is better than the current optimum
            if 2*length(Sub1.GL.Ki) < bestCost
                [Sub2.M, Sub2.color.Reg, Sub2.GL] = mpt_exHyperAdv(Sub2.P.Hi, Sub2.P.Ki, Sub2.PA, Sub2.color.Reg, Sub2.PAdom, OptGetGL);
                
                % cost is twice the inf-norm of the size of the HAs
                cost = 2*norm([length(Sub1.GL.Ki); length(Sub2.GL.Ki)], inf);
            else
                % update statistics:
                Sub2.GL.Ki = [];
                cost = NaN;
            end;
            
            % update statistics:
            % # of dividing hyperplane, # GL in Sub1, # GL in Sub2, cost
            stat(end+1,:) = [h, length(Sub1.GL.Ki), length(Sub2.GL.Ki), cost];
                
            fprintf('#hyperpl: %i --> %i %i (cost=%i)\n', length(GL.Ki), length(Sub1.GL.Ki), length(Sub2.GL.Ki), cost);  
            
            % update best (i.e. the so far smallest) cost
            if ~isempty(stat)
                minCost = min(stat(:,4));
            else
                minCost = inf;
            end;
            bestCost = min(bestCost, minCost);
        end;
    end;
    
    % determine hyperplane to branch on by evaluating stat
    minCost = min(stat(:,4));
    h = stat(find(stat(:,4)==minCost),1);
    h = h(1);
    H = GL.Hi(h,:); K = GL.Ki(h); 
end;

% get halfspaces
Psub1 = polytope(H, K);
Psub2 = polytope(-H, -K);

% divide problem into 2 subproblems
Sub1.PA = []; Sub2.PA = [];
ind1 = []; ind2 = [];
for i = 1:length(Prob.PA)
    P1 = Prob.PA(i) & Psub1;
    P2 = Prob.PA(i) & Psub2;
    if isfulldim(P1) 
        ind1(end+1) = i;
        Sub1.PA = horzcat(Sub1.PA, P1);
    end;
    if isfulldim(P2) 
        ind2(end+1) = i;
        Sub2.PA = horzcat(Sub2.PA, P2);
    end;
end;

% colors for subproblem
Sub1.color.Reg = Prob.color.Reg(ind1);
Sub2.color.Reg = Prob.color.Reg(ind2);

% update domains for subproblems
Sub1.PAdom = Prob.PAdom & Psub1;
Sub2.PAdom = Prob.PAdom & Psub2;  

return





function [Psub1, Psub2] = deriveSubproblemsOLD(Prob, GL, Opt)

% get hyperplane H*[x;u] = K to divide the problem into 2 subproblems
if Opt.divideSub == 0

    % approach 0: 
    % use a new hyperplane parallel to axis of the
    % coordinate with the longest edge of the bounding box
    % this can be done very quickly
    
    % bounding box of full problem
    BB.box = bounding_box(Prob.PAdom); 
    BB.vert = extreme(BB.box);
    BB.min = min(BB.vert);
    BB.max = max(BB.vert);
    BB.len = BB.max - BB.min;
    
    % get hyperplane parallel to axis of the
    % coordinate with the longest edge of the bounding box
    [dummy, i] = max(BB.len);
    dim = dimension(Prob.PA(1));
    H = zeros(1,dim); H(i) = 1;
    K = BB.min(i) + 0.5*BB.len(i);
    
else
    
    % approach 1:
    % go through all hyperplanes, choose the one that minimizes the cost
    % cost = abs(# GL in Sub1 - # GL in Sub2) + # GL in Sub1 + # GL in Sub2
    % this is computationally expensive (in particular for large problems)
    
    stat = [];
    
    % set of hyperplanes we can use to divide the problem into 2 subproblems
    hSet = 1:length(GL.Ki);
    hSet = hSet(randperm(length(hSet)));
    
    % restrict number of test hyperplanes according to size of HA
    if length(hSet) > 500, hSet = hSet(1:10); 
    elseif length(hSet) > 20, hSet = hSet(1:20); 
        %elseif length(hSet) > 40, hSet = hSet(1:round(length(hSet)/2));
    end;
    % set bound on maximal cost maxCost    
    maxCost = 1.2*length(GL.Ki); 

    % build statistics 'stat'
    bestCost = inf;
    while (length(hSet) > 0) & (bestCost > maxCost)
        % choose one hyperplane randomly and remove it from the set
        hSet = hSet(randperm(length(hSet)));
        h = hSet(1);
        hSet(find(hSet==h)) = [];
        
        fprintf('  testing hyperplane #%i: ', h);
        
        % get hyperplane
        H = GL.Hi(h,:); K = GL.Ki(h);
        
        % get halfspaces
        Psub1 = polytope(H, K);
        Psub2 = polytope(-H, -K);
        
        % divide problem into 2 subproblems
        Sub1.PA = []; Sub2.PA = [];
        ind1 = []; ind2 = [];
        for i = 1:length(Prob.PA)
            P1 = Prob.PA(i) & Psub1;
            P2 = Prob.PA(i) & Psub2;
            if isfulldim(P1) 
                ind1(end+1) = i;
                Sub1.PA = horzcat(Sub1.PA, P1);
            end;
            if isfulldim(P2) 
                ind2(end+1) = i;
                Sub2.PA = horzcat(Sub2.PA, P2);
            end;
        end;
        
        if isempty(Sub1.PA) | isempty(Sub2.PA)
            stat(end+1,:) = [h, 0, 0, inf];
            
            fprintf('hyperplane is not suitable (as it is bounding)\n');
        else
            % halfspace description of subproblems
            [Sub1.P.Hi, Sub1.P.Ki] = double(Sub1.PA);
            [Sub2.P.Hi, Sub2.P.Ki] = double(Sub2.PA);
            
            % colors for subproblem
            Sub1.color.Reg = Prob.color.Reg(ind1);
            Sub2.color.Reg = Prob.color.Reg(ind2);
            
            % update domains for subproblems
            Sub1.PAdom = Prob.PAdom & Psub1;
            Sub2.PAdom = Prob.PAdom & Psub2;  
            
            % get hyperplane arrangement (without the markings)
            OptGetGL = Opt;
            OptGetGL.maxHA = 0;
            OptGetGL.verbose = 0;
            [Sub1.M, Sub1.color.Reg, Sub1.GL] = mpt_exHyperAdv(Sub1.P.Hi, Sub1.P.Ki, Sub1.PA, Sub1.color.Reg, Sub1.PAdom, OptGetGL);
            [Sub2.M, Sub2.color.Reg, Sub2.GL] = mpt_exHyperAdv(Sub2.P.Hi, Sub2.P.Ki, Sub2.PA, Sub2.color.Reg, Sub2.PAdom, OptGetGL);
            
            % update statistics:
            % # of dividing hyperplane, # GL in Sub1, # GL in Sub2, cost
            % where cost = abs(# GL in Sub1 - # GL in Sub2) + # GL in Sub1 + # GL in Sub2
            % this cost is basically the infinity norm times 2
            cost = abs(length(Sub1.GL.Ki)-length(Sub2.GL.Ki)) + length(Sub1.GL.Ki) + length(Sub2.GL.Ki);
            stat(end+1,:) = [h, length(Sub1.GL.Ki), length(Sub2.GL.Ki), cost];
            
            fprintf('#hyperpl: %i --> %i %i (cost=%i)\n', length(GL.Ki), length(Sub1.GL.Ki), length(Sub2.GL.Ki), cost);
            
            % update best (i.e. the so far smallest) cost
            if ~isempty(stat)
                minCost = min(stat(:,4));
            else
                minCost = inf;
            end;
            bestCost = min(bestCost, minCost);
        end;
    end;
    
    % determine hyperplane to branch on by evaluating stat
    minCost = min(stat(:,4));
    h = stat(find(stat(:,4)==minCost),1);
    h = h(1);
    H = GL.Hi(h,:); K = GL.Ki(h); 
end;

% get halfspaces
Psub1 = polytope(H, K);
Psub2 = polytope(-H, -K);

% divide problem into 2 subproblems
Sub1.PA = []; Sub2.PA = [];
ind1 = []; ind2 = [];
for i = 1:length(Prob.PA)
    P1 = Prob.PA(i) & Psub1;
    P2 = Prob.PA(i) & Psub2;
    if isfulldim(P1) 
        ind1(end+1) = i;
        Sub1.PA = horzcat(Sub1.PA, P1);
    end;
    if isfulldim(P2) 
        ind2(end+1) = i;
        Sub2.PA = horzcat(Sub2.PA, P2);
    end;
end;

% colors for subproblem
Sub1.color.Reg = Prob.color.Reg(ind1);
Sub2.color.Reg = Prob.color.Reg(ind2);

% update domains for subproblems
Sub1.PAdom = Prob.PAdom & Psub1;
Sub2.PAdom = Prob.PAdom & Psub2;  

return







    % build statistics 'stat'
    bestCost = inf;
    while (length(hSet) > 0) & (bestCost > maxCost)
        % choose one hyperplane randomly and remove it from the set
        hSet = hSet(randperm(length(hSet)));
        h = hSet(1);
        hSet(find(hSet==h)) = [];
        
        fprintf('  testing hyperplane #%i: ', h);
        
        % indices of all hyperplanes except of the h-th
        I = 1:length(GL.Ki);
        I(h) = [];        
        
        % initialize position P
        % -1: i-th hyperplane on negative side of h-th hyperplane
        %  0: hyperlanes intersect on Prob.PAdom
        % +1: on positive side
        P = NaN*ones(1,length(GL.Ki)); P(h) = 0;

        % find a point on the h-th hyperplane inside the domain
        f = zeros(1,size(GL.Hi,2));
        [A, B]=double(Prob.PAdom);
        B = B*0.999;        % slightly reduce the domain by 0.1%
        Aeq = GL.Hi(h,:);
        Beq = GL.Ki(h);
        [xh,fval,lambda,exitflag,how] = mpt_solveLP(f,A,B,Aeq,Beq);
        
        if ~strcmp(how, 'ok'),
            fprintf('hyperplane is not suitable (as it is bounding or not intersecting domain)\n');
            stat(end+1,:) = [h, 0, 0, inf];
        else
            % check all hyperplanes in I
            for i = I
                f = zeros(1,size(GL.Hi,2));
                [A, B]=double(Prob.PAdom);
                B = B*0.99;
                Aeq = [GL.Hi(i,:); GL.Hi(h,:)];
                Beq = [GL.Ki(i); GL.Ki(h)];
                [xopt,fval,lambda,exitflag,how] = mpt_solveLP(f,A,B,Aeq,Beq);
                
                if strcmp(how, 'ok'),
                    % both hyperplanes intersect on the domain
                    P(i) = 0;
                else
                    % they do not intersect
                    if GL.Hi(i,:)*xh < GL.Ki(i)
                        P(i) = -1;
                    else
                        P(i) = +1;
                    end;
                end;
            end;
            
            numGL(1) = length(find(P==-1)) + length(find(P==0))-1;
            numGL(2) = length(find(P==+1)) + length(find(P==0))-1;
            
            % update statistics:
            % # of dividing hyperplane, # GL in Sub1, # GL in Sub2, cost
            % where cost = abs(# GL in Sub1 - # GL in Sub2) + # GL in Sub1 + # GL in Sub2
            % this cost is basically the infinity norm times 2
            stat(end+1,:) = [h, numGL(1), numGL(2), norm(numGL,inf)];
            
            fprintf('#hyperpl: %i --> %i %i (cost=%i)\n', length(GL.Ki), numGL(1), numGL(2), norm(numGL,inf));
            
            % update best (i.e. the so far smallest) cost
            minCost = min(stat(:,4));
            bestCost = min(bestCost, minCost);
        end;
    end;   
