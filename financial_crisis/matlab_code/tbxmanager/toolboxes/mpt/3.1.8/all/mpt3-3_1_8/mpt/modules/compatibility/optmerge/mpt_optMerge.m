%===============================================================================
%
% function [PAmer, colorMer] = mpt_optMerge(PA, Opt)
% 
% Title:        mpt_optMerge
%                                                             
% Project:      Optimal merging of a set of polyhedra 
%
% Input:        PA: polyhedral array
%               Opt: Options with the following fields
%                 .color: color of polyhedra
%                         this is an integer 1,2,... 
%                         the algorithm will merge polyhedra with the same color 
%                         if color is not specified, it is assumed to be 1
%                         for all polyhedra.
%                 .PAdom: the set of polyhedra (the so called domain) the problem 
%                         is defined in. If not given, the algorithm assumes the 
%                         domain is given by the hull (or envelope if the hull 
%                         computation fails) of PA.
%                 .PAcompl:  polyhedral array that is within the domain and not 
%                            in PA - it is the so called complement. If not specified, 
%                            the algorithm assumes PAcompl is empty. Then the 
%                            complement is filled up optimally by the merging algorithm.
%                            Refer also to the second remark below.
%                 .verbose: 0: silent (=default)
%                           1: verbose only important information
%                           2: verbose everything
%                 .algo: 0:        overlapping merging (boolean minimization)
%                        [1, inf]: optimal merging, with the given upper bound on nodes
%                                  (default is inf)
%
% Output:       PAmer: merged polyhedra array
%               colorMer: colors with the following fields
%                 .Reg:   color of each region
%                 .Table: regions with the same color 
%
% Authors:      Tobias Geyer <geyer@control.ee.ethz.ch>, Fabio D. Torrisi

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
%               2002.11.xx  1       initial version 
%               2002.11.xx  1.1     extension of hyperplanes of polyhedra 
%               2003.03.01  1.1.1   adapted to hys2pwa 1.1.2
%               2003.03.11  1.2     hyperplane extension
%               2003.03.14  1.2.1   generalize to PWA feedback laws
%               2003.05.22  1.2.2   bug fixed for m_black
%               2003.06.22  1.2.3   works now also for non-convex domains
%               2003.07.02  1.2.4   plot intermediate state of merging 
%               2003.07.03  1.3     restructured
%               2003.08.04  1.3.1   minimal representation of polyhedra
%               2003.08.23  1.4     associate colours to the polyhedra
%               2004.02.10  1.5     optional overlapping merging (Espresso),
%                                   upgraded to cddmex
%               2004.04.02  2.0     rewritten for MPT
%               2004.04.23  2.1     numerically more stable,
%                                   accepts now also overlapping polyhedra,
%                                   uses hyparr to compute markings
%               2004.06.15  2.2     requires domain
%               2005.07.20  2.3     bug fix and improved help
%
% Requires:     mpt_exHyperAdv,
%               hyparr,
%               merge5,
%               plotPolyMark, plotPaCol,
%               writeEspressoOnOff, readEspressoOne,
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


function [PAmer, colorMer] = mpt_optMerge(PA, Opt)

% options
% --------------------------------------------------------------------

if nargin < 2,
    Opt = [];
end;
if nargin < 1, 
    disp('mpt_optMerge ver 2.3   Tobias Geyer 2002-2005');
    return;
end;
if ~isfield(Opt, 'color'), 
    color.Reg = ones(1,length(PA));
else
    color.Reg = Opt.color(:)';
end
%color.Table = Opt.Table;
if ~isfield(Opt, 'algo'), 
    Opt.algo = 0;
end;
if ~isfield(Opt, 'verbose'), 
    Opt.verbose = 1; 
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

if ~isfield(Opt, 'verbose'),
    Opt.verbose=0;
end

Opt.plot = 0; 
Opt.safeMode = 0;   
Opt.lpsolver = 3;
if ismac
    Opt.espressoCmd = which('espresso.mac');
elseif ispc
    Opt.espressoCmd = which('espresso.exe');
else
    Opt.espressoCmd = which('espresso.linux');
end

if isempty(Opt.espressoCmd)
    error('mpt_optMerge: couldn''t find espresso binary in Matlab path!');
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
if min(colorSorted) < 1 | max(colorSorted) > size(colorSorted)
    error('colors need to be 1, 2, ...')
end;    

% if the user has not specified the domain,
% we quickly compute here the hull or the envelope
if isempty(PAdom)
	disp('Computing convex hull...')
	PAdom = PolyUnion(PA).convexHull();
    disp('Warning: you have not specified the set of polyhedra (the so called');
    disp('domain) the problem is defined in. We assume the domain is given by');
    disp('the hull (or envelope if the hull computation fails).');
    fprintf('\n');
end;

% parts, that are within the domain and not in PA, are the so called complement
% Often, the complement is empty. If not specified, we assume it is empty.
% Then the complement is filled up optimally by the merging algorithm.
if isempty(PAcompl)
    disp('Warning: you have not specified the set difference between the set of');
    disp('polyhedra to be merged and the domain. We thus assume it is empty.');
    disp('Note that a non-empty not specified complement is filled up by the');
    disp('merging algorithm optimally!');
elseif any(PAcompl.isFullDim())
    % add the complement to the polyhedra array PA and associate an
    % auxiliary color to the complement
    colorSorted = sort(unique(color.Reg));
    maxCol = colorSorted(end);
    PA = [PA; PAcompl];
    color.Reg = [color.Reg, ones(1,length(PAcompl))*(maxCol+1)];
end;

% turn polyhedral array in halfspace representation
P.Hi = []; P.Ki = [];
for i = 1:length(PA)
    P.Hi{i} = PA(i).A;
    P.Ki{i} = PA(i).b;
end;

% start to measure the time
t_total = cputime;



% derive hyperplane arrangement with markings
% --------------------------------------------------------------------

%[M, color, GL, Dom] = mpt_exHyper(P.Hi, P.Ki, color, Opt);
[M, color, GL, Dom] = mpt_exHyperAdv(P.Hi, P.Ki, PA, color, PAdom, Opt);
cputime-t_total;

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
    disp(' ************************* merge ************************* ');
    disp(' ')    
else
    %disp('  merge ')
end;

% If we have initially added polyhedra to complement the convex hull of
% PA, we will not merge them here and they will not show up in the final
% set of polyhedra. However, their markings are in M and will serve as
% counterexamples.
if ~isempty(PAcompl) && any(PAcompl.isFullDim())
    maxCol = length(color.Table)-1;
else
    maxCol = length(color.Table);
end;
Col = [];
for i = 1:maxCol
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
        Pmer = Polyhedron([Hi; Dom.A], [Ki; Dom.B]).intersect(PAdom);
        
        % store the corresponding polyhedron 
        if Pmer.isFullDim() 
            % add the polyhedron
            PAmer = [PAmer; Pmer];
            
            % add the color
            colorMer.Reg(end+1) = c;
            colorMer.Table{c}(end+1) = length(PAmer);
        else
            warning('white polyhedron is (almost) empty'); 
        end;
        
    end;
        
end

% plot
if Opt.plot
    figure; plotPaCol(PAmer, colorMer.Reg)
end;


% some final outputs / statistics
% --------------------------------------------------------------------

nr_org = length(PA) - length(PAcompl);
nr_mer = length(PAmer);

if Opt.verbose>0,
    disp(' ')
    fprintf('total computation time: %0.2f s\n', cputime-t_total);
    fprintf('# regions of original PWA model:          %i\n', nr_org);
    fprintf('# regions after extension of hyperplanes: %i\n', size(M,2));
    fprintf('# regions after merging:                  %i\n', nr_mer);
    fprintf('--> reduction by %i regions or %0.2f percent\n', nr_org-nr_mer, (nr_org-nr_mer)/nr_org*100);
    disp(' ')
end

return
        















% subfunctions
% --------------------------------------------------------------------


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
        
        if Opt.algo > 0
            % optimal merging based on branch and bound
            % resulting in non-overlapping polyhedra
            [m_mer, nodes] = mpt_merge5(rem_ind*0, m_white, m_black, 0, 0, inf, Opt.algo);
        else
            % suboptimal merging based on boolean minimization
            % resulting in overlapping polyhedra
            nStr = tempname;
            infname = [nStr '_in.txt'];
            outfname = [nStr '_out.txt'];
            writeEspressoOnOff(infname, round(m_white), round(m_black));
            if ispc,
                eval(['! "' Opt.espressoCmd '" ' infname ' > ' outfname]);
            else
                eval(['!' Opt.espressoCmd ' ' infname ' > ' outfname]);
            end
            m_mer = readEspressoOne(outfname);
            delete(infname);
            delete(outfname);
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
