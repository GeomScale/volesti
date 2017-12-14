function probStruct=mpt_verifyProbStruct(probStruct,Options)
%MPT_VERIFYPROBSTRUCT Verifies the probStruct structure
%
% probStruct=mpt_verifyProbStruct(probStruct,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Verifies if the structure defining the optimization problem has all fields
% filled out.
%
% Internal function.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% probStruct       - Problem structure in the probStruct format
% Options.verbose  - level of verbosity
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% probStruct       - Verified problem structure in the probStruct format
%
% see also MPT_VERIFYSYSSTRUCT, MPT_CONTROL
%

% Copyright is with the following author(s):
%
% (C) 2003-2006 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%               kvasnica@control.ee.ethz.ch

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

narginchk(1, 2);

global mptOptions;

if ~isstruct(mptOptions),
    mpt_error;
end

if nargin<2
    Options=[];
end
Options = mpt_defaultOptions(Options, ...
    'verbose', mptOptions.verbose, ...
    'probstructname', inputname(1), ...
    'guierrors', 0, ...
    'useyalmip', 0 );

if isempty(Options.probstructname),
    Options.probstructname = 'probStruct';
end

psn = Options.probstructname;


if ~isfield(probStruct,'N'),
    if Options.guierrors,
        error('Prediction horizon must be specified!');
    else
        error(['Prediction horizon "' psn '.N" must be specified!']);
    end
end
if probStruct.N<1,
    if Options.guierrors,
        error('Prediction horizon must be greater zero!');
    else
        error(['Prediction horizon "' psn '.N" must be greater zero!']);
    end
end
if round(probStruct.N) ~= probStruct.N,
    if Options.guierrors,
        error('Prediction horizon must be an integer!');
    else
        error(['"' psn '.N" must be an integer !']);
    end
end

if ~isfield(probStruct,'Q')
    if Options.guierrors,
        error('Penalty on states must be specified!');
    else
        error(['"' psn '.Q" must be specified!']);
    end
end
if ~isfield(probStruct,'R')
    if Options.guierrors,
        error('Penalty on inputs must be specified!');
    else
        error(['"' psn '.R" must be specified!']);
    end
end
if ~isfield(probStruct,'norm')
    if Options.guierrors,
        error('Type of cost function must be specified!');
    else
        error(['"' psn '.norm" must be specified! Allowed values are 1, 2, Inf.']);
    end
end
if ~any(probStruct.norm==[1 2 Inf]),
    if Options.guierrors,
        error('Unknown type of cost function!');
    else
        error(['Unknown norm in "' psn '.norm"!']);
    end
end
if ~isfield(probStruct,'subopt_lev')
    if Options.verbose>0,
        disp(['Value of "' psn '.subopt_lev" not specified, assuming 0...']);
    end
    probStruct.subopt_lev=0;
end
if ~any(probStruct.subopt_lev==[0 1 2]),
    if Options.guierrors,
        error('Unknown type of solution!');
    else
        error(['"' psn '.subopt_lev" can only have values 0, 1 and 2!']);
    end
end
if ~isinf(probStruct.N) & probStruct.subopt_lev==1,
    disp(['WARNING: Prediction horizon "' psn '.N" should be set to Infinity for sub-optimality levels 1!']);
    probStruct.N = Inf;
end

%% check for proper dimension of probStruct.Q
if iscell(probStruct.Q),
    % probStruct.Q is a cell -- time-varying weights
    if length(probStruct.Q)~=probStruct.N,
        error(sprintf('"%s.Q" must be a cell array with %d elements!', psn, probStruct.N));
    end
    if probStruct.norm==2,
        for ii=1:length(probStruct.Q),
            if size(probStruct.Q{ii},1)~=size(probStruct.Q{ii},2),
                error(['All entries of "' psn '.Q" must be a square matrix!']);
            end
        end
    end
else
    % probStruct.Q is a matrix
    if probStruct.norm==2 & size(probStruct.Q,1)~=size(probStruct.Q,2),
        if Options.guierrors,
            error('Penalty on states must be a square matrix!');
        else
            error(['"' psn '.Q" must be a square matrix!']);
        end
    end
end

%% check for proper dimension of probStruct.Q
if isfield(probStruct,'Qy'),
    if iscell(probStruct.Qy),
        % probStruct.Qy is a cell -- time-varying weights
        if length(probStruct.Qy)~=probStruct.N,
            error(sprintf('"%s.Qy" must be a cell array with %d elements!', psn, probStruct.N));
        end
        if probStruct.norm==2,
            for ii=1:length(probStruct.Qy),
                if size(probStruct.Qy{ii},1)~=size(probStruct.Qy{ii},2),
                    error(['All entries of "' psn '.Qy" must be a square matrix!']);
                end
            end
        end
    else
        % probStruct.Qy is a matrix
        if probStruct.norm==2 & size(probStruct.Qy,1)~=size(probStruct.Qy,2),
            if Options.guierrors,
                error('Penalty on outputs must be a square matrix!');
            else
                error(['"' psn '.Qy" must be a square matrix!']);
            end
        end
    end
end


if iscell(probStruct.R),
    % probStruct.R is a cell -- time-varying weights
    if length(probStruct.R)~=probStruct.N,
        error(sprintf('"%s.R" must be a cell array with %d elements!', psn, probStruct.N));
    end
    if probStruct.norm==2,
        for ii=1:length(probStruct.R),
            if size(probStruct.R{ii},1)~=size(probStruct.R{ii},2),
                error(['All entries of "' psn '.R" must be a square matrix!']);
            end
        end
    end
else
    if probStruct.norm==2 & size(probStruct.R,1)~=size(probStruct.R,2),
        if Options.guierrors,
            error('Penalty on inputs must be a square matrix!');
        else
            error(['"' psn '.R" must be a square matrix!']);
        end
    end
end

if isfield(probStruct,'Rdu'),
    if probStruct.norm==2 & size(probStruct.Rdu,1)~=size(probStruct.Rdu,2),
        if Options.guierrors,
            error('Penalty on delta U must be a square matrix!');
        else
            error(['"' psn '.Rdu" must be a square matrix!']);
        end
    end
end

if probStruct.norm==2 & iscell(probStruct.Q),
    for ii=1:length(probStruct.Q),
        if probStruct.norm==2 & min(eig(probStruct.Q{ii}))<=-1e-10
            error(['All entries of "' psn '.Q" must be positive semi-definite matrices!']);
        end
    end
elseif probStruct.norm==2 & min(eig(probStruct.Q))<=-1e-10
    if Options.guierrors,
        error('Penalty on states must be a positive semi-definite matrix!');
    else
        error(['"' psn '.Q" must be a positive semi-definite matrix!']);
    end
end

if probStruct.norm==2 & iscell(probStruct.R),
    for ii=1:length(probStruct.R),
        if probStruct.norm==2 & min(eig(probStruct.R{ii}))<=-1e-10,
            error(['All entries of "' psn '.R" must be positive definite matrices!']);
        end
    end
elseif probStruct.norm==2 & min(eig(probStruct.R))<=-1e-10,
    if Options.guierrors,
        error('Penalty on inputs must be a positive definite matrix!');
    else
        error(['"' psn '.R" must be a positive definite matrix!']);
    end
end

if isfield(probStruct,'Rdu'),
    if probStruct.norm==2 & min(eig(probStruct.Rdu))<=-1e-10
        if Options.guierrors,
            error('Penalty on delta U must be a positive definite matrix!');
        else
            error(['"' psn '.Rdu" must be a positive definite matrix!']);
        end
    end
end
if isfield(probStruct,'Qy') & probStruct.norm==2,
    if iscell(probStruct.Qy),
        for ii=1:length(probStruct.Qy),
            if min(eig(probStruct.Qy{ii}))<=-1e-10
                error(['All entries of "' psn '.Qy" must be positive semi-definite matrices!']);
            end
        end
    else
        if min(eig(probStruct.Qy))<=-1e-10
            if Options.guierrors,
                error('Penalty on outputs must be a positive semi-definite matrix!');
            else
                error(['"' psn '.Qy" must be a positive semi-definite matrix!']);
            end
        end
    end
end

if isfield(probStruct,'x0bounds')
    if isfield(probStruct,'y0bounds'),
        error('"y0bounds" and "x0bounds" cannot be used simultaneously! Remove "x0bounds" from your problem structure');
    end
    if Options.verbose>0,
        disp('WARNING: parameter "x0bounds" is obsolete, use "y0bounds" instead.');
    end
    probStruct.y0bounds=probStruct.x0bounds;
end
if ~isfield(probStruct,'y0bounds')
    if Options.verbose > 1,
        disp('Constraints on y0 enabled');
    end
    probStruct.y0bounds=1;
end

if ~isfield(probStruct,'tracking')
    if Options.verbose>1,
        disp('Tracking set to OFF');
    end
    probStruct.tracking=0;
end
if probStruct.tracking,
    %probStruct.Tconstraint = 0;
end
if ~isfield(probStruct,'Tconstraint')
    if probStruct.norm==2 & Options.verbose>1,
        disp('No user-specified terminal constraint. Using stabilizing terminal set...');
    end
    probStruct.Tconstraint=1;
end
if probStruct.Tconstraint==0,
    if ~isfield(probStruct,'P_N') & Options.verbose >= 0,
        disp(['No terminal weight "' psn '.P_N" specified... closed-loop stability will not be guaranteed']);
    end
end
if probStruct.Tconstraint==2,
    if ~isfield(probStruct,'Tset')
        if Options.guierrors,
            error('Terminal set must be specified!');
        else
            error(['Terminal set "' psn '.Tset" must be specified if "' psn '.Tconstraint=2"!']);
        end
    end
    if ~isa(probStruct.Tset,'polytope')
        if Options.guierrors,
            error('Terminal set must be a POLYTOPE object!');
        else
            error(['"' psn '.Tset" must be a POLYTOPE object!']);
        end
    end
    if ~isfulldim(probStruct.Tset),
        if Options.guierrors,
            error('Terminal set must not be an empty polytope!');
        else
            error(['"' psn '.Tset" must not be an empty polytope!']);
        end
    end
end
if isfield(probStruct,'Tset') & isfulldim(probStruct.Tset) & probStruct.Tconstraint~=2,
    if probStruct.tracking==0,
        if Options.verbose>0,
            disp(['WARNING: Terminal set provided, you should set "' psn '.Tconstraint=2" to enforce it...']);
        end
        probStruct.Tconstraint=2;
    else
        disp(['WARNING: Terminal set cannot be used for tracking problems, ignoring "' psn '.Tset" ...']);
        probStruct.Tset = polytope;
        probStruct.Tconstraint = 0;
    end
end
if ~isfield(probStruct,'useSymmetry')
    if Options.verbose>1,
        disp('Symmetry will not be used');
    end
    probStruct.useSymmetry=0;
end
if ~isfield(probStruct,'feedback')
    if Options.verbose>1,
        disp('No Pre-Stabilization with Feedback...')
    end
    probStruct.feedback=0;
elseif probStruct.feedback==1
    if Options.verbose>0,
        disp('pre-Stabilization with feedback enabled...')
    end
end
if ~isfield(probStruct,'Qy')
else
    if iscell(probStruct.Qy),
        if probStruct.Tconstraint==1,
            disp(['Not possible to have stabilizing terminal set in output tracking case, setting "' psn '.Tconstraint=0"... '])
            probStruct.Tconstraint=0;
        end
    elseif(~isempty(probStruct.Qy) & ~any(probStruct.Qy)~=0)
        if(probStruct.Tconstraint==1)
            disp(['Not possible to have stabilizing terminal set in output tracking case, setting "' psn '.Tconstraint=0"... '])
            probStruct.Tconstraint=0;
        end
    end
end

if probStruct.feedback & probStruct.norm~=2,
    error('Pre-Stabilization is implemented only for 2-norm.')
end
if probStruct.Tconstraint~=2,
    probStruct.Tset=polytope;
end

if isfield(probStruct,'yref') & isfield(probStruct,'xref')
    error(['"' psn '.yref" must not be used in connection with ' psn '.xref!']);
end
    
if isfield(probStruct,'yref') & isfield(probStruct,'uref')
    error(['"' psn '.yref" must not be used in connection with ' psn '.uref!']);
end

if isfield(probStruct,'yref') & ~isfield(probStruct,'Qy')
    if Options.guierrors,
        error('Please define penalty on outputs for output regulation problems!');
    else
        error(['Please set "' psn '.Qy" for output regulation problems!']);
    end
end

if probStruct.tracking & probStruct.Tconstraint==2,
    if Options.guierrors,
        error('User-defined target set not supported for tracking problems!');
    else
        error(['User-defined target sets not supported for tracking problems. Please set "' psn '.Tconstraint=0".']);
    end
end

if probStruct.tracking & probStruct.subopt_lev>0
    error(['Tracking not supported for time-optimal and low-complexity strategies. Please set "' psn '.subopt_lev=0".']);
end

if probStruct.tracking & (isfield(probStruct,'yref')),
    error(['Tracking cannot be used in combination with "' psn '.yref". Either set "' psn '.tracking=0" or remove "' psn '.yref".']);
end

if probStruct.tracking & (isfield(probStruct,'xref') | isfield(probStruct,'uref')),
    error(['Tracking cannot be used in combination with "' psn '.xref". Either set "' psn '.tracking=0" or remove "' psn '.xref".']);
end

if probStruct.subopt_lev>0 & probStruct.Tconstraint==0,
    error(sprintf('"%s.Tconstraint" must be either 1 or 2 for "%s.subopt_lev=%d".', psn, psn, probStruct.subopt_lev));
end

if (probStruct.subopt_lev>0 | probStruct.N==Inf) & iscell(probStruct.Q) 
    if probStruct.subopt_lev>0,
        error(sprintf('Time-varying "%s.Q" not allowed for "%s.subopt_lev=%d"!', psn, psn, probStruct.subopt_lev));
    elseif isinf(probStruct.N)
        error(sprintf('Time-varying "%s.Q" not allowed for "%s.N=Inf"!', psn, psn));
    end
end

if (probStruct.subopt_lev>0 | probStruct.N==Inf) & iscell(probStruct.R) 
    if probStruct.subopt_lev>0,
        error(sprintf('Time-varying "%s.R" not allowed for "%s.subopt_lev=%d"!', psn, psn, probStruct.subopt_lev));
    elseif isinf(probStruct.N)
        error(sprintf('Time-varying "%s.R" not allowed for "%s.N=Inf"!', psn, psn));
    end
end

if isinf(probStruct.N) & probStruct.subopt_lev==2,
    if Options.guierrors,
        error('Prediction horizon must be finite for low-complexity solutions!');
    else
        error(['"' psn '.N" must be finite for "' psn '.subopt_lev=2"!']);
    end
end

if isinf(probStruct.N) & probStruct.subopt_lev==0,
    if isfulldim(probStruct.Tset),
        error(['"' psn '.Tset" not allowed for Infinite-time solution!']);
    elseif (probStruct.feedback | isfield(probStruct,'FBgain')),
        error('Feedback pre-stabilization not supported for infinite-time solutions!');
    elseif probStruct.Tconstraint~=1,
        error(['"' psn '.Tconstraint" must be 1 for infinite-time solutions!']);
    elseif (isfield(probStruct,'Qy') | isfield(probStruct,'yref')),
        error('Output tracking not allowed for infinite-time solutions!');
    elseif (isfield(probStruct,'xref') | isfield(probStruct,'uref')),
        error('Fixed state tracking not allowed for infinite-time solutions!');
    elseif probStruct.tracking > 0,
        error('Tracking not allowed for infinite-time solutions.');
    end
end

if probStruct.feedback,
    if probStruct.tracking > 0,
        error('Feedback pre-stabilization not supported for tracking problems.');
    end
end


if probStruct.Tconstraint==2 & isfield(probStruct,'yref'),
    error(['"' psn '.Tset" and "' psn '.yref" are contradicting objectives, please remove one of them.']);
end

if isfield(probStruct, 'Nc'),
    if ~Options.useyalmip,
        % mpt_constructMatrices cannot deal with certain problem formulations
        if probStruct.norm~=2 & probStruct.Nc ~= probStruct.N,
            error('Control horizon is only supported for quadratic cost function.');
        end
        if probStruct.subopt_lev > 0 & probStruct.Nc ~= probStruct.N,
            error('Control horizon must be equal to prediction horizon for low-complexity strategies!');
        end
        if probStruct.Nc > probStruct.N,
            error(['"' psn '.Nc" must be smaller than "' psn '.N" !']);
        end
        if probStruct.Nc < 1,
            error(['"' psn '.Nc" must be greater 0 !']);
        end
        %     if isfield(probStruct, 'inputblocking'),
        %         error(['"' psn '.Nc" cannot be used in combination with "' psn '.inputblocking" !']);
        %     end
        %     if isfield(probStruct, 'deltablocking'),
        %         error(['"' psn '.Nc" cannot be used in combination with "' psn '.deltablocking" !']);
        %     end
        if probStruct.Nc < probStruct.N
            probStruct.inputblocking = [ones(1, probStruct.Nc-1) probStruct.N-probStruct.Nc+1];
        end
    end
end

if isfield(probStruct, 'Qy') & probStruct.subopt_lev==1,
    error('Output weights not supported for minimum-time solutions.');
end

probStruct.verified=1;
