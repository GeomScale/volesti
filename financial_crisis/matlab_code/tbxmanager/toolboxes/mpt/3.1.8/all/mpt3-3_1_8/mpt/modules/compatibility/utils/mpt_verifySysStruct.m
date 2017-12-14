function sysStruct=mpt_verifySysStruct(sysStruct,Options)
%MPT_VERIFYSYSSTRUCT Verifies the sysStruct structure
%
% sysStruct=mpt_verifySysStruct(sysStruct,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Verifies if the structure defining the system dynamics and constraints has all
% fields filled out.
%
% Internal function.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% sysStruct        - System structure in the probStruct format
% Options.verbose  - level of verbosity 
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% sysStruct        - System structure in the sysStruct format
%
% see also MPT_VERIFYPROBSTRUCT, MPT_CONTROL
%

% Copyright is with the following author(s):
%
% (C) 2003 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch

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
    'sysstructname', inputname(1), ...
    'guierrors', 0, ...
    'ybounds_optional', 1);

if isempty(Options.sysstructname),
    Options.sysstructname = 'sysStruct';
end

ssn = Options.sysstructname;

if ~isstruct(sysStruct),
    error('mpt_verifySysStruct: Input argument must be a structure! See help mpt_sysStruct.');
end

if ~isfield(sysStruct,'A')
    if Options.guierrors,
        error('"A" must be defined!');
    else
        error(['"' ssn '.A" must be defined!']);
    end
end
if ~isfield(sysStruct,'B')
    if Options.guierrors,
        error('"B" must be defined!');
    else
        error(['"' ssn '.B" must be defined!']);
    end
end
if ~isfield(sysStruct,'C')
    if Options.guierrors,
        error('"C" must be defined!');
    else
        error(['"' ssn '.C" must be defined!']);
    end
end
if ~isfield(sysStruct,'D')
    if Options.guierrors,
        error('"D" must be defined!');
    else
        error(['"' ssn '.D" must be defined!']);
    end
end

if ~isfield(sysStruct, 'umin')
    if Options.guierrors,
        error('Lower bound on u(k) must be defined!');
    else
        error(['"' ssn '.umin" must be defined!']);
    end
end
if ~isfield(sysStruct, 'umax')
    if Options.guierrors,
        error('Upper bound on u(k) must be defined!');
    else
        error(['"' ssn '.umax" must be defined!']);
    end
end
if ~isfield(sysStruct, 'ymin') & Options.ybounds_optional==0,
    if Options.guierrors,
        error('Lower bound on y(k) must be defined!');
    else
        error(['"' ssn '.ymin" must be defined!']);
    end
end
if ~isfield(sysStruct, 'ymax') & Options.ybounds_optional==0,
    if Options.guierrors,
        error('Upper bound on y(k) must be defined!');
    else
        error(['"' ssn '.ymax" must be defined!']);
    end
end

if (~isfield(sysStruct, 'ymax') & ~isfield(sysStruct, 'ymin')) & ...
        (~isfield(sysStruct, 'xmax') & ~isfield(sysStruct, 'xmin')),
    error('Either output or state constraints must be provided.');
end

if iscell(sysStruct.A),
    sysStruct.type='PWA';
else
    sysStruct.type='LTI';
end

if strcmp(sysStruct.type,'LTI') | sysStruct.type==0,
    if size(sysStruct.C,1)~=size(sysStruct.D,1) 
        if Options.guierrors,
            error('"C" and "D" matrices have incompatible dimensions!');
        else
            error(['"' ssn '.C" and "' ssn '.D" have incompatible dimensions!']);
        end
    end
    if size(sysStruct.B,2)~=size(sysStruct.D,2)
        if Options.guierrors,
            error('"B" and "D" matrices have incompatible dimensions!');
        else
            error(['"' ssn '.B" and "' ssn '.D" have incompatible dimensions!']);
        end
    end
    if isfield(sysStruct,'f')
        if ~isempty(sysStruct.f) | ~all(sysStruct.f==0),
%             sysStruct.f=zeros(size(sysStruct.A,1),1);
%             disp('Nonzero affine terms not allowed for LTI systems! Setting them to zero...');
        end
        if size(sysStruct.f,1)~=size(sysStruct.A,1) | size(sysStruct.f,2)~=1,
            if Options.guierrors,
                error('Incompatible dimension of "A" and "f" matrices!');
            else
                error(sprintf('Incompatible dimension of the affine term "%s.f" ! It must be a %dx1 vector.', ...
                    ssn, size(sysStruct.A,1)));
            end
        end

    end
    if isfield(sysStruct,'g')
        if size(sysStruct.g,1)~=size(sysStruct.C,1) | size(sysStruct.g,2)~=1,
            if Options.guierrors,
                error('Incompatible dimension of "C" and "g" matrices!');
            else
                error(sprintf('Incompatible dimension of the affine term "%s.g" ! It must be a %dx1 vector.', ...
                    ssn, size(sysStruct.C,1)));
            end
        end
    end
    if isfield(sysStruct,'guardX') | isfield(sysStruct,'guardU') | isfield(sysStruct,'guardC')
        error('No guardlines allowed for LTI systems! Use "ymin", "ymax", "umin", "umax" to specify constraints');
    end
    if isfield(sysStruct,'Aunc'),
        if ~isfield(sysStruct,'Bunc'),
            error(['You must provide a pair "' ssn '.Aunc" and "' ssn '.Bunc" !']);
        end
        if ~iscell(sysStruct.Aunc),
            error(['"' ssn '.Aunc" must be a cell array!']);
        end
        if ~iscell(sysStruct.Bunc),
            error(['"' ssn '.Bunc" must be a cell array!']);
        end
        if length(sysStruct.Aunc)~=length(sysStruct.Bunc),
            error(['Uncertainty matrices "' ssn '.Aunc" and "' ssn '.Bunc" must have the same number of elements!']);
        end
        for ii=1:length(sysStruct.Aunc),
            if ~all(size(sysStruct.Aunc{ii})==size(sysStruct.A)) | ~all(size(sysStruct.Bunc{ii})==size(sysStruct.B))
                error(['Uncertainty matrices "' ssn '.Aunc" and "' ssn '.Bunc" must have the same dimension as A and B!']);
            end
        end
    end
    if iscell(sysStruct.A) & length(sysStruct.A)==1,
        fprintf('\nWARNING: The given PWA system consists only of one dynamics.\n');
        fprintf('Convert it to an LTI system by removing cell arrays to speed up computation.\n\n');
    end
    A=sysStruct.A;
    B=sysStruct.B;
    C=sysStruct.C;
else
    if ~iscell(sysStruct.A)
        error(['"' ssn '.A" must be a cell array!']);
    end
    if ~iscell(sysStruct.B)
        error(['"' ssn '.B" must be a cell array!']);
    end
    if ~iscell(sysStruct.C)
        error(['"' ssn '.C" must be a cell array!']);
    end
    if ~iscell(sysStruct.D)
        error(['"' ssn '.D" must be a cell array!']);
    end
    if ~isfield(sysStruct,'guardX')
        if Options.guierrors,
            error('Active region must be specified for PWA systems!');
        else
            error(['"' ssn '.guardX" must be defined for PWA systems!']);
        end
    end
    if ~isfield(sysStruct,'guardC')
        if Options.guierrors,
            error('Active region must be specified for PWA systems!');
        else
            error(['"' ssn '.guardC" must be defined for PWA systems!']);
        end
    end
    if length(sysStruct.A)~=length(sysStruct.B) | length(sysStruct.A)~=length(sysStruct.C) | length(sysStruct.A)~=length(sysStruct.D),
        error(['"' ssn '.A", "' ssn '.B", "' ssn '.C", and "' ssn '.D" must have same number of elemenets!']);
    end
    if length(sysStruct.A)~=length(sysStruct.guardX) | length(sysStruct.A)~=length(sysStruct.guardC),
        if Options.guierrors,
            error('Active region must be specified for each dynamics!');
        else
            error('Number of guardlines inconsistent with number of dynamics!');
        end    
    end
    if ~isfield(sysStruct,'guardC') | isempty(sysStruct.guardC),
        if Options.guierrors,
            error('Active region must not be empty!');
        else
            error(['Guardline "' ssn '.guardC" must be defined!']);
        end
    end
    if ~isfield(sysStruct,'guardX') | isempty(sysStruct.guardX),
        error(['Guardline "' ssn '.guardX" must be defined!']);
    end
    if ~isfield(sysStruct,'guardU') | isempty(sysStruct.guardU),
        for ii=1:length(sysStruct.A),
            sysStruct.guardU{ii} = zeros(size(sysStruct.guardX{ii},1),size(sysStruct.B{ii},2));
        end
    else
        if length(sysStruct.A)~=length(sysStruct.guardC),
            error('Number of guardlines inconsistent with number of dynamics!');
        end
    end
    if ~isfield(sysStruct,'f'),
        for ii=1:length(sysStruct.A),
            sysStruct.f{ii}=zeros(size(sysStruct.B{ii},1),1);
        end
        if Options.verbose>1,
            disp(['No affine "' ssn '.f" term specified, assuming zero.']);
        end
    end
    for ii=1:length(sysStruct.f),
        if size(sysStruct.f{ii},1)~=size(sysStruct.A{ii},1) | size(sysStruct.f{ii},2)~=1
            if Options.guierrors,
                error(sprintf('Incompatible dimension of "f" matrix for dynamics %d.', ii));
            else
                error(sprintf('Incompatible dimension of the affine term "%s.f{%d}" ! It must be a %dx1 vector.', ...
                    ssn, ii, size(sysStruct.A{ii},1)));
            end
        end
    end
    if isfield(sysStruct,'Aunc') | isfield(sysStruct,'Bunc'),
        error(['Parametric uncertainty not allowed for PWA systems! Remove fields ' ssn '.Aunc and ' ssn '.Bunc from the system structure.']);
    end
    if ~isfield(sysStruct,'g'),
        for ii=1:length(sysStruct.A),
            sysStruct.g{ii}=zeros(size(sysStruct.C{ii},1),1);
        end
        if Options.verbose>1,
            disp(['No affine "' ssn '.g" term specified, assuming zero.']);
        end
    end
    for ii=1:length(sysStruct.g),
        if  size(sysStruct.g{ii},1)~=size(sysStruct.C{ii},1) | size(sysStruct.g{ii},2)~=1,
            if Options.guierrors,
                error(sprintf('Incompatible dimension of "g" matrix for dynamics %d.', ii));
            else
                error(sprintf('Incompatible dimension of the affine term "%s.g{%d}" ! It must be a %dx1 vector.', ...
                    ssn, ii, size(sysStruct.C{ii},1)));
            end
        end
    end
    for ii=1:length(sysStruct.A),
        if size(sysStruct.A{ii},1)~=size(sysStruct.A{ii},2),
            if Options.guierrors,
                error(sprintf('"A" matrix in dynamics %d must be a square matrix!', ii));
            else
                error(sprintf('"%s.A{%d}" must be a square matrix!', ssn, ii));
            end
        end
        if size(sysStruct.A{ii},1)~=size(sysStruct.B{ii},1) 
            if Options.guierrors,
                error(sprintf('"A" and "B" matrices in dynamics %d have incompatible dimensions!', ii));
            else
                error(sprintf('"%s.A{%d}" and "%s.B{%d}" have incompatible dimensions!', ssn, ii, ssn, ii));
            end
        end
        if size(sysStruct.C{ii},1)~=size(sysStruct.D{ii},1) 
            if Options.guierrors,
                error(sprintf('"C" and "D" matrices in dynamics %d have incompatible dimensions!', ii));
            else
                error(sprintf('"%s.C{%d}" and "%s.D{%d}" have incompatible dimensions!', ssn, ii, ssn, ii));
            end
        end
        if size(sysStruct.B{ii},2)~=size(sysStruct.D{ii},2)
            if Options.guierrors,
                error(sprintf('"B" and "D" matrices in dynamics %d have incompatible dimensions!', ii));
            else
                error(sprintf('"%s.B{%d}" and "%s.D{%d}" have incompatible dimensions!', ssn, ii, ssn, ii));
            end
        end
    end
    A=sysStruct.A{1};
    B=sysStruct.B{1};
    C=sysStruct.C{1};
    D=sysStruct.D{1};
    for ii=1:length(sysStruct.A),
        if size(sysStruct.guardX{ii},1)~=size(sysStruct.guardU{ii},1) | ... 
            size(sysStruct.guardX{ii},1)~=size(sysStruct.guardC{ii},1) | ...
            size(sysStruct.guardX{ii},2)~=size(sysStruct.A{ii},2) | ...
            size(sysStruct.guardU{ii},2)~=size(sysStruct.B{ii},2),
        if Options.guierrors,
            error(sprintf('Incompatible dimensions of guardlines in dynamics %d!', ii));
        else
            error(sprintf('Incompatible dimensions of guardlines %s.guardX{%d}, %s.guardU{%d} and %s.guardC{%d}!', ...
                ssn, ii, ssn, ii, ssn, ii));
        end
        end
        P = polytope([sysStruct.guardX{ii} sysStruct.guardU{ii}],sysStruct.guardC{ii});
        [xcen, xrad] = chebyball(P);
        %[xcen, xrad] = polyinnerball([sysStruct.guardX{ii} sysStruct.guardU{ii}], sysStruct.guardC{ii});
        if ~isfield(sysStruct,'Uset') & xrad <= 0,
            if Options.guierrors,
                error(sprintf('Active region of dynamics %d contains equality constraints. Define inputs as boolean/integer first.', ii));
            else
                disp('Polytope in an extended XU space, P:={z | [guardX guardU]*z <= guardC} should be fully dimensional for all dynamics.');
                disp('Note: Use sysStruct.Uset to define boolean/integer inputs.');
                disp('Otherwise use standard techniques (e.g. Gauss reduction) to project the polytope P to a space where it is fully dimensional.');
                fprintf('\n');
                error('Polytope P:={z | [guardX guardU]*z <= guardC} is not fully dimensional for all dynamics!');
            end
        end
        
        if ~all(size(sysStruct.A{ii})==size(A)) | ~all(size(sysStruct.B{ii})==size(B)) | ~all(size(sysStruct.C{ii})==size(C)) | ~all(size(sysStruct.D{ii})==size(D))
            error('Dimensions of system matrices A,B,C,D must be identical in all PWA dynamics!');
        end
    end
end

% checks common for LTI and PWA systems:
sysStruct.umax = sysStruct.umax(:);
sysStruct.umin = sysStruct.umin(:);

if ~isfield(sysStruct,'dumax'),
    sysStruct.dumax=Inf*ones(size(sysStruct.umax));
end
if ~isfield(sysStruct, 'dumin'),
    sysStruct.dumin=-Inf*ones(size(sysStruct.umin));
end

if(size(sysStruct.dumax,2)>size(sysStruct.dumax,1) | size(sysStruct.dumin,2)>size(sysStruct.dumin,1))
    %disp('du constraints must be column vectors (transposing current entry)');
    sysStruct.dumax = sysStruct.dumax(:);
    sysStruct.dumin = sysStruct.dumin(:);
end
if ~isfield(sysStruct,'ymax') | ~isfield(sysStruct,'ymin'),
    if Options.ybounds_optional==1,
        sysStruct.ymax = repmat(Inf, size(C, 1), 1);
        sysStruct.ymin = repmat(-Inf, size(C, 1), 1);
    elseif Options.guierrors,
        error('Output constraints must be defined!');
    else
        error(['Output constraints must be defined! Please set "' ssn '.ymax" and "' ssn '.ymin".']);
    end
end
sysStruct.ymax = sysStruct.ymax(:);
sysStruct.ymin = sysStruct.ymin(:);

if size(sysStruct.ymax,1)~=size(C,1)
    if Options.guierrors,
        error(sprintf('Upper bound on y(k) must be a %dx1 vector! You provided a %dx1 vector.', size(C,1), length(sysStruct.ymax)));
    else
        error(sprintf('"%s.ymax" must be a %dx1 vector! You provided a %dx1 vector.', ...
            ssn, size(C,1), length(sysStruct.ymax)));
    end
end
if size(sysStruct.ymin,1)~=size(C,1),
    if Options.guierrors,
        error(sprintf('Lower bound on y(k) must be a %dx1 vector! You provided a %dx1 vector.', size(C,1), length(sysStruct.ymin)));
    else
        error(sprintf('"%s.ymin" must be a %dx1 vector! You provided a %dx1 vector.', ...
            ssn, size(C,1), length(sysStruct.ymin)));
    end
end
if size(sysStruct.umax,1)~=size(B,2)
    if Options.guierrors,
        error(sprintf('Upper bound on u(k) must be a %dx1 vector! You provided a %dx1 vector.', size(B,2), length(sysStruct.umax)));
    else
        error(sprintf('"%s.umax" must be a %dx1 vector! You provided a %dx1 vector.', ...
            ssn, size(B,2), length(sysStruct.umax)));
    end
end
if size(sysStruct.umin,1)~=size(B,2),
    if Options.guierrors,
        error(sprintf('Lower bound on u(k) must be a %dx1 vector! You provided a %dx1 vector.', size(B,2), length(sysStruct.umin)));
    else
        error(sprintf('"%s.umin" must be a %dx1 vector! You provided a %dx1 vector.', ...
            ssn, size(B,2), length(sysStruct.umin)));
    end
end

if (isfield(sysStruct,'dymax') & ~isfield(sysStruct,'dymin')) | ...
        (isfield(sysStruct,'dymin') & ~isfield(sysStruct,'dymax')),
    if Options.guierrors,
        error('You must define both lower and upper bound on y(k)-y(k+1)');
    else
        error(['Both "' ssn '.dymax" and "' ssn '.dymin" must be defined!']);
    end
end

if isfield(sysStruct,'dymax'),
    sysStruct.dymin = sysStruct.dymin(:);
    sysStruct.dymax = sysStruct.dymax(:);
    if length(sysStruct.dymin)~=size(C,1),
        if Options.guierrors,
            error(sprintf('Lower bound on y(k)-y(k+1) must be a %dx1 vector! You provided a %dx1 vector.', size(C,1), length(sysStruct.dymin)));
        else
            error(sprintf('"%s.dymin" must be a %dx1 vector! You provided a %dx1 vector.', ...
                ssn, size(C,1), length(sysStruct.dymin)));
        end
    end
    if length(sysStruct.dymax)~=size(C,1),
        if Options.guierrors,
            error(sprintf('Upper bound on y(k)-y(k+1) must be a %dx1 vector! You provided a %dx1 vector.', size(C,1), length(sysStruct.dymax)));
        else
            error(sprintf('"%s.dymax" must be a %dx1 vector! You provided a %dx1 vector.', ...
                ssn, size(C,1), length(sysStruct.dymax)));
        end
    end
end

if (isfield(sysStruct, 'xmax') & ~isfield(sysStruct, 'xmin')) | ...
        (isfield(sysStruct, 'xmin') & ~isfield(sysStruct, 'xmax'))
    if Options.guierrors,
        error('You must define both lower and upper bound on x(k).');
    else
        error(['Both "' ssn '.xmin" and "' ssn '.xmax" must be defined!']);
    end
end

if isfield(sysStruct, 'xmax')
    % we have state constraints, check if dimensions match
    sysStruct.xmax = sysStruct.xmax(:);
    sysStruct.xmin = sysStruct.xmin(:);
    if length(sysStruct.xmax)~=size(A,1)
        if Options.guierrors,
            error(sprintf('Upper bound on x(k) must be a %dx1 vector! You provided a %dx1 vector.', size(A,1), length(sysStruct.xmax)));
        else
            error(sprintf('"%s.xmax" must be a %dx1 vector! You provided a %dx1 vector.', ...
                ssn, size(A,1), length(sysStruct.xmax)));
        end
    end
    if length(sysStruct.xmin)~=size(A,1)
        if Options.guierrors,
            error(sprintf('Lower bound on x(k) must be a %dx1 vector! You provided a %dx1 vector.', size(A,1), length(sysStruct.xmin)));
        else
            error(sprintf('"%s.xmin" must be a %dx1 vector! You provided a %dx1 vector.', ...
                ssn, size(A,1), length(sysStruct.xmin)));
        end
    end
end

if ~isfield(sysStruct,'noise'),
    if Options.verbose>1,
        disp('No noise specified, assuming zero');
    end
    sysStruct.noise=polytope;
end
if ~isfield(sysStruct,'Pbnd'),
    % unitbox() is more efficient than constructing the polytope directly
    sysStruct.Pbnd = unitbox(size(A, 2), mptOptions.infbox);
else
    if ~isa(sysStruct.Pbnd,'polytope'),
        error(['"' ssn '.Pbnd" must be a polytope object!']);
    end
end
if dimension(sysStruct.Pbnd)~=size(A,2),
    error(['"' ssn '.Pbnd" must be a polytope in ' num2str(size(A,2)) 'D !']);
end
if isfield(sysStruct,'Uset'),
    if ~iscell(sysStruct.Uset),
        if size(B,2)~=1,
            error(['"' ssn '.Uset" must be a cell array with ' num2str(size(B,2)) ' elements!']);
        else
            Uset = {sysStruct.Uset};
            sysStruct.Uset = Uset;
        end
    else
        if length(sysStruct.Uset)~=size(B,2),
            error(['"' ssn '.Uset" must be a cell array with ' num2str(size(B,2)) ' elements!']);
        end
    end
end

nx = size(A,1);
nu = size(B,2);
ny = size(C,1);

if isfield(sysStruct,'StateName'),
    if ~iscell(sysStruct.StateName) & nx~=1,
        error(['"' ssn '.StateName" must be a cell array with ' num2str(nx) ' elements!']);
    elseif ~iscell(sysStruct.StateName) & nx==1,
        sysStruct.StateName = {sysStruct.StateName};
    end
end

if isfield(sysStruct,'InputName'),
    if ~iscell(sysStruct.InputName) & nu~=1,
        error(['"' ssn '.InputName" must be a cell array with ' num2str(nu) ' elements!']);
    elseif ~iscell(sysStruct.InputName) & nu==1,
        sysStruct.InputName = { sysStruct.InputName };
    end
end

if isfield(sysStruct,'OutputName'),
    if ~iscell(sysStruct.OutputName) & ny~=1,
        error(['"' ssn '.OutputName" must be a cell array with ' num2str(ny) ' elements!']);
    elseif ~iscell(sysStruct.OutputName) & ny==1,
        sysStruct.OutputName = { sysStruct.OutputName };
    end
end
    
if any(sysStruct.umax - sysStruct.umin <= 0),
    if Options.guierrors,
        error('Lower bound on u(k) must be strictly smaller than the upper bound!');
    else
        error(['"' ssn '.umin" must be strictly smaller than "' ssn '.umax"!']);
    end
end

if any(sysStruct.ymax - sysStruct.ymin <= 0),
    if Options.guierrors,
        error('Lower bound on y(k) must be strictly smaller than the upper bound!');
    else
        error(['"' ssn '.ymin" must be strictly smaller than "' ssn '.ymax"!']);
    end
end

if isfield(sysStruct, 'xmin'),
    if any(sysStruct.xmax - sysStruct.xmin <= 0),
        if Options.guierrors,
            error('Lower bound on x(k) must be strictly smaller than the upper bound!');
        else
            error(['"' ssn '.xmin" must be strictly smaller than "' ssn '.xmax"!']);
        end
    end
end

if ~isfield(sysStruct, 'Ts'),
    if Options.verbose > 1,
        fprintf('WARNING: No sampling time provided, assuming %s.Ts=1\n', ssn);
    end
    sysStruct.Ts = 1;
end

if length(sysStruct.ymax(:)) ~= length(sysStruct.ymin(:)),
    error(sprintf('"%s.ymax" and "%s.ymin" must be of same length.', ssn, ssn));
end
if length(sysStruct.umax(:)) ~= length(sysStruct.umin(:)),
    error(sprintf('"%s.umax" and "%s.umin" must be of same length.', ssn, ssn));
end
if length(sysStruct.umax(:)) ~= length(sysStruct.dumax(:)),
    error(sprintf('"%s.umax" and "%s.dumax" must be of same length.', ssn, ssn));
end
if length(sysStruct.umin(:)) ~= length(sysStruct.dumin(:)),
    error(sprintf('"%s.umin" and "%s.dumin" must be of same length.', ssn, ssn));
end
if isfield(sysStruct, 'xmax'),
    if length(sysStruct.xmax(:)) ~= length(sysStruct.xmin(:)),
        error(sprintf('"%s.xmax" and "%s.xmin" must be of same length.', ssn, ssn));
    end
end
if isfield(sysStruct, 'dymax'),
    if length(sysStruct.dymax(:)) ~= length(sysStruct.dymin(:)),
        error(sprintf('"%s.dymax" and "%s.dymin" must be of same length.', ssn, ssn));
    end
end

% check if noise has correct dimension
if isfield(sysStruct, 'noise') && isfulldim(sysStruct.noise),
	if dimension(sysStruct.noise) ~= nx
		error(sprintf('"%s.noise" has wrong dimension.', ssn));
	end
end

sysStruct.verified=1;
