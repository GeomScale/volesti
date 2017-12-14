function [nx,nu,ny,ndyn,nbool,ubool,intInfo] = mpt_sysStructInfo(sysStruct)
%MPT_SYSSTRUCTINFO Returns information about system structure
%
% [nx,nu,ny,ndyn,nbool,ubool,intInfo] = mpt_sysStructInfo(sysStruct)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Returns number of states, inputs, outputs and number of dynamics contained in
% a given system structure
%
% When run with no output arguments, displays textual information about a given
% system, i.e.
%   mpt_sysStructInfo(sysStruct)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% sysStruct  - system structure describing an LTI system
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% nx       - number of states
% nu       - number of control inputs
% ny       - number of outputs
% ndyn     - number of dynamics
% nbool    - number of boolean inputs
% ubool    - indexes of integer (or boolean) inputs
% intInfo  - structure with information about overlapping dynamics
%   .Xintersect - cell array, "Xintersect{i}" is a vector of indices of dynamics
%                 which intersect with dynamics "i" in the X-space
%   .PdynX      - cell array, "PdynX{i}" is a polytope which defines dynamics
%                 "i" in  the X-space
%

% Copyright is with the following author(s):
%
% (C) 2003-2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
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

narginchk(1, 1);

if ~isstruct(sysStruct),
    error('mpt_sysStructInfo: Input argument must be a sysStruct structure!');
end

if ~isfield(sysStruct,'verified')
    verOpt.verbose=0;
    sysStruct=mpt_verifySysStruct(sysStruct,verOpt);
end

ispwa = iscell(sysStruct.A);
if ispwa,
    % PWA system
    ndyn = length(sysStruct.A);
    nx = size(sysStruct.A{end},2);
    nu = size(sysStruct.B{end},2);
    ny = size(sysStruct.C{end},1);
else
    % LTI system
    ndyn = 1;
    nx = size(sysStruct.A,2);
    nu = size(sysStruct.B,2);
    ny = size(sysStruct.C,1);
end
nbool = 0;
ubool = [];

if nargout==0,
    % display textual information
    if ispwa,
        systype = sprintf('PWA system (%d dynamics)', ndyn);
    else
        systype = 'LTI system';
    end
    plurx = [num2str(nx) ' state' repmat('s', 1, (nx>1))];
    pluru = [num2str(nu) ' input' repmat('s', 1, (nu>1))];
    plury = [num2str(ny) ' output' repmat('s', 1, (ny>1))];   
    fprintf('\n%s, %s, %s, %s\n', systype, plurx, pluru, plury);
    if ispwa,
        % print information about each dynamics
        for idyn = 1:ndyn,
            gC = sysStruct.guardC{idyn};
            nc = length(gC);
            if isfield(sysStruct, 'guardX'),
                gX = sysStruct.guardX{idyn};
            else
                gX = zeros(nc, nx);
            end
            if isfield(sysStruct, 'guardU'),
                gU = sysStruct.guardU{idyn};
            else
                gU = zeros(nc, nu);
            end
            fprintf('\nGuards of dynamics %d:\n', idyn);
            for ii = 1:nc,
                fprintf('%s\n', char(mptaffexpr(nx, nu, gX(ii, :), gU(ii, :), gC(ii, :))))
            end
        end
    end
    fprintf('\n');
    clear nx
end
    
if nargout <= 4,
    return
end

if isfield(sysStruct,'Uset')
    % boolean inputs present
    if iscell(sysStruct.Uset),
        for ii=1:length(sysStruct.Uset),
            if any(isinf(sysStruct.Uset{ii})) | any(isinf(-sysStruct.Uset{ii})),
                % this input is continuous
            else
                nbool = nbool+1;
                ubool = [ubool; ii];
            end
        end
    else
        if any(isinf(sysStruct.Uset)) | any(isinf(-sysStruct.Uset)),
            nbool = 0;
        else
            nbool = 1;
            ubool = 1;
        end
    end
else
    nbool = 0;
end

if nargout <= 6
    return
end

emptypoly = polytope;

Pdyn = {};
PdynX = {};
if iscell(sysStruct.A),
    Xintersect = {};
    for ii=1:ndyn,
        Xintersect{ii} = [];
        gX = sysStruct.guardX{ii};
        gC = sysStruct.guardC{ii};
        gU = sysStruct.guardU{ii};
        nonzerorows = [];
        for jj=1:size(gX,1)
            if any(gX(jj,:)~=0),
                nonzerorows = [nonzerorows; jj];
            end
        end
        gXnz = gX(nonzerorows,:);
        gCnz = gC(nonzerorows,:);
        gUnz = gU(nonzerorows,:);
        pureX_gX{ii} = gXnz;
        pureX_gU{ii} = gUnz;
        pureX_gC{ii} = gCnz;
        
        if isempty(gXnz),
            P = emptypoly;
        else
            P = polytope([gXnz gUnz],gCnz);
        end
        
        % project it down to X space
        if ~isfulldim(P) | isempty(gXnz)
            PX = sysStruct.Pbnd;
            PdynX{ii} = emptypoly;
        else
            try
                TT = evalc('PX = projection(P,1:nx) & sysStruct.Pbnd;');
            catch
                PX = sysStruct.Pbnd;
            end
            PdynX{ii} = PX;
        end
        Pdyn{ii} = PX;
        
    end
    for ii=1:ndyn,
        for jj=1:ndyn
            if ii==jj,
                Xintersect{ii} = [Xintersect{ii} ii];
            elseif dointersect(Pdyn{ii}, Pdyn{jj}),
                Xintersect{ii} = [Xintersect{ii} jj];
            end
        end
    end
else
    Xintersect = 1;
end


if ~iscell(sysStruct.A),
    intInfo = [];
    return
end

Dyn = {};
for ii=1:ndyn,
    dyns = Xintersect{ii};
    D = ii;
    for jj=1:length(dyns),
        if dyns(jj)~=ii,
            D = [D Xintersect{dyns(jj)}];
        end
    end
    Dyn{ii} = unique(D);
end
Davail = ones(length(Dyn),1);
for ii=1:length(Dyn),
    for jj=1:length(Dyn),
        if ii==jj | Davail(jj)==0, continue, end
        if isempty(setdiff(Dyn{ii},Dyn{jj})),
            % Dyn{ii} is a subset of Dyn{jj}
            Davail(ii)=0;
            break
        end
    end
end
ovl_dyns = find(Davail==1);
Dyns = {Dyn{ovl_dyns}};
dynamics_links = [];
for dyn=1:ndyn,
    for jj=1:length(Dyns),
        if any(Dyns{jj} == dyn)
            dynamics_links = [dynamics_links; dyn jj];
            break
        end
    end
end

intInfo.Xintersect = Xintersect;
intInfo.dyns_stack = Dyns;
intInfo.dyns_links = dynamics_links;
intInfo.stacks = length(Dyns);
intInfo.Pdyn = Pdyn;
intInfo.PdynX = PdynX;
