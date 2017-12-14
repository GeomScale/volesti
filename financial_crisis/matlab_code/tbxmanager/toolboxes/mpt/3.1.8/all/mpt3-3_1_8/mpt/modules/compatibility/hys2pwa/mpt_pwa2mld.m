function S = mpt_pwa2mld(sysStruct)
%MPT_PWA2MLD Converts a PWA system into the MLD representation
%
% S = mpt_pwa2mld(sysStruct)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Converts a PWA system to an equivalent MLD description:
%
%        x(k+1) = A x + B1 u + B2 d + B3 z
%         y(k)  = C x + D1 u + D2 d + D3 z
%  E2 d + E3 z <= E1 u + E4 x + E5
%
% The conversion is done by using one binary variable per one region of the PWA
% system.
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% sysStruct - System structure
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% S         - Structure containing matrices of the MLD model
%

% Copyright is with the following author(s):
%
% (C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
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

if ~isfield(sysStruct, 'verified'),
    sysStruct = mpt_verifySysStruct(sysStruct);
end

if ~iscell(sysStruct.A),
    error('System must be piecewise-affine!');
end

if ~isfield(sysStruct, 'xmax') | ~isfield(sysStruct, 'xmin')
    error('State constraints must be provided!');
end

[nx,nu,ny,ndyn,nbool,ubool] = mpt_sysStructInfo(sysStruct);

E1 = []; E2 = []; E3 = []; E4 = []; E5 = [];
A = []; B1 = []; B2 = []; B3 = [];
C = []; D1 = []; D2 = []; D3 = [];

nz = ndyn*nx;
ndelta = ndyn;

xumin = [sysStruct.xmin; sysStruct.umin];
xumax = [sysStruct.xmax; sysStruct.umax];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Si*x + Ri*u - Ti <= Mstar_i * (1 - d_i)
for idyn = 1:ndyn,
    Si = sysStruct.guardX{idyn};
    Ri = sysStruct.guardU{idyn};
    Ti = sysStruct.guardC{idyn};
    nc = size(Si,1);
    
    [mstar,m] = derivebounds([Si Ri], -Ti, xumin, xumax);
    Mstar = repmat(max(mstar), nc, 1);
    
    e4 = -Si;
    e1 = -Ri;
    e5 = Ti + Mstar;
    e2 = Mstar;
    e3 = zeros(nc, nz);
    
    E1 = [E1; e1];
    if idyn==1,
        E2 = e2;
    else
        E2 = [E2 zeros(size(E2, 1),1); zeros(nc, idyn-1) e2];
    end
    E3 = [E3; e3];
    E4 = [E4; e4];
    E5 = [E5; e5];
end

%%%%%%%%%%%%%%
% \sum d_i = 1
e2 = [ones(1, ndelta); -ones(1, ndelta)];
e3 = zeros(2, nz);
e4 = zeros(2, nx);
e5 = [1; -1];
e1 = zeros(2, nu);

E1 = [E1; e1];
E2 = [E2; e2];
E3 = [E3; e3];
E4 = [E4; e4];
E5 = [E5; e5];


%%%%%%%%%%%%%%%%%%%
% x(t+1) = \sum z_i
% zi <= M * di
% zi >= m * di
% zi <= Ai * x + Bi * u + fi - m + m * di
% zi >= Ai * x + Bi * u + fi - M + M * di
[A, B1, B2, B3, aE1, aE2, aE3, aE4, aE5] = defineaux(sysStruct);
E1 = [E1; aE1];
E2 = [E2; aE2];
E3 = [E3; aE3];
E4 = [E4; aE4];
E5 = [E5; aE5];



%%%%%%%%%%%%%%%%%%%%%%
% y(t) = C*x + D*u + g

% first detect if output equation is identical for all dynamics
allsame = 1;
for idyn = 2:ndyn,
    if any(any(sysStruct.C{idyn}~=sysStruct.C{1})) | ...
            any(any(sysStruct.D{idyn}~=sysStruct.D{1})) | ...
            any(any(sysStruct.g{idyn}~=sysStruct.g{1}))
        allsame = 0;
        break
    end
end

if allsame,
    % special case - output equation identical for all dynamics
    % but no affine "g" term is allowed
    if all(sysStruct.g{1}==0),
        C = sysStruct.C{1};
        D1 = sysStruct.D{1};
        D2 = zeros(ny, ndelta);
        D3 = zeros(ny, nz);
    else
        % non-zero "g" terms must be handled specially
        allsame = 0;
    end
end
if ~allsame
    % different output equations in different dynamics and/or non-zero "g" term

    % NOTE! NOTE! NOTE!
    % we can actually do better than this. currently we assign one additional
    % auxiliary variable per one dynamics, ignoring the fact that several
    % dynamics can share the same output equation. if it is the case, we can
    % reduce the number of auxiliary variables and just work with their
    % combination.
    
    %%%%%%%%%%%%%%%%%%%
    % y(t) = \sum z_j
    % zi <= M * di
    % zi >= m * di
    % zi <= Ci * x + Di * u + gi - m + m * di
    % zi >= Ci * x + Di * u + gi - M + M * di
    [C, D1, D2, D3, aE1, aE2, aE3, aE4, aE5] = defineaux(sysStruct, 'Yequation');
    
    E1 = [E1; aE1];
    E2 = [E2; aE2];
    E4 = [E4; aE4];
    E5 = [E5; aE5];
    
    % introduce additional auxiliary "z" variables which describe the output
    % equation
    E3 = [E3 zeros(size(E3, 1), ny*ndyn); zeros(size(aE3, 1), nz) aE3];
    
    % we introduced 'ny*ndyn' additional auxiliary variables, increase the counter
    nz = nz + ny*ndyn;
end


S.A = A;
S.B1 = B1;
S.B2 = B2;
S.B3 = B3;
S.C = C;
S.D1 = D1;
S.D2 = D2;
S.D3 = D3;
S.E1 = E1;
S.E2 = E2;
S.E3 = E3;
S.E4 = E4;
S.E5 = E5;
S.name = 'mpt_pwa2mld_autoconversion';
S.MLDisvalid = 1;
S.MLDstructver = 2;
S.MLDsymtable = 0;
S.MLDrowinfo = 0;

S.rowinfo = [];
S.symtable = [];

S.ne = size(E1, 1);
S.nxr = nx;  % check this!
S.nxb = 0;   % check this!
S.nx = nx;
S.nur = nu-nbool;
S.nub = nbool;
S.nu = nu;
S.nyr = ny;  % check this!
S.nyb = 0;   % check this!
S.ny = ny;
S.nd = ndelta;
S.nz = nz;
S.ul = sysStruct.umin;
S.uu = sysStruct.umax;
S.xl = sysStruct.xmin;
S.xu = sysStruct.xmax;
S.yl = sysStruct.ymin;
S.yu = sysStruct.ymax;
S.dl = -inf * ones(S.nd, 1);
S.du = +inf * ones(S.nd, 1);
S.zl = -inf * ones(S.nz, 1);
S.zu = +inf * ones(S.nz, 1);
return

%---------------------------------------------------
function [A, B1, B2, B3, E1, E2, E3, E4, E5] = defineaux(sysStruct, yeq),

[nx,nu,ny,ndyn,nbool,ubool] = mpt_sysStructInfo(sysStruct);
nxt = nx;
if nargin>1,
    % call this function with 2 input arguments to define auxiliary variables
    % for the output equation
    Ai = sysStruct.C;
    Bi = sysStruct.D;
    fi = sysStruct.g;
    
    % this is important to set since the code below was originally written to
    % handle only the state-update equation (hence it uses 'nx')
    nx = ny;
else
    Ai = sysStruct.A;
    Bi = sysStruct.B;
    fi = sysStruct.f;
end

E1 = []; E2 = []; E3 = []; E4 = []; E5 = [];
A = []; B1 = []; B2 = []; B3 = [];
C = []; D1 = []; D2 = []; D3 = [];


nz = ndyn*nx;
ndelta = ndyn;

xumin = [sysStruct.xmin; sysStruct.umin];
xumax = [sysStruct.xmax; sysStruct.umax];

%%%%%%%%%%%%%%%%%%%%%%
% estimate "M" and "m"
Ms = [];
ms = [];
for idyn = 1:ndyn,
    [M,m] = derivebounds([Ai{idyn} Bi{idyn}], fi{idyn}, xumin, xumax);
    Ms = [Ms; max(M)];
    ms = [ms; min(m)];
end
Ms = max(Ms);
ms = min(ms);

%%%%%%%%%%%%%%%%%%%
% x(t+1) = \sum z_i
A = zeros(nx);
for idyn = 1 : ndyn,
    B3 = [B3 eye(nx)];
end
B2 = zeros(nx, ndelta);
B1 = zeros(nx, nu);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% zi <= M * di
% zi >= m * di
% zi <= Ai * x + Bi * u + fi - m + m * di
% zi >= Ai * x + Bi * u + fi - M + M * di
for idyn = 1:ndyn,
    M = repmat(Ms, nx, 1);
    m = repmat(ms, nx, 1);
    aE1 = []; aE2 = []; aE3 = []; aE4 = []; aE5 = [];
    
    % zi <= M * di
    e3 = eye(nx);
    e2 = -M;
    e1 = zeros(nx, nu);
    e4 = zeros(nx, nxt);
    e5 = zeros(nx, 1);
    aE1 = [aE1; e1]; aE2 = [aE2; e2]; aE3 = [aE3; e3]; aE4 = [aE4; e4]; aE5 = [aE5; e5];
   
    % zi >= m * di
    e3 = -eye(nx);
    e2 = m;
    e1 = zeros(nx, nu);
    e4 = zeros(nx, nxt);
    e5 = zeros(nx, 1);
    aE1 = [aE1; e1]; aE2 = [aE2; e2]; aE3 = [aE3; e3]; aE4 = [aE4; e4]; aE5 = [aE5; e5];
    
    % zi <= Ai * x + Bi * u + fi - m + m * di
    e3 = eye(nx);
    e2 = -m;
    e1 = Bi{idyn};
    e4 = Ai{idyn};
    e5 = fi{idyn} - m;
    aE1 = [aE1; e1]; aE2 = [aE2; e2]; aE3 = [aE3; e3]; aE4 = [aE4; e4]; aE5 = [aE5; e5];
        
    % zi >= Ai * x + Bi * u + fi - M + M * di
    e3 = -eye(nx);
    e2 = M;
    e1 = -Bi{idyn};
    e4 = -Ai{idyn};
    e5 = -fi{idyn} + M;
    aE1 = [aE1; e1]; aE2 = [aE2; e2]; aE3 = [aE3; e3]; aE4 = [aE4; e4]; aE5 = [aE5; e5];
    
    E1 = [E1; aE1];
    E2 = [E2; [zeros(nx*4, idyn-1) aE2 zeros(nx*4, ndyn-idyn)]];
    E3 = [E3; [zeros(nx*4, nx*(idyn-1)) aE3 zeros(nx*4, nx*(ndyn-idyn))]];
    E4 = [E4; aE4];
    E5 = [E5; aE5];
end
return

%---------------------------------------------------
function [M,m] = derivebounds(A, b, lb, ub)
% derive min/max bounds on each row of
% A*z + b
% lb <= z <= ub

M = repmat(-inf,size(b,1),1);
m = -M;
lower = lb;
upper = ub;
for i = 1:length(b)
    ind = A(i,:)>0;
    i1 = find(ind);
    i2 = find(~ind);
    a = A(i,[i1 i2]);
    M(i) = a*[upper(i1);lower(i2)];
    m(i) = a*[lower(i1);upper(i2)];
end
M(isnan(M)) = inf;
m(isnan(m)) = -inf;
M(isinf(M)) = 1e4; M = M+b+0.1;
m(isinf(m)) = -1e4;m = m+b-0.1;

return
