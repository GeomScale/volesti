function [xc,R,lambda] = chebyball_f(H,K,lpsolver)
%CHEBYBALL_F Computes center and radius of the largest ball inscribed in a polytope
%
% [xc,R,lambda] = chebyball_f(H,K,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Returns center and radius of a largest ball that can be inscribed in a polytope P
%
% NOTE: internal routine. Takes matrices H and K as input, not a polytope!
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% H,K               - Matrices defining the polytope, i.e. Hx<=K
% Options.lpsolver  - LP solver to use (see help mpt_solveLP)
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------------------------
% xc                - Center of Chebyshev's ball
% R                 - Radius of Chebyshev's ball
% lambda            - Set of lagrangian multipliers
%
% see also CHEBYBALL
%

% Copyright is with the following author(s):
%
% (C) 2003 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2003 Mato Baotic, Automatic Control Laboratory, ETH Zurich,
%          baotic@control.ee.ethz.ch

% ---------------------------------------------------------------------------
% Legal note:
%          This library is free software; you can redistribute it and/or
%          modify it under the terms of the GNU Lesser General Public
%          License as published by the Free Software Foundation; either
%          version 2.1 of the License, or (at your option) any later version.
%
%          This library is distributed in the hope that it will be useful,
%          but WITHOUT ANY WARRANTY; without even the implied warranty of
%          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%          Lesser General Public License for more details.
% 
%          You should have received a copy of the GNU Lesser General Public
%          License along with this library; if not, write to the 
%          Free Software Foundation, Inc., 
%          59 Temple Place, Suite 330, 
%          Boston, MA  02111-1307  USA
%
% ---------------------------------------------------------------------------

global mptOptions

if nargin==3,
    Options.lpsolver = lpsolver;
else
    Options.lpsolver = mptOptions.lpsolver;
end

% If any boundary is -Inf polytope P is empty
%--------------------------------------------
if any(K==-Inf)
    [nc,nx]=size(H);
    xc=zeros(nx,1);
    R=-Inf;
    lambda=zeros(nc,1);
    return;
end

% Remove rows with Inf boundaries
%--------------------------------
ii=(K==Inf);
H(ii,:)=[];
K(ii)=[];

[nc,nx]=size(H);

if nc==0
    xc=zeros(nx,1);
    R=Inf;
    lambda=zeros(nc,1);
    return;
end

% use 'rescue' function - resolve an LP automatically if it is infeasible
% with default solver
[xopt,fval,lambda,exitflag,how]=mpt_solveLPs([zeros(1,nx) 1],[H, -sqrt(sum(H.*H,2))],...
    K, [], [], [zeros(nx,1); 1000], Options.lpsolver);

if ~strcmp(how,'ok'),
    % maybe there is a numerical problem, thus we normalize H and K
    rel_tol=mptOptions.rel_tol;    % relative tolerance
    abs_tol=mptOptions.abs_tol;    % absolute tolerance
    %[nc,nx]=size(H);
    Anorm=sqrt(sum(H .* H,2));
    ii=find(Anorm<=abs_tol | Anorm<=rel_tol*K);
    nii=length(ii);
    if nii>0
        Anorm(ii)=1;
        H(ii,:)=repmat([1 zeros(1,nx-1)],nii,1);
        % decide if some constraint is always true
        jj=(K(ii)>=-abs_tol);
        K(ii(jj))=Inf;
        % or is it always false
        K(ii(~jj))=-Inf;
    end
    temp=1./Anorm;
    %H=H .* repmat(temp,1,nx);
    H = H .* temp(:,ones(1,nx));   
    % replaces H = H .* repmat(temp,1,nx); - call to repmat is rather slow
    K=K .* temp;
    
    % If any boundary is -Inf polytope P is empty
    %--------------------------------------------
    if any(K==-Inf)
        xc=zeros(nx,1);
        R=-Inf;
        lambda=zeros(nc,1);
        return;
    end
    
    % Remove rows with Inf boundaries
    %--------------------------------
    ii=(K==Inf);
    H(ii,:)=[];
    K(ii)=[];
    
    if size(H,1)==0
        xc=zeros(nx,1);
        R=Inf;
        lambda=zeros(nc,1);
        return;
    end
    
    x0 = [zeros(nx,1); 1000];         % hard-coded initial conditions
    
    % use 'rescue' function - resolve an LP automatically if it is infeasible
    % with default solver
    [xopt,fval,lambda,exitflag,how]=mpt_solveLPs([zeros(1,nx) 1],[H, -sqrt(sum(H.*H,2))],...
        K,[],[],x0,Options.lpsolver);
end

xc=xopt(1:nx); % Center of the ball
R=-xopt(nx+1); % Radius of the ball
