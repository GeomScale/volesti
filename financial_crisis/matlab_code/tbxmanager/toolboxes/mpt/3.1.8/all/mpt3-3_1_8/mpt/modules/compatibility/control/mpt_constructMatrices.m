function [G,W,E,H,F,Y,Cf,Cx,Cc,symmetric,bndA,bndb,Pinvset]=mpt_constructMatrices(sysStruct,probStruct,Options,setHorizon)
%MPT_CONSTRUCTMATRICES Constructs matrices for the finite time constrained optimal control problem
%
% [G,W,E,H,F,Y,Cf,Cx,Cc,symmetric,bndA,bndb,Pinvset]=
%    mpt_constructMatrices(sysStruct,probStruct,Options,setHorizon)
%
% [Matrices]=mpt_constructMatrices(sysStruct,probStruct,Options,setHorizon)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Constructs cost and constraint matrices for the finite time constrained 
% optimal control problems for linear and PWA systems as a function of the 
% prediction horizon "horizon".
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% sysStruct        - System structure in the sysStruct format
%                    Dynamics:       x(k+1)=Ax(k)+Bu(k)+f
%                                    y(k)=Cx(k)+Du(k)
%                    Constraints:    ymin  <=    y(k)     <= ymax    
%                                    umin  <=    u(k)     <= umax
%                                    dumin <= u(k)-u(k+1) <=dumax  
%                    Uncertainty:    Either polytopic ("Aunc", "Bunc") or 
%				     additive "noise". All constraints will be 
%				     enforced for the uncertain system.
%		     
%		     Consult the MPT manual for additional details
%
% probStruct       - Problem structure in the probStruct format
%                    NORM=2
%                    R,Q - Objective:    min J=x'(t+N_y|t) P x(t+N_y|t) +
%                     \sum_{k=0}^N_y-1  x'(t+k|t) Q x(t+k|t) + u'(t+k) R u(t+k)
%                    NORM=1
%                    R,Q - Objective:    min J=  ||P x(t+N_y|t)||_1 + 
%                     \sum_{k=0}^N_y-1  ||Q x(t+k|t)||_1 + ||R u(t+k)||_1
%                    NORM=Inf
%                    R,Q - Objective:    min J=  ||P x(t+N_y|t)||_Inf + 
%		             \sum_{k=0}^N_y-1   ||Q x(t+k|t)||_Inf + ||R u(t+k)||_Inf
%
%                  - horizon          
%		             Prediction horizon N; How many time steps are considered.
%
%                  - Tset
%                    Terminal set constraint: For the final state x_N the 
%		             following must hold: Tset.H * x_N <= Tset.K; If the system
%           	     is subject to uncertainty, the terminal set constraint is 
%		             automatically "robustified";
%
%                  - Qy
%                    If a weight Qy is specified, then the cost will be on the output (y(k)'Qy(k)) 
%                    and not on the state (x(k)'Qx(k)) for the optimization problem at hand 
%                    (see Problem structure above).
%
%
% Options.includeLQRset
%           If set to 1 and probStruct.Tconstraint==1, the LQR target set around
%           the origin will be computed and added to the matrices
% Options.noNoiseOnTset 
%		    If set to 1 the set constraint is not robustified with respect
%		    to additive noise. This is only relevant if robust controllers 
%		    for PWA systems are computed. For those problems the Tset is
%		    "robustified" outside of this function. Default 0;
% Options.lpsolver 
%	            Solver for LPs (see help mpt_solveLP for details)
% Options.verbose   
%  		    Level of verbosity (see help mpt_init for details)
% Options.abs_tol    
%		    absolute tolerance
%
% setHorizon	   
%		    Time step at which the terminal set constraint Tset is
%		    enforced. This is equal to the prediction horizon by default.
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
%
% NORM=2
% ------
%
% G,W,E,H,F,Y,Cf,Cx,Cc  - matrices of the problem, i.e. 
%
%       J=min_U  (0.5 U' H U + (x(0)' F + Cf) U + x(0) Y x(0) + Cx x + Cc)
%       G U <= W + E x(0)
%
%  symmetric   0/1 if constraints are symmetric
%
% Note: The elements Cf, Cx, Cc will be zero if the affine dynamics 
%	"f" are zero (x+=Ax+Bu+f).
%
%
% NORM=1 / Infinity
% ------
% G,W,E,H,F  - matrices of the problem, i.e. 
%
%       J=min_U   H U + Fx
%       G U <= W + E x(0)
%
%  symmetric   0/1 if constraints are symmetric
%
%  Note: in order to compute these problems slack variables epsilon are introduced.
%        Therefore, the optimizer U will actually be [U epsilon]. The first 
%	     (no. Inputs * prediction Horzion) elemenst correspond to the input
%            

% Copyright is with the following author(s):
%
% (C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
% (C) 2004 Pascal Grieder, Automatic Control Laboratory, ETH Zurich,
% (C) 2003 Pascal Grieder, Automatic Control Laboratory, ETH Zurich,
%          grieder@control.ee.ethz.ch

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

mpt_obsoleteFunction;

narginchk(2, 4);

global mptOptions;
if ~isstruct(mptOptions),
    mpt_error;
end
if nargin<3,
    Options=[];
end
if ~isfield(Options,'debug_level'),
    Options.debug_level=mptOptions.debug_level;
end
if ~isfield(Options,'verbose'),
    Options.verbose=mptOptions.verbose;
end
if ~isfield(Options,'abs_tol'),
    Options.abs_tol=mptOptions.abs_tol;
end
if ~isfield(Options,'lpsolver'),
    Options.lpsolver=mptOptions.lpsolver;
end
if ~isfield(Options,'includeLQRset')
    Options.includeLQRset=1;
end
if ~isfield(sysStruct,'verified'),
    sysStruct=mpt_verifySysStruct(sysStruct);
end
if ~isfield(probStruct,'verified'),
    probStruct=mpt_verifyProbStruct(probStruct);
end
if ~isfield(Options,'pwa_index'),
    Options.pwa_index=1;  % which dynamics to extract
end
if ~isfield(Options,'noConstraintReduction')
    Options.noConstraintReduction=0;
end
if ~isfield(Options,'noReduce'),
    Options.noReduce=0;
end
if ~isfield(Options, 'guierrors')
    Options.guierrors = 0;
end
if ~isfield(Options, 'infsetctr'),
    Options.infsetctr = 100;
end

ycost=1; %assume output cost
if ~isfield(probStruct,'Qy'),
    ycost=0;    %no output cost
else
    if iscell(probStruct.Qy),
        for ii=1:length(probStruct.Qy),
            if ~isempty(probStruct.Qy{ii}) & ~any(probStruct.Qy{ii})~=0
                if(probStruct.tracking & size(sysStruct.Cy,1)~=length(probStruct.Qy{ii}))
                    fprintf('\n\n')
                    disp('If you wish to punish the ouput, the dimension of the Qy-weight matrix must be adjusted!');
                    error('probStruct.Qy and output y have incompatible dimensions');
                end
            end
        end
    else
        if ~isempty(probStruct.Qy) & ~any(probStruct.Qy)~=0
            if(probStruct.tracking & size(sysStruct.Cy,1)~=length(probStruct.Qy))
                fprintf('\n\n')
                disp('If you wish to punish the ouput, the dimension of the Qy-weight matrix must be adjusted!');
                error('probStruct.Qy and output y have incompatible dimensions');
            end
        end
    end
end

if ~isa(probStruct.Tset,'polytope'),
    fprintf('\n\n')
    error('mpt_construcMatrices: terminal set MUST be a polytope object!');
end

if(~isfulldim(probStruct.Tset))
    if Options.verbose>1,
        disp('mpt_constructMatrices: Warning, Pfinal either undefined or not fully dimensional!');
    end
    Pfinal=polytope;
else
    Pfinal=probStruct.Tset;
end

horizon = probStruct.N;

if(nargin<4 | isempty(setHorizon))
    setHorizon=horizon; %add set constraint as terminal set constraint
end

if isfield(sysStruct, 'noise') && isfulldim(sysStruct.noise)
    noise=1;
    if(isfield(Options,'noNoiseOnTset') & Options.noNoiseOnTset)
        noiseOnTset=0;  %do not take minkowski difference on target set
    else
        noiseOnTset=1;
    end
else 
    noiseOnTset=0;
    noise=0;
end
if(isfield(sysStruct,'Aunc') & ~isempty(sysStruct.Aunc))
    unc_dyn=1;
    max_uncert=length(sysStruct.Aunc);
else
    unc_dyn=0;
    max_uncert=1;
end

constraints_x0=probStruct.y0bounds;       %state constraints on the state x0
constraints_xN=1;                         %state constraints on the state xN

%this automatically removes redundant constraints in the formulation
%this needs to be switched off for the infinite horizon algorithm for LTI systems
if (sysStruct.type==0 | strcmp(sysStruct.type,'LTI')) & probStruct.N==inf & ~isfield(sysStruct,'guardC'),
    reduce_constraints=0;   
else
    reduce_constraints=1;
end
if Options.noConstraintReduction,
    reduce_constraints=0;
end
lpsolver=Options.lpsolver;
DEBUG=Options.debug_level;

if strcmp(sysStruct.type,'PWA') | sysStruct.type==1 | iscell(sysStruct.A),
    if(iscell(sysStruct.A))
        sysStruct.A=sysStruct.A{Options.pwa_index};
        sysStruct.B=sysStruct.B{Options.pwa_index};
        sysStruct.C=sysStruct.C{Options.pwa_index};
        sysStruct.D=sysStruct.D{Options.pwa_index};
        sysStruct.guardX=sysStruct.guardX{Options.pwa_index};
        sysStruct.guardU=sysStruct.guardU{Options.pwa_index};
        sysStruct.guardC=sysStruct.guardC{Options.pwa_index};
        sysStruct.f=sysStruct.f{Options.pwa_index};
        sysStruct.g=sysStruct.g{Options.pwa_index};
    else
        fprintf('\n\n')
        error('System is defined as PWA but contains only one dynamic')
    end
end

[A,B,C,D,Q,R,ymin,ymax,umin,umax,dumin,dumax,bndA,bndb]=mpt_evalSystem(sysStruct,probStruct);
if(ycost)
    Q=probStruct.Qy;
    probStruct.Q = Q;
end

if(size(umax,2)>size(umax,1) | size(umin,2)>size(umin,1))
    fprintf('\n\n')
    error('mpt_contructMatrices: u constraints must be column vectors (transpose current entry)')
end
if(size(dumax,2)>size(dumax,1) | size(dumin,2)>size(dumin,1))
    fprintf('\n\n')
    error('mpt_contructMatrices: du constraints must be column vectors (transpose current entry)')
end
if(size(ymax,2)>size(ymax,1) | size(ymin,2)>size(ymin,1))
    fprintf('\n\n')
    error('mpt_contructMatrices: y constraints must be column vectors (transpose current entry)')
end

if Options.includeLQRset & (~isfulldim(probStruct.Tset) & probStruct.Tconstraint==1) &...
        probStruct.norm==2 & probStruct.tracking==0 & ycost==0,
    %compute invariant set for the LQR controller

    tolerance=1e-6;
    %COMPUTE TARGET SET / TERMINAL SET
    %get LQR feedback law
    if iscell(probStruct.Q),
        QQ = probStruct.Q{end};
    else
        QQ = probStruct.Q;
    end
    if iscell(probStruct.R),
        RR = probStruct.R{end};
    else
        RR = probStruct.R;
    end
    [FLQR,PLQR] = dlqr(sysStruct.A,sysStruct.B,QQ,RR);
    FLQR=-FLQR;
    GLQR=zeros(length(sysStruct.umax),1);
    
    % Fx >= u_min
    % Fx <= u_max
    % y  >= y_min
    % y  <= y_max
    % du >= du_min
    % du <= du_max
    totA = [bndA; FLQR; -FLQR; sysStruct.C; -sysStruct.C; ...
            FLQR-FLQR*(sysStruct.A+sysStruct.B*FLQR); -(FLQR-FLQR*(sysStruct.A+sysStruct.B*FLQR))]; 
    totb = [bndb*1e2; sysStruct.umax; -sysStruct.umin; sysStruct.ymax; -sysStruct.ymin; sysStruct.dumax; -sysStruct.dumin];
    
    if isfield(sysStruct,'dymax'),
        totA = [totA; sysStruct.C - sysStruct.C*(sysStruct.A + sysStruct.B*FLQR); ...
            -(sysStruct.C - sysStruct.C*(sysStruct.A + sysStruct.B*FLQR))];
        totb = [totb; sysStruct.dymax; -sysStruct.dymin];
    end
    if isfield(sysStruct,'xmin'),
        nx = length(sysStruct.xmin);
        totA = [totA; eye(nx); -eye(nx)];
        totb = [totb; sysStruct.xmax; -sysStruct.xmin];
    end
    
    X = polytope(totA, totb);
    if(~isfulldim(X))
        fprintf('\n\n')
        error('mpt_constructMatrices: No invariant target set found for LQR controller! Try setting Tconstraint=0 instead (no stability guarantee).')
    end
    if(~isfield(sysStruct,'Aunc') | isempty(sysStruct.Aunc))
        %no polytopic uncertainty
        A_CL{1}=sysStruct.A+sysStruct.B*FLQR;              %closed loop dynamics
    else
        %polytopic uncertainty
        for i=1:length(sysStruct.Aunc)
            A_CL{i}=sysStruct.Aunc{i}+sysStruct.Bunc{i}*FLQR; %closed loop dynamics
        end
    end
    
    [PinvSet,tstar,fd] = mpt_infset(A_CL,X,Options.infsetctr,sysStruct.noise,Options);
    
    Pfinal = PinvSet & sysStruct.Pbnd;
    
    probStruct.Tconstraint = 2;    % we just computed LQR invariant set, hence we enforce it as a taget set
    if(~isfulldim(Pfinal))
        fprintf('\n\n')
        error('No full dimensional invariant target set found, Stability cannot be proven. Try setting Tconstraint=0 instead (no stability guarantee). ')
    end
    probStruct.Tset = Pfinal;
end

nu=size(B,2);
nx=size(A,1);
ny=size(C,1);
nuH=nu*horizon;         %degrees of freedom = horizon  * number of inputs

if(~isfield(sysStruct,'f'))
    if isfield(probStruct,'xref'),
        f = sysStruct.A*probStruct.xref + sysStruct.B*probStruct.uref - probStruct.xref;
    else
        f=zeros(nx,1);
    end
else
    f=sysStruct.f;
end

% xn = x - xref  ->  x = xn + xref
% un = u - uref  ->  u = un + uref
% y = Cx + Du + g
% yref = C xref + D uref + g
% y-yref = Cx + Du + g - yref
% yn = C(xn+xref) + D(un+uref) + (g - yref)
% yn = C xn + D un + (g - yref + C xref + D uref)

if ~isfield(sysStruct,'g'),
    gg = zeros(ny,1);
else
    gg = sysStruct.g;
end


% y <= ymax
% y - yref <= ymax - yref
% yn <= ymax - C xref - D uref - g

% y >= ymin
% y - yref >= ymin - yref
% yn >= ymin - yref
% yn >= ymin - C xref - D uref - g 

% y = Cx + Du + g
% yn = C xn + D un + g
% yn = C (x - xref) + D (u - uref) + g
% yn = C x + D u + (g - C xref - D uref)

if isfield(probStruct,'xref'),
    ymax = ymax - C*probStruct.xref - D*probStruct.uref - gg;
    ymin = ymin - C*probStruct.xref - D*probStruct.uref - gg;
    if isfield(sysStruct, 'xmin')
        sysStruct.xmin = sysStruct.xmin - probStruct.xref;
        sysStruct.xmax = sysStruct.xmax - probStruct.xref;
    end
end


if isfield(probStruct,'P_N'),
    P=probStruct.P_N;
    if(~ycost & any(size(P,2)~=nx))
        if Options.guierrors,
            error('Dimension of penalty on final state must be equal to dimension of A matrix!');
        else
            fprintf('\n\n')
            error('mpt_contructMatrices: Dimension of probStruct.P_N must be equal to dimension of A matrix!');
        end
    end
    if(ycost & length(probStruct.P_N)~=ny)
        if Options.guierrors,
            error('Penalty on final state must be a square matrix of the dimension of the output!')
        else
            fprintf('\n\n')
            error('The terminal weight probStruct.P_N must be a square matrix of the dimension of the output y!')
        end
    end
elseif probStruct.norm==2 & (probStruct.tracking==0 | isfield(probStruct,'xref')),
    if ycost,
        P=Q;
    else
        [KLQ,P]=dlqr(A,B,Q,R); % Solution of Riccati equation
    end
elseif probStruct.norm==2 & probStruct.tracking>0 & ~isfield(probStruct,'xref') & ~ycost,
    if ycost,
        nyd = sysStruct.dims.ny;
        nx = sysStruct.dims.nx;
        nu = sysStruct.dims.nu;
        [KLQ,P]=dlqr(A(1:nx,1:nx),B(1:nx,1:nu),Q(1:nx,1:nx),R); % Solution of Riccati equation
        P = Q(1:nx,1:nx);
        %P= [P zeros(nx,nu) -P;zeros(nu,nx) R zeros(nu,nx); zeros(nyd,2*nx+nu)]; %terminal cost
        %KLQ=[-KLQ zeros(nu,nx+nu)];
        nu=size(B,2);
        nx=size(A,1);
    else
        nu=size(B,2);
        if probStruct.tracking==1,
            nx=(size(A,2)-nu)/2;
        else
            nx = size(A,2)/2;
        end
        % if you get an error on this line, most probably you already augmented
        % system and problem matrices to deal with tracking. In that case, use
        % Options.autoTracking = 0 to switch off automatic tracking augmentation.
        [KLQ,P]=dlqr(A(1:nx,1:nx),B(1:nx,1:nu),Q(1:nx,1:nx),R); % Solution of Riccati equation
        if probStruct.tracking==1,
            P= [P zeros(nx,nu) -P;zeros(nu,nx) R zeros(nu,nx); -P zeros(nx,nu) P]; %terminal cost
        else
            P = [P -P; -P P];
        end
            
        KLQ=[-KLQ zeros(nu,nx+nu)];
        nu=size(B,2);
        nx=size(A,1);
    end
elseif probStruct.tracking>0 & ~isfield(probStruct,'xref') & ~ycost,
    if ycost,
        nyd = sysStruct.dims.ny;
        nx = sysStruct.dims.nx;
        nu = sysStruct.dims.nu;
        P=probStruct.Q(1:nx,1:nx);
        %P= [P zeros(nx,nu) -P;zeros(nu,nx) R zeros(nu,nx); zeros(nyd,2*nx+nu)]; %terminal cost
        nu=size(B,2);
        nx=size(A,1);
    else
        nu=size(B,2);
        if probStruct.tracking==1,
            nx=(size(A,2)-nu)/2;
        else
            nx = size(A,2)/2;
        end
        
        % if you get an error on this line, most probably you already augmented
        % system and problem matrices to deal with tracking. In that case, use
        % Options.autoTracking = 0 to switch off automatic tracking augmentation.
        %disp('tracking for 1/Inf norm');
        P=probStruct.Q(1:nx,1:nx);
        
        if probStruct.tracking==1,
            P= [P zeros(nx,nu) -P;zeros(nu,nx) R zeros(nu,nx); -P zeros(nx,nu) P]; %terminal cost
        else
            if probStruct.norm==2,
                P = [P -P; -P P];
            else
                P = [P -P];
            end
        end
        nu=size(B,2);
        nx=size(A,1);
    end
else
    P = Q;
end

% if ycost 
%     if(probStruct.tracking)
%         P = Q;
%     else
%         P=C'*P*C;
%     end
% end

if(isfield(probStruct,'feedback') & probStruct.feedback==1)
    %system is being pre-stabilized by linear feedback FB (i.e., x^+=(A+B*FB)x+Bu+f);
    %of course, constraints are imposed on the "new input" u_new=FB*x+Bu
    if isfield(probStruct,'FBgain')
        % if the feedback matrix is provided, use it
        FB=probStruct.FBgain;%update dynamic matrices in case pre-stabilization is used
        %A=A+B*FB;
    else
        % if feedback gain matrix is not give, use LQR
        [FB,S,E] = dlqr(sysStruct.A,sysStruct.B,probStruct.Q,probStruct.R);
        [FB,S,E] = dlqr(A,B,Q,R);
        FB=-FB;
    end
    A=A+B*FB; %update dynamics
    C=C+D*FB;
else
    FB=zeros(nu,nx);
end
probStruct.FBgain = FB;

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% BUILDING THE CONSTRAINT MATRICES:
%  The original constraints      g u(k) <= w + e x(k), forall k=0...N-1
%  are transformed into          G U <= W + E x(0),  where U=[u_0 ... u_{N-1}]
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% ymin  <=    y(k)     <= ymax     (y(k)=Cx(k)+Du(k) ; x(k+1)=Ax(k)+Bu(k))
%% umin  <=    u(k)     <= umax
%% dumin <= u(k)-u(k+1) <=dumax      


%---------------------------------------------------------------------
%%first transform min-max constraints into form  g u(k) <= w + e x(k)
g=[];
w=[];
e=[];

% G U <= W + E x

% y = C x + D u + g
%  y <=  ymax   ->   D u <= (ymax-g) - C x
% -y <= -ymin   ->  -D u <= -ymin+g  + C x

if isfield(sysStruct, 'xmin')
    % x(k+1) = A x + B u + f
    % B u = -A x + f
    %
    gx = [0*B; 0*B];
    wx = [sysStruct.xmax; -sysStruct.xmin];
    ex = [-eye(nx); eye(nx)];
    
    g = [g; gx];
    w = [w; wx];
    e = [e; ex];

    constraints_x0 = 1;
end

g=[g; D;-D];               %output constraints
w=[w;ymax-gg;-ymin+gg];    %output constraints
e=[e;-C;C];      %output constraints

g=[g; eye(nu);-eye(nu)];   %input constraints u=FB*x+u_c
w=[w; umax;-umin];         %input constraints
e=[e; -FB;FB];             %input constraints

if(constraints_x0) %Are there constraints on the initial state?
    %add state and input constraints for the first step
    G=[g zeros(size(g,1),nuH-size(g,2))];  %augment g matrix to full horizon
    W=w;
    E=e;       
else
    %add input constraints only for the first step
    Gt=[eye(nu);-eye(nu)];
    G=[Gt zeros(2*nu,nuH-nu)];
    W=[umax;-umin];
    E=[-FB;FB];  
end

Anom=A;
Bnom=B;
for unc_ctr=1:max_uncert
    if(unc_dyn)
        A=sysStruct.Aunc{unc_ctr};
        B=sysStruct.Bunc{unc_ctr};
    end

    %---------------------------------------------------------------------
    %%Now construct the full constraint matrices by iterating over the horizon
    Ai=eye(nx);
    fsum=zeros(nx,1);
    GU=zeros(size(g,1),nuH);     %vector which will store [e*[A^{i-1}B ... AB B] g] in the following iteration
    vecX=zeros(nx,nuH);          %vector which will store [A^{i-1}B ... AB B] in the following iteration
    noiseVal=zeros(length(w),1);
    for i=1:horizon
        
        % patch for systems with parametric uncertainty
        % we need to consider all combinations of dynamic matrices, i.e.
        % A1*A1, A1*A2, A2*A2, ...
        
        if unc_dyn
            Ai_stack = sub_allAuncCombs(sysStruct.Aunc, i-1, nx);
        else
            Ai_stack = {};
            Ai_stack{1} = Ai;
        end

        for iunc = 1:length(Ai_stack),
            Ai = Ai_stack{iunc};
            if unc_dyn & length(Ai_stack)>1 & iunc>1,
                % restore GU which is modified if (constraints_xN & all(D==0))
                GU = GUstore;
            end
            
            dymaxadded=0;

            fsum=fsum+Ai*f;

            %Set: x(i)=AAi*x(0)+vecX*U
            AAi=A*Ai;               %=A^i used for faster computation
            vecX=[Ai*B vecX];       %store [A^{i-1}B ... AB B]
            vecX(:,nuH+1:end)=[];   %store [A^{i-1}B ... AB B]

            %Set: u(i)=vecU*U
            vecU=[zeros(nu,(i-1)*nu) eye(nu) zeros(nu,(horizon-i)*nu)];


            GU=[-e*Ai*B GU(:,1:(i-1)*nu) g zeros(size(g,1),nuH-(i+1)*nu)];       %store [A^{i-1}B ... AB B]
            GU(:,(nuH+1):end)=[];                                                %remove elements that were shifted too far
            GUstore = GU;

            if(noise)
                if isfield(sysStruct, 'xmin'),
                    tmp=subtractNoise(A^(i-1),e(1:2*nx+2*ny,:),sysStruct.noise,Options); %worst case noise on state
                    noiseVal=noiseVal+[tmp;zeros(length(w)-2*nx-2*size(C,1),1)];         %no noise on input signal
                else
                    tmp=subtractNoise(A^(i-1),e(1:2*ny,:),sysStruct.noise,Options); %worst case noise on state
                    noiseVal=noiseVal+[tmp;zeros(length(w)-2*size(C,1),1)];         %no noise on input signal
                end

            end
            
            if(i<horizon)
                %FBtemp=[zeros(size(g,1)-2*nu,nuH); FB*vecX;-FB*vecX];
                G=[G;GU];            %Add input and state constraints for step i
                W=[W;w+e*fsum+noiseVal];    %Add input and state constraints for step i
                %FBtemp=[zeros(size(g,1)-2*nu,nx); -FB*AAi;+FB*AAi];
                E=[E;e*AAi];          %Add input and state constraints for step i

                %ADD SLEW CONSTRAINTS
                dMat=[zeros(nu,(i-1)*nu) eye(nu) -eye(nu) zeros(nu,nuH-(i+1)*nu)];  %dMat*U= u(k-1)-u(k)
                vecXt=[vecX(:,(nu+1):end) zeros(nx,nu)]; %Shift vecX such that:   x(i-1)=Ai*x(0)+vecXt*U
                %x(k-1)-x(k)= (Ai*x(0)+vecXt*U) - (AAi*x(0)+vecX*U)
                %y(k-1)-y(k)= C (x(k-1)-x(k)) + D(u(k-1)-u(k))

                G=[G;dMat+FB*(vecXt-vecX);-dMat-FB*(vecXt-vecX)];           %slew constraints;  u(k-1)-u(k)<=umax, etc.
                W=[W;dumax;-dumin];         %slew constraints
                E=[E;FB*(AAi-Ai);-FB*(AAi-Ai)];       %slew constraints

                if isfield(sysStruct,'dymax'),
                    dMat=[zeros(nu,(i-1)*nu) eye(nu) -eye(nu) zeros(nu,nuH-(i+1)*nu)];  %dMat*U= u(k-1)-u(k)
                    vecXt=[vecX(:,(nu+1):end) zeros(nx,nu)]; %Shift vecX such that:   x(i-1)=Ai*x(0)+vecXt*U
                    G = [G; C*(vecXt - vecX) + D*dMat + D*FB*(vecXt-vecX);...
                        -C*(vecXt-vecX) - D*dMat - D*FB*(vecXt-vecX)];
                    W = [W; -sysStruct.dymin; sysStruct.dymax];
                    E = [E; C*(AAi-Ai); -C*(AAi-Ai)];
                    dymaxadded=1;
                end

            elseif(constraints_xN & all(all(D==0)))
                %extract the state constraints for the final state x_N (there is no input u_N)
                GU=GU(1:(end-2*nu),:);          %Add state constraints for step horizon (=xN)
                G=[G;GU];                       %Add state constraints for step horizon (=xN)
                if 1 & isfield(sysStruct, 'xmax'),
                    %Add state constraints for step horizon (=xN)
                    W=[W; sysStruct.xmax-fsum+noiseVal(1:nx); -sysStruct.xmin+fsum+noiseVal(nx+1:2*nx)];
                    W=[W; ymax-C*fsum+noiseVal(2*nx+1:2*nx+ny);-ymin+C*fsum+noiseVal(2*nx+ny+1:2*nx+2*ny)];
                else
                    %Add state constraints for step horizon (=xN)
                    W=[W;ymax-C*fsum+noiseVal(1:ny);-ymin+C*fsum+noiseVal(ny+1:2*ny)];     
                end
                tmp=e(1:(end-2*nu),:);
                E=[E;tmp*AAi];                    %Add state constraints for step horizon (=xN)
            end
            if(i==setHorizon & isfulldim(Pfinal)) %Add extra set constraint
                [Hfinal, Kfinal]=double(Pfinal);
                G=[G;Hfinal*vecX(1:nx,:)];
                if(noiseOnTset)
                    noiseValTset=zeros(length(Kfinal),1);
                    for tt=1:i
                        noiseValTset=noiseValTset+subtractNoise(A^(tt-1),Hfinal,sysStruct.noise,Options);      %worst case noise on state
                    end
                    W=[W;Kfinal-Hfinal*fsum+noiseValTset];
                else
                    W=[W;Kfinal-Hfinal*fsum];
                end
                E=[E;-Hfinal*AAi];
            end

            if 0 & isfield(sysStruct, 'xmax')
                % de-activated, state constraints are already included in 'g',
                % 'w' and 'e'
                I2nx = [eye(nx); -eye(nx)];
                G = [G; I2nx*vecX(1:nx,:)];
                W = [W; wx - I2nx*fsum];
                E = [E; -I2nx*AAi];
            end
            
            if ~dymaxadded & isfield(sysStruct,'dymax'),
                %dMat=[zeros(nu,(i-1)*nu) eye(nu) -eye(nu) zeros(nu,nuH-(i+1)*nu)];  %dMat*U= u(k-1)-u(k)
                vecXt=[vecX(:,(nu+1):end) zeros(nx,nu)]; %Shift vecX such that:   x(i-1)=Ai*x(0)+vecXt*U
                G = [G; C*(vecXt - vecX) + D*FB*(vecXt-vecX);...
                    -C*(vecXt-vecX) - D*FB*(vecXt-vecX)];
                W = [W; -sysStruct.dymin; sysStruct.dymax];
                E = [E; C*(AAi-Ai); -C*(AAi-Ai)];
            end
      
            Ai=A*Ai;                %=A^{i-1} (in next iteration) used for faster computation
        end
    end
end%unc_ctr


if(isfield(sysStruct,'guardX'))
    %add constraint where dynamics are defined: GU <= W + Ex
    G=[G;sysStruct.guardU zeros(size(sysStruct.guardU,1),nuH-nu)];
    W=[W;sysStruct.guardC];
    E=[E;-sysStruct.guardX-sysStruct.guardU*FB];
end

% if isfield(sysStruct, 'xmin')
%     % x(k+1) = A x + B u + f
%     % B u = -A x + f
%     % 
%     G = [G; B zeros(nx,nuH-nu); -B zeros(nx,nuH-nu)];
%     W = [W; sysStruct.xmax; -sysStruct.xmin];
%     E = [E; -eye(nx); eye(nx)];
% end


if(~ycost)
    %no cost on output; adjust matrices
    D=zeros(nx,nu);
    D=D*0;
    C=eye(nx);
    gg = zeros(nx,1);
elseif(probStruct.tracking)
    D=sysStruct.Dy;
    C=sysStruct.Cy;
    if iscell(D),
        D = D{Options.pwa_index};
    end
    if iscell(C),
        C = C{Options.pwa_index};
    end
end

if(ycost & probStruct.tracking==1)
    gg=gg(1:size(C,1));
end

if ycost %& isfield(probStruct,'yref'),
    if isfield(sysStruct,'Cy'),
        C = sysStruct.Cy;
        D = sysStruct.Dy;
        if iscell(D),
            D = D{Options.pwa_index};
        end
        if iscell(C),
            C = C{Options.pwa_index};
        end
        gg = gg(1:size(C,1));
        if size(Q,1)>size(C,1),
            Q = Q(1:size(C,1),1:size(C,1));
        end
        if size(P,1)>size(C,1),
            if probStruct.norm==2,
                P = P(1:size(C,1),1:size(C,1));
            else
                P = P(:, 1:size(C,1));
            end
        end
    end
    if isfield(probStruct, 'yref')
        gg = gg - probStruct.yref;
    end
end

A = Anom;
B = Bnom;
Porig = P;

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% BUILDING THE COST FUNCTION:
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if(probStruct.norm==2)
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    %  The original cost function    J= (sum_{i=0}^N x_i' Q x_i + u_i' R u_i) + x_N' P x_N,
    %  is transformed into           J= U' H U + x' F U + CfU + x' Y x,  where U=[u_0 ... u_{N-1}]
    %                                x(k+1)=Ax(k)+Bu(k)+f
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    vecX=zeros(nx,nuH);     %vector which will store [A^{i-1}B ... AB B] in the following iteration
                            %J=x'Yx + x'F U + U' H U 
    H=zeros(nuH,nuH);       %this is the H matrix
    F=zeros(nx,nuH);        %this is the F matrix
    Y=zeros(nx,nx);         %this is the Y matrix
    Cf=zeros(1,nuH);
    Cx=zeros(1,nx);
    Cc=0;
    
    Ai=eye(nx);             %used for faster computation
    fA=zeros(1,nx);
    Asum=zeros(nx,nx);
    
    for i=1:horizon         
        Asum=Asum+Ai;           
        %Set: x(i)=AAi*x(i)+vecX*U+Asum*f
        %     y(i)=C  *x(i)+D*vecX*U+gg 
        AAi=A*Ai;               %=A^i used for faster computation
        vecX=[Ai*B vecX];       %store [A^{i-1}B ... AB B] 
        vecX(:,nuH+1:end)=[];   %store [A^{i-1}B ... AB B] 
   
        vecU=[zeros(nu,(i-1)*nu) eye(nu) zeros(nu,(horizon-i)*nu)];   %u(i)=vecU*U
        
        if iscell(probStruct.R),%in case of time varying cost matrices
            R = probStruct.R{i};
        end 
        if iscell(probStruct.Q),%in case of time varying cost matrices
            Q = probStruct.Q{i};
        end
        if iscell(Porig),
            P = Porig{i};
        else
            P = Porig;
        end
        
        %Sum J=J+u(i)'Ru(i) below  
        %   For Pre-Stabilization: Shift Xvec such that: x(i-1) = Ai*x(0)+vecXt*U
        %   Therefore: U_real(i-1) = (vecU+FB*vecXt)*U + FB*Ai * x(0)
        vecXt=[vecX(:,(nu+1):end) zeros(nx,nu)];  
        H=H+(vecU+FB*vecXt)'*R*(vecU+FB*vecXt);
        F=F+2*(FB*Ai)'*R*(vecU+FB*vecXt);
        Y=Y+Ai'*FB'*R*FB*Ai;
        
        if(i==1)
            %Since the state vector X=[x_0 ... x_N] vector contains
            %one element more than the U=[u_0 ... u_{N-1}] vector,
            %we add one more element to the matrices here...
            H=H+(vecU'*D')*Q*(D*vecU);                      %build H
            F=F+2*(Ai'*C'*Q*(D*vecU));                      %build F
            Cf=Cf+2*(gg')*Q*(D*vecU);                       %build Cf*U
            Cx=Cx+2*(f'*Asum'*C'+gg')*Q*C*Ai;               %build Cx*x
            Cc=Cc+(f'*Asum'*C'+gg')*Q*(C*Asum*f+gg);        %build Cc
            Y=Y+Ai'*C'*Q*C*Ai; 
        end
        
        if(i<horizon)%cost function J= U' H U + x' F U + FfU + x' Y x
            %     x(i)=AAi*x(0)+vecX*U+Asum*f and  y(i)=C*x(i)+D*vecX*U+gg 
            %     Therefore:
            %     y(i)=(C*AAi)*x(0) + (C*vecX+D*vecX)*U + (C*Asum*f+gg)
            
            %vecU is only used for direct feedthrough (D not 0), hence
            %the time index must be shifted by 1
            vecU=[zeros(nu,i*nu) eye(nu) zeros(nu,(horizon-i-1)*nu)];   %u(i)=vecU*U
            
            %Sum J=J+y(i)'Qy(i) below
            H=H+(vecX'*C'+vecU'*D')*Q*(C*vecX+D*vecU);  %build H
            F=F+2*(AAi'*C'*Q*(C*vecX+D*vecU));           %build F 
            Y=Y+AAi'*C'*Q*C*AAi;                              %build Y
            Cf=Cf+2*(f'*Asum'*C'+ gg')*Q*(C*vecX+D*vecU);   %build Cf*U
            Cx=Cx+2*(f'*Asum'*C'+gg')*Q*C*AAi;              %build Cx*x
            Cc=Cc+(f'*Asum'*C'+gg')*Q*(C*Asum*f+gg);        %build Cc
        elseif(all(all(D==0)) | all(all(P==0)))
            %     x(i)=AAi*x(0)+vecX*U+Asum*f and  y(i)=C*x(i)+D*vecX*U+gg 
            %     Therefore:
            %     y(i)=(C*AAi)*x(0) + (C*vecX+D*vecX)*U + (C*Asum*f+gg)
            
            %Sum J=J+y(i)'Py(i) below
            H=H+(vecX'*C'+vecU'*D')*P*(C*vecX+D*vecU);   %build H    
            F=F+2*AAi'*C'*P*(C*vecX+D*vecU);             %build F 
            Y=Y+AAi'*C'*P*C*AAi;                            %build Y
            Cf=Cf+2*(f'*Asum'*C'+gg')*P*C*vecX;             %build Cf*U
            Cx=Cx+2*(f'*Asum'*C'+gg')*P*C*AAi;              %build Cx*x
            Cc=Cc+(f'*Asum'*C'+gg')*P*(C*Asum*f+gg);        %build Cc
        else
            fprintf('\n\n')
            disp('PROBLEM STATEMENT MAKES NO SENSE')
            error('Terminal weight ''probStruct.P_N'' must be zero if ''sysStruct.D'' is non-zero!')
        end
        
        Ai=A*Ai;                %=A^{i-1} (in next iteration) used for faster computation
    end

    %eigH=max(abs(eig(H)));        %scale the matrices
    eigH=1;
    H=2*H/eigH;                     %scale the matrices (multiply with 2 because H is automatically multiplied with 0.5 by quadprog)
    F=F/eigH;                     %scale the matrices
    Y=Y/eigH;                     %scale the matrices
    Cf=Cf/eigH;                   %scale the matrices
    Cx=Cx/eigH;                   %scale the matrices
    Cc=Cc/eigH;                   %scale the matrices
    
elseif(probStruct.norm==1)
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    %  The original cost function    J= (sum_{i=0}^N  |Q x_i|_1 +  |R u_i|_1) +  |P x_N|_1,
    %  is transformed into           J=  H [U epsilon] + F x0                          
    %                                
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if iscell(Q),
        xcost = size(Q{1}, 1);
    else
        xcost=size(Q,1);
    end
    if iscell(R),
        ucost = size(R{1}, 1);
    else
        ucost=size(R,1);
    end
    if(~ycost)
        xNcost=size(P,1);
    else
        %no terminal weight on output possible, because "terminal input" is unknown
        xNcost=0;
    end
        
    Ai=eye(nx);             %used for faster computation
    
    %    G=[G zeros(size(G,1),(xcost+ucost)*horizon)]; %extend state space to introduce slack variables
     G=[G zeros(size(G,1),(xcost+ucost)*horizon+xNcost)]; %extend state space to introduce slack variables (add cost for x0 as well)
    
    nuH=nu*horizon;         %degrees of freedom = horizon  * number of inputs
    vecX=zeros(nx,nuH+(xcost+ucost)*horizon+xNcost);     %vector which will store [A^{i-1}B ... AB B zeros(nx,length(epsilon))]
                                                         %in the following iteration
    Asum=zeros(nx,nx);

    
    %eps_xN(i)=vecExN*U (eps_xN(i)=>P*x(i) & eps_xN(i)=>-P*x(i))
    vecExN=[zeros(xNcost,horizon*xcost+nuH) eye(xNcost) zeros(xNcost,horizon*ucost)]; %eps_xN(i)=vecExN*U
    
    for i=1:horizon         
        if iscell(probStruct.R),
            R = probStruct.R{i};
        end 
        if iscell(probStruct.Q),
            Q = probStruct.Q{i};
        end
        AAi=A*Ai;               %=A^i used for faster computation
        Asum=Asum+Ai;
        %x(i)=AAi*x+vecX*U+Asum*f
        %u(i)=vecU*U
        %eps_x(i)=vecEx*U   (eps_x(i)=>Q*x(i) & eps_x(i)=>-Q*x(i))
        %eps_u(i)=vecEu*U   (eps_u(i)=>R*u(i) & eps_u(i)=>-R*u(i))
        vecX=[Ai*B vecX];       %store [A^{i-1}B ... AB B [bunch of zeros till its full length]] 
        vecX(:,(nuH+(xcost+ucost)*horizon+xNcost+1):end)=[];           %store [A^{i-1}B ... AB B] 
        vecU=[zeros(nu,(i-1)*nu) eye(nu) zeros(nu,(horizon-i)*nu+(xcost+ucost)*horizon+xNcost)];   %u(i)=vecU*U
        %         vecEx=[zeros(xcost,(i-1)*xcost+nuH) eye(xcost) zeros(xcost,(horizon-i)*xcost) zeros(xcost,ucost*horizon)]; %eps_x(i)=vecEx*U
        %         vecEu=[zeros(ucost,horizon*xcost+nuH) zeros(ucost,ucost*(i-1)) eye(ucost) zeros(ucost,ucost*(horizon-i))]; %eps_u(i)=vecEu*U
        vecEx=[zeros(xcost,i*xcost+nuH) eye(xcost) zeros(xcost,(horizon-i-1)*xcost+xNcost) zeros(xcost,ucost*horizon)]; %eps_x(i)=vecEx*U
        vecEu=[zeros(ucost,horizon*xcost+xNcost+nuH) zeros(ucost,ucost*(i-1)) eye(ucost) zeros(ucost,ucost*(horizon-i))]; %eps_u(i)=vecEu*U
        
        if(i==1)
            %add cost for state x0/y0
            vecEx0=[zeros(xcost,nuH) eye(xcost) zeros(xcost,(horizon-1)*xcost+xNcost) zeros(xcost,ucost*horizon)];
            Matrices.vecEx0 = vecEx0;
            G=[G;-vecEx0+Q*(D*vecU)];
            W=[W;-Q*gg];
            E=[E;-Q*C];
            G=[G;-vecEx0-Q*(D*vecU)];
            W=[W;Q*gg];
            E=[E;Q*C];
        end
        if(i<horizon)
            if(ycost)
                %shift current input by one for output tracking
                vecU=[zeros(nu,i*nu) eye(nu) zeros(nu,(horizon-i-1)*nu+(xcost+ucost)*horizon+xNcost)];   %u(i)=vecU*U
            end
            %Q*(A^i*x+[A^{i-1}B ... AB B]*x)<=eps1
            %G=[G;Q*vecX zeros(xcost,(i-1)*xcost) -eye(xcost) zeros(xcost,(horizon-i)*xcost) zeros(xcost,ucost*horizon)];
            G=[G;Q*C*vecX-vecEx+Q*(D*vecU)];
            W=[W;-Q*gg-Q*C*Asum*f];
            E=[E;-Q*C*AAi];
            %-eps1<=Q*(A^i*x+[A^{i-1}B ... AB B]*x)
            G=[G;-Q*C*vecX-vecEx-Q*(D*vecU)];
            W=[W;Q*gg+Q*C*Asum*f];
            E=[E;Q*C*AAi];
            if(ycost)
                %shift current input back
                 vecU=[zeros(nu,(i-1)*nu) eye(nu) zeros(nu,(horizon-i)*nu+(xcost+ucost)*horizon+xNcost)];   %u(i)=vecU*U
            end
        elseif ~ycost
            %P*(A^i*x+[A^{i-1}B ... AB B]*x)<=eps1
            G=[G;P*C*vecX-vecExN+P*D*vecU];
            W=[W;-P*gg-P*C*Asum*f];
            E=[E;-P*C*AAi];
            %-eps1<=P*(A^i*x+[A^{i-1}B ... AB B]*x)
            G=[G;-P*C*vecX-vecExN-P*D*vecU];
            W=[W;P*gg+P*C*Asum*f];
            E=[E;P*C*AAi];
        end
        %R*u<=eps2
        G=[G;R*vecU-vecEu];
        W=[W;zeros(ucost,1)];
        E=[E;zeros(ucost,size(E,2))];
        %-eps<=R*u2
        G=[G;-R*vecU-vecEu];
        W=[W;zeros(ucost,1)];
        E=[E;zeros(ucost,size(E,2))];
        
        Ai=A*Ai;                %=A^{i-1} (in next iteration) used for faster computation
    end
    
    %New cost function (cost is actually defined by the new constraints)
    H=[zeros(1,nuH) ones(1,(xcost+ucost)*horizon+xNcost)];       %this is the H matrix
    F=zeros(1,nx);  %for linear objectives there is no F part, cost is expressed through slacks
    Y=[];
    Cx=[];
    Cf=[];
    Cc=[];
elseif(probStruct.norm==Inf)
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    %  The original cost function    J= (sum_{i=0}^N  |Q x_i|_inf +  |R u_i|_inf +  |P x_N|_inf,
    %  is transformed into           J=  H [U epsilon] + F x0                          
    %                                
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    xcost=size(Q,1);
    ucost=size(R,1);
    if(~ycost)
        xNcost=size(P,1);   %only one slack needed for inf norm
        yTer=0;             %binary variable which is used to crop the slack variable vector. 
                            %read additional info below...
    else
        %no terminal weight on output possible, because "terminal input" is unknown
        yTer=1;             %binary variable; if output is being punished, there is no terminal weight
                            %i.e. there is one slack variable less
        xNcost=0;           %dimension of terminal state weight
    end
    
    Ai=eye(nx);             %used for faster computation
    
    G=[G zeros(size(G,1),2*horizon+~yTer)]; %extend state space to introduce slack variables
                                      %one slack for state and one for inputs
    
    nuH=nu*horizon;         %degrees of freedom = horizon  * number of inputs
    vecX=zeros(nx,nuH+2*horizon+~yTer);     %vector which will store [A^{i-1}B ... AB B] in the following iteration
    Asum=zeros(nx,nx);

    %eps_xN(i)=vecExN*U   (eps_xN(i)=>P*x(i) & eps_xN(i)=>-P*x(i))
    vecExN=[zeros(xNcost,nuH+horizon) ones(xNcost,1) zeros(xNcost,horizon)]; %eps_xN(i)=vecExN*U
    
    for i=1:horizon         
        if iscell(probStruct.R),
            R = probStruct.R{i};
        end 
        if iscell(probStruct.Q),
            Q = probStruct.Q{i};
        end
        AAi=A*Ai;               %=A^i used for faster computation
        Asum=Asum+Ai;
        %x(i)=AAi*x+vecX*U
        %u(i)=vecU*U
        %eps_x(i)=vecEx*U   (eps_x(i)=>Q*x(i) & eps_x(i)=>-Q*x(i))
        %eps_u(i)=vecEu*U   (eps_u(i)=>R*u(i) & eps_u(i)=>-R*u(i))
        vecX=[Ai*B vecX];       %store [A^{i-1}B ... AB B] 
        vecX(:,(nuH+2*horizon+~yTer+1):end)=[];           %store [A^{i-1}B ... AB B] 
        vecU=[zeros(nu,(i-1)*nu) eye(nu) zeros(nu,(horizon-i)*nu+2*horizon+~yTer)];   %u(i)=vecU*U
        vecEx=[zeros(xcost,nuH+i) ones(xcost,1) zeros(xcost,horizon-i-yTer) zeros(xcost,horizon)];     %eps_x(i)=vecEx*U
        vecEu=[zeros(ucost,nuH+horizon+~yTer) zeros(ucost,i-1) ones(ucost,1) zeros(ucost,horizon-i)]; %eps_u(i)=vecEu*U

        
        if(i==1)
            %add cost for state x0
            vecEx0=[zeros(xcost,nuH) ones(xcost,1) zeros(xcost,horizon-yTer) zeros(xcost,horizon)]; 
            
            G=[G;-vecEx0+Q*(D*vecU)];
            W=[W;-Q*gg];
            E=[E;-Q*C];
            G=[G;-vecEx0-Q*(D*vecU)];
            W=[W;Q*gg];
            E=[E;Q*C];
        end
        if(i<horizon)
            if(ycost)
                %shift current input by one for output tracking
                vecU=[zeros(nu,i*nu) eye(nu) zeros(nu,(horizon-i-1)*nu+2*horizon+~yTer)];   %u(i)=vecU*U
            end
            %Q*(A^i*x+[A^{i-1}B ... AB B]*x)<=eps1
            G=[G;Q*C*vecX-vecEx+Q*(D*vecU)];
            W=[W;-Q*gg-Q*C*Asum*f];
            E=[E;-Q*C*AAi];
            %-eps1<=Q*(A^i*x+[A^{i-1}B ... AB B]*x)
            G=[G;-Q*C*vecX-vecEx-Q*(D*vecU)];
            W=[W;Q*gg+Q*C*Asum*f];
            E=[E;Q*C*AAi];
            if(ycost)
                %shift current input back
                vecU=[zeros(nu,(i-1)*nu) eye(nu) zeros(nu,(horizon-i)*nu+2*horizon+~yTer)];   %u(i)=vecU*U
            end
        elseif ~ycost
            %Q*(A^i*x+[A^{i-1}B ... AB B]*x)<=eps1
            G=[G;P*C*vecX-vecExN+P*D*vecU];
            W=[W;-P*gg-P*C*Asum*f];
            E=[E;-P*C*AAi];
            %-eps1<=Q*(A^i*x+[A^{i-1}B ... AB B]*x)
            G=[G;-P*C*vecX-vecExN-P*D*vecU];
            W=[W;P*gg+P*C*Asum*f];
            E=[E;P*C*AAi];
        end
        %R*u<=eps2
        G=[G;R*vecU-vecEu];
        W=[W;zeros(ucost,1)];
        E=[E;zeros(ucost,size(E,2))];
        %-eps<=R*u2
        G=[G;-R*vecU-vecEu];
        W=[W;zeros(ucost,1)];
        E=[E;zeros(ucost,size(E,2))];
        
        Ai=A*Ai;                %=A^{i-1} (in next iteration) used for faster computation
    end
    
    %New cost function (cost is actually defined by the new constraints)
    H=[zeros(1,nuH) ones(1,2*horizon+~yTer)];       %this is the H matrix
    F=zeros(1,nx); %for linear objectives there is no F part, cost is expressed through slacks
    Y=[];
    Cx=[];
    Cf=[];
    Cc=[];
    
% NEW FEATURE, NOT DOCUMENTED, ALLOWS SPECIFYING 1 NORM ON THE STATE AND INPUT + INF NORM ON THE LAST STATE
% USED WHEN WE WANT TO PASS COST EXPRESSION FOR COST TO GO AS A MAX EXPRESSION (USING INF NORM IS NATURAL
% IN SUCH A CASE)
elseif(probStruct.norm==[1 1 inf])
    disp('WARNING: In mpt_constructMatrices. You are using MIXED norm expression [1 1 inf] !!!');
    disp('This feature is experimental / not documented, and it may be removed in the future!!')
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    %  The original cost function    J= (sum_{i=0}^N  |Q x_i|_1 +  |R u_i|_1 +  |P x_N|_inf,
    %  is transformed into           J=  H [U epsilon] + F x0                          
    %                                
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    xcost=size(Q,1);
    ucost=size(R,1);
    xNcost=size(P,1);
    Ai=eye(nx);             %used for faster computation
    
    %    G=[G zeros(size(G,1),(xcost+ucost)*horizon)]; %extend state space to introduce slack variables
     G=[G zeros(size(G,1),(xcost+ucost)*horizon+1)]; %extend state space to introduce slack variables (add cost for x0 as well)
    
    nuH=nu*horizon;         %degrees of freedom = horizon  * number of inputs
    vecX=zeros(nx,nuH+(xcost+ucost)*horizon+1);     %vector which will store [A^{i-1}B ... AB B zeros(nx,length(epsilon))]
                                                         %in the following iteration
    Asum=zeros(nx,nx);

    
    %eps_xN(i)=vecExN*U (eps_xN(i)=>P*x(i) & eps_xN(i)=>-P*x(i))
    vecExN=[zeros(xNcost,horizon*xcost+nuH) ones(xNcost,1) zeros(xNcost,horizon*ucost)]; %eps_xN(i)=vecExN*U
    
    for i=1:horizon         
        AAi=A*Ai;               %=A^i used for faster computation
        Asum=Asum+Ai;
        %x(i)=AAi*x+vecX*U+Asum*f
        %u(i)=vecU*U
        %eps_x(i)=vecEx*U   (eps_x(i)=>Q*x(i) & eps_x(i)=>-Q*x(i))
        %eps_u(i)=vecEu*U   (eps_u(i)=>R*u(i) & eps_u(i)=>-R*u(i))
        vecX=[Ai*B vecX];       %store [A^{i-1}B ... AB B [bunch of zeros till its full length]] 
        vecX(:,(nuH+(xcost+ucost)*horizon+1+1):end)=[];           %store [A^{i-1}B ... AB B] 
        vecU=[zeros(nu,(i-1)*nu) eye(nu) zeros(nu,(horizon-i)*nu+(xcost+ucost)*horizon+1)];   %u(i)=vecU*U
        %         vecEx=[zeros(xcost,(i-1)*xcost+nuH) eye(xcost) zeros(xcost,(horizon-i)*xcost) zeros(xcost,ucost*horizon)]; %eps_x(i)=vecEx*U
        %         vecEu=[zeros(ucost,horizon*xcost+nuH) zeros(ucost,ucost*(i-1)) eye(ucost) zeros(ucost,ucost*(horizon-i))]; %eps_u(i)=vecEu*U
        vecEx=[zeros(xcost,i*xcost+nuH) eye(xcost) zeros(xcost,(horizon-i-1)*xcost+1) zeros(xcost,ucost*horizon)]; %eps_x(i)=vecEx*U
        vecEu=[zeros(ucost,horizon*xcost+1+nuH) zeros(ucost,ucost*(i-1)) eye(ucost) zeros(ucost,ucost*(horizon-i))]; %eps_u(i)=vecEu*U

        
        if(i==1)
            %add cost for state x0/y0
            vecEx0=[zeros(xcost,nuH) eye(xcost) zeros(xcost,(horizon-1)*xcost+1) zeros(xcost,ucost*horizon)];
            
            G=[G;-vecEx0+Q*D*vecU];
            W=[W;-Q*gg];
            E=[E;-Q*C];
            G=[G;-vecEx0-Q*D*vecU];
            W=[W;Q*gg];
            E=[E;Q*C];
        end
        if(i<horizon)
            %Q*(A^i*x+[A^{i-1}B ... AB B]*x)<=eps1
            %G=[G;Q*vecX zeros(xcost,(i-1)*xcost) -eye(xcost) zeros(xcost,(horizon-i)*xcost) zeros(xcost,ucost*horizon)];
            G=[G;Q*C*vecX-vecEx+Q*D*vecU];
            W=[W;-Q*gg-Q*C*Asum*f];
            E=[E;-Q*C*AAi];
            %-eps1<=Q*(A^i*x+[A^{i-1}B ... AB B]*x)
            G=[G;-Q*C*vecX-vecEx-Q*D*vecU];
            W=[W;Q*gg+Q*C*Asum*f];
            E=[E;Q*C*AAi];
        elseif ~ycost
            %P*(A^i*x+[A^{i-1}B ... AB B]*x)<=eps1
            G=[G;P*C*vecX-vecExN+P*D*vecU];
            W=[W;-P*gg-P*C*Asum*f];
            E=[E;-P*C*AAi];
            %-eps1<=P*(A^i*x+[A^{i-1}B ... AB B]*x)
            G=[G;-P*C*vecX-vecExN-P*D*vecU];
            W=[W;P*gg+P*C*Asum*f];
            E=[E;P*C*AAi];
        end
        %R*u<=eps2
        G=[G;R*vecU-vecEu];
        W=[W;zeros(ucost,1)];
        E=[E;zeros(ucost,size(E,2))];
        %-eps<=R*u2
        G=[G;-R*vecU-vecEu];
        W=[W;zeros(ucost,1)];
        E=[E;zeros(ucost,size(E,2))];
        
        Ai=A*Ai;                %=A^{i-1} (in next iteration) used for faster computation
    end
    
    %New cost function (cost is actually defined by the new constraints)
    H=[zeros(1,nuH) ones(1,(xcost+ucost)*horizon+1)];       %this is the H matrix
    F=zeros(1,nx);  %for linear objectives there is no F part, cost is expressed through slacks
    Y=[];
    Cx=[];
    Cf=[];
    Cc=[];
    
    
else
    fprintf('\n\n')
    error('No norm specified for cost objective');
end




%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%   Post-Processing of the Constraints 
%   (e.g. removal of redundant rows, reordering of constraints, etc.)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%Check if problem is symmetric / reorder constraints
%If the mp-solvers want to take advantage of symmetry, the constraints need to be ordered in a 
%predefined manner, i.e. they must alternate symmetrically. This is enforced by the lines below. 
GWE=[G W E];
tG=[];
tW=[];
tE=[];

foundmatch=0;
symmetric=probStruct.useSymmetry;
while(~isempty(GWE) & probStruct.useSymmetry)
    tG(end+1,:)   =GWE(1,1:nuH);
    tW(end+1,:)   =GWE(1,nuH+1);
    tE(end+1,:)   =GWE(1,(nuH+2):end);
    
    foundmatch=0;
    for i=2:size(GWE,1)
        normVec1=GWE(i,:);
        normVec1(find(isinf(normVec1)))=[];
        normVec2=[-tG(end,:) tW(end) -tE(end,:)];
        normVec2(find(isinf(normVec2)))=[];

        if(norm(normVec1-normVec2)<Options.abs_tol*size(GWE,2))
            foundmatch=1;
            tG(end+1,:)=GWE(i,1:nuH);
            tW(end+1,:)=GWE(i,nuH+1);
            tE(end+1,:)=GWE(i,(nuH+2):end);
            GWE([1 i]',:)=[];
            break
        end
    end
    if(foundmatch==0)
        if Options.verbose>1,
            symmetric=0;
            disp('mpt_constructMatrices: Constraints are not symmetric')
        end
        break
    end     
end
if(probStruct.useSymmetry)
    G=tG;
    W=tW;
    E=tE;
end
if Options.verbose>1,
    if(foundmatch==1)
        disp('mpt_constructMatrices: Constraints are symmetric')
    end
end

if ~isempty(bndA),
    % include Pbnd
    E=[E; -bndA];
    W=[W; bndb];
    G=[G; zeros(size(bndA,1),size(G,2))];
    Matrices.PbndIncluded = 1;
end

% Take out limits which are Inf
aux=find(isinf(W));
G(aux,:)=[];
W(aux,:)=[];
E(aux,:)=[];

% Take out possible constraints where the matrix [G K] has null rows
aux=find(all(([G E]==zeros(size([G E])))')'); % Rows which are all 0
G(aux,:)=[];
W(aux,:)=[];
E(aux,:)=[];

constraints_reduced = 0;
if(reduce_constraints)
    %remove all redundant constraints
    if Options.noReduce,
        GEW=polytope([G -E],W,0,2);
    else
        GEW=polytope([G -E],W);
    end
    if(~isfulldim(GEW)) 
        if Options.verbose>0,
            disp('mpt_constructMatrices: Problem is infeasible.')
        end
        G=zeros(1,nuH);
        E=zeros(1,nx);
        W=-Inf;
        S=E;
        H=-1;
        F=0;
        Hinv=1;
        Pinvset = polytope;
        symmetric = 0;
        if nargout==1,
            Matrices.G = G;
            Matrices.W = W;
            Matrices.E = E;
            Matrices.H = H;
            Matrices.F = F;
            Matrices.Y = Y;
            Matrices.Cf = Cf;
            Matrices.Cx = Cx;
            Matrices.Cc = Cc;
            Matrices.symmetric = symmetric;
            Matrices.bndA = bndA;
            Matrices.bndb = bndb;
            Matrices.Pinvset = Pinvset;
            Matrices.GEW = GEW;
            Matrices.constraints_reduced = 1;
            G = Matrices;
        end
        return
    end
    [GE,W]=double(GEW);
    G=GE(:,1:(end-nx));
    E=-GE(:,(end-nx+1):end);
    constraints_reduced = 1;
end


if(probStruct.norm==2)
    H=(H+H')/2;      %make sure hessian is symmetric
    Hinv=inv(H);
    S=E+G*Hinv*F';
    
    if(cond(H)>1e8)
        disp('mpt_constructMatrices:')
        disp(' +++++++++++++    WARNING     WARNING     WARNING     WARNING     +++++++++++++')
        disp('               THE CONDITION NUMBER OF THE H MATRIX IS VERY LARGE              ')
        disp('              THIS MAY LEAD TO NUMERICAL ERRORS AND INCONSISTENCIES        ')
        disp(' CHECK THE VALIDITY OF THE OBTAINED RESULTS AFTER THE ALGORITHM HAS COMPLETED')
    end

else
    %set to empty for 1 / Inf norm
    Hinv=[];
    S=[];
end



Pinvset = probStruct.Tset;

%final error check
if(probStruct.norm==2)
    if(size(H,2)~=size(G,2) | size(E,2)~=size(F,1) | size(E,1) ~= size(W,1) | size(F,1)~=size(E,2) | size(Y,2)~=size(E,2))
        error('mpt_constructMatrices: error in function, incompatible dimension of output matrices.')
    end
else
    if(size(H,2)~=size(G,2) | size(E,2)~=size(F,2) | size(E,1) ~= size(W,1))
        error('mpt_constructMatrices: error in function, incompatible dimension of output matrices.')
    end
end


if nargout==1,
    Matrices.G = G;
    Matrices.W = W;
    Matrices.E = E;
    Matrices.H = H;
    Matrices.F = F;
    Matrices.Y = Y;
    Matrices.Cf = Cf;
    Matrices.Cx = Cx;
    Matrices.Cc = Cc;
    Matrices.symmetric = symmetric;
    Matrices.bndA = bndA;
    Matrices.bndb = bndb;
    Matrices.Pinvset = Pinvset;
    Matrices.constraints_reduced = constraints_reduced;
    G = Matrices;
end


return



%--------------------------------------------------------
function noiseVal=subtractNoise(Apow,e,noise,Options)
   
polynoise = isa(noise, 'polytope');
noiseVal=zeros(size(e,1),1);
eA=e*Apow;
if ~polynoise,
    % NOTE! noise is represented column-wise if it is in V-representation
    noiseVal = min(eA * noise,[],2);
else
    [Hnoise,Knoise]=double(noise);
    for ii=1:size(e,1)
        % minimize in each direction
        [xopt,fval,lambda,exitflag,how]=mpt_solveLPi(eA(ii,:),Hnoise,Knoise,[],[],[],Options.lpsolver);
        if(~strcmp(how,'ok'))
            fprintf('\n\n')
            error(' mpt_contructMatrices: Not possible to robustify constraints');
        end
        noiseVal(ii)=eA(ii,:)*xopt;        
    end
end


%--------------------------------------------------------
function Ai_stack = sub_allAuncCombs(Aunc, N, nx)

nunc = length(Aunc);

if N==0,
    % special case, return identity matrix
    Ai_stack{1} = eye(nx);
    return
end


indices = [];
for ii=1:nunc,
    indices = [indices repmat(ii, 1, N)];
end

% compute all possible combinations of matrices
% e.g.
% A1, A2 for N=1
% A1*A1, A1*A2, A2*A2 for N==2
% A1*A1*A1, A1*A1*A2, A1*A1*A3, A1*A2*A2, A1*A2*A3
% etc.

unique_indices = unique(nchoosek(indices, N), 'rows');
nrows = size(unique_indices, 1);
Ai_stack = cell(1, nrows);
for ir = 1:nrows,
    onerow = unique_indices(ir, :);
    Ai_stack{ir} = Aunc{onerow(1)};
    for ic = 2:length(onerow),
        Ai_stack{ir} = Ai_stack{ir} * Aunc{onerow(ic)};
    end
end

%--------------------------------------------------------
function [A,B,C,D,Q,R,ymin,ymax,umin,umax,dumin,dumax,bndA,bndb]=mpt_evalSystem(sysStruct,probStruct)
%MPT_EVALSYSTEM Extracts data from sysStruct and probStruct structures
%
% [A,B,C,D,Q,R,ymin,ymax,umin,umax,dumin,dumax,bndA,bndb]=mpt_evalSystem(sysStruct,probStruct)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Extracts informations from the system and problem definition structures
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% sysStruct      - System structure in the sysStruct format
% probStruct     - Problem structure in the probStruct format
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% A,B,C,D        - System dynamics matrices
% Q,R            - Weighting matrices in the cost index
% ymin, ymax     - min/max constraints on the output
% umin, umax     - min/max constraints on the control input
% dumin, dumax   - min/max constraints on the slew rate of control input
% bndA, bndb     - region of exploration bndA*x<=bndb
%

% Copyright is with the following author(s):
%
% (C) 2004 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2003 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2003 Pascal Grieder, Automatic Control Laboratory, ETH Zurich,
%          grieder@control.ee.ethz.ch

%
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

narginchk(2, 2);

if ~isfield(sysStruct,'verified'),
    sysStruct=mpt_verifySysStruct(sysStruct);
end

if ~isfield(probStruct,'verified'),
    probStruct=mpt_verifyProbStruct(probStruct);
end

if iscell(sysStruct.A),
    A=sysStruct.A{1};
    B=sysStruct.B{1};
    C=sysStruct.C{1};
    D=sysStruct.D{1};
else
    A=sysStruct.A;
    B=sysStruct.B;
    C=sysStruct.C;
    D=sysStruct.D;
end
    
ymin=sysStruct.ymin;
ymax=sysStruct.ymax;
umin=sysStruct.umin;
umax=sysStruct.umax;
dumin=sysStruct.dumin;
dumax=sysStruct.dumax;

if iscell(probStruct.Q),
    Q = probStruct.Q{end};
else
    Q=probStruct.Q;
end
if iscell(probStruct.R),
    R = probStruct.R{end};
else
    R=probStruct.R;
end

if(isfield(sysStruct,'Pbnd'))
    [bndA, bndb]=double(sysStruct.Pbnd);
    Pbnd = sysStruct.Pbnd;
elseif isfield(sysStruct,'bndA'),
    bndA=sysStruct.bndA;
    bndb=sysStruct.bndb;
else
    bndA=[];
    bndb=[];
end

if(~isfield(sysStruct,'polyUncert'))
    polyUncert=0;
else
    polyUncert=sysStruct.polyUncert;
end
if(~isfield(sysStruct,'addUncert'))
    addUncert=0;
else
    addUncert=sysStruct.addUncert;
end


%------------------------------------------------------------
function [PinvSet,tstar,fd] = mpt_infset(A_CL, X, varargin)

sys = LTISystem('A', A_CL{1});
sys.setDomain('xu', toPolyhedron(X));
PinvSet = polytope(sys.invariantSet());
tstar = 0;
fd = [];

