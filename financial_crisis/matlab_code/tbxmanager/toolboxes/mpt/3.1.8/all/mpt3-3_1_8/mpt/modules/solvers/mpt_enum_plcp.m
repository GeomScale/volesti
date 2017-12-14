function sol = mpt_enum_plcp(problem)
%
% enumeration based PLCP solver for parametric linear complementarity
% problems (PLCP) of the form
% 
%     find w, z
%      s.t.    w - M*z = q + Q*th
%              w(th)'*z(th) = 0
%              w(th) >= 0
%              z(th) >= 0
%              Ath*th <= bth
%
%  where w, z are the decision variables and th is the parametric variable
%  The data of the problem are
%        - matrix M (which does not have to be positive semidefinite)
%        - vector q
%        - matrix Q
%        - matrix Ath (constraints on the parameters th)
%        - vector bth (constraints on the parameters th)
%
% SYNTAX:
%   U = mpt_enum_plcp(problem)
%
% INPUT:
%       problem - PLCP formulation of parametric optimization problem
%                  given as "Opt" class of MPT3
%
% OUTPUT:
%       U - PolyUnion object representing the explicit solution 
%

% AUTHOR: Colin N. Jones 2013, EPFL, colin.jones@epfl.ch
%         Revised in 2014 by M. Herceg, ETH, herceg@control.ee.ethz.ch
%

% global options
global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

% Extract the input data
M = problem.M;
Q = problem.Q;
q = problem.q;

n = size(M,1);
A = full([eye(n) -M]);

if ~isempty(problem.recover)
    problem.Data.P = speye(size(problem.recover.uX,1));
    % reorder variables is this came from Yalmip
    if ~isempty(problem.varOrder)
        problem.Data.P = problem.Data.P(problem.varOrder.requested_variables,:);
    end
end

% bigM = 1e6;
% % Generate a YALMIP optimizer object for testing feasibility
% w = sdpvar(n,1); wBnd = sdpvar(n,1);
% z = sdpvar(n,1); zBnd = sdpvar(n,1);
% x = sdpvar(size(Q,2),1);
% con = [w - M*z == Q*x + q] + [(1-wBnd)*bigM >= w >= 0] + [(1-zBnd)*bigM >= z >= 0];
% con = con + [ problem.Ath*x <= problem.bth ];
% if ~isempty(problem.Ae)
%     con = con + [ problem.Ae*z == problem.be ];
% end
% % If wBnd(i) == 0, then w(i) = 0, same for z
% s = sdpsettings;
% opt = optimizer(con, 0, s, [wBnd;zBnd], [w;z;x]);

% if isempty(problem.Ae)
%     Ae = [eye(n) -M -Q];
%     be = q;
% else
%     Ae = [eye(n) -M -Q;
%             zeros(problem.me,n) problem.Ae zeros(problem.me,problem.d)];
%     be = [q; problem.be];
% end
% 

S.A = [-M -Q;  -eye(n) zeros(n,problem.d)];
S.b = [q; zeros(n,1)];
S.f = zeros(size(S.A, 2), 1);
S.lb = [];
S.ub = [];
S.quicklp = true; % performance!
I = eye(n);


%%
B      = zeros(n,1); % Feasible bases
INFEAS = zeros(n,0); % Infeasible bases

%% 
fprintf('%-10s%-10s%-12s%-10s\n', 'Level', 'Feasible', 'Infeasible', 'LPs');

% prepare flags
flag = MPTOPTIONS.OK;
how = 'ok';

% iterate over equations
tStart = clock; tic;
numLPs = 0;
for iii = 1:n
    
    if toc > MPTOPTIONS.report_period
        fprintf('%2i / %2i', iii, n);
        fprintf('%10i', size(B,2));
        fprintf('%12i', size(INFEAS,2));
        fprintf('%10.2e', numLPs);
        fprintf('\n');
        tic;
    end
    
    % For each basis on the stack, generate the children for level iii
    TMP = B; TMP(iii,:) = 1; B(iii,:) = -1;
    B = [TMP B];
    
    % Bases to be tested:
    Izero  = [B<0;B>0]; % Variables fixed to zero
    Ibasis = [B>0;B<0]; % Variables in the basis
    
    % First check - are they full rank?
    rankTest = ones(1,size(Izero,2));
    for i = 1:size(Izero,2)
        if rank(A(:,Ibasis(:,i))) < iii
            rankTest(i) = 0;
        end
    end
       indTest = find(rankTest);
    
    % Test if each basis has any feasible children
%    [res0, lp_infeas] = opt{Izero(:,indTest)};
    lp_infeas = zeros(1,size(indTest,2));
%    lp_infeasn = lp_infeas;
    for jj=1:size(indTest,2)
        % MPT call
        wbnd = Izero(1:n,indTest(jj));
        zbnd = Izero(n+1:2*n,indTest(jj));
        nz = nnz(zbnd);

        if isempty(problem.Ae)
            % no equality constraints on binary variables
            S.Ae = [-M(wbnd,:) -Q(wbnd,:);
                I(zbnd,:) zeros(nz,problem.d)];
            S.be = [q(wbnd); zeros(nz,1)];
        else
            % append equality constraints on binary variables if present
            S.Ae = [-M(wbnd,:) -Q(wbnd,:);
                I(zbnd,:) zeros(nz,problem.d);
                problem.Ae zeros(problem.me,problem.d)];
            S.be = [q(wbnd); zeros(nz,1); problem.be];
        end
        

        res = mpt_solve(S);
        
        lp_infeas(jj) = (res.exitflag==MPTOPTIONS.INFEASIBLE);
           
    end
    
        
    numLPs = numLPs + sum(rankTest);
    rankTest(indTest(lp_infeas==1)) = 0;
    
    % Store the infeasible solutions
    INFEAS = [INFEAS B(:,rankTest==0)];
    B(:,rankTest==0) = [];
    
    % break in case of maximum LPs
    if numLPs>=MPTOPTIONS.modules.solvers.enum_plcp.maxLPs
        disp('Maximum LPs achieved. Interrupting ...');
        flag = -3;
        how = 'maximum LPs achieved';
        break;
    end
    
    % break in case of maximum regions
    if size(B,2)>=MPTOPTIONS.modules.solvers.enum_plcp.maxregions
        disp('Maximum number of regions achieved. Interrupting ...');
        flag = -2;
        how = 'maximum number of regions achieved';
        break;
    end
    
    
end


% final statistics
fprintf('%2i / %2i', iii, n);
fprintf('%10i', size(B,2));
fprintf('%12i', size(INFEAS,2));
fprintf('%10.2e', numLPs);
fprintf('\n');

% construct regions
nreg = size(B,2);
R = [];
for i=1:nreg
    basis = [B(:,i)>0; B(:,i)<0];
    region = getRegion(basis, problem);
    if ~isempty(region)
        if region.isFullDim && region.isBounded
            R = [R, region];
        end
    end   
end

% output
if ~isempty(R)
    U = PolyUnion('Set',R,'FullDim',true);
else
    U = PolyUnion;
    flag = MPTOPTIONS.INFEASIBLE;
    how = 'infeasible';
end

% return variable
sol.xopt = U;
sol.exitflag = flag;
sol.how = how;
sol.stats.solveTime = etime(clock, tStart);
sol.stats.numLPs = numLPs;

fprintf('mpt_enum_plcp: %d regions\n',U.Num);

if MPTOPTIONS.verbose >= 1
    % show statistics
    fprintf('Computation time : %.2f seconds\n', time);
end


% store internal data
U.setInternal('feasible_bases',B);
U.setInternal('infeasible_bases',INFEAS);



end

function region = getRegion(basis,problem)
% construct region for a given basis
global MPTOPTIONS

A = [eye(problem.n) -problem.M];

H = A(:,basis)\[-problem.Q problem.q];
if any(any(isnan(H)))
    region = [];
    return;
end
Hn = [H; problem.Ath problem.bth];
zero_rows = sum(abs(Hn),2)<MPTOPTIONS.zero_tol;
Hn(zero_rows,:)=[];

% H-represenation of a region
region = Polyhedron(Hn(:, 1:end-1), Hn(:, end));
region.normalize;

% add function handles
x = zeros(problem.n*2,problem.d+1); % Solution
H(:,1:problem.d) = -H(:,1:problem.d);
x(basis,:) = H;

% set function handles for w and z variables
Lz = AffFunction(x(problem.n+1:2*problem.n,1:end-1),x(problem.n+1:2*problem.n,end));
Lw = AffFunction(x(1:problem.n,1:end-1),x(1:problem.n,end));
region.addFunction(Lz,'z');
region.addFunction(Lw,'w');

% if LCP was created from LP/QP, compute also primal, dual variables and
% the value function
if ~isempty(problem.recover)
    TT = problem.recover.uX*x + problem.recover.uTh;
    
    % compute primal variables
    Lprimal = problem.Data.P*TT;
    region.addFunction(AffFunction(Lprimal(:,1:end-1),Lprimal(:,end)),'primal');
    
    % compute dual variables for inequalities
    %Ldual = problem.recover.lambdaX*x + problem.recover.lambdaTh;
    Lineq = problem.recover.lambda.ineqlin.lambdaX*x + problem.recover.lambda.ineqlin.lambdaTh;
    region.addFunction(AffFunction(Lineq(:,1:end-1),Lineq(:,end)),'dual-ineqlin');
    
    % dual variables for equalities
    Leq = problem.recover.lambda.eqlin.lambdaX*x + problem.recover.lambda.eqlin.lambdaTh;
    region.addFunction(AffFunction(Leq(:,1:end-1),Leq(:,end)),'dual-eqlin');
    
    % dual variables for lower bounds
    Llb = problem.recover.lambda.lower.lambdaX*x + problem.recover.lambda.lower.lambdaTh;
    region.addFunction(AffFunction(Llb(:,1:end-1),Llb(:,end)),'dual-lower');
    
    % dual variables for upper bounds
    Lub = problem.recover.lambda.upper.lambdaX*x + problem.recover.lambda.upper.lambdaTh;
    region.addFunction(AffFunction(Lub(:,1:end-1),Lub(:,end)),'dual-upper');
    
    % compute the objective value
    Y = TT(:,1:end-1);
    T = TT(:,end);
    
    if ~isempty(problem.Internal.H)
        qt = 0.5*Y'*problem.Internal.H*Y + problem.Internal.pF'*Y + problem.Internal.Y;
        lt = T'*problem.Internal.H*Y + T'*problem.Internal.pF + problem.Internal.f'*Y + problem.Internal.C;
        at = 0.5*T'*problem.Internal.H*T + problem.Internal.f'*T + problem.Internal.c;
        region.addFunction(QuadFunction(qt,lt,at),'obj');
    else
        qt = problem.Internal.pF'*Y + problem.Internal.Y;
        lt = T'*problem.Internal.pF + problem.Internal.f'*Y + problem.Internal.C;
        at = problem.Internal.f'*T + problem.Internal.c;
        if any(any(qt))
            region.addFunction(QuadFunction(qt,lt,at),'obj');
        else
            region.addFunction(AffFunction(lt,at),'obj');            
        end
    end
end

% add basis to Internal
region.setInternal('basis',basis);

% check complementarity conditions
% Lw  = Fw*th + gw
% Lz  = Fz*th + gz
% w'*z = 0
% (Fw*th + gw)'*(Fz*th +gz) = 0
% th'*(Fw'*Fz)*th + (gz'*Fw + gw'*Fz)*th + gw'*gz = 0

if Lw.g'*Lz.g > MPTOPTIONS.abs_tol    
    region = [];    
else
    if norm(Lz.g'*Lw.F + Lw.g'*Lz.F,Inf)>MPTOPTIONS.abs_tol
        region = [];
    else
       if norm(Lw.F'*Lz.F,Inf)>MPTOPTIONS.abs_tol
           region = [];
       end
    end
end

 
end
