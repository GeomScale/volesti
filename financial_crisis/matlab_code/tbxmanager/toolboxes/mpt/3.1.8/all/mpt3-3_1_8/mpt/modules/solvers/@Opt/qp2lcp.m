function S = qp2lcp(S,reduce_flag)
%
%  ELIMINATEEQUATIONS: Transforms LP/QP/MPLP/MPQP to LPC/PLCP 
%  ===========================================================
%  
%  
%  SYNTAX
%  ------
%     
%      problem.qp2lcp
%      qp2lcp(problem)
%    
%  
%  DESCRIPTION
%  -----------
%     Transformation of LP, QP, MPLP, and MPQP to LCP/PLCP formulation. Consider
%  the following MPQP format that is accepted by Opt class: 
%                                   1  T             T              
%                             min   - x Hx+(Ftheta+f) x         (1) 
%                                   2                               
%                            s.t.   Ax <= b + Btheta            (2) 
%                                                                   
%                                   A x = b  + Etheta           (3) 
%                                    e     e                        
%                                                                   
%                                   A     theta = b             (4) 
%                                    theta         theta            
%                                                               (5) 
%     which contains minequality constrains and m_e  equality constraints and
%  constraints on the parameter theta. If there are lower and upper bounds on the
%  variables x  present, i.e. lb and ub, these can be merged to inequalities (??).
%  This format is not appropriate for transformation to LCP form because it
%  contains equality constraints and the inequalities are not nonnegative. To get
%  appropriate LCP representation, the equality constrains of the problem (??)-(??)
%   are removed using eliminateEquations method of the Opt class. The intermediate
%  form of the optimization problem is given as 
%                               1    T                 T               
%                         min   - x   Hx  +(Ftheta + f) x          (6) 
%                               2  Nc   Nc               Nc            
%                                                                      
%                        s.t.   Ax   <= b + Btheta                 (7) 
%                                 Nc                                   
%                                                                      
%                               A     theta = b                    (8) 
%                                theta         theta                   
%                                                                  (9) 
%     with x_Nc  as the non-basic variables that map to x  affinely. Problem
%  (??)-(??)  can be transformed effectively when considering the rank of the
%  matrix A. If the rank of matrix A  is less than the number of inequalities in
%  (??), then vector x_N  can be expressed as a difference of two positive numbers,
%  i.e. 
%                                            +      -           
%                                x    =   x    - x         (10) 
%                                 Nc       Nc     Nc            
%                                  +                            
%                               x     >=  0                (11) 
%                                Nc                             
%                                  -                            
%                               x     >=  0                (12) 
%                                Nc                             
%     Using the substitution 
%                                    (    +  )        
%                                    ( x     )        
%                                    (  Nc   )        
%                                v = (    -  )    (13)
%                                    ( x     )        
%                                    (  Nc   )        
%     and putting back to (??)-(??)  we get 
%                     1  T (        )    ((       )      (       ))T              
%               min   - v  ( H  -H  )v + (( F -F  )theta ( f -f  ))  v       (14) 
%                     2    ( -H H   )    ((       )      (       ))               
%                     (       )                                                   
%              s.t.   ( -A A  )v >= -b -Btheta                               (15) 
%                     (       )                                                   
%                     v >= 0                                                 (16) 
%     which is a form suitable for PLCP formulation. If the rank of matrix A  is
%  greater or equal than the number of inequalities in (??), we can factorize the
%  matrix A  rowwise 
%                                        (     )
%                                        ( A   )
%                                        (  B  )
%                                    A = (     )
%                                        ( A   )
%                                        (  N  )
%    where B, N  are index sets corresponding to rows from which submatrices are
%  built. The factored system can be written as 
%                                        -b                       
%                              -A x  =       -Btheta + y     (17) 
%                                B         B                      
%                                        -b                       
%                              -A x >=       -Btheta         (18) 
%                                N         N                      
%                                 y  >=  0                   (19) 
%     where the matrix A_B  form by rows in the set B  must be invertible. Using
%  this substitution the system (??)(??)  can be rewritten in variable y 
%                                 1  T                 T              
%                           min   - y  H y + (Ftheta+f)  y       (20) 
%                                 2                                   
%                          s.t.   A >= b + Btheta                (21) 
%                                 y >= 0                         (22) 
%     which is suitable for PLCP formulation where 
%                                     -T   -1                      
%                             H  =  A   HA                    (23) 
%                                    B    B                        
%                                     A -T   -1      -T            
%                             F  =  -(    HA   B  -A   F      (24) 
%                                      B    B   B   B              
%                                      -T   -1     -T              
%                             f  =  -A   HA   b -A   f        (25) 
%                                     B    B   B  B                
%                                       -1                         
%                             A  =  A A                       (26) 
%                                    N B                           
%                                       -1                         
%                             b  =  A A   b -b                (27) 
%                                    N B   B  N                    
%                                       -1                         
%                             B  =  A A   B -b  -B            (28) 
%                                    N B   B  N   N                
%     The corresponding PLCP can be written as follows: 
%                               w - Mz  =   q + Qtheta      (29) 
%                                    w  >=  0               (30) 
%                                    z  >=  0               (31) 
%                                   T                            
%                                  w z  =   0               (32) 
%                                                                
%     where the problem data are built from MPQP (??)-(??) 
%                                        (     T  )          
%                                        ( H -A   )          
%                                  M  =  (        )     (33) 
%                                        ( A  0   )          
%                                        (     )             
%                                  q  =  ( f   )        (34) 
%                                        ( -b  )             
%                                        (     )             
%                                  Q  =  ( F   )        (35) 
%                                        ( -B  )             
%     Original solution to MPQP problem (??)-(??)  can be obtained by affine map
%  from the variables w(theta)  and z(theta)  to x. The matrices of the backward
%  map are stored inside recover property of the Opt class as follows 
%                                  (    )     (        )
%                            x = uX( w  )+ uTh( theta  )
%                                  ( z  )     (   1    )
%     If the problem was formulated using YALMIP, it is possible that some
%  variables are in the different order. The original order of variables is stored
%  in problem.varOrder.requested_variables and the map to original variables is
%  given by 
%                                  (    )          (        )
%                       x = primalX( w  )+ primalTh( theta  )
%                                  ( z  )          (   1    )
%  
%  
%  INPUT
%  -----
%     
%        
%          problem LP/QP/MPLP/MPQP optimization problem     
%                  given as Opt class.                      
%                  Class: Opt                               
%                    
%  
%  
%  SEE ALSO
%  --------
%     solve,  mpt_call_lcp
%  

%  AUTHOR(s)
%  ---------
%     
%    
%   (c) 2010-2013  Colin Neil Jones: EPF Lausanne
%   mailto:colin.jones@epfl.ch 
%     
%    
%   (c) 2010-2013  Martin Herceg: ETH Zurich
%   mailto:herceg@control.ee.ethz.ch 
%  
%  

%  LICENSE
%  -------
%    
%    This program is free software; you can redistribute it and/or modify it under
%  the terms of the GNU General Public License as published by the Free Software
%  Foundation; either version 2.1 of the License, or (at your option) any later
%  version.
%    This program is distributed in the hope that it will be useful, but WITHOUT
%  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
%  FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%    You should have received a copy of the GNU General Public License along with
%  this library; if not, write to the  Free Software Foundation, Inc.,  59 Temple
%  Place, Suite 330,  Boston, MA 02111-1307 USA
%  -------------------------------------------------------------------------------
%    
%      This document was translated from LaTeX by HeVeA (1).
%  ---------------------------------------
%    
%    
%   (1) http://hevea.inria.fr/index.html
 
 
global MPTOPTIONS

if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

% if no arguments, do redundancy elimination automatically
if nargin<=1
    reduce_flag = true;
else
    % check input argument if provided
    validate_logicalscalar(reduce_flag);
end

% deal with arrays
if numel(S)>1
    for i=1:numel(S)
        S(i).qp2lcp(reduce_flag);
    end
    return
end

% for LCP do nothing
if strcmpi(S.problem_type,'LCP')
    return;
end

% store the original objective function
S.Internal.H = S.H;
S.Internal.pF = S.pF;
S.Internal.f = S.f;
S.Internal.Y = S.Y;
S.Internal.C = S.C;
S.Internal.c = S.c;

% store the original constraints (used in Opt/feasibleSet)
S.Internal.constraints.A = S.A;
S.Internal.constraints.b = S.b;
S.Internal.constraints.pB = S.pB;
S.Internal.constraints.Ae = S.Ae;
S.Internal.constraints.be = S.be;
S.Internal.constraints.pE = S.pE;

if reduce_flag
%     % Remove redundancies to reduce the problem size if possible
%     Hp = [-S.pB S.A S.b;
%         S.Ath zeros(length(S.bth),S.n) S.bth;    
%         zeros(S.n,S.d) eye(S.n) S.ub;
%         zeros(S.n,S.d) -eye(S.n) -S.lb];
%     % remove Inf bounds if any
%     Hp(isinf(Hp(:,end)),:)=[];
% 
%     P = Polyhedron('H',Hp);
%     if isEmptySet(P),
%         error('Opt:qp2lcp', 'Feasible set is empty.');
%     end
%     hull = P.minHRep();
%     pB = -hull.H(:,1:S.d);
%     A = hull.H(:,S.d+1:S.n+S.d);
%     b  =  hull.H(:,end);

    % take internally stored values
    if S.isParametric       
        % the redundancy elimination was done by validation, take the
        % internal values
        pB = S.Internal.pB;
        A = S.Internal.A;
        b = S.Internal.b;
    else
        % perform the elimination, don't throw error if empty
        P = Polyhedron('A',S.Internal.A,'b',S.Internal.b,'Ae',S.Ae,'be',S.be);
        [P, hull] = P.minHRep();
        A = P.A;
        b = P.b;
        pB = S.Internal.pB;

        % update internal values
        S.Internal.A = A;
        S.Internal.b = b;
        
        % indices of constraints to be extracted
        ie = 1:S.Internal.m-numel(S.Internal.removed_rows.ineqlin);
        ilb = S.Internal.m-numel(S.Internal.removed_rows.ineqlin)+1:S.Internal.m-numel(S.Internal.removed_rows.ineqlin)+S.Internal.n-numel(S.Internal.removed_rows.lower);
        iub = S.Internal.m-numel(S.Internal.removed_rows.ineqlin)+S.Internal.n-numel(S.Internal.removed_rows.lower)+1:...
            S.Internal.m-numel(S.Internal.removed_rows.ineqlin)+S.Internal.n-numel(S.Internal.removed_rows.lower)+S.Internal.n-numel(S.Internal.removed_rows.upper);
        % store removed rows
        S.Internal.removed_rows.ineqlin = [S.Internal.removed_rows.ineqlin; find(hull.I(ie))];
        S.Internal.removed_rows.lower = [S.Internal.removed_rows.lower; find(hull.I(ilb))];
        S.Internal.removed_rows.upper = [S.Internal.removed_rows.upper; find(hull.I(iub))];
        
    end
    
else
    % take merged constraints 
    A = S.Internal.A;
    b = S.Internal.b;
    pB = S.Internal.pB;
end

% put all constraints to A, b, pB
S.A = A;
S.b = b;
S.pB = pB;
S.lb = [];
S.ub = [];
S.m = size(A,1);

% Remove equality constraints
equations = false;
if S.me > 0    
    % extract variables for recovering multipliers
    Ae = S.Ae;
    Am = A;
    if isempty(S.H)
        Hm = zeros(S.n);
    else
        Hm = S.H;
    end
    fm = S.f;    
    Fm = S.pF;
    
    % check if the rank of the factored matrix is equal to zero - only
    % in case the elimination cannot be executed    
    if isempty(S.vartype)
        S.vartype = repmat('C',1,S.n);
    end
    if rank(S.Ae(:,S.vartype=='C'),MPTOPTIONS.abs_tol)~=0
        equations = true;

        % eliminate equations
        S.eliminateEquations;
    end    
end

% store remaining equality constraints
if S.me>0
   eqlin.Ae = S.Ae;
   eqlin.be = S.be;
   eqlin.pE = S.pE;
end

% identify integer/binary variables
if ~isempty(S.vartype)
    if any(S.vartype=='S') || any(S.vartype=='N')
        error('qp2lcp: Semicontinuous (S) or semiinteger (N) variables are not supported for transformation to LCP.');
    end
    
    % find continuous variables
    ind_b = find(S.vartype=='B');
    if isempty(ind_b)
        ind_b = zeros(0,1);
    end
    ind_i = find(S.vartype=='I');
    ind_c = find(S.vartype=='C');
else
    % index sets in case no vartype is specified
    ind_c = 1:S.n;
    ind_i = zeros(0,1);
    ind_b = zeros(0,1);
end

% number of continuous variables
nc = numel(ind_c);

% Convert the problem to the standard form:
% min 0.5 u'*H*u + (F*th + f)'*u
% s.t. A*u <= b + B*th

if isempty(S.H)
    H = zeros(S.n);
else
    H = S.H;
end

% extract integer variables and map into inequality constraints
% alpha1*xB + alpha2*yB + beta*xN <= gamma + delta*th
if isempty(ind_i)
    % no integers
    % inequality constraints
    alpha1 = S.A(:,ind_b);
    alpha2 = zeros(S.m,0);
    beta = S.A(:,ind_c);
    gamma = S.b;
    if S.isParametric
        delta = S.pB;
    end
    % no change in the cost function
    f = S.f;
    if S.isParametric
        pF = S.pF;
    end
    
    % number of binary and continuous variables
    nxb = numel(ind_b);
    nyb = 0;

else
    % integer map xI = T*yB + t
    T = S.Internal.T;
    t = S.Internal.t;
    % inequality constraints
    alpha1 = S.A(:,ind_b);
    alpha2 = S.A(:,ind_i)*T;
    beta = S.A(:,ind_c);
    gamma = S.b-S.A(:,ind_i)*t;
    if S.isParametric
        delta = S.pB;
    end
    
    % number of binary and continuous variables
    nxb = numel(ind_b);
    nyb = size(T,2);    
    
    % cost function
    if ~isempty(S.H)
        % QP
        H = [S.H(ind_b,ind_b) S.H(ind_b,ind_i)*T S.H(ind_b,ind_c);
            T'*S.H(ind_i,ind_b) T'*S.H(ind_i,ind_i)*T T'*S.H(ind_i,ind_c);
            S.H(ind_c,ind_b) S.H(ind_c,ind_i)*T S.H(ind_c,ind_c)];
        f = [S.f(ind_b) + 0.5*(S.H(ind_i,ind_b)' + S.H(ind_b,ind_i))*t;
            T'*S.f(ind_i) + T'*S.H(ind_i,ind_i)*t;
            S.f(ind_c) + 0.5*(S.H(ind_i,ind_c)' + S.H(ind_c,ind_i))*t];
    else
        % LP
        H = zeros(nxb+nyb+nc);
        f = [S.f(ind_b);
            T'*S.f(ind_i);
            S.f(ind_c)];        
    end
    if S.isParametric
        pF = [S.pF(ind_b,:);
            T'*S.pF(ind_i,:);
            S.pF(ind_c,:)];
    end
end

% number of binary variables
nb = nxb + nyb;

% adjust the index sets if there are integers
ind_bn = 1:nb;
ind_cn = nb+1:nb+nc;

% get the rank
r = rank(beta,MPTOPTIONS.abs_tol);

if r<nc  
    % express xR as a difference of two positive numbers, i.e. xR = x+ - x-
       
    % new objective function
    Hnew = [ H(ind_bn,ind_bn)  H(ind_bn,ind_cn) -H(ind_bn,ind_cn);
             H(ind_cn,ind_bn)  H(ind_cn,ind_cn) -H(ind_cn,ind_cn);
            -H(ind_cn,ind_bn) -H(ind_cn,ind_cn)  H(ind_cn,ind_cn)];
    fnew = [f(ind_bn); f(ind_cn); -f(ind_cn)];
    if S.isParametric
        pFnew = [pF(ind_bn,:); pF(ind_cn,:); -pF(ind_cn,:)];
    end
    ind_bnew = 1:nb;
    ind_cnew = nb+1:nb+2*nc;
    
    % new constraints AnewB*[xB;yB] + AnewR*xR <= bnew + Bnew*th
    AnewB = [alpha1 alpha2];
    AnewR = [beta -beta];
    bnew = gamma;
    if S.isParametric
        Bnew = delta;
    end
    
    % modify equality constraint (if any) s.t. Ae*z = be + pE*th
    % [Ae(:,xb)*xb + Ae(:,x+)*x+ + Ae(:,x-)*x- + 0*lam ] = be + pE*th 
    if S.me>0
        if isempty(ind_i)
            S.Ae = [eqlin.Ae(:,ind_b) eqlin.Ae(:,ind_c) -eqlin.Ae(:,ind_c) zeros(S.me,S.m)];
        else            
            S.Ae = [eqlin.Ae(:,ind_b) eqlin.Ae(:,ind_i)*T eqlin.Ae(:,ind_c) -eqlin.Ae(:,ind_c) zeros(S.me,S.m)];
            S.be = eqlin.be - eqlin.Ae(:,ind_i)*t;
        end
    end
    
    % create M,q for LCP
    M = [-eye(nb), zeros(nb,2*nc), zeros(nb,S.m);
         0.5*(Hnew(ind_bnew,ind_cnew)'+Hnew(ind_cnew,ind_bnew)), Hnew(ind_cnew,ind_cnew), AnewR';
         -AnewB, -AnewR, zeros(S.m)];
    q = [ones(nb,1); fnew(ind_cnew); bnew];
    if S.isParametric
        Q = [zeros(nb,S.d); pFnew(ind_cnew,:); Bnew];
    else
        Q = [];
    end

    % Build mapping from the lcp solution back to the original variables
    % xR = z(nb+1:nb+nc) - z(nb+nc+1:nb+2*nc)
    %     
    % w  = [vB; vR+; vR-; s], z = [[xB;yB]; x+; x-; lam]
    %
    % x = [0 I -I]*[w] + [0]*[th]
    %              [z]       [1 ]
    % u = uX * [w;z] + uTh * [th;1]
    recover.uX = zeros(S.n,2*(nb+2*nc+S.m));
    recover.uX(ind_b,nb+2*nc+S.m+1:nb+2*nc+S.m+nxb) = eye(nxb);
    if ~isempty(ind_i)
        recover.uX(ind_i,nb+2*nc+S.m+nxb+1:nb+2*nc+S.m+nxb+nyb) = T;
    end
    recover.uX(ind_c,nb+2*nc+S.m+nb+1:nb+2*nc+S.m+nb+nc) = eye(nc);
    recover.uX(ind_c,nb+2*nc+S.m+nb+nc+1:nb+2*nc+S.m+nb+2*nc) = -eye(nc);
    recover.uTh = zeros(S.n,S.d+1);
    if ~isempty(ind_i)
        recover.uTh(ind_i,S.d+1) = t;
    end

    % Lagrange multipliers are:
    % lambda_ineq = z(nb+2*nc+1:end);
    % lambda_ineq = z(nb+2*nc+1:end) + [0]*[th;1]
    % lambda_ineq = lambdaX * [w;z] + lambdaTh *[th;1];
    nlam = S.m;
    lambdaX = zeros(S.m,2*(nb+2*nc+S.m));
    lambdaX(:,nb+2*nc+S.m+nb+2*nc+1:nb+2*nc+S.m+nb+2*nc+S.m) = eye(S.m);
    lambdaTh = zeros(S.m,S.d+1);

    
    % multipliers for the original inequalities
    kept_rows.ineq = setdiff(1:S.Internal.m,S.Internal.removed_rows.ineqlin);
    recover.lambda.ineqlin.lambdaX = zeros(S.Internal.m,2*(nb+2*nc+S.m));
    recover.lambda.ineqlin.lambdaX(kept_rows.ineq,nb+2*nc+S.m+nb+2*nc+1:nb+2*nc+S.m+nb+2*nc+S.Internal.m) = eye(numel(kept_rows.ineq),S.Internal.m);
    recover.lambda.ineqlin.lambdaTh = zeros(S.Internal.m,S.d+1);
    
    % multipliers for the original lower bound
    kept_rows.lb = setdiff(1:S.Internal.n,S.Internal.removed_rows.lower);    
    recover.lambda.lower.lambdaX = zeros(S.Internal.n,2*(nb+2*nc+S.m));        
    recover.lambda.lower.lambdaX(kept_rows.lb,nb+2*nc+S.m+nb+2*nc+S.Internal.m+1:...
        nb+2*nc+S.m+nb+2*nc+S.Internal.m+numel(kept_rows.lb)) = eye(numel(kept_rows.lb));
    recover.lambda.lower.lambdaTh = zeros(S.Internal.n,S.d+1);

    % multipliers for the original upper bound
    kept_rows.ub = setdiff(1:S.Internal.n,S.Internal.removed_rows.upper);
    recover.lambda.upper.lambdaX = zeros(S.Internal.n,2*(nb+2*nc+S.m));        
    recover.lambda.upper.lambdaX(kept_rows.ub,nb+2*nc+S.m+nb+2*nc+S.Internal.m+numel(kept_rows.lb)+1:...
        nb+2*nc+S.m+nb+2*nc+S.Internal.m+numel(kept_rows.lb)+numel(kept_rows.ub)) = eye(numel(kept_rows.ub));
    recover.lambda.upper.lambdaTh = zeros(S.Internal.n,S.d+1);

    
%     recover.lambdaX = zeros(S.m,2*nd);
%     recover.lambdaX(:,2*S.n+1:2*S.n+S.m) = eye(S.m);
%     recover.lambdaTh = zeros(S.m,S.d+1); 

    
%     % A'*lam = -H*x -[F f]*[th;1]
%     % A'*lam = -H*(uX*[w;z]+uTh*[th;1])-[F f]*[th;1]
%     % lam =  (-A'\H*uX)*[w;z] -A'\(H*uTh-[F f])*[th;1]
%     
%     recover.lambdaX = -S.A'\H*recover.uX;
%     recover.lambdaTh = -S.A'\(H*recover.uTh+[S.pF S.f]);

else

    % factorize beta to get
    %  alpha1(P,:)*xB + alpha2(P,:)*yB + beta(P,:)*xR = gamma(P) + delta(P,:)*th - y  % must be invertible mapping
    %  alpha1(Q,:)*xB + alpha2(Q,:)*yB + beta(Q,:)*xR <= gamma(Q) + delta(Q,:)*th
    %          y >= 0
    [L,U,p] = lu(sparse(beta),'vector');
    Pv = p(1:nc); % basic variables
    Qv = setdiff(1:S.m,Pv); % non-basic variables
    nlam = numel(Qv);
    
    % substitute
    % use factorized solution to compute inv(beta(P,:))
    ibetal = linsolve(full(L(1:nc,:)),eye(nc),struct('LT',true));
    ibetaP = linsolve(full(U),ibetal,struct('UT',true));
    %ibetaP = inv(beta(P,:));
    
    % form new objective function
    % xR = C1*xB + C2*yB + C3*y + D1 + D2*th
    C1 = -ibetaP*alpha1(Pv,:);
    C2 = -ibetaP*alpha2(Pv,:);
    C3 = -ibetaP;
    D1 = ibetaP*gamma(Pv);
    if S.isParametric
       D2 = ibetaP*delta(Pv,:); 
    end
    
    % new objective function
    if isempty(ind_i)
        % no integers
        if ~isempty(S.H)
            % QP
            Hnew = [S.H(ind_b,ind_b) + S.H(ind_b,ind_c)*C1 + C1'*S.H(ind_c,ind_b) + C1'*S.H(ind_c,ind_c)*C1,   S.H(ind_b,ind_c)*C3 + C1'*S.H(ind_c,ind_c)*C3;
                C3'*S.H(ind_c,ind_b) + C3'*S.H(ind_c,ind_c)*C1, C3'*S.H(ind_c,ind_c)*C3];
            fnew = [0.5*(S.H(ind_b,ind_c)+S.H(ind_c,ind_b)')*D1 + C1'*S.H(ind_c,ind_c)*D1 + S.f(ind_b) + C1'*S.f(ind_c);
                C3'*S.H(ind_c,ind_c)*D1 + C3'*S.f(ind_c)];
        else
            % LP
            Hnew = zeros(nb+nc);
            fnew = [S.f(ind_b) + C1'*S.f(ind_c);
                C3'*S.f(ind_c)];
        end
    else
        if ~isempty(S.H)
            % QP
            Hnew = [S.H(ind_b,ind_b) + S.H(ind_b,ind_c)*C1 + C1'*S.H(ind_c,ind_b) + C1'*S.H(ind_c,ind_c)*C1, ...
                S.H(ind_b,ind_i)*T + S.H(ind_b,ind_c)*C2 + C1'*S.H(ind_c,ind_i)*T + C1'*S.H(ind_c,ind_c)*C2, ...
                S.H(ind_b,ind_c)*C3 + C1'*S.H(ind_c,ind_c)*C3;
                T'*S.H(ind_i,ind_b) + T'*S.H(ind_i,ind_c)*C1 + C2'*S.H(ind_c,ind_b) + C2'*S.H(ind_c,ind_c)*C1, ...
                T'*S.H(ind_i,ind_i)*T + T'*S.H(ind_i,ind_c)*C2 + C2'*S.H(ind_c,ind_i)*T + C2'*S.H(ind_c,ind_c)*C2, ...
                T'*S.H(ind_i,ind_c)*C3 + C2'*S.H(ind_c,ind_c)*C3;
                C3'*S.H(ind_c,ind_b) + C3'*S.H(ind_c,ind_c)*C1, ...
                C3'*S.H(ind_c,ind_i)*T + C3'*S.H(ind_c,ind_c)*C2, ...
                C3'*S.H(ind_c,ind_c)*C3];
            fnew = [0.5*(S.H(ind_b,ind_c)+S.H(ind_c,ind_b)')*D1 + C1'*S.H(ind_c,ind_c)*D1 + S.f(ind_b) + C1'*S.f(ind_c) + ...
                0.5*(S.H(ind_i,ind_b)'+S.H(ind_b,ind_i))*t + 0.5*C1'*(S.H(ind_i,ind_c)'+S.H(ind_c,ind_i))*t;
                0.5*T'*(S.H(ind_i,ind_c)+S.H(ind_c,ind_i)')*D1 + C2'*S.H(ind_c,ind_c)*D1 + T'*S.f(ind_i) + C2'*S.f(ind_c) + ...
                T'*S.H(ind_i,ind_i)*t + 0.5*C2'*(S.H(ind_i,ind_c)'+S.H(ind_c,ind_i))*t;
                C3'*S.H(ind_c,ind_c)*D1 + C3'*S.f(ind_c) + 0.5*C3'*(S.H(ind_i,ind_c)'+S.H(ind_c,ind_i))*t];
        else
           % LP 
           Hnew = zeros(nb+nc);
           fnew = [S.f(ind_b) + C1'*S.f(ind_c);
               T'*S.f(ind_i) + C2'*S.f(ind_c);
               C3'*S.f(ind_c)];
        end
    end
    if S.isParametric
        if isempty(ind_i)
            % no integers
            if ~isempty(S.H)
                % QP
                pFnew = [ 0.5*(S.H(ind_b,ind_c)+S.H(ind_c,ind_b)')*D2 + C1'*S.H(ind_c,ind_c)*D2 + S.pF(ind_b,:) + C1'*S.pF(ind_c,:);
                    C3'*S.H(ind_c,ind_c)*D2 + C3'*S.pF(ind_c,:) ];
            else
                % LP
                pFnew = [ S.pF(ind_b,:) + C1'*S.pF(ind_c,:);
                    C3'*S.pF(ind_c,:) ];
            end
        else
            if ~isempty(S.H)
                % QP                
                pFnew = [ 0.5*(S.H(ind_b,ind_c)+S.H(ind_c,ind_b)')*D2 + C1'*S.H(ind_c,ind_c)*D2 + S.pF(ind_b,:) + C1'*S.pF(ind_c,:);
                    0.5*T'*(S.H(ind_i,ind_c)+S.H(ind_c,ind_i)')*D2 + C2'*S.H(ind_c,ind_c)*D2 + T'*S.pF(ind_i,:) + C2'*S.pF(ind_c,:);
                    C3'*S.H(ind_c,ind_c)*D2 + C3'*S.pF(ind_c,:) ];
            else
                % LP
                pFnew = [ S.pF(ind_b,:) + C1'*S.pF(ind_c,:);
                    T'*S.pF(ind_i,:) + C2'*S.pF(ind_c,:);
                    C3'*S.pF(ind_c,:) ];
            end
        end
    end
    
    
    % new constraints AnewB*[xB;yB] + AnewR*y <= bnew + Bnew*th
    alpha = [alpha1, alpha2];
    AnewB = alpha(Qv,:)-beta(Qv,:)*ibetaP*alpha(Pv,:);
    AnewR = -beta(Qv,:)*ibetaP;
    bnew  = gamma(Qv)-beta(Qv,:)*ibetaP*gamma(Pv);
    if S.isParametric
        pBnew = delta(Qv,:)-beta(Qv,:)*ibetaP*delta(Pv,:);
    end
    if S.me>0
        if isempty(ind_i)
            S.Ae = [eqlin.Ae(:,ind_b)+eqlin.Ae(:,ind_c)*C1, eqlin.Ae(:,ind_c)*C3, zeros(S.me,nlam)];
            S.be = eqlin.be -eqlin.Ae(:,ind_c)*D1;
        else            
            S.Ae = [eqlin.Ae(:,ind_b)+eqlin.Ae(:,ind_c)*C1, eqlin.Ae(:,ind_i)*T+eqlin.Ae(:,ind_c)*C2, eqlin.Ae(:,ind_c)*C3, zeros(S.me,nlam)];
            S.be = eqlin.be -eqlin.Ae(:,ind_i)*t -eqlin.Ae(:,ind_c)*D1;
        end
        if S.isParametric
            S.pE = eqlin.pE - eqlin.Ae(:,ind_c)*D2;
        end        
    end    
    
    % create M,q for LCP
    M = [-eye(nb) zeros(nb,nc+nlam);
        0.5*(Hnew(ind_bn,ind_cn)'+Hnew(ind_cn,ind_bn)) Hnew(ind_cn,ind_cn) AnewR';
        -AnewB -AnewR zeros(nlam)];
    q = [ones(nb,1); fnew(ind_cn); bnew];
    if S.isParametric
        Q = [zeros(nb,S.d); pFnew(ind_cn,:); pBnew];
    else
        Q = [];
    end
    
    
    % Build mapping from the lcp solution back to the original variables
    % xR = C1*xB + C2*yB + C3*y + D1 + D2*th
    % 
    % w  = [vB; vR; s], z = [[xB;yB]; y; lam]
    %
    % xR = [0 C1 C2 C3]*[w] + [D2 D1]*[th]
    %                   [z]           [1 ]
    % xB = [0 I  0  0 ]*[w] + [0  0 ]*[th]
    %                   [z]           [1 ]
    % yB = [0 0  I  0 ]*[w] + [0  0 ]*[th]
    %                   [z]           [1 ]
    
    % u = uX * [w;z] + uTh * [th;1]
    recover.uX = zeros(S.n, 2*(nb+nc+nlam));
    recover.uX(ind_b, nb+nc+nlam+1:nb+nc+nlam+nxb) = eye(nxb);
    if ~isempty(ind_i)
        recover.uX(ind_i, nb+nc+nlam+nxb+1:nb+nc+nlam+nxb+nyb) = T;
    end
    recover.uX(ind_c, nb+nc+nlam+1:nb+nc+nlam+nb+nc) = [C1 C2 C3];
   
    recover.uTh = zeros(S.n, S.d+1);
    if ~isempty(ind_i)
        recover.uTh(ind_i,S.d+1) = t;
    end
    if S.isParametric
        recover.uTh(ind_c,:) = [D2 D1];
    else
        recover.uTh(ind_c,:) = D1;
    end

    % Lagrange multipliers
    %  lambda_P = w(nb+1:nb+nc)  - corresponds to vR of vector w
    %  lambda_Q = z(nb+nc+1:nb+nc+nlam) - corresponds to lam in vector z
    %  lambda = lambdaX * [w;z] + lambdaTh * [th;1]

    % A'*lam = -H*x -[F f]*[th;1]
    % A'*lam = -H*(uX*[w;z]+uTh*[th;1])-[F f]*[th;1]
    % lam =  (-A'\H*uX)*[w;z] -A'\(H*uTh-[F f])*[th;1]

    lambdaX = zeros(S.m,2*(nb+nc+nlam));
    lambdaX(Pv,nb+1:nb+nc) = eye(length(Pv),nc);
    lambdaX(Qv,(nb+nc+nlam)+nb+nc+1:2*(nb+nc+nlam)) = eye(length(Qv),nlam);
    lambdaTh = zeros(S.m, S.d+1);

    % multipliers for the original inequalities    
    kept_rows.ineq = setdiff(1:S.Internal.m,S.Internal.removed_rows.ineqlin);    
    recover.lambda.ineqlin.lambdaX = zeros(S.Internal.m,2*(nb+nc+nlam));
    recover.lambda.ineqlin.lambdaX(kept_rows.ineq,:) = lambdaX(1:numel(kept_rows.ineq),:);
    recover.lambda.ineqlin.lambdaTh = zeros(S.Internal.m,S.d+1);
    
    % multipliers for the original lower bound
    kept_rows.lb = setdiff(1:S.Internal.n,S.Internal.removed_rows.lower);    
    recover.lambda.lower.lambdaX = zeros(S.Internal.n,2*(nb+nc+nlam));        
    recover.lambda.lower.lambdaX(kept_rows.lb,:) = ...
        lambdaX(numel(kept_rows.ineq)+1:numel(kept_rows.ineq)+length(kept_rows.lb),:);
    recover.lambda.lower.lambdaTh = zeros(S.Internal.n,S.d+1);

    % multipliers for the original upper bound
    kept_rows.ub = setdiff(1:S.Internal.n,S.Internal.removed_rows.upper);
    recover.lambda.upper.lambdaX = zeros(S.Internal.n,2*(nb+nc+nlam));        
    recover.lambda.upper.lambdaX(kept_rows.ub,:) = ...
        lambdaX(numel(kept_rows.ineq)+length(kept_rows.lb)+1:numel(kept_rows.ineq)+length(kept_rows.lb)+length(kept_rows.ub),:);
    recover.lambda.upper.lambdaTh = zeros(S.Internal.n,S.d+1);

        
           
end

if equations
    % Map from inequalities to equalities
    % u = uX * x + uTh * [th;1]
    % origU = recover.Y * u + recover.th * [th;1]
    % origU = recover.Y * [uX * x + uTh * [th;1]] + recover.th * [th;1]
    % origU = (recover.Y * uX) * x + (recover.Y * uTh + recover.th) * [th;1]
    recover.uX  = S.recover.Y*recover.uX;
    recover.uTh = S.recover.th + S.recover.Y*recover.uTh;
    
    % Lagrange multipliers for equalities
    % ni = -Ae'\H*x  -Ae'\A'*lam -Ae'\[F f]*[th; 1]
    % ni = -Ae'\H*(uX*[w;z]+uTh) -Ae'\'*(uLam*[w;z]+uTh)-Ae'\[F f]*[th;1]
    % ni = [(-Ae'\H)*Ux-(Ae'\A')*uLam]*[w;z] + 
    %      [(-Ae'\H)*uTh -(Ae'\A')*lamTh -(Ae'\[F f])]*[th; 1]

    AeH = -Ae'\Hm;
    AeA = -Ae'\Am';
    
    niX = AeH*recover.uX + AeA*lambdaX;
    if S.isParametric
        niTh = AeH*recover.uTh + AeA*lambdaTh -Ae'\[Fm fm];
    else
        niTh = AeH*recover.uTh + AeA*lambdaTh -Ae'\fm;
    end
    
    % lagrange multipliers for the original system
    kept_rows.eq = setdiff(1:S.Internal.me,S.Internal.removed_rows.eqlin);
    if r<nc
        recover.lambda.eqlin.lambdaX = zeros(S.Internal.me,2*(nb+2*nc+nlam));
        recover.lambda.eqlin.lambdaX(kept_rows.eq,:) = niX;
    else
        recover.lambda.eqlin.lambdaX = zeros(S.Internal.me,2*(nb+nc+nlam));
        recover.lambda.eqlin.lambdaX(kept_rows.eq,:) = niX;
    end
    recover.lambda.eqlin.lambdaTh = zeros(S.Internal.me,S.d+1);
    recover.lambda.eqlin.lambdaTh(kept_rows.eq,:) = niTh;
    
    
%     recover.lambdaX = [recover.lambdaX; niX];
%     recover.lambdaTh = [recover.lambdaTh; niTh];
else
    if r<nc
        recover.lambda.eqlin.lambdaX = zeros(0,2*(nb+2*nc+nlam));
    else
        recover.lambda.eqlin.lambdaX = zeros(0,2*(nb+nc+nlam));
    end
    recover.lambda.eqlin.lambdaTh = zeros(0,S.d+1);
end

% reorder variables is this came from Yalmip
if ~isempty(S.varOrder)
     Px = speye(size(recover.uX,1));
     Px = Px(S.varOrder.requested_variables,:);
     recover.primalX = Px*recover.uX;
     recover.primalTh= Px*recover.uTh;
end


% clear unnecessary fields
S.A  = []; S.b  = []; S.pB = [];
S.H  = []; S.f  = []; S.pF = [];
S.Y  = []; S.C  = []; S.c  = [];
S.lb = []; S.ub = [];

S.n  = size(M,1); % Problem dimension
S.m  = 0; % Number of inequalities
S.me = numel(S.be); % Number of equality constraints
% Number of parameters
if ~isempty(Q)
    S.d = size(Q,2);
else
    S.d  = 0; 
end
S.solver = '';
S.problem_type = ''; 
S.isParametric = [];

% set up LCP
S.M = M;
S.q = q;
if ~isempty(Q)
    S.Q = Q;
end
S.recover = recover;

% indicate which variables correspond to binaries
if ~isempty(S.vartype)
    if r<nc
        S.vartype = [repmat('B',1,nb), repmat('C',1,2*nc+nlam)];
    else
        S.vartype = [repmat('B',1,nb), repmat('C',1,nc+nlam)];
    end
end

% validate data
S.validate;

end
