function R = mpt_call_lcp(S)
%
%  MPT_CALL_LCP: A gateway function to LCP solver (without errorchecks) 
%  =====================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      R = mpt_call_lcp(S)
%    
%  
%  DESCRIPTION
%  -----------
%     The function implements calls to LCP solver based on the optimization problem
%  to be solved. Supported problems are LCP, LP, and QP. If the problem type is
%  LCP, then LCP mex-function is called directly. Otherwise a transformation to LCP
%  problem is performed. Assume that QP/LP is written in a form 
%                                        1  T    T             
%                                  min   - x Hx+f x        (1) 
%                                        2                     
%                                 s.t.   Ax <= b           (2) 
%                                                              
%                                        A x = b           (3) 
%                                         e     e              
%     We want to transform it to 
%                                       1  T      T            
%                                 min   - x H x+f  x       (4) 
%                                       2    f   f             
%                                                              
%                                s.t.   A x >= b           (5) 
%                                        f      f              
%                                       x >= 0             (6) 
%     which corresponds to the following LCP 
%                                   w - Mz  =   q       (7) 
%                                        w  >=  0       (8) 
%                                        z  >=  0       (9) 
%                                       T                   
%                                      w z  =   0      (10) 
%                                                           
%     where 
%                                       (       T  )          
%                                       ( H  -A    )          
%                                       (  f   f   )          
%                                 M  =  (          )     (11) 
%                                       ( A   0    )          
%                                       (  f       )          
%                                       (      )              
%                                       ( f    )              
%                                       (  f   )              
%                                 q  =  (      )         (12) 
%                                       ( -b   )              
%                                       (   f  )              
%     If LP or QP contains equality constraints, these are removed first. It is
%  required that the system of linear equations A_ex=b_e  is consistent, i.e. no
%  linearly dependent rows are found and the number of equalities is strictly less
%  than number of variables. The principle is based on factorizing equality
%  constraints A_ex=b_e  in basic  x_Bc and non-basic variables x_Nc, i.e. 
%                                    (              )
%                               A  = ( A     A      )
%                                e   (  e,Bc  e,Nc  )
%    which gives 
%                                                      
%                             A    x   + A    x   = b  
%                              e,Bc Bc    e,Nc Nc    e 
%     where the index sets Bc, Nc denote the columns from which factored system is
%  built. The factored submatrix A_e,Bc  must be invertible in order to express
%  basic variables as a function of non-basic variables, i.e. 
%                                   -1                -1      
%                       x   = -A      A    x   + A      b     
%                        Bc     e,Bc   e,Nc Nc    e,Bc   e,Bc 
%     With the substitution 
%                                           -1      
%                                 C = -A      A     
%                                       e,Bc   e,Nc 
%     and 
%                                          -1      
%                                 D = A      b     
%                                      e,Bc   e,Bc 
%    the relation between basic and non-basic variables is simplified to 
%                                                      
%                              x   = Cx   + D      (13)
%                               Bc     Nc              
%     The above QP/LP problem (??)-(??)  can be expressed only in non-basic
%  variables x_Nc  as follows: 
%                                  1    T      T                    
%                            min   - x   Hx  +f x   + c        (14) 
%                                  2  Nc   Nc    Nc                 
%                                                                   
%                           s.t.   Ax   <= b                   (15) 
%                                    Nc                             
%     where 
%                          T           T                                    
%                   H  =  C H     C + C H      + H     C + H           (16) 
%                            Bc,Bc       Bc,Nc    Nc,Bc     Nc,Nc           
%                          T           T            T       T               
%                   f  =  D H     C + D H      + f   C + f             (17) 
%                            Bc,Bc       Nc,Bc    Bc      Nc                
%                         1 T             T                                 
%                   c  =  -D H     D + f   D                           (18) 
%                         2   Bc,Bc     Bc                                  
%                                                                           
%                   A  =  A  C+A                                       (19) 
%                          Bc   Nc                                          
%                                                                           
%                   b  =  b - A  D                                     (20) 
%                              Bc                                           
%     Original solution to QP problem (??)-(??)  can be obtained via relation (??).
%  Problem (??)-(??)  can be transformed to (??)-(??)  effectively when considering
%  the rank of the matrix A. If the rank of matrix A  is less than the number of
%  inequalities in (??)  vector x_N  can be expressed as a difference of two
%  positive numbers, i.e. 
%                                            +      -           
%                                x    =   x    - x         (21) 
%                                 Nc       Nc     Nc            
%                                  +                            
%                               x     >=  0                (22) 
%                                Nc                             
%                                  -                            
%                               x     >=  0                (23) 
%                                Nc                             
%                                                          (24) 
%     With the substitution 
%                                    (    +  )        
%                                    ( x     )        
%                                    (  Nc   )        
%                                v = (    -  )    (25)
%                                    ( x     )        
%                                    (  Nc   )        
%     and putting back to (??)-(??)  we get 
%                           1  T (        )    (  T   T  )                 
%                     min   - v  ( H  -H  )v + ( f  -f   )v + c       (26) 
%                           2    ( -H H   )    (         )                 
%                           (       )                                      
%                    s.t.   ( -A A  )v >= -b                          (27) 
%                           (       )                                      
%                           v >= 0                                    (28) 
%     which is a form equivalent to (??)-(??). If the rank of matrix A  is greater
%  or equal than the number of inequalities in (??), we can factorize the matrix A 
%  rowwise 
%                                        (     )
%                                        ( A   )
%                                        (  B  )
%                                    A = (     )
%                                        ( A   )
%                                        (  N  )
%    where B, N  are index sets corresponding to rows from which submatrices are
%  built. The factored system (??)  can be written as 
%                                          >= -b             
%                                   -A x                (29) 
%                                     B Nc      B            
%                                          >= -b             
%                                   -A x                (30) 
%                                     N Nc      N            
%     where the matrix A_B  form by rows in the set B  must be invertible. Using
%  the substitution 
%                                     = -b              
%                              -A x         + y     (31)
%                                B Nc     B             
%     the system (??)(??)  can be rewritten in variable y 
%                     1  T (  -T   -1 )    (   -T   -1b    -T  )                 
%               min   - y  (A   HA    )y + (-A   HA     -A   f )y + c       (32) 
%                     2    ( B    B   )    (  B    B   B  B    )                 
%                         -1         -1b -b                                      
%              s.t.   A A   y >= A A                                        (33) 
%                      N B        N B   B  N                                     
%                     y >= 0                                                (34) 
%     which is equivalent to (??)-(??).
%  
%  INPUT
%  -----
%     
%        
%          S              Structure of the Opt class.              
%                         Class: struct                            
%          S.H            Quadratic part of the objective          
%                         function.                                
%                         Class: double                            
%                         Default: []                              
%          S.f            Linear part of the objective function.   
%                         Class: double                            
%          S.A            Linear part of the inequality            
%                         constraints Ax <= b.                     
%                         Class: double                            
%          S.b            Right hand side of the inequality        
%                         constraints Ax <= b.                     
%                         Class: double                            
%          S.Ae           Linear part of the equality constraints  
%                         A_ex=b_e.                                
%                         Class: double                            
%                         Default: []                              
%          S.be           Right hand side of the equality          
%                         constraints A_ex=b_e.                    
%                         Class: double                            
%                         Default: []                              
%          S.lb           Lower bound for the variables x >= lb.   
%                         Class: double                            
%                         Default: []                              
%          S.ub           Upper bound for the variables x <= ub.   
%                         Class: double                            
%                         Default: []                              
%          S.M            Data matrix for linear-complementarity   
%                         problem w-Mx=q.                          
%                         Class: double                            
%          S.q            Right hand side vector for               
%                         linear-complementarity problem w-Mx=q.   
%                         Class: double                            
%          S.n            Problem dimension (number of variables). 
%                                                                  
%                         Class: double                            
%          S.m            Number of inequalities in Ax <= b.       
%                         Class: double                            
%          S.me           Number of equalities in A_ex=b_e.        
%                         Class: double                            
%          S.problem_type A string specifying the problem to be    
%                         solved (only LP, QP and LCP problems are 
%                         allowed).                                
%                         Class: char                              
%          S.routine      An integer specifying which subroutine   
%                         to use for factorization. For details,   
%                         see help for LCP solver.                 
%                                                                  
%                            - The number 0 indicates to use fast  
%                            recursive factorization based on      
%                            rank-1 updates from LUMOD package.    
%                            - The number 1 indicates to use LU    
%                            factorization based on DGESV routine  
%                            of LAPACK.                            
%                            - The number 2 indicates to use QR    
%                            factorization based on DGELS routine  
%                            of LAPACK.                            
%                                                                  
%                         Class: double                            
%          S.test         Call (false) or not to call (true) MPT   
%                         global settings.                         
%                         Class: logical                           
%                         Default: false                           
%                           
%  
%  
%  OUTPUT
%  ------
%     
%        
%          R          Result structure                         
%                     Class: struct                            
%          R.xopt     The optimal values for variable z.       
%                     Class: double                            
%          R.lambda   The optimal values for variable w.       
%                     Class: double                            
%          R.obj      If the LCP problem was created by        
%                     conversion from LP/QP, this value        
%                     represent the optimal cost of the        
%                     appropriate LP/QP.                       
%                     Class: double                            
%          R.exitflag An integer value that informs if the     
%                     result was feasible (1), or otherwise    
%                     (different from 1).                      
%                     Class: double                            
%          R.how      A string that informs if the result was  
%                     feasible ('ok'), or if any problem       
%                     appeared through optimization.           
%                     Class: char                              
%                       
%  
%  
%  SEE ALSO
%  --------
%     mpt_solve
%  

%  AUTHOR(s)
%  ---------
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

if isempty(MPTOPTIONS) && ~S.test
    MPTOPTIONS = mptopt;
end

if ~any(strcmpi(S.problem_type,{'LP','QP','LCP'}))
    error('mpt_call_lcp: LCP solver does not solve %s problems!',S.problem_type);
end

% convert sparse matrices to full matrices
if issparse(S.H)
    S.H = full(S.H);
end
if issparse(S.f)
    S.f = full(S.f);
end
if issparse(S.A)
    S.A = full(S.A);
end
if issparse(S.b)
    S.b = full(S.b);
end
if issparse(S.Ae)
    S.Ae = full(S.Ae);
end
if issparse(S.be)
    S.be = full(S.be);
end

%% for LP, QP we need to transform to LCP form
if any(strcmpi(S.problem_type,{'LP','QP'}))
            
    % if no H is present
    if isempty(S.H)
        S.H = zeros(S.n);
    end
    H = S.H;
    f = S.f;
    
    % merge inequality constraints
    A=S.A;
    b=S.b;
    % detect Inf boundaries
	if ~S.test && (~isempty(S.lb) || ~isempty(S.ub))
        ilb = (S.lb==-Inf) | (S.lb<=-MPTOPTIONS.infbound);
        iub = (S.ub==Inf)  | (S.ub>=MPTOPTIONS.infbound);
		% store kept rows
		kept_rows.lb = find(~ilb);
		kept_rows.ub = find(~iub);
		if any(~ilb)
			% put ones at the positions where there is lb/ub
			Alb = zeros(nnz(~ilb),S.n);
			Alb(:,~ilb) = -eye(nnz(~ilb));
			A = [A; Alb];
			b = [b; -S.lb(~ilb)];
		end
		if any(~iub)
			Aub = zeros(nnz(~iub),S.n);
			Aub(:,~iub) = eye(nnz(~iub));
			A = [A; Aub];
			b = [b; S.ub(~iub)];
		end
	else
		kept_rows.lb = [];
		kept_rows.ub = [];
	end
    % store merged constraints
    Am = A;
    bm = b;
    
    % remove equality constraints if any
    kept_rows.eq = 1:S.me;
    Ae = S.Ae;
    be = S.be;
    me = S.me;
    if me>0
        % factorize Ae to get
        %  Ae(Br,Bc)*x(Bc) + Ae(Br,Nc)*x(Nc) = be(Br) % must be invertible mapping
            
        % check rank of Ae
        if S.test
            re = rank(Ae);
        else
            re = rank(Ae,MPTOPTIONS.rel_tol);
        end
        % check consistency
        if S.test
            rc = rank([Ae be]);
        else
            rc = rank([Ae be],MPTOPTIONS.rel_tol);
		end

		if re>=S.n
			% overdetermined system, no degrees of freedom
			R.xopt = S.Ae\S.be;
			
			% check constraints
			eqc = ( norm(S.Ae*R.xopt-S.be,Inf) <= MPTOPTIONS.rel_tol );
			if ~isempty(S.b)
				ineqc = ( S.A*R.xopt <= S.b + MPTOPTIONS.rel_tol );
			else
				ineqc = true;
			end
			
			% multipliers
			% H*x + f + A'*lambda_ineq + Ae'*lambda_eq = 0
			R.lambda.ineqlin = zeros(S.m,1);
			R.lambda.eqlin = -transpose(S.Ae)\(H*R.xopt + S.f);
			R.lambda.lower = zeros(S.n,1);
			R.lambda.upper = zeros(S.n,1);
			
			if eqc && all(ineqc)
				% feasible solution
				if S.test
					R.exitflag = 1;
				else
					R.exitflag = MPTOPTIONS.OK;
				end
				R.how = 'ok';
			else
				% infeasible
				if S.test
					R.exitflag = 2;
				else
					R.exitflag = MPTOPTIONS.INFEASIBLE;
				end
				R.how = 'infeasible';
			end
			R.obj = 0.5*R.xopt'*H*R.xopt + S.f'*R.xopt + S.c;
			return
		end
		
        % for underdetermined system check linearly dependent rows
        if re<me || rc<me
            
            % if the right hand side is not linearly dependent, infeasible solution
            if rc>re
                R.xopt = zeros(S.n,1);
                R.obj = 0;                
                R.lambda.ineqlin = [];
                R.lambda.eqlin = [];
                R.lambda.lower = [];
                R.lambda.upper = [];
                if S.test
                    R.exitflag = 2;
                else
                    R.exitflag = MPTOPTIONS.INFEASIBLE;
                end
                R.how = 'infeasible';
                return
            end
            
            while me ~= re
                % find linearly dependent rows
                [~,~,p] = lu(sparse([Ae be]),'vector');
                rd = p(re+1:end);
                
                % remove linearly dependent rows
                Ae(rd,:) = [];
                be(rd) = [];
                me = me-length(rd);
                kept_rows.eq(rd) = [];
                
                if S.test
                    re = rank(full([Ae be]));
                else
                    re = rank(full([Ae be]),MPTOPTIONS.rel_tol);
                end
            end
            
        end
        
        if ~isempty(Ae)
            
            [Le,Ue,pe,qe] = lu(sparse(Ae),'vector');
            
            if S.test
                rk = rank(full(Ue(:,1:re)));
            else
                rk = rank(full(Ue(:,1:re)),MPTOPTIONS.abs_tol);
            end
            if rk~=re
                % if invertibility is not achieved try reduced echelon
                % elimination
                
                % find linear independent columns of Ae
                if S.test
                    [~,jb]=rref(Ae,1e-4);
                else
                    [~,jb]=rref(Ae,MPTOPTIONS.abs_tol);
                end
                [Le,Ue,pe] = lu(sparse(Ae(:,jb)),'vector');
                % update column selection
                qe = [jb, setdiff(1:S.n,jb)];
                
                % here check rank with default tolerance ->
                % test_opt_eliminateEquations_12_pass
                rk = rank(full(Ue(:,1:re)));
                if rk~=re
                    % all possibilities failed, don't know what to do more, for now just throw error
                    error('mpt_call_lcp: Could not find invertible submatrix for removing equalities.');
                end
                
            end
            Br = pe(1:re); Bc = qe(1:re);
            Nr = pe(re+1:end); Nc = qe(re+1:end);
            % Nr must be empty -> otherwise there are depentent rows in Ae
            % which must be removed
            
            % substitute x(Bc) = C*x(Nc) + D
            Aebn = Ae(Br,Nc);
            %iAebb = inv(Ae(Br,Bc));
            beb = be(Br);
            
            % use factorized solution to compute C
            % C = -Ae(Br,Bc)\Aebn;
            Cl = -linsolve(full(Le),Aebn,struct('LT',true));
            C = linsolve(full(Ue(:,1:re)),Cl,struct('UT',true));
            
            % use factorized solution to compute D
            % D = Ae(Br,Bc)\beb;
            Dl = linsolve(full(Le),beb,struct('LT',true));
            D = linsolve(full(Ue(:,1:re)),Dl,struct('UT',true));
            
            Abc = A(:,Bc); Anc = A(:,Nc);
            
            % modify inequality constraints
            %A = -Abc*iAebb*Aebn + Anc;
            A = Abc*C + Anc;
            %b = b - Abc*iAebb*beb;
            if ~isempty(b)
                b = b - Abc*D;
            else
                b = [];
            end
            
            % modify cost
            %H = S.H(Nc,Nc) + Aebn'*iAebb'*S.H(Bc,Bc)*iAebb*Aebn - Aebn'*iAebb'*S.H(Bc,Nc) - S.H(Nc,Bc)*iAebb*Aebn;
            H = S.H(Nc,Nc) + C'*S.H(Bc,Bc)*C + C'*S.H(Bc,Nc) + S.H(Nc,Bc)*C;
            %f = S.H(Nc,Bc)*iAebb*beb - Aebn'*iAebb'*S.f(Bc) + S.f(Nc) - Aebn'*iAebb'*S.H(Bc,Bc)*iAebb*beb;
            f = S.H(Nc,Bc)*D + C'*S.f(Bc) + S.f(Nc) + C'*S.H(Bc,Bc)*D;
            
        else
           % matrix Ae has been completely eliminated 
           S.me = 0;
           S.Ae = zeros(0,S.n);
           S.be = [];
        end
    end
    
    % actual dimensions
    if S.test
        r = rank(A);
    else
        r = rank(A,MPTOPTIONS.abs_tol);
    end
    %m = size(A,1);
    n = size(A,2);
    if r<n
        % express x as a difference of two positive numbers, i.e. x = x+ - x-
        
        % new objective function
        Hnew = [H -H; -H H];
        fnew = [f; -f];
        
        % new constraints Anew*y >= bnew
        Anew = [-A A];
        bnew = -b;
        
        % catch error when allocating on extra large dimensions
        try
            % create M,q for LCP
            M = [Hnew -Anew'; Anew zeros(size(Anew,1))];
            q = [fnew; -bnew];
            
            % solve LCP
            % dimension of LCP to solve is "2*size(H,1)+size(A,1)"
            if S.test
                [z,w,basis,exfl] = lcp(M, q);
            else
                if ~isempty(S.routine)
                    lcpopt = MPTOPTIONS.modules.solvers.lcp;
                    lcpopt.routine = S.routine;
                    [z,w,basis,exfl] = lcp(M, q, lcpopt );
                else
                    [z,w,basis,exfl] = lcp(M, q, MPTOPTIONS.modules.solvers.lcp );
                end
            end
        catch M1
            % display what happened wrong
            disp(M1.message);            
            z = zeros(size(H,1)+size(Anew,1),1);
            w = z;
            % return error status 
            exfl = -4;                      
        end
        
        % recover solution
        R.xopt = z(1:n) - z(n+1:2*n);

        % recover mulpliers
        lambda_ineq = z(2*n+1:end);
        R.lambda.ineqlin = lambda_ineq(1:S.m);
        if ~isempty(S.lb)
            R.lambda.lower = zeros(S.n,1);
            R.lambda.lower(kept_rows.lb) = lambda_ineq(S.m+1:S.m+numel(kept_rows.lb));
        else
            R.lambda.lower = zeros(S.n,1);
        end
        if ~isempty(S.ub) && isempty(S.lb)
            R.lambda.upper = zeros(S.n,1);
            R.lambda.upper(kept_rows.ub) = lambda_ineq(S.m+1:S.m+numel(kept_rows.ub));
        elseif ~isempty(S.ub) && ~isempty(S.lb)
            R.lambda.upper = zeros(S.n,1);
            R.lambda.upper(kept_rows.ub) = lambda_ineq(S.m+numel(kept_rows.lb)+1:S.m+numel(kept_rows.lb)+numel(kept_rows.ub));
        else
            R.lambda.upper = zeros(S.n,1);
        end
        
        
    else        
        % factorize A to get
        %  -A(B,:)*x = -b(B) + y  % must be invertible mapping
        %  -A(N,:)*x >= -b(N)
        %          y >= 0
        [L,U,p] = lu(sparse(-A),'vector');
        B = p(1:n);
        N = p(n+1:end);
        
        % substitute
        % use factorized solution to compute inv(A(B,:))
        iAbl = -linsolve(full(L(1:n,:)),eye(n),struct('LT',true));
        iAb = linsolve(full(U),iAbl,struct('UT',true));
        %iAb = inv(A(B,:)); 
        bb = b(B);
        An = A(N,:);
        bn = b(N);
        
        % form new objective function
        Hnew = iAb'*H*iAb;
        fnew = -iAb'*H*iAb*bb - iAb'*f;
        
        % new constraints Anew* y>= bnew
        Anew = An*iAb;
        if isempty(bn)
            bnew = [];
        else
            bnew = Anew*bb - bn;
        end
        
        % catch error when allocating on extra large dimensions
        try
            % create M,q for LCP
            M = [Hnew -Anew'; Anew zeros(size(Anew,1))];
            q = [fnew; -bnew];
                    
            % solve LCP
            % dimension of LCP to solve is "size(H,1)+size(Anew,1)"
            if S.test
                [z,w,basis,exfl] = lcp(M, q);
            else
                if ~isempty(S.routine)
                    lcpopt = MPTOPTIONS.modules.solvers.lcp;
                    lcpopt.routine = S.routine;
                    [z,w,basis,exfl] = lcp(M, q, lcpopt );
                else
                    [z,w,basis,exfl] = lcp(M, q, MPTOPTIONS.modules.solvers.lcp );
                end
            end
        catch M2
            % display what happened wrong
            disp(M2.message);
            z = zeros(size(H,1)+size(Anew,1),1);
            w = z;
            % return error status 
            exfl = -4;
        end
        
        % recover variables
        xopt = z(1:n);
        lambda_ineq = zeros(length(bm),1);
        lambda_ineq(B) = w(1:n);
        lambda_ineq(N) = z(n+1:end);        
                
        % recover solution
        R.xopt = iAb*(bb - xopt);
        % recover mulpliers
        R.lambda.ineqlin = lambda_ineq(1:S.m);
        if ~isempty(S.lb)
            R.lambda.lower = zeros(S.n,1);
            R.lambda.lower(kept_rows.lb) = lambda_ineq(S.m+1:S.m+numel(kept_rows.lb));
        else
            R.lambda.lower = zeros(S.n,1);
        end
        if ~isempty(S.ub) && isempty(S.lb)
            R.lambda.upper = zeros(S.n,1);
            R.lambda.upper(kept_rows.ub) = lambda_ineq(S.m+1:S.m+numel(kept_rows.ub));
        elseif ~isempty(S.ub) && ~isempty(S.lb)
            R.lambda.upper = zeros(S.n,1);
            R.lambda.upper(kept_rows.ub) = lambda_ineq(S.m+numel(kept_rows.lb)+1:S.m+numel(kept_rows.lb)+numel(kept_rows.ub));
        else
            R.lambda.upper = zeros(S.n,1);
        end

        
    end
    
    
    % if equalities were present, map back to original variables
    if S.me>0
        xopt = zeros(S.n,1);
        xopt(Nc) = R.xopt;
        %xopt(Bc) = iAebb*(-Aebn*R.xopt + beb);
        xopt(Bc) = C*R.xopt + D;
        R.xopt = xopt;
        % solve overdetermined system to get multipliers for equalities
        % H*x + f + Am'*lambda_ineq + Ae'*lambda_eq = 0
        lambda_eq = zeros(S.me,1);
        lambda_eq(kept_rows.eq) = -Ae'\(S.H*R.xopt + S.f + Am'*lambda_ineq);
        
        % extend multipliers
        R.lambda.eqlin = lambda_eq;
    else
        R.lambda.eqlin = []; 
    end
            
    
else
    %% solve LCP directly
    if S.test
        [z,w,basis,exfl] = lcp(S.M, S.q);
    else        
        [z,w,basis,exfl] = lcp(S.M, S.q, MPTOPTIONS.modules.solvers.lcp );
    end

    R.xopt = z;
    R.lambda = w;
    
end

% this call of LCP solver is disabled, because the size of the problem is
% 2*n+m which is larger than with constraints factorization

%     if strcmpi(S.problem_type,'QP')
%         % if the problem is QP with H>0, then we can call LCP directly
%         v = eig(S.H);
%         if ~( any(abs(imag(v)) > MPTOPTIONS.abs_tol) || any(real(v) < -MPTOPTIONS.abs_tol) )
%             % Hessian is positive semidefinite
%             
%             % merge constraints
%             A = [-S.A; eye(S.n); -eye(S.n)];
%             b = [-S.b; S.lb; -S.ub];
%             
%             
%             if S.me==0
%                 % no equalities
%                 M = A*(S.H\A'); 
%                 q = -A*(S.H\S.f)-b;
%                 
%                 [z,w,basis,exfl] = lcp(full(M), full(q), MPTOPTIONS.modules.solvers.lcp );
%                 R.xopt = S.H\(A'*z - S.f);
%                 R.lamda = z(1:S.m);
%                 
%             else
%                 % equalities
%                 M = [A zeros(S.m+2*S.n,S.me)] * ( [S.H S.Ae'; -S.Ae zeros(S.me)] \ [A'; zeros(S.me,S.m+2*S.n)] );
%                 q = -[A zeros(S.m+2*S.n,S.me)] * ( [S.H S.Ae'; -S.Ae zeros(S.me)] \ [S.f; S.be] ) - b;
%                 
%                 [z,w,basis,exfl] = lcp(full(M), full(q), MPTOPTIONS.modules.solvers.lcp );
%                 xx = [S.H S.Ae'; -S.Ae zeros(S.me)] \ [A'*z - S.f; -S.be];
%                 R.xopt = xx(1:S.n);
%                 % merging multipliers on ineqality and equality constraints
%                 R.lambda = [z(1:S.m); xx(S.n+1:end)];
%                 
%             end
%             

% recalculate the objective function for LP, QP
R.obj = [];
if any(strcmpi(S.problem_type,{'LP','QP'}))
    if ~isempty(R.xopt) && ~isempty(S.H)
        R.obj = 0.5*R.xopt'*S.H*R.xopt + S.f'*R.xopt;
    elseif ~isempty(R.xopt) && ~isempty(S.f)
        R.obj = S.f'*R.xopt;
    end
end

switch exfl
    case 1
        R.how = 'ok';
        if S.test
            R.exitflag = 1;
        else
            R.exitflag = MPTOPTIONS.OK;
        end
    case -1
        R.how = 'infeasible';
        if S.test
            R.exitflag = 2;
        else
            R.exitflag = MPTOPTIONS.INFEASIBLE;
        end
    case -2
        R.how = 'unbounded';
        if S.test
            R.exitflag = 3;
        else
            R.exitflag = MPTOPTIONS.UNBOUNDED;
        end
    case -3
        R.how = 'preterminated (due to time limit or maximum pivot limit)';
        if S.test
            R.exitflag = -1;
        else
            R.exitflag = MPTOPTIONS.ERROR;
        end
    case -4
        R.how = 'other (numerical) error';
        if S.test
            R.exitflag = -1;
        else
            R.exitflag = MPTOPTIONS.ERROR;
        end
end

end
