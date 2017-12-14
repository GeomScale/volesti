function opt = validate(opt)
%
% Validation function for Opt class
%
%
global MPTOPTIONS

if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

% deal with arrays
if numel(opt)>1
    for i=1:numel(opt)
        opt(i).validate;
    end
    return
end


%% Validate LCP inputs
if ~isempty(opt.M)
    opt.n = size(opt.M,1);
    opt.d = max([size(opt.Q,2), size(opt.Ath,2)]);
    
    if size(opt.M,1) ~= size(opt.M,2)
        error('Opt: M must be square');
    end
    if ~isempty(opt.Q) 
        if size(opt.Q,1) ~= opt.n || size(opt.Q,2) ~= opt.d
            error('Opt: Q must be %i x %d\n', opt.n, opt.d);
        end
    end
    if isempty(opt.q)
        error('Opt: Vector q must be provided for LCP.');
    end
    opt.q = opt.q(:);
    if length(opt.q) ~= opt.n
        error('Opt: q must be %i x 1\n', opt.n);
    end
    
    if ~isempty(opt.Ath)
        if size(opt.Ath,2) ~= opt.d
            error('Opt: Ath must have %i columns', opt.d);
        end
    else
        opt.Ath = zeros(0, opt.d);
    end
    if ~isempty(opt.bth)
        opt.bth = opt.bth(:);
        if length(opt.bth) ~= size(opt.Ath,1)
            error('Opt: bth must have length %i', size(opt.Ath,1))
        end
    else
        opt.bth = zeros(0,1);
    end

    
    if opt.d > 0
        opt.isParametric = true;
%         if isempty(opt.Ath) || isempty(opt.bth)
%             error('Opt: Bounds on the parameters must be provided.');
%         end
%         Pbnd = Polyhedron(opt.Ath,opt.bth);
%         if ~Pbnd.isBounded
%             error('Opt: Parameter space is not bounded, please provide bounds on all parameters.');
%         end
    else
        opt.isParametric = false;
    end

    % set as LCP problem
    opt.problem_type = 'LCP';
else
%% Validate LP/QP/MILP/MIQP inputs
    % Get sizes
    opt.n  = max([size(opt.H,1),length(opt.f),size(opt.A,2),size(opt.Ae,2),length(opt.lb),length(opt.ub)]);
    opt.d = max([size(opt.pF,2) size(opt.pB,2) size(opt.pE,2) size(opt.Ath,2)]);
    opt.m  = size(opt.A,1);
    opt.me = size(opt.Ae,1);
    
    % Validate parametric inputs
    if ~isempty(opt.pF) || ~isempty(opt.pB) || ~isempty(opt.pE)        
        opt.isParametric = true;
    else
        opt.isParametric = false;
    end   
    
    if isempty(opt.A) && isempty(opt.Ae) && isempty(opt.lb) && isempty(opt.ub) && ~opt.isParametric
        if isempty(opt.H) && isempty(opt.f)
            error('Opt: Empty problem');
        elseif isempty(opt.H) && ~isempty(opt.f)
            error('Opt: Unconstrained problem.');
        end
    end
    
    % validate vartype
    if ~isempty(opt.vartype)
        if isnumeric(opt.vartype)
            % convert to char if it is numeric
            opt.vartype = char(opt.vartype);
        end
        if ~isvector(opt.vartype) || ~ischar(opt.vartype)
            error('The argument must be a vector of strings.');
        end
        if length(opt.vartype)~=opt.n
            error('The field "vartype" must be a vector of dimension %i,',opt.n);
        end
        % checking if string is correct
        for i=1:length(opt.vartype)
            if ~any(strcmpi(opt.vartype(i),{'C','I','B','S','N'}))
                %C-continuous, I-integer, B-binary, S-semicontinuous, N-semiinteger
                error('Given string does not belong to group "C-continuous, I-integer, B-binary, S-semicontinuouos, N-semiinteger.');
            end
        end
        % must be a row vector
        if size(opt.vartype,1)>size(opt.vartype,2)
            opt.vartype = opt.vartype';
        end
    end

    
    % if all elements of hessian zeros, make H empty
    if all(all(abs(opt.H)<MPTOPTIONS.zero_tol))
        opt.H = [];
    end
    
    % find any binary/integer/semiinteger variables
    c = regexpi(opt.vartype,'(B|I|N)');
    
    if ~isempty(opt.f)  
        opt.f = opt.f(:);
        if length(opt.f) ~= opt.n
            error('Opt: Linear cost function f must be %i x 1\n', opt.n); 
        end
        if isempty(c)
            opt.problem_type = 'LP';
        else
            opt.problem_type = 'MILP';
        end
    else
        opt.f = zeros(opt.n,1);
        if isempty(opt.H)
            if isempty(c)
                opt.problem_type = 'LP';
            else
                opt.problem_type = 'MILP';
            end            
        end
    end
    
    if ~isempty(opt.H)
        if size(opt.H,1) ~= size(opt.H,2) || size(opt.H,1) ~= opt.n
            error('Opt: Hessian H must be %i x %i\n', opt.n, opt.n);
        end
        v = eig(opt.H);
        if any(real(v) < -MPTOPTIONS.abs_tol)
            %if any(abs(imag(v)) > MPTOPTIONS.abs_tol) ||any(real(v) < -MPTOPTIONS.abs_tol)
            % test on imaginary parts of eigenvalues caused matrix to loose
            % positive-semidefiniteness by elimitating equations
            % from an optimalization problem
            error('Opt: Hessian is not positive semidefinite');
        end
        % Make hessian symmetric
        opt.H = 0.5 * (opt.H + opt.H');
        %opt.type = opt.OPT_QP;
        if isempty(c)
            opt.problem_type = 'QP';
        else
            opt.problem_type = 'MIQP';
        end
    else
        opt.H = [];
    end
    
    if ~isempty(opt.A)
        if size(opt.A,2) ~= opt.n
            error('Opt: A must have %i columns\n', opt.n);
        end
        opt.b = opt.b(:);
        if length(opt.b) ~= opt.m,
            error('Opt: b must be %i x 1\n', opt.m); 
        end
    else
        opt.A = zeros(opt.m, opt.n);
    end
    if isempty(opt.b)
       opt.b = zeros(opt.m, 1);
    end
    
    opt.b = opt.b(:);
    if length(opt.b) ~= opt.m
        error('Opt: RHS b must be %i x 1\n', opt.m); 
    end
    
    if ~isempty(opt.Ae)
        if size(opt.Ae,2) ~= opt.n
            error('Opt: Ae must have %i columns\n', opt.n);
        end
        opt.be = opt.be(:);
        if length(opt.be) ~= opt.me
            error('Opt: be must be %i x 1\n', opt.me);
        end
    else
        opt.Ae = zeros(opt.me, opt.n);
    end
    if isempty(opt.be)
        opt.be = zeros(opt.me, 1);
    end
    
    if ~isempty(opt.lb)
        opt.lb = opt.lb(:);
        if length(opt.lb) ~= opt.n
            error('Opt: lb must be %i x 1\n', opt.n);
        end
    else
        opt.lb = -inf(opt.n,1);
    end
    
    if ~isempty(opt.ub)
        opt.ub = opt.ub(:);
        if length(opt.ub) ~= opt.n
            error('Opt: ub must be %i x 1\n', opt.n);
        end
    else
        opt.ub = inf(opt.n,1);
    end
    
    % check if lb is actually lower that ub
    for i=1:opt.n
        if opt.lb(i)>opt.ub(i)
            error('Lower bound %i is higher than its upper bound.',i);
        end
    end
    
    if ~isempty(opt.c)
        if numel(opt.c)~=1
            error('Opt: The constant term "D" must be scalar.');
        end
    else        
        opt.c = 0;
    end
        
    if opt.isParametric
        if ~isempty(opt.pF)
            if size(opt.pF,2) ~= opt.d || size(opt.pF,1) ~= opt.n
                error('Opt: pF must be %i x %i', opt.n, opt.d);
            end
        else
            opt.pF = zeros(opt.n,opt.d);
        end
        
        if ~isempty(opt.pB)
            if size(opt.pB,2) ~= opt.d || size(opt.pB,1) ~= opt.m
                error('Opt: pB must be %i x %i', opt.m, opt.d);
            end
        else
            opt.pB = zeros(opt.m,opt.d);
        end
        
        if ~isempty(opt.pE)
            if size(opt.pE,2) ~= opt.d || size(opt.pE,1) ~= opt.me
                error('Opt: pE must be %i x %i', opt.me, opt.d);
            end
        else
            opt.pE = zeros(opt.me,opt.d);
        end
        
        if ~isempty(opt.Y)
            if size(opt.Y,1) ~= opt.d || size(opt.Y,2) ~=opt.d
                error('Opt: The matrix "Y" must be square of dimension %d.',opt.d);
            end
        else
            opt.Y = zeros(opt.d);
        end
        if ~isempty(opt.C)
            if size(opt.C,1) ~= 1 || size(opt.C,2) ~=opt.d
                error('Opt: The vector "C" must of dimension 1 x %d.',opt.d);
            end
        else
            opt.C = zeros(1,opt.d);
        end
        
        
    end        
    
    if ~isempty(opt.Ath)
        if size(opt.Ath,2) ~= opt.d
            error('Opt: Ath must have %i columns', opt.d);
        end
    else
        opt.Ath = zeros(0, opt.d);
    end
    if ~isempty(opt.bth)
        opt.bth = opt.bth(:);
        if length(opt.bth) ~= size(opt.Ath,1)
            error('Opt: bth must have length %i', size(opt.Ath,1))
        end
    else
        opt.bth = zeros(0,1);
    end

%     if ~opt.isParametric
%         % remove zero rows from inequality constraints
%         zr = all(abs(opt.A)<MPTOPTIONS.zero_tol, 2);
%         opt.A(zr,:) = [];
%         opt.b(zr) = [];
%         opt.m = opt.m-nnz(zr);
%         
%         % remove zero rows from equality constraints
%         zre = all(abs(opt.Ae)<MPTOPTIONS.zero_tol, 2);
%         opt.Ae(zre,:) = [];
%         opt.be(zre) = [];
%         opt.me = opt.me-nnz(zre);
%     end   
    
    % check consistency of equality constraints
    % check rank of equality constraints
    re = rank(full(opt.Ae));
    
    % overdetermined system, no degrees of freedom
    if re>opt.n
        error('Opt: Overdetermined system, no degrees of freedom.');
    end
    
    % Try to avoid Inf terms here because it makes numerical troubles
    ilb = (opt.lb==-Inf) | (opt.lb<-MPTOPTIONS.infbound);
    iub = (opt.ub==Inf)  | (opt.ub>MPTOPTIONS.infbound);
    Hn = [-opt.pB opt.A opt.b];
    
    % merge constraints with lb/ub
    if any(~ilb)
        % put ones at the positions where there is lb/ub
        Zlb = [zeros(nnz(~ilb),opt.d+opt.n), -opt.lb(~ilb)];
        Zlb(:,[false(opt.d,1); ~ilb; false]) = -eye(nnz(~ilb));
        Hn = [Hn; Zlb];
    end
    if any(~iub)
        Zub = [zeros(nnz(~iub),opt.d+opt.n), opt.ub(~iub)];
        Zub(:,[false(opt.d,1); ~iub; false]) = eye(nnz(~iub));
        Hn = [Hn; Zub];
    end
    Hn = [Hn; opt.Ath zeros(size(opt.Ath,1),opt.n) opt.bth];
    
    % store internally merged constraints % used in eliminateEquations and qp2lcp
    opt.Internal.pB = -Hn(:,1:opt.d);
    opt.Internal.A = Hn(:,opt.d+1:opt.n+opt.d);
    opt.Internal.b  =  Hn(:,end);

    % also store dimensions and cost function
    opt.Internal.me = opt.me; % store the number of original equalities
    opt.Internal.m = opt.m; % store the number of original inequalities
    opt.Internal.n = opt.n; % store the number of original variables
    
    % store removed rows for proper postprocessing
    opt.Internal.removed_rows.lower = find(ilb);
    opt.Internal.removed_rows.upper = find(iub);
    opt.Internal.removed_rows.ineqlin = [];
    opt.Internal.removed_rows.eqlin = [];
    
    % Preprocess the parametric inputs
    if opt.isParametric
        % Change to full-matrix representation
        vars = {'H','f','A','b','Ae','be','lb','ub','pF','pB','pE','Ath','bth','Y','C','c','M','Q','q'};
        for i=1:length(vars)
            if issparse(opt.(vars{i}))
                opt.(vars{i}) = full(opt.(vars{i}));
            end
        end
        
        % Compute tight upper and lower bounds on all variables and
        % parameters
        %P = Polyhedron('H', Hn, 'He', [-opt.pE opt.Ae opt.be]);
        P = Polyhedron(Hn(:, 1:end-1), Hn(:, end));
		
		if P.isEmptySet
			% exit quickly if problem is infeasible
			opt.Internal.pB = -P.H(:,1:opt.d);
			opt.Internal.A = P.H(:,opt.d+1:opt.n+opt.d);
			opt.Internal.b  =  P.H(:,end);
			validate_solvernames(opt);
			return
		end
		
		[P, hull] = P.minHRep();
		P.outerApprox;
		lb = P.Internal.lb;
		ub = P.Internal.ub;
        
        % Parameter bounds
        if any(lb(1:opt.d)==Inf) || any(ub(1:opt.d)==-Inf),
            % P was computed without considering equalities, retry with
            % equalities
            P = Polyhedron('H', Hn, 'He', [-opt.pE opt.Ae opt.be]);
            [P, hull] = P.minHRep();
            P.outerApprox;
            lb = P.Internal.lb;
            ub = P.Internal.ub;
            % do not throw error if infeasible
            % if any(lb(1:opt.d)==Inf) || any(ub(1:opt.d)==-Inf)
            %    error('Opt: Parameter space is empty'); 
            % end
        end
        plb = lb(1:opt.d); iplb = (plb==-Inf) | (lb(1:opt.d)<-MPTOPTIONS.infbound);
        pub = ub(1:opt.d); ipub = (pub==Inf) | (ub(1:opt.d)>MPTOPTIONS.infbound);
        if any(iplb) || any(ipub)
            %warning('Opt:validate','Parameter space is unbounded');
            % replace INF values with infbound
            plb(iplb) = -MPTOPTIONS.infbound;
            pub(ipub) = MPTOPTIONS.infbound;            
        end

        % overwrited internally irredundant representation with parametric
        % data
        opt.Internal.pB = -hull.H(:,1:opt.d); % merged inequality constraints on decision and parametric variables A*x <= b + pB*th
        opt.Internal.A = hull.H(:,opt.d+1:opt.n+opt.d);
        opt.Internal.b  =  hull.H(:,end);
        opt.Internal.plb = plb; % extracted lower bound on parameters
        opt.Internal.pub = pub; % extracted upper bound on parameters
        opt.Internal.lb = lb(opt.d+1:opt.d+opt.n); % extracted lower bound on decision variables
        opt.Internal.ub = ub(opt.d+1:opt.d+opt.n); % extracted upper bound on decision variables        
        opt.Internal.removed_rows.ineqlin = find(hull.I(1:opt.m));
        remained.lb = setdiff(transpose(1:opt.n),opt.Internal.removed_rows.lower); % indices that have remained after elimination of -Inf lower bound
        remained.ub = setdiff(transpose(1:opt.n),opt.Internal.removed_rows.upper); % indices that have remained after elimination of Inf upper bound
        % append the indices that have been eliminated from the remaining lb/ub constraints
        opt.Internal.removed_rows.lower = [opt.Internal.removed_rows.lower; remained.lb(hull.I(opt.m+1:opt.m+nnz(~ilb)))];
        opt.Internal.removed_rows.upper = [opt.Internal.removed_rows.upper; remained.ub(hull.I(opt.m+nnz(~ilb)+1:opt.m+nnz(~ilb)+nnz(~iub)))];
        
        
        % Inequalities of P that only involve the parameter
        H = hull.H;
        I = sum(abs(H(:,opt.d+1:end-1)) > MPTOPTIONS.zero_tol, 2) == 0;
        Hth = [H(I,1:opt.d) H(I,end)];
        hull.H(I,:) = [];
        
        if any(lb(1:opt.d)==Inf) || any(ub(1:opt.d)==-Inf)
            opt.Ath = zeros(0,opt.d);
            opt.bth = [];        
        else
            thHull = Polyhedron('H',[opt.Ath opt.bth;Hth],'lb',plb,'ub',pub).minHRep();
            opt.Ath = thHull.H(:,1:end-1);
            opt.bth = thHull.H(:,end);
        end
        
%         % Variable bounds
%         opt.lb = lb(opt.d+1:end);
%         opt.ub = ub(opt.d+1:end);
                        
%         % Remove inequalities that only involve one variable (these are now in lb/ub)
%         H = hull.H;
%         I = sum(abs(H(:,1:end-1)) > MPTOPTIONS.zero_tol, 2) == 1;
%         H(I,:) = [];
%         opt.A  =  H(:,opt.d+1:end-1);
%         opt.pB = -H(:,1:opt.d);
%         opt.b  =  H(:,end);
%         opt.Internal.removed_rows.ineqlin = [opt.Internal.removed_rows.ineqlin; find(hull.I)];
        
%         opt.Ae =  hull.He(:,opt.d+1:end-1);
%         opt.pE = -hull.He(:,1:opt.d);
%         opt.be =  hull.He(:,end);

    else
        % Remaining tests for non-parametric problems
        
                        
        % for underdetermined system check linearly dependent rows
        if re<opt.me
            % check rank consistency
            rc = rank(full([opt.Ae opt.be]));
            
            % if the right hand side is not linearly dependent, infeasible solution
            if rc>re
                error('Opt: Equality constraints are not consistent.')
            end
%             % otherwise remove linearly dependent equalities
%             disp('Opt: Removing linearly dependent equality constraints.');
%             
%             % find linearly dependent rows
%             [~,~,pp] = lu(sparse(opt.Ae),'vector');
%             rd = pp(re+1:end);
%             
%             % remove linearly dependent rows
%             opt.Ae(rd,:) = [];
%             opt.be(rd) = [];
%             opt.me = opt.me-length(rd);
        end
        
    end
    
    % find affine map between integer and binary variables
    % x_I = T*y_B + t
    % where x_I are integers and y_B are binary variables and T, t are
    % matrices
    ind_i = find(opt.vartype=='I');
    ni = numel(ind_i);
    if ni>0
        % truncate zeros
        if ~opt.isParametric
            P = Polyhedron('H', Hn, 'He', [opt.Ae opt.be]);
            P.outerApprox;
            opt.Internal.lb = P.Internal.lb;
            opt.Internal.ub = P.Internal.ub;
        end
        lbi = opt.Internal.lb(ind_i);
        ubi = opt.Internal.ub(ind_i);
        izlbi = abs(lbi) < MPTOPTIONS.zero_tol;
        izubi = abs(ubi) < MPTOPTIONS.zero_tol;
        lbi(izlbi) = zeros(nnz(izlbi),1);
        ubi(izubi) = zeros(nnz(izubi),1);
        
        % lower and upper bounds on integer variables
        lbi = sign(lbi).*floor(abs(lbi));
        ubi = sign(ubi).*floor(abs(ubi));
        % correct inf values
        % using infbound/2 value because intmax/2 or infbound gives
        % large values in the transformation matrix which causes
        % problems for mixed-integer solvers (see
        % test_opt_integer_02_pass)
        lbi(isinf(lbi)) = -MPTOPTIONS.infbound/2;
        ubi(isinf(ubi)) = MPTOPTIONS.infbound/2;
        
        % find negative lb that is the offset
        t = zeros(ni,1);
        indt = (lbi<0);
        t(indt) = lbi(indt);
        
        % how many binary variables are needed to represent each
        % integer?
        imax = ubi-t;
        nb = zeros(ni,1);
        for i=1:ni
            if imax(i)>intmax
                if MPTOPTIONS.verbose>1
                    fprintf('The computed upper bound on the integer variable %i exceeds the limit %i. Using the upper limit instead.\n',i,intmax);
                end
                imax(i) = intmax;
            end
            nb(i) = numel(dec2bin(imax(i)));
        end
        % number of auxiliary binary variables to introduce
        sb = sum(nb);
        
        % prepare the matrix T
        T = zeros(ni,sb);
        stmp = 0;
        for i=1:ni
            T(i,stmp+1:stmp+nb(i)) = power(2,nb(i)-1:-1:0);
            stmp = stmp+nb(i);
        end
        
        % store internally
        opt.Internal.T = T;
        opt.Internal.t = t;
    end

    
end

% validate solvernames
validate_solvernames(opt);

end
