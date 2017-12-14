function opt = setMPT26Data(opt, mat)
%
%  Convert from old MPT2.6 format
%

global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

% convert sparse matrices to full
f = fieldnames(mat);
for i = 1:length(f)
	if issparse(mat.(f{i}))
		mat.(f{i}) = full(mat.(f{i}));
	end
end

opt.A  = mat.G;
opt.b  = mat.W;
opt.pB = mat.E;

if isfield(mat,'Y')
    opt.Y = mat.Y;
end
if isfield(mat,'Cx')
    opt.C = mat.Cx;
end
if isfield(mat,'Cc')
    opt.c = mat.Cc;
end

if isfield(mat, 'Aeq')
    opt.Ae = mat.Aeq;
    opt.be = mat.beq;
end
if isfield(mat, 'Beq')
    opt.pE = -mat.Beq;
end

n = size(opt.A,2);
opt.Ath = [];
opt.bth = [];
if isfield(mat, 'lb')
    opt.lb = mat.lb(1:n); % First n variables are the decision vars
    mat.lb(1:n) = [];
    d = length(mat.lb); % Last d variables are the parameters
    opt.Ath = [opt.Ath;-eye(d)];
    opt.bth = [opt.bth;-mat.lb];
end
if isfield(mat, 'ub')
    opt.ub = mat.ub(1:n);
    mat.ub(1:n) = [];
    d = length(mat.ub); % Last d variables are the parameters
    opt.Ath = [opt.Ath;eye(d)];
    opt.bth = [opt.bth;mat.ub];
end

if ~isempty(mat.bndA)
    opt.Ath = [opt.Ath;mat.bndA];
    opt.bth = [opt.bth;mat.bndb];
end
if isfield(mat, 'requested_variables')
    opt.varOrder.requested_variables = mat.requested_variables;
end
if isfield(mat, 'param_var')
    opt.varOrder.param_var = mat.param_var;
end
if isfield(mat, 'free_var')
    opt.varOrder.free_var = mat.free_var;
end
if isfield(mat, 'binary_var_index')
    opt.varOrder.binary_var_index = mat.binary_var_index;
    opt.vartype = 'C'*ones(n,1);
    opt.vartype(opt.varOrder.binary_var_index) = 'B';
end

if ~isfield(mat,'qp')
    if size(mat.H,1) == size(mat.H,2)
        mat.qp = true;
        % if all elements of hessian zeros ->LP
        if all(all(abs(mat.H)<MPTOPTIONS.zero_tol))
            mat.qp = false;
        end
    else
        mat.qp = false;
    end
end

if mat.qp
    opt.H  = mat.H;
    opt.pF = mat.F';
    opt.f  = mat.Cf';
else
    d = 0;
    if isfield(opt, 'pB'),  d = max([d size(opt.pB,2)]); end
    if isfield(opt, 'pE'),  d = max([d size(opt.pE,2)]); end
    if isfield(opt, 'Ath'), d = max([d size(opt.Ath,2)]); end
    
    opt.H  = [];
    if isfield(mat, 'D')
        % parameterized cost th'*D*x (issue #117)
        opt.pF = mat.D;
    else
        opt.pF = zeros(n,d);
    end
    opt.f  = mat.H';
end

% assign solver in uppercase if present
if isfield(mat,'solver')
    opt.solver = upper(strtrim(mat.solver));
end


end
