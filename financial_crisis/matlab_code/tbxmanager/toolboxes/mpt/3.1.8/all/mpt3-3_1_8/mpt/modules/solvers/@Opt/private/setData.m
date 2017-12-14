function opt = setData(opt, varargin)
%
% Set the data for Opt constructor
%

ip = inputParser;
ip.KeepUnmatched = true;
ip.CaseSensitive = true;
ip.addOptional('P', [], @validate_polyhedron);
ip.addOptional('solver',  '', @ischar);
ip.addOptional('vartype', '', @validate_vartype);
ip.addParamValue('A',  [], @validate_realmatrix);
ip.addParamValue('b',  [], @validate_realvector);
ip.addParamValue('pB', [], @validate_realmatrix);
ip.addParamValue('Ae', [], @validate_realmatrix);
ip.addParamValue('be', [], @validate_realvector);
ip.addParamValue('pE', [], @validate_realmatrix);
ip.addParamValue('H',  [], @validate_realmatrix);
ip.addParamValue('f',  [], @validate_realvector);
ip.addParamValue('pF', [], @validate_realmatrix);
ip.addParamValue('lb', [], @validate_realinfvector);
ip.addParamValue('ub', [], @validate_realinfvector);
ip.addParamValue('Ath', [], @validate_realmatrix);
ip.addParamValue('bth', [], @validate_realvector);
ip.addParamValue('M',  [], @validate_realmatrix);
ip.addParamValue('q',  [], @validate_realvector);
ip.addParamValue('Q',  [], @validate_realmatrix);
ip.addParamValue('Y',  [], @validate_realmatrix);
ip.addParamValue('C',  [], @validate_realvector);
ip.addParamValue('c',  [], @validate_realvector);

ip.parse(varargin{:});
p = ip.Results;

% assign parsed outputs
opt.A  = p.A;  opt.b  = p.b;  opt.pB = p.pB;
opt.Ae = p.Ae; opt.be = p.be; opt.pE = p.pE;
opt.H  = p.H;  opt.f  = p.f;  opt.pF = p.pF;
opt.lb = p.lb; opt.ub = p.ub;
opt.Ath = p.Ath; opt.bth = p.bth;
opt.Y = p.Y; opt.C = p.C; opt.c = p.c;
opt.M  = p.M;  opt.q  = p.q;  opt.Q  = p.Q;

% if polyhedron is given, replace equalities and inequalitites
% with the one stored in a polyhedron
if ~builtin('isempty',p.P)
    if ~isempty(p.P.A)
        opt.A  = p.P.A;  opt.b  = p.P.b;
    end
    if ~isempty(p.P.Ae)
        opt.Ae = p.P.Ae; opt.be = p.P.be;
    end
end

% assign solver in uppercase
if ~isempty(p.solver)
    opt.solver = upper(strtrim(p.solver));
end

% assign vartype in uppercase
if ~isempty(p.vartype)
    opt.vartype = upper(p.vartype);
end

end
