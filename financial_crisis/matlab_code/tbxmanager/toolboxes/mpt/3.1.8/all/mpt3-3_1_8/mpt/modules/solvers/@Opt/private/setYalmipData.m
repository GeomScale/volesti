function opt = setYalmipData(opt, con, obj, th, u)
%
%         Set data from a YALMIP object
%
% Convert a yalmip problem min obj s.t. con to the form of a
% parametric optimization problem
%
% th == parametric variables
% u  == optimization variables of interest

global MPTOPTIONS

if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end


if nargin < 5, error('Opt: Syntax : opt.setYalmipData(con, obj, th, u)'); end

% call solvesdp with options.pureexport=true to get the parsed model
options = sdpsettings;
options.pureexport = true;
interfacedata = solvesdp(con, obj, options, th, u);

% Convert from Yalmip => mpt2.6 => mpt3 format
mat = yalmip2mpt(interfacedata);
mat.lb(mat.lb==Inf) = MPTOPTIONS.infbound;
mat.lb(mat.lb==-Inf) = -MPTOPTIONS.infbound;
mat.ub(mat.ub==Inf) = MPTOPTIONS.infbound;
mat.ub(mat.ub==-Inf) = -MPTOPTIONS.infbound;
opt = opt.setMPT26Data(mat);

end
