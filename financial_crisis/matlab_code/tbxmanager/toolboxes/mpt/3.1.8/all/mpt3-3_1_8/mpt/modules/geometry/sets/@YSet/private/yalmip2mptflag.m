function sol = yalmip2mptflag(diagnostic)
%
% Convert from a YALMIP optimization error type to an MPT type
%

global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

narginchk(1, 1);

if ~isa(diagnostic,'struct')
    error('Input argument must be structure - diagnostic info about the solution from YALMIP.');
end

% convert to MPT flags
switch diagnostic.problem
    case {0,4}
        sol.exitflag = MPTOPTIONS.OK;        
    case 1, 
        sol.exitflag = MPTOPTIONS.INFEASIBLE;        
    case 2,    
        sol.exitflag = MPTOPTIONS.UNBOUNDED;        
    otherwise, 
        sol.exitflag = MPTOPTIONS.ERROR;
end

% put exit status from the solver
if isfield(diagnostic,'info')
    sol.how = diagnostic.info;
elseif isfield(diagnostic,'infostr')
    sol.how = diagnostic.infostr;
else
    sol.how = yalmiperror(diagnostic.problem);
end

end
