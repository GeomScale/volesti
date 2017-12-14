function opt=validate_solvernames(opt)
%% validate solvernames
% if "solver" is empty, take the first from the group, otherwise check if
% the solver is in the list

global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS=mptopt;
end
%disp('solver_validation')
s = MPTOPTIONS.solvers_list;
if ~opt.isParametric
    switch opt.problem_type
        case {'LP'}
            if isempty(opt.solver)
                opt.solver = MPTOPTIONS.lpsolver;
            elseif ~any(strcmpi(opt.solver,s.LP))
                error('Given solver is not in the list of LP solvers.');
            end
        case 'QP'
            if isempty(opt.solver)
                opt.solver = MPTOPTIONS.qpsolver;
            elseif ~any(strcmpi(opt.solver,s.QP))
                error('Given solver is not in the list of QP solvers.');
            end
        case 'LCP'
            if isempty(opt.solver)
                opt.solver = MPTOPTIONS.lcpsolver;
            elseif ~any(strcmpi(opt.solver,s.LCP))
                error('Given solver is not in the list of LCP solvers.');
            end
        case {'MILP'}
            if isempty(opt.solver)
                opt.solver = MPTOPTIONS.milpsolver;
            elseif ~any(strcmpi(opt.solver,s.MILP))
                error('Given solver is not in the list of MILP solvers.');
            end
        case 'MIQP'
            if isempty(opt.solver)
                opt.solver = MPTOPTIONS.miqpsolver;
            elseif ~any(strcmpi(opt.solver,s.MIQP))
                error('Given solver is not in the list of MIQP solvers.');
            end
        otherwise
            error('Unknown type of optimization problem.');
    end
else
    switch opt.problem_type
        case {'LP','MILP'}
            if isempty(opt.solver)
                opt.solver = MPTOPTIONS.plpsolver;
            elseif ~any(strcmpi(opt.solver,s.parametric.LP))
                error('Given solver is not in the list of parametric LP solvers.');
            end
        case {'QP','MIQP'}
            if isempty(opt.solver)
                opt.solver = MPTOPTIONS.pqpsolver;
            elseif ~any(strcmpi(opt.solver,s.parametric.QP))
                error('Given solver is not in the list of parametric QP solvers.');
            end            
        case {'LCP'}
            if isempty(opt.solver)
                opt.solver = MPTOPTIONS.plcpsolver;
            elseif ~any(strcmpi(opt.solver,s.parametric.LCP))
                error('Given solver is not in the list of parametric LCP solvers.');
            end            
        otherwise
            error('Unknown type of optimization problem.');        
    end    
end

% throw error if there's no solver for given problem
if isempty(opt.solver)
    error('For %s problems there is no solver assigned.',opt.problem_type);
end

end