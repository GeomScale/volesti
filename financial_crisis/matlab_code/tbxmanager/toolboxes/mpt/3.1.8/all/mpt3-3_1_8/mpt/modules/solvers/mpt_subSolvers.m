function y = mpt_subSolvers(varargin)
%% detect solvers

global MPTOPTIONS
persistent MPT_SOLVERS_LIST

if ~isempty(MPTOPTIONS)
    s = MPTOPTIONS.solvers_list;
else
    %if the list of solvers has been saved, use that one
    if ispref('MPT','MPTOPTIONS')
        v = getpref('MPT','MPTOPTIONS');
        s = v.solvers_list;
    elseif ~isempty(MPT_SOLVERS_LIST)
        s = MPT_SOLVERS_LIST;
    else
        s = mpt_detect_solvers;
        MPT_SOLVERS_LIST = s;
    end
end

% return the first solvers in the list
if nargin==0
    y = s;
elseif nargin==1
    switch upper(varargin{1})
        case 'LP'
            y = s.LP{1};
        case 'QP'
            y = s.QP{1};
        case 'MILP'
            if isempty(s.MILP)
                y = '';
            else
                y = s.MILP{1};
            end
        case 'MIQP'            
            if isempty(s.MIQP)
                y = '';
            else
                y = s.MIQP{1};
            end
        case 'LCP'
			if isempty(s.LCP)
				y = '';
			else
				y = s.LCP{1};
			end
        case 'PLP'
            y = s.parametric.LP{1};
        case 'PQP'
            y = s.parametric.QP{1};
        case 'PLCP'
            if isempty(s.parametric.LCP)
                y = '';
			else
                y = s.parametric.LCP{1};
            end
        otherwise
            error('Unrecognized string.');
    end
else
    y = '';
end

end
