function ret = mpt_call_enum_plcp(opt)
%
% a gateway routine to parametric ENUMPLCP solver
%

global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

if strcmpi(opt.problem_type,'LCP')
    % LCP
    
    % direct call
    ret = mpt_enum_plcp(opt);
   
elseif any(strcmpi(opt.problem_type,{'LP','QP','MILP','MIQP'}))
    % LP or QP
    
    % indirect call
    % 1. Convert to LCP format
    lc = opt.copy;
    lc.qp2lcp;
    
    % 2. Solve the LCP
    ret = mpt_enum_plcp(lc);    
    
else
    error('mpt_call_enumplcp: ENUMPLCP solver does not solve %s problems.',opt.problem_type);    
end