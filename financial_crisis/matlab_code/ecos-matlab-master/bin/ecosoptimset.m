function opts = ecosoptimset(varargin)
% Create default options struct for ECOS solver.
%
% OPTIONS = ECOSOPTIMSET returns a pre-initialized options struct for ECOS.
%
% OPTIONS = ECOSOPTIMSET(PARAMETER1,VALUE1,...,PARAMETERN,VALUEN) returns an
% inialized options struct for ECOS with the parameters PARAMETER set to
% VALUE. Multiple parameters have to be defined in PARAMETER-VALUE pairs.
% PARAMETER must be a string of one of the following keywords:
%
% 
%    For the core ECOS solver:
%       VERBOSE: set print level to 0 (off), 1 (solve status), 2 (each IPM iteration)
%       ABSTOL, FEASTOL, RELTOL: numerical values of stopping criteria
%       ABSTOL_INACC, FEASTOL_INACC, RELTOL_INACC: reduced precision stopping criteria
%       MAXIT: maximum number of iterations
% 
%    For BOOLEAN or INTEGER PROGRAMMING (the module is (C) Han Wang, Stanford University):
%      BOOL_VARS_IDX: index array of boolean variables (1-based indexing)
%      INT_VARS_IDX:  index array of integer variables (1-based indexing)
%      MI_MAXIT: maximum number of SOCP solves (branch&bound sub-problems)
%      MI_ABS_GAP_TOL: required absolute gap between lower and upper bound
%      MI_REL_GAP_TOL: required relative gap between lower and upper bound
%                      w.r.t. upper bound
%
% (C) A. Domahidi, ETH Zurich & embotech GmbH, Zurich, Switzerland, 2012-15.


if( mod(nargin,2) ~= 0)
    error('A parameter list with string-value pairs is expected as argument')
end

% default entries
opts.feastol = eps^(1/2);
opts.reltol = eps^(1/2);
opts.abstol = eps^(1/2);
opts.maxit = 50;
opts.verbose = 2;

% set entries as requested by user
for i = 1:2:nargin
    assert(ischar(varargin{i}),'Parameter %d is not a string');
    opts.(lower(varargin{i})) = varargin{i+1};
end

