function display(opt)
%
%  DISPLAY: Overload display for Opt class. 
%  =========================================
%  
%  
%  SYNTAX
%  ------
%     
%      display(problem)
%      problem.display()
%    
%  
%  DESCRIPTION
%  -----------
%     Default display for Opt class.
%  
%  INPUT
%  -----
%     
%        
%          problem Object of the Opt class that defines     
%                  optimization problem.                    
%                  Class: Opt                               
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
%   (c) 2010-2013  Colin Neil Jones: EPF Lausanne
%   mailto:colin.jones@epfl.ch 
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

if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

if numel(opt) > 1
    fprintf('Array of %i Opt objecs.\n', numel(opt));
    return
end
if numel(opt)==0
    fprintf('Empty problem.\n');
    return;
end

if opt.n==0 && ~opt.isParametric
    fprintf('Empty problem.\n');
    return;
end

fprintf('-------------------------------------------------\n');
try
    %opt = opt.validate;
    switch opt.problem_type
        case 'FEAS', Dtype = 'Feasibility problem';
        case 'QP', Dtype = 'Quadratic program';
        case 'LP', Dtype = 'Linear program';
        case 'LCP', Dtype = 'Linear complementarity problem';
        case 'MILP', Dtype = 'Mixed Integer Linear Problem';
        case 'MIQP',  Dtype = 'Mixed Integer Quadratic Problem';
        otherwise, Dtype = 'unknown';
    end
    if opt.isParametric
        %switch opt.PQPSOLVER
        %    case opt.PQPSOLVER_F2F,   Dsolver = 'Facet-to-facet solver';
        %    case opt.PQPSOLVER_MPT26, Dsolver = 'MPT 2.6 solver';
        %    case opt.PQPSOLVER_PLCP,  Dsolver = 'Parametric LCP';
        %    otherwise, Dsolver = [];
        %end
        Dtype = sprintf('Parametric %s%s', lower(Dtype));
    end
    
    fprintf('%s\n', Dtype);
    
    fprintf('\tNum variables:              %3i\n', opt.n);
    fprintf('\tNum inequality constraints: %3i\n', opt.m);
    fprintf('\tNum equality constraints:   %3i\n', opt.me);
    if any(~isinf(opt.lb))
        fprintf('\tNum lower bounds            %3i\n', sum(~isinf(opt.lb)));
    end
    if any(~isinf(opt.ub))
        fprintf('\tNum upper bounds            %3i\n', sum(~isinf(opt.ub)));
    end
    if opt.isParametric
        fprintf('\tNum parameters:             %3i\n', opt.d);
    end
    
    fprintf('\tSolver:                     %3s\n', opt.solver);
catch err
    fprintf('Invalid optimization problem.\n\n');
    fprintf('\tError : %s\n', err.message);
end
fprintf('-------------------------------------------------\n');

end
