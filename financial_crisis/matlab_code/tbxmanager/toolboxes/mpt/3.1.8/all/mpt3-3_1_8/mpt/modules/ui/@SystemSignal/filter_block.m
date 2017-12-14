function filter = filter_block(varargin)
%
%  BLOCK: Adds a move blocking constraint 
%  =======================================
%  
%  
%  SYNTAX
%  ------
%     
%      
%    
%  
%  DESCRIPTION
%  -----------
%     Adding this filter to an MPC setup will modify the constraints in such a way
%  that the differences between several consecutive optimization variables are
%  equal to zero.
%    In the most common scenario, adding this filter to a signal representing
%  control inputs will add a move-blocking constraint, which is equivalent to
%  setting a control horizon.
%    To enable this filter, call model.u.with('block') (note that you can add this
%  constrain to any type of signals, e.g., to state and output signals as well).
%    Once the filter is enabled, parameters of the blocking scheme can be specified
%  in the model.u.block.from and model.u.block.to parameters. Setting these values
%  to non-zero integers will add the constraint u_from = u_from+1 = ??? = u_to-1 =
%  u_to.
%    The filter can be removed by calling model.x.without('block').
%  

%  AUTHOR(s)
%  ---------
%     
%    
%   (c) 2003-2013  Michal Kvasnica: STU Bratislava
%   mailto:michal.kvasnica@stuba.sk 
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
 
 
filter = FilterSetup;
filter.addField('from', [], @isnumeric);
filter.addField('to', [], @isnumeric);

% the filter impacts the following calls:
filter.callback('constraints') = @on_constraints;

end

%------------------------------------------------
function out = on_constraints(obj, varargin)
% called when constructing constraints

out = [];
for k = obj.block.from:obj.block.to-1
    out = out + [ obj.var(:, k) == obj.var(:, k+1) ];
end

end
