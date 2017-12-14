classdef ConvexSetInterface < handle
  %%
  % ConvexSetInterface
  %
  % Represents a convex set in R^d.
  %
  % Getting the basic functionality of all operations given in this class
  % requires that only two additional functions be specified -
  % extreme and project.
  %
  
  %%
  properties(Abstract = true, SetAccess = protected)
%   properties(SetAccess = private)
    Dim; % Dimension of representation
  end
  
  % These should be made abstract, but matlab seems to have trouble
  % counting the number of output arguments in that case...
   methods(Abstract = true)
%  methods
    
    sol = extreme(obj, y)
    %
    % [e, flag] = extreme(obj, x)
    %
    % Compute an extreme point of this set in the direction x.
    %
    % Parameters:
    % x - Vector of length d
    %
    % Returns:
    % e    - An extreme point of C from the face , or NULL if empty.
    % flag -  1  if infeasible
    %         2  if unbounded
    %         12 if infeasible or unbounded (can't tell which)
    %         0  if ok
    %
    
    %%
    sol = project(obj, y)
    %
    % Compute the projection of the point x onto this set.
    % Computes the closest point in this set to the point x:
    %
    % Parameters:
    % x - Vector of length d
    %
    % Returns:
    % p    - Closest point to x contained in this set, or NULL if empty.
    % flag -  1  if infeasible
    %         2  if unbounded
    %         12 if infeasible or unbounded (can't tell which)
    %         0  if ok
    %

    %%
    sol = shoot(obj, r, x0)
    %
    % Compute max alpha s.t. x*alpha in Set
    %
    % Parameters:
    % x - Vector of length d
    %
    % Returns:
    % alpha - Maximum value of alphs s.t. x*alpha is in set, or NaN if
    % empty
    %
  end
end

