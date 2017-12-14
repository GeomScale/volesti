function S = mpt_ui_options
%
% Option settings for "ui" module.
%


% maximal number of iterations for LTISystem.invariantSet()
S.invariantSet.maxIterations = 100;

% maximal prediction horizon for minimum-time controllers
S.EMinTimeController.maxIterations = 100;

% settings for search tree export
S.EMPCSTController.roundcoef = 1e12;

end
