function h = plot_internal(obj, options)
% Plots a single convex set
%
% This is an internal helper called from ConvexSet/plot. No error checks
% are performed. Implicitly assumes the object is a singleton. Both input
% arguments are required. "options" must be a structure provided by
% ConvexSet/plot.

global MPTOPTIONS

if obj.isEmptySet()
	% do not plot empty sets
	h = [];
	return
end
% if any(~obj.isBounded)
% 	error('Can only plot bounded sets.');
% end

if obj.Dim == 1
	% compute extremal point in both directions
	X = [-1; 1];
	
elseif obj.Dim == 2
    % Grid the circle
    th = linspace(0, 2*pi, options.grid+1)';
    th(end) = [];
    
    X = [sin(th) cos(th)];
else
    [X Y Z] = sphere(options.grid);
    X = [X(:) Y(:) Z(:)];
end

% Compute an extreme point in each direction
E = []; % Extreme points
R = []; % Rays
tic; first = true;
for i=1:size(X,1)
    if toc > MPTOPTIONS.report_period
        tic
        if first, fprintf('Plotting...\n'); first=false; end
        fprintf('%i of %i\n', i, size(X,1));
    end
    ret = obj.extreme(X(i,:));
    switch ret.exitflag
        case MPTOPTIONS.OK
            E(end+1,:) = ret.x';
        case MPTOPTIONS.UNBOUNDED
            R(end+1,:) = ret.x';
        case MPTOPTIONS.ERROR
            error('Cannot tell if object is unbounded or infeasible. Try different solver.');
        otherwise
            error('A problem occured inside the solver. Try different solver.');
    end
end

% Create a polyhedron and plot it
P = Polyhedron('V', E, 'R', R);
h = plot_internal(P, options);

end
