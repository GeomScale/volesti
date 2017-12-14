function h = fplot_internal(obj, function_name, options)
%
% Plot function over a single convex set
%
% This is an internal helper called from ConvexSet/fplot. No error checks
% are performed. Implicitly assumes the object is a singleton.

% global options
global MPTOPTIONS
if isempty(MPTOPTIONS)
	MPTOPTIONS = mptopt;
end

fun = obj.getFunction(function_name);

% if "obj" is Polyhedron, then Polyhedron/fplot_internal() should have been
% called automatically by Matlab

% for any other convex set do gridding

if obj.Dim ==1
	B = obj.outerApprox;
	% lb and ub only
	X = [B.Internal.lb;
		B.Internal.ub];
elseif obj.Dim == 2
	% Grid the circle
	th = linspace(0,2*pi, options.grid+1)';
	th(end) = [];
	
	X = [sin(th) cos(th)];
	
elseif obj.Dim == 3
	[X Y Z] = sphere(options.grid);
	X = [X(:) Y(:) Z(:)];
end

% Compute an extreme point in each direction
E = []; % Extreme points
R = []; % Rays
tic; 
first = true;
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

% Create a polyhedron out of given extreme points
P = Polyhedron('V', E, 'R', R);
% Copy functions
P.copyFunctionsFrom(obj);
% Plot
h = P.fplot_internal(function_name, options);

end
