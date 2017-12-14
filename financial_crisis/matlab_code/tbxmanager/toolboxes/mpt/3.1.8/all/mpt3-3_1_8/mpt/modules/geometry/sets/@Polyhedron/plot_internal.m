function h = plot_internal(P, options)
% Plots a single polyhedron.
%
% This is an internal helper called from ConvexSet/plot. No error checks
% are performed. Implicitly assumes the object is a singleton. Both input
% arguments are required. "options" must be a structure provided by
% ConvexSet/plot.

global MPTOPTIONS

h = [];

if P.isEmptySet() || P.isFullSpace()
	% do not plot empty sets and R^n
	return
end

% If color is a letter, rather than a vector then convert
if ischar(options.color) || isempty(options.color)
	clr = charToColor(options.color, options.colormap, options.array_index, options.colororder);
else
	clr = options.color;
end

if options.wire,
	options.alpha = 0;
end

% If we're trying to plot an affine hull, then create a bounding polytope and re-plot
if size(P.He,1) > 0 && size(P.H,1) == 0
	% some objects can be already plotted, take current axis
	xlim = get(gca,'Xlim');
	lb = xlim(1);
	ub = xlim(2);
	if P.Dim >= 2
		ylim = get(gca,'Ylim');
		lb = [lb;ylim(1)];
		ub = [ub;ylim(2)];
	end
	if P.Dim >= 3
		zlim = get(gca,'Zlim');
		lb = [lb;zlim(1)];
		ub = [ub;zlim(2)];
	end
	Q = Polyhedron('lb',lb,'ub',ub,'He',P.He);
	% in case Q is empty, put larger bounds and repeat eventually 100-times
	k=0;
	d = abs(ub-lb);
	while Q.isEmptySet
		lb = lb-d/2;
		ub = ub+d/2;
		Q = Polyhedron('lb',lb,'ub',ub,'He',P.He);
		k = k+1;
		if k>100
			break;
		end
	end
	h = plot_internal(Q, options);
	return
end

% Compute the incidence map
% (also computes irredundant V-rep and H-rep)
iMap = P.incidenceMap;

V  = P.V;
R  = P.R;
H  = P.H;
He = P.He;

% if the Polyhedon consist of just one point change the default marker to "."
if size(He,1)>=P.Dim
	if strcmpi(options.marker,'none')
		options.marker='.';
	end
end

% Deal with silly cases. Dimension 0 and 1
if size(V,2) == 1
	if size(R,1) > 0 % Zero-dimensional or 1D ray
		if strcmp(options.marker, 'none')
			options.marker = '.';
			options.markerSize = max([15 options.markerSize]);
		end
		h = pplot(V, '.', 'marker', options.marker, 'markersize', options.markerSize, 'color', clr);
		
		if size(R,1) == 1 % One dimensional ray
			h2 = pplot([V;V+R/norm(R)],'-','linestyle',options.linestyle,'linewidth',options.linewidth);
			h = [h; h2];
		end
	else % 1D or 0D polytope
		h = pplot(V, '.', 'marker', options.marker, 'markersize', options.markerSize, 'color', clr, ...
			'linestyle',options.linestyle,'linewidth',options.linewidth);
	end
	return
end

% Lower-dimensional polyhedron
lowerDim = false;
if size(He,1) > 0
	H = He(1,:);
	lowerDim = true;
end

% Lift the 2D plot to 3D
if size(V,2) == 2
	V = [V zeros(size(V,1),1)];
	R = [R zeros(size(R,1),1)];
	
	% There's only one 'facet'
	H = [0 0 1 0];
	lowerDim = true;
end

newV = zeros(0,3);
if size(R,1) > 0
	
	%   % Choose a hyperplane that intersects all rays as close to the unit
	%   % circle as possible
	%   R = normalize(R);
	%   sol = Opt('H',2*(R'*R),'f',-2*sum(R,1),'A',-R,'b',zeros(size(R,1),1)).solve;
	% %   sol = mptSolve('H',2*(R'*R),'f',-2*sum(R,1),'A',-R,'b',zeros(size(R,1),1));
	%   if sol.exitflag ~= MPTOPTIONS.OK,
	%     warning('Solver error bounding the recession cone. This is likely because the polyhedron has a non-empty lineality space.');
	%     h = -1;
	%     return
	%   end
	
	if size(R,1) == 1
		sol.x = R';
	else
		sol.x = mean(normalize(R))';
	end
	sol.x = sol.x / norm(sol.x);
	if any(R*sol.x < -MPTOPTIONS.abs_tol),
		error('Could not find ray in strict interior of the recession cone');
	end
	
	% Move the hyperplane so that it's well outside the bounded portion of
	% the polyhedron
	PV   = Polyhedron('V',V);
	supp = PV.support(sol.x);
	B    = PV.outerApprox;
	sz   = max(B.Internal.ub - B.Internal.lb);
	supp = supp + max([2*sz 1]);
	% Our hyperplane is now sol.x'*y == supp
	
	% Add a ray to each vertex that shares at least two facets with the ray
	incVR = iMap.incVH*iMap.incRH' >= P.Dim - size(He,1) - 1;
	
	% newV are aritificial vertices generated from the rays
	nNewV = sum(incVR(:));
	newV = zeros(nNewV,3);
	k = 1;
	newLen = zeros(nNewV,1);
	for i=1:size(V,1)
		for j=1:size(R,1)
			if incVR(i,j)
				% Add vertex at x = V(i) + alpha*R(j)
				% Choose alpha s.t. x is in the hyperplane sol.x'*x == supp
				alpha = (supp-sol.x'*V(i,:)') / (sol.x'*R(j,:)');
				newV(k,:) = V(i,:) + alpha*R(j,:);
				newLen(k) = norm(alpha*R(j,:));
				k = k + 1;
			end
		end
	end
	newLen = mean(newLen);
end

% Data points
X = [newV;V];

% Compute incidence map
inc = abs(H*[X -ones(size(X,1),1)]')' < MPTOPTIONS.abs_tol;

% Convert from incidence map to Matlab face matrix
[I,J] = find(inc);
Faces = zeros(size(H,1),max([sum(inc,1) size(newV,1)]));
for i=1:size(H,1)
	v = I(J==i);
	
	% This inequality is weakly redundant
	% This happens e.q. when plotting a cone, since we add an artificial facet
	%   if isempty(v), remove(end+1) = i; continue; end
	if isempty(v),
		Faces(i,:) = NaN*ones(1,size(Faces,2));
		continue
	end
	
	% Order vertices for plot
	ord = orderForPlot(X(v,:),H(i,1:end-1));
	v = v(ord);
	
	Faces(i,:) = [v' NaN*ones(1,size(Faces,2)-length(v))];
end

% Plot!
set(gcf,'renderer','opengl');
h(1) = patch('Vertices', X, ...
	'Faces', Faces, ...
	'FaceColor', clr,...
	'EdgeColor', options.edgecolor, ...
	'EdgeAlpha', options.edgealpha, ...
	'FaceAlpha', options.alpha,...
	'LineStyle', options.linestyle,...
	'LineWidth', options.linewidth,...
	'Marker', options.marker,...
	'MarkerSize', options.markerSize);

% Plot the 'endcap' in white
if size(R,1) > 0
	if lowerDim == false
		ord = orderForPlot(newV,sol.x);
		h(2) = patch('Vertices', newV, ...
			'Faces', ord, ...
			'FaceColor', [0 0 0],...
			'FaceAlpha', 0.75);
	else
		
		% newV contains two elements
		h(2) = pplot(newV, 'w-', 'linewidth', options.linewidth*2);
		
	end
	
	% Add an arrow on the unbounded facet
	if size(newV,1)>1
		x0 = mean(newV)';
	else
		x0 = newV';
	end
	x1 = x0 + 0.5*newLen*sol.x;
	
	%   alpha = (1.5*supp-sol.x'*x0) / (sol.x'*sol.x);
	%   x1 = x0 + alpha*sol.x;
	hl(1) = line([x0(1);x1(1)],[x0(2);x1(2)],[x0(3);x1(3)]);
	set(hl, 'linewidth', options.linewidth*1.5, 'color', 'k');
	hl(2) = plot3(x1(1),x1(2),x1(3),'k.','markersize',options.linewidth*15);
	h = [h(:); hl(:)];
end

if isfield(options, 'showindex') && options.showindex && P.Dim==2 && P.isBounded()
	xc = P.chebyCenter.x;
	text(xc(1), xc(2), num2str(options.array_index));
end

smoothLines(h);

end
