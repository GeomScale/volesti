function h = fplot_internal(obj, function_name, options)
%
% Plot function over a single polyhedron
%
% This is an internal helper called from ConvexSet/fplot. No error checks
% are performed. Implicitly assumes the object is a single polyhedron.

fun = obj.getFunction(function_name);

h = [];

% if the Polyhedon consist of just one point change the default marker to "."
if size(obj.He,1)>=obj.Dim
	options.marker = '.';
else
	options.marker = 'none';
end
options.markerSize = 6;

% get color as RGB vector
if ischar(options.color) || isempty(options.color)
	options.color = charToColor(options.color, options.colormap, ...
		options.array_index, options.colororder);
end

if obj.Dim<2
	% plot 1D sets
	obj.outerApprox();
	x = linspace(obj.Internal.lb, obj.Internal.ub, options.grid)';
	y = zeros(size(x));
	for j=1:numel(x)
		v = fun.feval(x(j));
		y(j) = v(options.position);
	end
	
	hl = line(x, y, 'Color', options.color, ...
		'LineStyle', options.linestyle, ...
		'LineWidth', options.linewidth, ...
		'Marker', options.marker);
	h = [h; hl(:)];
	
else
	% plot 2D sets
	if isa(fun, 'AffFunction') && ~options.contour
		obj.minVRep();
		X = obj.V;
		V = zeros(size(X,1), 1);
		
		% TODO: exploit vectorization once Function/feval supports it
		% V = fun.feval(V');
		for j=1:size(X,1)
			t = fun.feval(X(j, :)');
			V(j, 1) = t(options.position);
		end
		[~, ~, vv]= svd([X V]);
		n = vv(:, 3);
		I = orderForPlot([X V], n);
		
		hp = patch('Vertices', [X V], 'Faces', I(:)', ...
			'LineStyle', options.linestyle, ...
			'LineWidth', options.linewidth, ...
			'FaceAlpha', options.alpha,...
			'EdgeColor', options.edgecolor,...
			'FaceLighting', 'phong',...
			'AmbientStrength', 0.7,...
			'FaceColor', options.color,...
			'Marker', options.marker);
		h = [h; hp(:)];
		
	else
		[X,Y] = obj.meshGrid(options.grid);
		
		V = NaN*X;
		for j = 1:size(X,1)
			for k = 1:size(X,2)
				x = [X(j,k);Y(j,k)];
				if isnan(x(1)) || isnan(x(2)), continue; end
				t = fun.feval(x);
				if isempty(t) || any(isnan(t)), continue; end
				V(j,k) = t(options.position);
			end
		end

		if options.showgrid
			edgealpha = 1;
			AlphaData = ones(size(V));
		else
			% plot all edges with alpha=0 to hide them
			AlphaData = zeros(size(V));
			edgealpha = 'flat';
		end
		
		% Plot the function
		if options.alpha == 0 && strcmp(options.linestyle, 'none')
			hs = [];
		else
			if options.contour
				surf_fun = @surfc;
			else
				surf_fun = @surf;
			end
			hs = surf_fun(X,Y,V,'facecolor', options.color,...
				'linestyle', options.linestyle, ...
				'linewidth', options.linewidth, ...
				'facealpha', options.alpha,...
				'facelighting', 'phong',...
				'AmbientStrength', 0.7,...
				'AlphaData', AlphaData, ...
				'edgealpha', edgealpha, ...
				'Marker', options.marker);
			
			if ~options.showgrid && obj.isFullDim()
				% plot the outline
				
				% detect outer boundaries
				B = zeros(size(X));
				for i = 1:size(B, 2)
					p_nan = find(~isnan(X(:, i)));
					if ~isempty(p_nan)
						B(p_nan(1), i) = 1;
						B(p_nan(end), i) = 1;
					end
				end
				for i = 1:size(B, 1)
					p_nan = find(~isnan(X(i, :)));
					if ~isempty(p_nan)
						B(i, p_nan(1)) = 1;
						B(i, p_nan(end)) = 1;
					end
				end

				% prepare points
				P = [];
				[irow, icol] = find(B);
				for i = 1:numel(irow)
					px = X(irow(i), icol(i));
					py = Y(irow(i), icol(i));
					pz = V(irow(i), icol(i));
					P = [P; px py pz];
				end
				
				% sort vertices in a cyclic way
				xc = obj.chebyCenter.x;
				ang = angle((P(:, 1)-xc(1))+(P(:, 2)-xc(2))*sqrt(-1));
				[~, ind] = sort(ang);
				x1 = P(ind, 1);
				x2 = P(ind, 2);
				x3 = P(ind, 3);
				
				% plot the outline
				hb = line([x1; x1(1)], [x2; x2(1)], [x3; x3(1)]);
				set(hb, 'Color', options.edgecolor, ...
					'linestyle', options.linestyle, ...
					'linewidth', options.linewidth);
			end
		end
		h = [h; hs(:)];
		
	end
end

% plot the polyhedron if required
if options.show_set
	hv = obj.plot_internal(options);
	h = [h; hv(:)];
end

smoothLines(h);

end
