classdef polytope
	% Compatibility class for polytopes
	
	properties(SetAccess = private, Hidden)
		P % internal Polyhedron representation of the polytope
	end
	
	methods
		
		function obj = polytope(arg1, arg2, varargin)
			% Constructor
			
			if nargin==1 && isa(arg1, 'Polyhedron') && length(arg1)==1
				% import from a Polyhedron object
				% (this case is the most frequent one)
				obj.P = arg1;
				if ~arg1.irredundantHRep
					obj.P.minHRep();
				end
				return
				
			elseif nargin>=2
				% import Hrep from A, b
				idx = isinf(arg2);
				arg1(idx, :) = [];
				arg2(idx) = [];
				obj.P = Polyhedron(arg1, arg2);
				if ~obj.P.isFullDim
					obj = polytope;
				end
				
			elseif nargin==1 && ...
                    (isa(arg1, 'Polyhedron') || isa(arg1, 'PolyUnion'))
				% import a polyhedral array
                if isa(arg1, 'PolyUnion')
                    arg1 = arg1.Set;
                end
				n = numel(arg1);
				if n==0
					obj.P = Polyhedron;
				elseif n==1
					obj.P = arg1;
					if ~arg1.irredundantHRep
						obj.P.minHRep();
					end
				else
					obj.P = arg1;
					for i = 1:n
						if ~arg1(i).irredundantHRep
							obj.P(i).minHRep();
						end
					end
				end
				return
				
			elseif nargin==0
				% an empty polytope
				obj.P = Polyhedron;
				return

			elseif nargin==1
				% import Vrep
				if isa(arg1, 'double')
					obj.P = Polyhedron('V', arg1);
					obj.P.minHRep();
				else
					error('Unrecognized type of first input argument.');
				end
			end
			
			% normalize the H-representation
			if nargin<=2 || varargin{1}==0
				obj.P.normalize();
			end
			
			% eliminate redundant constraints
			if nargin<=3 || varargin{2}==0
				for i = 1:numel(obj.P)
					obj.P(i).minHRep();
				end
			end
		end
		
		function out = toPolyhedron(obj)
			% Convert the object to an instance of @Polyhedron
			
			% use the copy constructor
			out = Polyhedron(obj.P);
		end
		
		function out = display(obj)
			% Display method
			
			obj.P.display();
			fprintf('\nThe @polytope class will be deprecated, use @Polyhedron instead.\n');
		end
		
		function out = end(obj, varargin)
			
			out = length(obj.P);
		end
		
		function [obj, keptrows] = reduce(obj, varargin)
			% Eliminates redundant constraints
			
			if obj.P.irredundantHRep
				% already in minal representation
				keptrows = 1:size(obj.P.A, 1);
			else
				% recompute
				[~, details] = obj.P.minHRep();
				keptrows = find(details.I==0);
			end
		end
		
		function [out, kept] = reduceunion(obj, varargin)
			% Reduces union of polytopes
			
			U = PolyUnion(obj.P);
			kept = U.reduce();
			out = polytope(U.Set);
		end
		
		function out = volume(obj, varargin)
			% Returns volume of the polytope
			
			out = obj.P.volume();
		end
		
		function out = merge(obj, varargin)
			% Greedy merging of union of polytopes
			
			out = polytope(PolyUnion(obj.P).merge());
		end
		
		function out = isminrep(obj)
		
			out = obj.P.irredundantHRep;
		end
		
		function [A, b] = double(obj)
			% Return the H-representation
			
			H = obj.P.H;
			A = H(:, 1:end-1);
			b = H(:, end);
			if nargout<=1
				A = H;
			end
		end
		
		function h = plot(varargin)
			% Polytope plotter
            %
            % Syntax:
            %   plot(P)
            %   plot(P, 'y') -- plot the polytope in yellow color
            %   plot(P, [0.9 0.9 0.9]) -- use RGB color code
            %   plot(P, options) -- use a structure of options (see below)
            %   plot(P, 'y', options)
            %
            % Note: To plot a colored wireframe, use
            %   plot(P, struct('wire', 1, 'edgecolor', 'b'))
            %
            % Supported options:
            %       'alpha': transparency (1=opaque, 0=complete transparency)
            %                double, default=1
            %        'wire': if true, the set is plotted in wire frame
            %                logical, default=false
            %   'linestyle': style of the set's border
            %                string, default='-'
            %   'linewidth': width of the set's border
            %                double, default = 1
            %   'edgecolor': color of edges
            %                string, default='k'
            %   'edgealpha': transparency of edges (1=opaque, 0=complete transparency)
            %                double, default=1
            %      'marker': marker of the plot
            %                string, default = 'none'
            %  'markerSize': size of the marker
            %                double, default=6
            %    'colormap': color map to use
            %                string or a function handle, default='mpt'
			
            h = [];
            
            % split input arguments into objects and corresponding options
            [objects, options] = parsePlotOptions('polytope', varargin{:});
            if numel(objects)==0
                % no objects to plot
                return
            end
            
            % plot each object separately, hence we need to hold on
            prevHold = ishold;
            if ~ishold,
                newplot;
            end
            hold('on');
            
            % plot each polytope
            for i = 1:length(objects)
                
                % determine options for each object
                args = {};
                for j = 1:numel(options{i})
                    if isstruct(options{i}{j})
                        % convert structure into a set of {'param', value}
                        % pairs
                        f = fieldnames(options{i}{j});
                        for k = 1:numel(f)
                            args{end+1} = f{k};
                            args{end+1} = options{i}{j}.(f{k});
                        end
                    else
                        % color needs to be prefixed
                        args{end+1} = 'color';
                        args{end+1} = options{i}{j};
                    end
                end
                
                h = [h; objects{i}.P.plot(args{:})];
            end
            
            if ~prevHold
                hold('off');
            end
			if nargout==0
				clear h
			end
		end
		
		function out = regiondiff(obj, arg, varargin)
			% Set difference
			
			out = polytope(obj.P \ arg.P);
		end
		
		function out = mldivide(obj, arg, varargin)
			% Set difference
			
			out = polytope(obj.P \ arg.P);
		end
		
		function [B, lb, ub] = bounding_box(obj, Options, varargin)
			% Outer box approximation of the polytope

			if nargin<2
				Options = [];
			end
			
			% let obj.P.Internal.lb/ub be computed
			obj.P.outerApprox;

			% even for arrays we return a single bounding box, so let us
			% get the tightest bounds first
			d = obj.P(1).Dim;
			lb = Inf(d, 1);
			ub = -Inf(d, 1);
			for i = 1:length(obj.P)
				lb = min(lb, obj.P(i).Internal.lb);
				ub = max(ub, obj.P(i).Internal.ub);
			end
			
			if isfield(Options, 'noPolyOutput') && Options.noPolyOutput
				% we are only interested in the bounds
				B = obj;
			else
				% return the bounding box as a polytope
				B = polytope([eye(d); -eye(d)], [ub; -lb]);
			end
			
			if nargin>1 && isfield(Options, 'Voutput') && Options.Voutput
				if isfield(Options, 'bboxvertices') && Options.bboxvertices
					B = extreme(B);
				else
					B = [lb ub];
				end
			end
		end
		
		function out = projection(obj, dims, varargin)
			
			out = polytope(obj.P.projection(dims));
		end

		function [x, r] = chebyball(obj, varargin)
			% Computes center and radius of the Chyshev ball
			
			a = obj.P.chebyCenter();
			x = a.x;
			r = a.r;
		end
		
		function out = dimension(obj)
			% Dimension
			
			out = obj.P.Dim;
		end
		
		function out = isfulldim(obj)
			% Returns true if the polytope is fully dimensional
			
			out = any(obj.P.isFullDim());
		end
		
		function out = isbounded(obj)
			% Returns true if the polytope is bounded
			
			out = any(obj.P.isBounded());
		end
		
		function [V, R, obj] = extreme(obj, varargin)
			% Computes extremal vertices of a polytope
			
			obj.P.computeVRep();
			V = obj.P.V;
			R = obj.P.R;
		end
		
		function [x, r] = facetcircle(obj, idx, varargin)
			% Compute Chebyshev's center of a given face

			global MPTOPTIONS

			if ~obj.P.irredundantHRep
				obj.P.minHRep();
			end
			if nargin < 2
				idx = 1:size(obj.P.H, 1);
			end
			nidx = length(idx);
			x = zeros(obj.P.Dim, nidx);
			r = zeros(1, nidx);
			for i = 1:nidx
				z = obj.P.chebyCenter(idx(i));
				% using F = obj.P.getFacet(idx(i)); z = F.chebyCenter; is
				% more natural, but much slower
				if z.exitflag~=MPTOPTIONS.OK
					% use a different method if chebyCenter failed
					try
						q = obj.P.facetInteriorPoints;
						z.x = q(idx(i), :)';
						z.r = 1;
					catch
						fprintf('Couldn''t compute point on a facet.\n');
						z.x = NaN(obj.P.Dim, 1);
						z.r = -Inf;
					end
				end
				x(:, i) = z.x;
				r(i) = z.r;
			end
		end
		
		function P = intersect(obj, Q)
			% Intersection of two polytopes
			
			if obj.P.isEmptySet() || Q.P.isEmptySet()
				P = obj;
			else
				P = polytope(obj.P.intersect(Q.P));
			end
		end
		
		function out = dointersect(P, Q)
			% Returns true if two polytopes intersect
			
			out = isfulldim(intersect(P, Q));
		end
		
		function P = and(obj, Q)
			% Intersection of two polytopes
			
			P = intersect(obj, Q);
		end
		
		function P = plus(obj, arg, varargin)
			% Minkowski addition of polytopes
			
			if isa(arg, 'polytope')
				arg2 = arg.P;
			else
				arg2 = arg;
			end
			if isa(obj, 'polytope')
				arg1 = obj.P;
			else
				arg1 = obj;
			end
			P = polytope(arg1 + arg2);
		end
		
		function P = minus(obj, arg, varargin)
			% Pontryagin difference
			
			P = polytope(obj.P - arg.P);
		end
		
		function [isin, inwhich, closest] = isinside(obj, x, varargin)
			% Determines whether the polytope contains a given point
			
			if numel(obj.P)==1 && isempty(obj.P.H)
				% MPT2 does not check for dimension of 'x' if the polytope
				% is empty and just returns false
				isin = false;
				inwhich = [];
				closest = 0;
			else
				if nargout>2
					[isin, inwhich, closest] = obj.P.isInside(x, varargin{:});
				else
					[isin, inwhich] = obj.P.isInside(x, varargin{:});
				end
			end
		end
		
		function obj = horzcat(obj, varargin)
			% Horizontal concatenation
			
			if isa(obj, 'double') && isempty(obj)
				P = [];
			elseif length(obj.P)==1 && isempty(obj.P.H)
				P = [];
			else
				P = obj.P;
			end
			for i = 1:length(varargin)
				Q = varargin{i};
				if ~isempty(Q.P(1).H)
					P = [P, Q.P];
				end
			end
			if isempty(P)
				obj = polytope;
			else
				obj.P = P;
			end
			if ~isa(obj, 'polytope')
				obj = polytope(obj.P);
			end
		end
		
		function obj = vertcat(obj, varargin)
			% Vertical concatenation
			
			obj = horzcat(obj, varargin{:});
		end
		
		function [out, how] = union(obj, varargin)
			% Convex union of an array of polytopes
			
			U = PolyUnion(obj.P);
			how = U.isConvex();
			if how
				out = polytope(U.convexHull());
			else
				out = obj;
			end
		end
		
		function out = isconvex(obj, varargin)
			% Determines whether an array of polytopes is covnex
			
			out = PolyUnion(obj.P).isConvex();
		end
		
		function out = length(obj)
			% Returns number of elements of a polytope array
			
			out = length(obj.P);
		end
		
		function out = eq(obj, arg)
			% Returns true if two polytopes are identical
			
			out = (obj.P == arg.P);
		end
		
		function out = ne(obj, arg)
			% Returns true if two polytopes are not equal

			out = ~(obj.P == arg.P);
		end

		function out = le(obj, arg, varargin)
			% Non-strict subsetness check
			
			out = PolyUnion(obj.P) <= PolyUnion(arg.P);
		end

		function out = lt(obj, arg, varargin)
			% Strict subsetness check
			
			out = PolyUnion(obj.P) < PolyUnion(arg.P);
		end

		function out = ge(obj, arg, varargin)
			% Non-strict supersetness check
			
			out = PolyUnion(obj.P) >= PolyUnion(arg.P);
		end

		function out = gt(obj, arg, varargin)
			% Strict supersetness check
			
			out = PolyUnion(obj.P) > PolyUnion(arg.P);
		end
		
		function out = subsref(obj, X)
			% Indexing of polytope arrays
			
			out = polytope(subsref(obj.P, X));
		end
		
		function out = subsasgn(obj, X, Q)

			if isa(Q, 'polytope')
				arg = Q.P;
			else
				arg = Q;
			end
			if isempty(obj)
				obj = polytope;
			end
			out = polytope(subsasgn(obj.P, X, arg));
		end
		
		function obj = set(obj, varargin)
			% Internal setter function
			
		end
		
		function out = get(obj, what)
			% Return polytope data
			
			switch lower(what)
				case 'h',
					out = obj.P.A;
				case 'k'
					out = obj.P.b;
				otherwise
					error('Unrecognized property "%s".');
			end
		end
		
		function [P, keptrows, feasible] = domain(obj, A, f, Q, varargin)
			% Computes polytope that is mapped to an another polytope using affine map

			if nargin<3
				f = zeros(size(A, 1), 1);
			end
			P = polytope(obj.P.invAffineMap(A, f));
			if nargin>2
				P = P.intersect(Q);
			end
			keptrows = [];
			feasible = ~P.P.isEmptySet();
		end
		
		function R = range(obj, A, f, varargin)
			% Forward affine transformation
			
			R = polytope(obj.P.affineMap(A));
			if nargin>2
				R = R + f;
			end
		end
		
		function R = mtimes(arg1, arg2, varargin)
			% Multiplication of polytopes
			
			if isa(arg1, 'polytope')
				arg1 = arg1.P;
			end
			if isa(arg2, 'polytope')
				arg2 = arg2.P;
			end
			R = polytope(arg1*arg2);
		end
		
		function [dist, x, y] = distance(obj, arg, varargin)
			% Distance between two polytopes
			
			if isa(arg, 'polytope')
				arg = arg.P;
			end
			z = obj.P.distance(arg);
			dist = z.dist;
			x = z.x;
			y = z.y;
		end
		
		function out = nconstr(obj)
			% Returns number of inequality constraints defining the polytope
			
			out = size([obj.P.H; obj.P.He], 1);
		end
		
		function out = slice(obj, varargin)
			% Slice of a polytope
			
			out = polytope(obj.P.slice(varargin{:}));
		end
		
		function out = triangulate(obj, varargin)
			% Triangulation of a polytope
			
			out = polytope(obj.P.triangulate());
		end
			
		function varargout = pelemfun(fname, obj, varargin)
			% Apply a function to each element of a polytope array
			
			nout = nargout;
			if length(obj.P)==1
				out = cell(1, nout);
				[out{:}] = feval(fname, obj, varargin{:});
				for i = 1:nout
					varargout{i}{1} = out{i};
				end
			else
				lenP = length(obj.P);
				varargout = cell(1, nout);
				for i = 1:nout
					varargout{i} = cell(1, lenP);
				end
				for i = 1:lenP
					out = cell(1, nout);
					[out{:}] = feval(fname, polytope(obj.P(i)), varargin{:});
					for j = 1:nout
						varargout{j}{i} = out{j};
					end
				end
			end
		end
		
		function [nx, ny] = size(obj, a)
			% Size method behaving as in MTP2
			
			nx = 1;
			ny = length(obj.P);
			if nargin==2 && a==2
				nx = ny;
			end
			if nargout==1
				nx = [nx ny];
			end
		end
		
		function H = hull(obj)
			% Convex hull
			
			PU = PolyUnion('Set', toPolyhedron(obj));
			H = polytope(PU.convexHull);
        end
        
        function E = envelope(obj)
            % Convex envelope
            
            PU = PolyUnion(toPolyhedron(obj));
            E = polytope(PU.envelope());
        end
        
        function Q = uminus(obj)
            % Unary minus
            
            Q = polytope(-obj.P);
        end
	end
	
	methods(Static)
		function out = loadobj(obj)
			% Loads polytopes
			
			if ~isfield(obj, 'Array')
				% unitialized object
				out = polytope;
			elseif isempty(obj.Array)
				out = polytope(obj.H, obj.K);
			else
				out = [];
				for i = 1:length(obj.Array)
					out = [out obj.Array{i}];
				end
				fprintf('\nThe @polytope class will be deprecated, use @Polyhedron instead.\n\n');
			end
		end
	end
end
