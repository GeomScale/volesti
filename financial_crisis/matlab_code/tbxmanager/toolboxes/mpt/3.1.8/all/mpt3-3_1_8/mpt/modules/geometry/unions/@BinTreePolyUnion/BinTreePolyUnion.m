classdef BinTreePolyUnion < PolyUnion
	% Class representing memory-optimized binary search trees

	properties (SetAccess=protected)
		Tree
	end
	
	methods
		
		function obj = BinTreePolyUnion(in)
			% Constructs a memory-optimized binary search tree
			%
			% Syntax:
			%   T = BinTreePolyUnion(P)
			%
			% Inputs:
			%   P: single PolyUnion object (all sets fully dimensional)
			% Outputs:
			%   T: BinTreePolyUnion object
			%
			% The linear representation of the binary search tree is
			% available in "T.Tree" as a matrix with the following format:
			%   T.Tree(i, :) = [hp, index_left, index_right)
			%     hp: hyperplane to test
			%     index_left: next row to visit if hp*[x; -1]<=0
			%     index_right: next row to visit if hp*[x; -1]>0
			% A negative index represents a leaf node, in which case its
			% absolute value denotes index of the region which contains the
			% query point "x".
			%
			% Additional information about the tree (depth, number of
			% nodes, structure-based representation of the tree, etc.) is
			% available in "T.Internal.BinaryTree".
			%
			% Point location in binary trees is performed by
			%   [isin, inwhich] = T.contains(x)
			%
			% Note that the tree discards "outer" boundaries of regions.
			% Hence T.contains(x) always returns true even if the point
			% lies outside of the convex hull.
			%
			% Note that the tree is constructed in a recursive fashion.
			% Since Matlab imposes an upper limit on the number of
			% recursive calls, in rare circumstances an error can be
			% trigger. In that case simply increase the recursion limit and
			% re-run the tree construction.
			
			if nargin==0
				% empty constructor (used in copying)
				return
			end
			
			%% validation
			narginchk(1, 1);
			if ~isa(in, 'PolyUnion') || numel(in)~=1
				error('Input must be a single PolyUnion object.');
			elseif in.Num<1
				error('The input must be a non-empty PolyUnion.');
			end
			% all sets must be fully dimensional
			if ~all([in.Set.isFullDim()])
				error('All sets must be fully dimensional.');
			end

			%% set data
			% all sets must be in minimal H-representation
			in.Set.minHRep();
			% always operate on a copy of the input polyunion
			obj.Set = in.Set.copy();
			% normalize (work with the copy)
			obj.Set.normalize();
			% set internal properties
			obj.Internal = in.Internal;
			obj.Data = in.Data;
			obj.Dim = in.Dim;
			
			%% compute the tree
			obj.constructTree();
			
			% TODO: we should explicitly disallow adding new sets to the
			% polyunion, since doing so would require re-construction of
			% the tree. Same applies to any other method which modifies the
			% underlying set.
			
			% TODO: optionally include boundaries of the convex hull
			
		end
		
		function [isin, inwhich] = contains(obj, x)
			% Point location using binary search trees

			%% validation and error checks
			narginchk(2, 2);
			% use obj.forEach(@(u) u.contains(x)) to evaluate arrays
			error(obj.rejectArray());
			isin = false;
			inwhich = [];
			if numel(obj)==0 || ( numel(obj)==1 && obj.Num==0 )
				% empty object
				return
			end
			error(validate_vector(x, obj.Dim, 'point'));
			
			%% search
			idx = 1;
			while true
				hpx = obj.Tree(idx, 1:obj.Dim+1)*[x; -1];
				if hpx<=0
					idx = obj.Tree(idx, end-1);
				else
					idx = obj.Tree(idx, end);
				end
				if idx<0
					% reached leaf node
					isin = true;
					inwhich = -idx;
					return
				elseif idx==0
					% infeasible
					return
				end
			end
			
		end
		
		function toC(obj, function_name, file_name)
			% Exports the binary tree to a C-code
			%
			%   tree.toC('function', 'output')
			
			global MPTOPTIONS
			
			narginchk(2, 3);
            error(obj.rejectArray());
            
            if ~ischar(function_name)
                error('The function name must be given as a string.');
            end
            if nargin<3,
                file_name = 'mpt_getInput';
            else
                if isempty(file_name)
                    file_name = 'mpt_getInput';
                end
                if ~ischar(file_name)
                    error('The file name must be given as a string.');
                end
            end
            
            % get the short name if the full path is provided
            [path_to_file,short_name] = fileparts(file_name);
            if ~isempty(regexp(short_name,'[^a-zA-Z0-9_]', 'once'))
                error('The file name must contain only alphanumerical characters including underscore "_".');
            end
            % correct the file name with the proper extenstion ".c"
            if isempty(path_to_file)
                file_name = [short_name,'.c'];
                mex_name = [short_name,'_mex.c'];
            else
                file_name = [path_to_file, filesep, short_name,'.c'];
                mex_name = [path_to_file, filesep, short_name,'_mex.c'];
            end
            
			% is the request function present?
			if ~obj.hasFunction(function_name)
				error('No such function "%s" in the object.', function_name);
			end
			% support only for PWA/PWQ functions
            if ~(isa(obj.Set(1).Functions(function_name), 'AffFunction') || ...
                    isa(obj.Set(1).Functions(function_name), 'QuadFunction') )
                error('Only quadratic and affine functions are supported.');
            end
            isquadratic = isa(obj.Set(1).Functions(function_name), 'QuadFunction');
			
            % single or double precision to export?
            precision = MPTOPTIONS.modules.geometry.unions.BinTreePolyUnion.toC.precision;
            precision = strtrim(lower(precision));
            if isempty(precision)
                precision = 'double';
            elseif ~isequal(precision,'single') && ~isequal(precision,'double')
                error('The specified precision in the option can be either "single" or "double".');
            end
            if isequal(precision,'single')
                precision = 'float';
            end
            
			fun = obj.Set(1).Functions(function_name);

            %% write to a file
            outfid = fopen(file_name, 'w');
            if outfid < 0,
                error('Cannot open file "%s" for writing!', file_name);
            end
            
			fprintf(outfid,'/*  Identifies a control law associated to a given state X using a binary search tree.\n\n');
            fprintf(outfid,'  Usage:\n  long region = %s( %s *X, %s *U)\n\n',short_name, precision, precision);
                
            % header    
            header={''
                '   where X is a vector of dimension MPT_DOMAIN and U is a vector of dimension '
                '   MPT_RANGE for PWA functions. The output variable "region" indicates index of'
                '   a region where the point X is located. If "region" index is smaller than 1 '
                '   (region < 1), there is no control law associated to a given state.'
                ''
                ''
                '   Please note that all code in this file is provided under the terms of the'
                '   GNU General Public License, which implies that if you include it directly'
                '   into your commercial application, you will need to comply with the license.'
                '   If you feel this is not a good solution for you or your company, feel free '
                '   to contact me at michal.kvasnica@stuba.sk, I can re-license this specific '
                '   piece of code to you free of charge.'
                ''
                ''
                '   Copyright (C) 2006-2013 by Michal Kvasnica (michal.kvasnica@stuba.sk) '
                '   Revised in 2013 by Martin Herceg, (herceg@control.ee.ethz.ch)'
                ''
                ''
                '  This program is free software; you can redistribute it and/or modify'
                '  it under the terms of the GNU General Public License as published by'
                '  the Free Software Foundation; either version 2 of the License, or'
                '  (at your option) any later version.'
                ''
                '  This program is distributed in the hope that it will be useful,'
                '  but WITHOUT ANY WARRANTY; without even the implied warranty of'
                '  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the'
                '  GNU General Public License for more details.'
                ''
                '  You should have received a copy of the GNU General Public License'
                '  along with this program; if not, write to the Free Software'
                '  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.'
                '*/ '
                ''};

            for i=1:numel(header)
                fprintf(outfid,[header{i},'\n']);
            end
            
            fprintf(outfid,'/* Generated on %s by MPT %s */ \n\n', datestr(now), MPTOPTIONS.version);
		
			% data of the search tree
			ST = obj.Tree'; ST = ST(:);
            if isquadratic
                out = ''; 
            else
                out = char(sprintf('#define MPT_RANGE %d', fun.R));
            end
			out = char(out, sprintf('#define MPT_DOMAIN %d', fun.D));
			out = char(out, sprintf('static %s MPT_ST[] = {', precision));
			temp = '';
			for i = 1:length(ST),
                if isequal(precision,'float')
                    temp = [temp sprintf('%.7e,\t', ST(i))];
                else
                    temp = [temp sprintf('%.14e,\t', ST(i))];
                end
				if mod(i, 5)==0 || i==length(ST),
					out = char(out, temp);
					temp = '';
				end
			end
			out = char(out, '};');
			
			% data of the function

            % quadratic term
            if isquadratic
                out = char(out, sprintf('static %s MPT_H[] = {', precision));
                for i = 1:obj.Num
                    H = obj.Set(i).Functions(function_name).H;
                    for j = 1:fun.D
                        h = H(j, :);
                        temp = '';
                        for k = 1:length(h),
                            if isequal(precision,'float')
                                temp = [temp sprintf('%.7e,\t', h(k))];
                            else
                                temp = [temp sprintf('%.14e,\t', h(k))];
                            end
                        end
                        out = char(out, temp);
                    end
                end
                out = char(out, '};');
            end
            
			% linear term
			out = char(out, sprintf('static %s MPT_F[] = {', precision));
			for i = 1:obj.Num
				F = obj.Set(i).Functions(function_name).F;
				for j = 1:fun.R
					f = F(j, :);
					temp = '';
					for k = 1:length(f),
                        if isequal(precision,'float')
                            temp = [temp sprintf('%.7e,\t', f(k))];
                        else
                            temp = [temp sprintf('%.14e,\t', f(k))];
                        end
					end
					out = char(out, temp);
				end
			end
			out = char(out, '};');
			
			% constant term
			out = char(out, sprintf('static %s MPT_G[] = {', precision));
			for i = 1:obj.Num,
				g = obj.Set(i).Functions(function_name).g;
				for j = 1:fun.R,
					f = g(j, :);
                    temp = '';
                    for k = 1:length(f),
                        if isequal(precision,'float')
                            temp = [temp sprintf('%.7e,\t', f(k))];
                        else
                            temp = [temp sprintf('%.14e,\t', f(k))];
                        end
                    end
					out = char(out, temp);
				end
			end
			out = char(out, '};');

			% convert into a single string, add line breaks
			out_nl = [];
            for i = 1:size(out, 1),
                out_nl = [out_nl sprintf('%s\n', deblank(out(i, :)))];
            end

            % write the search tree
            fprintf(outfid, out_nl);

            % write the function for evaluation of the binary tree
            fprintf(outfid,'static long %s(const %s *X, %s *U)\n', short_name, precision, precision);
            
            footer = {''
                '{'
                '    long node = 1, row;'
                };
            if isquadratic
                footer = [footer; {
                    '    int ix, jx;'
                    sprintf('    %s hx, sx, k;',precision)
                    ''
                    }];                
            else
                footer = [footer; {
                    '    int ix, iu;'
                    sprintf('    %s hx, k;',precision)
                    ''
                    '    /* initialize U to zero*/'
                    '    for (iu=0; iu<MPT_RANGE; iu++) {'
                    '        U[iu] = 0;'
                    '    }'
                    ''
                    }];
            end
            footer = [footer; {
                '    /* find region which contains the state x0 */'
                '    while (node > 0) {'
                '        hx = 0;'
                '        row = (node-1)*(MPT_DOMAIN+3);'
                '        for (ix=0; ix<MPT_DOMAIN; ix++) {'
                '            hx += MPT_ST[row+ix]*X[ix];'
                '        }'
                '        k = MPT_ST[row+MPT_DOMAIN];'
                ''
                '        if ((hx - k) <= 0) {'
                '            /* x0 on the negative side of the hyperplane */'
                '            node = (long)MPT_ST[row+MPT_DOMAIN+1];'
                '        } else {'
                '            /* x0 on positive side the hyperplane */'
                '            node = (long)MPT_ST[row+MPT_DOMAIN+2];'
                '        }'
                '    }'
                ''
                '    node = -node;'
                ''
                }];
            if isquadratic
                footer = [footer; {
                    '    /* compute control action associated to state x0 */'
                    '    U[0] = 0;'
                    '    for (ix=0; ix<MPT_DOMAIN; ix++) {'
                    '        sx = 0;'
                    '        for (jx=0; jx<MPT_DOMAIN; jx++) {'
                    '             sx += MPT_H[(node-1)*MPT_DOMAIN*MPT_DOMAIN + ix*MPT_DOMAIN + jx]*X[jx];'
                    '        }'
                    '        U[0] += sx*X[ix];'
                    '    }'
                    '    for (ix=0; ix<MPT_DOMAIN; ix++) {'
                    '        U[0] += MPT_F[(node-1)*MPT_DOMAIN + ix]*X[ix];'
                    '    }'
                    '    U[0] += MPT_G[node-1];'                    
                }];
            else
                footer = [footer; {
                    '    /* compute control action associated to state x0 */'
                    '    for (iu=0; iu<MPT_RANGE; iu++) {'
                    '        for (ix=0; ix<MPT_DOMAIN; ix++) {'
                    '            U[iu] += MPT_F[(node-1)*MPT_DOMAIN*MPT_RANGE + iu*MPT_DOMAIN + ix]*X[ix];'
                    '        }'
                    '        U[iu] += MPT_G[(node-1)*MPT_RANGE + iu];'
                    '    }'
                    }];
            end
            footer = [footer; {
                ''
                '    return node;'
                '}'
                }];
            
            % write the footer
            for i=1:numel(footer)
                fprintf(outfid,[footer{i},'\n']);
            end

			fclose(outfid);
            fprintf('Output written to "%s".\n', file_name);
            
            %% write mpt_getInput_mex.c
            
            fdn = fopen(mex_name, 'w');
            if fdn<0,
                error('Cannot open file "%s" for writing!',mex_name);
            end
            
            % write the mex-file for evaluation under Matlab
            cmex={''
                '/*'
                '  Autogenerated C-mex file for evalution of explicit controllers.'
                ''
                '  This file is to be compiled under Matlab using the syntax'
                ''
                sprintf('   mex -largeArrayDims %s.c',[short_name,'_mex'])
                ''
                '  Usage in Matlab:'
                ''
                sprintf('   u = %s(x)', [short_name,'_mex'])
                ''
                '  where x is the vector of double for which to evaluate PWA/PWQ function.'
                ''
                ''
                '  Copyright 2013 by Martin Herceg, Automatic Control Laboratory,'
                '  ETH Zurich, herceg@control.ee.ethz.ch'
                ''
                ''
                '  This program is free software; you can redistribute it and/or modify'
                '  it under the terms of the GNU General Public License as published by'
                '  the Free Software Foundation; either version 2 of the License, or'
                '  (at your option) any later version.'
                ''
                '  This program is distributed in the hope that it will be useful,'
                '  but WITHOUT ANY WARRANTY; without even the implied warranty of'
                '  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the'
                '  GNU General Public License for more details.'
                ''
                '  You should have received a copy of the GNU General Public License'
                '  along with this program; if not, write to the Free Software'
                '  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.'
                '*/ '
                ''
                ''
                sprintf('/* Generated on %s by MPT %s */ \n\n',datestr(now), MPTOPTIONS.version)
                ''
                ''
                '#include <math.h>'
                '#include "mex.h"'
                sprintf('#include "%s.c"',short_name)
                ''
                '#ifndef MPT_RANGE'
                '#define MPT_RANGE 1'
                '#endif'
                ''
                ''
                'void mexFunction( int nlhs, mxArray *plhs[], '
                '		  int nrhs, const mxArray*prhs[] )'
                ''
                '{ '
                '    int j;'
                '    double *xin, *uout;'
                '    char msg[70];'
                '    long region;'
                ''
                '    mwSize M,N,D; '
                };
            if isequal(precision,'float')
                cmex = [cmex; '    float x[MPT_DOMAIN], u[MPT_RANGE];'];
            end
            cmex = [cmex; {
                ''
                '    /* Check for proper number of arguments */'
                ''
                '    if (nrhs != 1) { '
                '	    mexErrMsgTxt("One input arguments required."); '
                '    } else if (nlhs > 1) {'
                '	    mexErrMsgTxt("Too many output arguments."); '
                '    } '
                '    if (mxIsEmpty(prhs[0]))'
                '        mexErrMsgTxt("The argument must not be empty.");'
                '	if ( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) )'
                '        mexErrMsgTxt("The argument must be real and of type DOUBLE.");'
                '    if ( mxGetNumberOfDimensions(prhs[0])!=2 )'
                '		mexErrMsgTxt("The argument must be a vector"); '
                ''
                '    M = mxGetM(prhs[0]);'
                '    N = mxGetN(prhs[0]);'
                '    if ( (M!=1) && (N!=1) )'
                '		mexErrMsgTxt("The argument must be a vector.");'
                ''
                '    /* dimension */'
                '    D = mxGetNumberOfElements(prhs[0]);'
                '    if (D!=MPT_DOMAIN) {'
                '        sprintf(msg, "The size of the input argument must be %%d x 1.",MPT_DOMAIN);'
                '        mexErrMsgTxt(msg);'
                '    }'
                ''
                '    /* get and verify input data */'
                '    xin = mxGetPr(prhs[0]);   '
                '    for (j=0; j<D; j++) {'
                '        if ( mxIsNaN(xin[j]) || mxIsInf(xin[j]) )'
                '            mexErrMsgTxt("No ''NaN'' or ''Inf'' terms are allowed.");'
                '    }'
                ''
                '    /* create a matrix for the output */ '
                '    plhs[0] = mxCreateDoubleMatrix( MPT_RANGE, 1, mxREAL);    '
                '    uout = mxGetPr(plhs[0]);'
                ''
                }];
            if isequal(precision,'float')
                cmex = [cmex; {
                    '    /* recast double as float */'
                    '    for (j=0; j<MPT_DOMAIN; j++){'
                    '        x[j] = (float)xin[j];'
                    '    }'
                    ''
                    '    /* Do the evaluation in a subroutine */'
                    sprintf('    region = %s(x, u);', short_name)
                    ''
                    }];
            else
                cmex = [cmex; {
                    ''
                    '    /* Do the evaluation in a subroutine */'
                    sprintf('    region = %s(xin, uout);', short_name)
                    ''
                    }];
            end
            if isequal(precision,'float')
                cmex = [cmex; {
                    '    /* cast as double for output */'
                    '    for (j=0; j<MPT_RANGE; j++){'
                    '        uout[j] = (double)u[j];'
                    '    }'
                    ''}];
            end
            cmex = [cmex; {
                ''
                '    return;'
                ''
                '}'
                }];
            
            % write the core
            for i=1:numel(cmex)
                fprintf(fdn,[cmex{i},'\n']);
            end
            
            fclose(fdn);
            fprintf('C-mex function written to "%s".\n', mex_name);
            
            
		end
		
		function display(obj)
			% Display method for BinTreePolyUnion objects
			
			% call super-class'es display() method first
			obj.display@PolyUnion();
			fprintf('Memory-optimized binary tree, depth: %d, no. of nodes: %d\n', ...
				obj.Internal.BinaryTree.depth, ...
				obj.Internal.BinaryTree.n_nodes);
		end
		
	end
	
	methods (Access=private)

		function constructTree(obj)
			% Constructs the binary search tree
			
			global MPTOPTIONS
			r_coef = MPTOPTIONS.modules.geometry.unions.BinTreePolyUnion.round_places;
			start_time = clock;
			
			% Find all unique half-spaces (remember the sets are already in
			% minimal, normalized H-representation)
			%
			% TODO: if the union is convex, we can discard outer boundaries
			% as they will never be optimal separating hyperplanes. Just do
			% not forget to add them back afterward.
			H = cat(1, obj.Set.H);
			
			% find unique "a" from a'*x<=0 and a'*x>=0
			H = round(H*10^r_coef)/10^r_coef;
			nH_orig = size(H, 1);
			H = unique(H, 'rows');
			[i, j] = ismember(H, -H, 'rows');
			idx1 = find(i==0);
			idx2 = j(j > (1:length(j))');
			H = H(unique([idx1; idx2]), :);
			nH = size(H, 1);
			fprintf('Found %d unique hyperplanes (out of %d)\n', nH, nH_orig);
			
			% Determine position of each region w.r.t. each half-space
			fprintf('Determining position of %d regions w.r.t. unique hyperplanes...\n', ...
				obj.Num);
			Neg = false(nH, obj.Num);
			Pos = false(nH, obj.Num);
            % tolerance for checking full dimensionality of intersections
            tol = 10*MPTOPTIONS.rel_tol;
			tic
			for i = 1:nH
				if toc>MPTOPTIONS.report_period
					fprintf('Progress: %d/%d\n', i, nH);
					tic
				end
				[Neg(i, :), Pos(i, :)] = obj.getPosition(H(i, :), tol);
            end
			obj.Internal.BinaryTree.runtime.H = etime(clock, start_time);
			fprintf('...done in %.1f seconds.\n', ...
				obj.Internal.BinaryTree.runtime.H);
			
			% Discard half-spaces which are satisfied by all regions
			keep = find(sum(Pos, 2)>=1);
			fprintf('Discarding %d outer boundaries.\n', nH-length(keep));
			H = H(keep, :);
			Neg = Neg(keep, :);
			Pos = Pos(keep, :);
			fprintf('Considering %d candidates for separating hyperplanes.\n', size(H, 1));
			
			% Store the data in the object
			obj.Internal.BinaryTree.H = H;
			obj.Internal.BinaryTree.Neg = Neg;
			obj.Internal.BinaryTree.Pos = Pos;
			obj.Internal.BinaryTree.n_nodes = 0;
			obj.Internal.BinaryTree.depth = 1;
			% number of open nodes (only for progress reports)
			obj.Internal.BinaryTree.n_open = 0;
			
			% linear representation of the tree
			info = obj.emptyNode();
			obj.Internal.BinaryTree.Linear = info([]);
			
			% construct the tree
			fprintf('Constructing the tree...\n');
			tic
			start_time = clock;
			root = obj.emptyNode();
			root.hp = zeros(1, obj.Dim+1);
			root.depth = 0;
			root.parent = 0;
			root.index = 0;
			obj.Internal.BinaryTree.Tree = obj.createNode(root, 1:obj.Num);
			obj.Internal.BinaryTree.runtime.tree = etime(clock, start_time);
			fprintf('...done in %.1f seconds.\n', ...
				obj.Internal.BinaryTree.runtime.tree);

			if isempty(H)
				% trivial case: no candidates
				obj.Tree = [zeros(1, obj.Dim+1), 1, 1];
			else
				% put the tree into a matrix
				obj.Tree = zeros(length(obj.Internal.BinaryTree.Linear), ...
					obj.Dim+1+2);
				for i = 1:length(obj.Internal.BinaryTree.Linear)
					idx = obj.Internal.BinaryTree.Linear(i).index;
					obj.Tree(idx, :) = [ obj.Internal.BinaryTree.Linear(i).hp, ...
						obj.Internal.BinaryTree.Linear(i).left, ...
						obj.Internal.BinaryTree.Linear(i).right ];
				end
			end
			
			fprintf('Depth: %d, no. of nodes: %d\n', ...
				obj.Internal.BinaryTree.depth, ...
				obj.Internal.BinaryTree.n_nodes);
		end

		function [negative, positive] = getPosition(obj, hp, tol)
			% Determines location of regions w.r.t. hyperplane hp*[x; -1]=0:
			%  negative(i)=1 if the i-th region has full-dimensional
			%                intersection with {x | hp*[x; -1] <= 0}
			%  positive(i)=1 if the i-th region has full-dimensional
			%                intersection with {x | hp*[x; -1] >= 0}
			
			% TODO: exploit nested half-spaces, i.e., if { x | a*x<=b }
			% \subset { x | c*x<=d }, then any region on the negative side
			% of "a, b" will also be on the negative side of "c, d".
			negative = false(1, obj.Num);
			positive = false(1, obj.Num);
			for i = 1:obj.Num
                [s, r] = fast_isFullDim([obj.Set(i).H; hp]);
                if s && -r.obj > tol
					% full dimensional intersection with {x | hp*[x;-1]<=0} 
                    % with the chebyradius greater than tol
                    negative(i) = true;
                end
                [s, r] = fast_isFullDim([obj.Set(i).H; -hp]);
                if s && -r.obj > tol
                    % full dimensional intersection with {x | hp*[x;-1]>=0}
                    % with the chebyradius greater than tol
                    positive(i) = true;
                end
			end
			
		end
		
		function new = createNode(obj, node, Idx)
			% Creates a node of the binary search tree
			
			global MPTOPTIONS
			
			if isempty(Idx)
				% no more regions, terminate this branch
				new = [];
				return
			end

			obj.Internal.BinaryTree.n_open = obj.Internal.BinaryTree.n_open+1;
			
			% report from time to time
			if toc>MPTOPTIONS.report_period
				tic
				fprintf('Progress: depth: %d, no. of nodes: %d, open: %d\n', ...
					obj.Internal.BinaryTree.depth, ...
					obj.Internal.BinaryTree.n_nodes, ...
					obj.Internal.BinaryTree.n_open);
			end
			
			% Select the best hyperplane, i.e., the one which minimizes the
			% maximal number of splitted regions
			D = Inf(size(obj.Internal.BinaryTree.Neg, 1), 1);
			for i = setdiff(1:size(D, 1), abs(node.hps))
				% regions on the negative side of each hyperplane
				NegIdx = intersect(find(obj.Internal.BinaryTree.Neg(i, :)), Idx);
				PosIdx = intersect(find(obj.Internal.BinaryTree.Pos(i, :)), Idx);
				% fitness of the i-th hyperplane
                if isempty(NegIdx) || isempty(PosIdx)
                    % exclude hyperplanes which do not split the space
                    D(i) = Inf;
                else
                    D(i) = max(numel(NegIdx), numel(PosIdx));
                end
			end
			% which hyperplane is the best cut?
			[~, BestIdx] = min(D(:, 1));

			% Find which regions are on the left and on the right of the best
			% hyperplane
			NegIdx = intersect(find(obj.Internal.BinaryTree.Neg(BestIdx, :)), Idx);
			PosIdx = intersect(find(obj.Internal.BinaryTree.Pos(BestIdx, :)), Idx);
			if MPTOPTIONS.verbose>0
				fprintf('\n Idx: %s\n', mat2str(Idx));
				fprintf('-Idx: %s\n', mat2str(NegIdx));
				fprintf('+Idx: %s\n', mat2str(PosIdx));
			end

			% create a new node
			new = obj.emptyNode();
			new.hp = obj.Internal.BinaryTree.H(BestIdx, :);
			new.hps = node.hps;
			new.parent = node.index;
			new.depth = node.depth+1;
			new.index = obj.Internal.BinaryTree.n_nodes+1;

			% update tree status
			obj.Internal.BinaryTree.n_nodes = obj.Internal.BinaryTree.n_nodes+1;
			obj.Internal.BinaryTree.depth = max(obj.Internal.BinaryTree.depth, new.depth);

			% which regions are on the negative side of the hyperplane?
			if numel(NegIdx)==1
				% single region = leaf node
				new.left = NegIdx;
			elseif numel(NegIdx)>1
				% multiple regions = split this node
				new_neg = new;
				new_neg.hps = [new_neg.hps; -BestIdx];
				new.left = obj.createNode(new_neg, NegIdx);
			else
				% empty on the left
				new.left = 0;
			end
			
			% which regions are on the positive side of the hyperplane?
			if numel(PosIdx)==1
				% single region = leaf node
				new.right = PosIdx;
			elseif numel(PosIdx)>1
				% multiple regions = split this node
				new_pos = new;
				new_pos.hps = [new_pos.hps; BestIdx];
				new.right = obj.createNode(new_pos, PosIdx);
			else
				% empty on the right
				new.right = 0;
			end
			
			% update the linear structure of the tree
			linear = obj.emptyNode();
			linear.index = new.index;
			linear.hp = new.hp;
			linear.hps = new.hps;
			linear.depth = new.depth;
			linear.parent = new.parent;
			if isstruct(new.left)
				linear.left = new.left.index;
			else
				linear.left = -new.left;
			end
			if isstruct(new.right)
				linear.right = new.right.index;
			else
				linear.right = -new.right;
			end
			obj.Internal.BinaryTree.Linear(end+1) = linear;
			obj.Internal.BinaryTree.n_open = obj.Internal.BinaryTree.n_open-1;
			
		end
	end
	
	methods (Static, Hidden)
		
		function N = emptyNode()
			% Creates an empty node structure
			
			N = struct('hp', [], 'hps', [], ...
				'depth', 0, 'parent', 0, 'index', 0, ...
				'left', [], 'right', []);
		end
	end
	
end

