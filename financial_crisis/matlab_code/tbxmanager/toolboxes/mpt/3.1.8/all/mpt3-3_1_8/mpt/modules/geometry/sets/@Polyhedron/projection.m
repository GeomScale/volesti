function q = projection(P, dims, method, solver)
%
%  PROJECTION: Compute the projection of the Polyhedron. 
%  ======================================================
%  
%  
%  SYNTAX
%  ------
%     
%      p = P.projection(dims, method)
%      p = projection(P, dims, method)
%      p = projection(P, dims, method, solver)
%    
%  
%  DESCRIPTION
%  -----------
%     Computes the polyhedron p 
%                           p = {x | x = y(dims), y in P} 
%  
%  
%  INPUT
%  -----
%     
%        
%          P      Polyhedron in any format                 
%                 Class: Polyhedron                        
%          dims   Dimensions upon which to project         
%                 Class: double                            
%          method Sets the method used to compute the      
%                 projection. If omitted, then the method  
%                 is chosen based on the properties of P   
%                 and dims.                                
%                 Class: char                              
%                 Allowed values:                          
%                                                          
%                   fourier  : Use fourier elimination.    
%                    Good if projecting off a small number 
%                    of dimensions. Projection will be     
%                    highly redundant (use q.minHRep() to  
%                    get irredundant form).                
%                   ifourier  : Fourier elimination with   
%                    intermediate redundancy elimination.  
%                    Dimensions are projected off one by   
%                    one, followed by removal of redundant 
%                    constraints. The result is still      
%                    redundant, though.                    
%                   mplp  : Use mplp algorithm to compute  
%                    projection. Good if the projection is 
%                    not very degenerate (dimensions of    
%                    projected facets are equal to the     
%                    faces that they were projected from)  
%                   vrep  : Compute vertex representation  
%                    and then project                      
%                                                          
%          solver If "mplp" method is selected, then this  
%                 argument determines which solver to use  
%                 for solving MPLP problem. By default the 
%                 first parametric solver is selected from 
%                 the list of available solvers.           
%                 Class: char                              
%                   
%  
%  
%  OUTPUT
%  ------
%     
%        
%          p Projection of P  onto the dimensions     
%            dims                                     
%            Class: Polyhedron                        
%              
%  
%  
%  SEE ALSO
%  --------
%     project
%  

%  AUTHOR(s)
%  ---------
%     
%    
%   (c) 2010-2013  Colin Neil Jones: EPF Lausanne
%   mailto:colin.jones@epfl.ch 
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

% Choose fourier elimination if we're projecting off less than this number
% of dimensions
DELTA_DIMS_CHOOSE_FOURIER = 3;

if nargin < 3, method = ''; solver = ''; end
if nargin < 4, solver = ''; end

if isempty(method)
	% read default settings from global options. Note that the default can
	% still be empty, in which case we perform an automatic selection.
	method = MPTOPTIONS.modules.geometry.sets.Polyhedron.projection.method;
end

% deal with arrays
no = numel(P);
if no>1
	q(size(P)) = Polyhedron;
	for i=1:no
		q(i) = P(i).projection(dims, method, solver);
	end
	return;
end

% check dimensions
for i=1:numel(dims)
	validate_dimension(dims(i));
end
% no dimensions or no projection, copy the polyhedron
if isempty(dims) || isequal(dims, 1:P.Dim)
	q = P.copy();
	return;
end
if any(dims>P.Dim)
	mpt_kblink(2); % provide the ID of the corresponding KB article
	error('Cannot compute projection on higher dimension than %i.', P.Dim);
end

if P.isEmptySet,
	q = Polyhedron.emptySet(length(dims));
	return;
elseif P.isFullSpace()
	q = Polyhedron.fullSpace(length(dims));
	return
end

% Compute projection of the polyhedron

if P.hasVRep
	% Easy way - V-rep
	q = Polyhedron('V', P.V(:,dims), 'R', P.R(:,dims));
else
	
	% Hard way - inequality rep
	if isempty(method) % No method specified
		% choose based on properties of P
		if P.Dim <= 4
			% this call can be problematic due to computations of vertices in CDD,
			%method = 'vrep';
			method = 'fourier';
		elseif P.Dim - length(dims) <= DELTA_DIMS_CHOOSE_FOURIER
			method = 'fourier';
        elseif P.isFullDim()
            % must have a fully dimensional polyhedron
			method = 'mplp';
        else
            method = 'ifourier';
		end
		
	end
	
	%%
	if size(P.He,1) > 0
		% Use the affineMap function, which will convert into a projection of a
		% full-dimensional polyhedron and call back to here
		T = zeros(length(dims),P.Dim);
		for i=1:length(dims)
			T(i,dims(i)) = 1;
		end
		q = P.affineMap(T, method);
	else
		% shift polyhedron to origin for better numerics
		xc = P.chebyCenter;
		if xc.exitflag == MPTOPTIONS.OK
			Pn = P - xc.x;
		else
			Pn = Polyhedron(P);
		end
		
		% We can now safely assume that He is empty
		switch lower(method)
			
			case 'vrep'
				% Enumerate vertices and project
				Pn.minVRep();
				q = projection(Pn, dims);
				q.minHRep(); % return H-rep, because the entry was in H
				
			case 'fourier'
				Pn.minHRep(); % Redundancy elimination
				
				% we need to respect ordering of dimensions (issue #101)
				dims_eliminate = setdiff(1:Pn.Dim, dims);
				% reorder columns such that dimensions to keep are first
				H = Pn.H(:, [dims(:); dims_eliminate(:); Pn.Dim+1]');
				Hn = fourier(H, 1:length(dims));

				if isempty(Hn)
					% polyhedron is full R^P.dim
					q = Polyhedron.fullSpace(numel(dims));
				else
					q = Polyhedron(Hn(:, 1:end-1), Hn(:, end));
				end
				
			case 'ifourier'
				% fourier with intermediate redundancy elimination
				
				Pn.minHRep(); % Redundancy elimination
				
				% we need to respect ordering of dimensions (issue #101)
				dims_eliminate = setdiff(1:Pn.Dim, dims);
				% reorder columns such that dimensions to keep are first
				H = Pn.H(:, [dims(:); dims_eliminate(:); Pn.Dim+1]');
				
				% eliminate one dimension at a time
                t = tic;
				for i = 1:length(dims_eliminate)
                    if MPTOPTIONS.verbose > 0 && toc(t) > MPTOPTIONS.report_period
                        fprintf('progress: %d/%d\n', i, length(dims_eliminate));
                        t = tic;
                    end
					Hn = fourier(H, 1:Pn.Dim-i);
					if isempty(Hn)
						% polyhedron is full R^P.dim
						q = Polyhedron.fullSpace(numel(dims));
						break
					else
						% remove redundant constraints...
						q = Polyhedron(Hn(:, 1:end-1), Hn(:, end));
						if i < length(dims_eliminate)
							% ... but only if we are not finished yet
							q.minHRep();
						end
						H = q.H;
					end
				end
				
			case 'mplp'
				% On Polyhedral Projection and Parametric Programming
				%  by C.N. Jones, E.C. Kerrigan and J.M. Maciejowski
				%  J Optim Theory Appl (2008) 138: 207?220
				%  DOI 10.1007/s10957-008-9384-4
				
				% Shift P to contain the origin
				%         ip = P.interiorPoint;
				%         if ~ip.isStrict, error('Polyhedron should have a non-empty interior at this point in the code'); end
				if any(Pn.b < -MPTOPTIONS.abs_tol)
					error('HAVE NOT IMPLEMENTED PROJECTION FOR POLYHEDRON THAT DO NOT CONTAIN THE ORIGIN YET');
				end
				
				D = Pn.A; C = D(:,dims); D(:,dims) = [];
				d = size(C,2); k = size(D,2); n = size(C,1);
				
				% Setup parametric problem whose cost is the projection
				%dat = Opt;
				%dat.A  = [D -P.b];
				%dat.b  = zeros(n,1);
				%dat.pB = -C;
				%dat.f  = [zeros(k,1);1];
				%dat.Ath = [eye(d);-eye(d)];
				%dat.bth = ones(2*d,1);
				if isempty(solver)
					dat = Opt('A',[D -Pn.b],'b', zeros(n,1), 'pB',-C,'f',[zeros(k,1);1], 'Ath',[eye(d); -eye(d)], 'bth', ones(2*d,1));
				else
					% in case MPQP solver is provided, solve MPLP instead
					if strcmpi(solver,'MPQP')
						dat = Opt('A',[D -Pn.b],'b', zeros(n,1), 'pB',-C,'f',[zeros(k,1);1], 'Ath',[eye(d); -eye(d)], 'bth', ones(2*d,1),'solver', 'MPLP');
					else
						dat = Opt('A',[D -Pn.b],'b', zeros(n,1), 'pB',-C,'f',[zeros(k,1);1], 'Ath',[eye(d); -eye(d)], 'bth', ones(2*d,1),'solver', solver);
					end
                end
				
                if MPTOPTIONS.verbose > 0
                    % verbose output enabled
                    sol = dat.solve();
                else
                    % keep silent on verbose<=0
                    evalc('sol = dat.solve;');
                end
				
				% The normals of the facets is the cost function
				A = zeros(sol.xopt.Num,dat.d);
				for i=1:sol.xopt.Num
					objfun = sol.xopt.Set(i).getFunction('obj');
					A(i,:) = objfun.F;
				end
				q = Polyhedron(A,ones(size(A,1),1));
				
			otherwise
				error('Supported methods are "vrep", "fourier", "ifourier", and "mplp".');
		end
		
		% shift back to original coordinates
		if xc.exitflag == MPTOPTIONS.OK
			q = q+xc.x(dims);
		end
		
	end
end

end
