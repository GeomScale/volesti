function [obj, sol] = minHRep(obj)
%
%  MINHREP: Compute an irredundant H-representation of a polyhedron. 
%  ==================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      P.minHRep()
%      [P, sol] = P.minHRep()
%    
%  
%  DESCRIPTION
%  -----------
%     Computes an irredundant H-representation of the polyhedron: 
%                                         ---                   
%                    P = {x  |  Ax <= b } | | {x  |  A  x = b  }
%                                         | |         e      e  
%     Notes: 
%    
%     - If an H-representation is already known, then this function does redundancy
%     elimination. 
%     - Calling P.minHRep() will store the irredundant H-representation in P. 
%     - For empty polyhedron P the result remains the same (no redundancy
%     elimination). 
%  
%  
%  INPUT
%  -----
%     
%        
%          P Polyhedron in any format                 
%            Class: Polyhedron                        
%              
%  
%  
%  OUTPUT
%  ------
%     
%        
%          sol    Irredundant H-representation of P        
%                 Class: struct                            
%          sol.H  Matrix of inequalities                   
%                                [      ]                  
%                              H [  x   ]<= 0              
%                                [  -1  ]                  
%                                                          
%                 Class: double matrix                     
%          sol.He Matrix of equalities                     
%                                 [      ]                 
%                              H  [  x   ]= 0              
%                               e [  -1  ]                 
%                                                          
%                 Class: double matrix                     
%                   
%  
%  
%  SEE ALSO
%  --------
%     minVRep
%  

%  AUTHOR(s)
%  ---------
%     
%    
%   (c) 2010-2013  Colin Neil Jones: EPF Lausanne
%   mailto:colin.jones@epfl.ch 
%     
%    
%   (c) 2010-2013  Martin Herceg: ETH Zurich
%   mailto:herceg@control.ee.ethz.ch 
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
 
 
global MPTOPTIONS
if isempty(MPTOPTIONS)
	MPTOPTIONS = mptopt;
end

if numel(obj)>1
	% element-wise operation for arrays
	if nargout==0
		obj.forEach(@minHRep);
	elseif nargout==1
		obj = obj.forEach(@minHRep);
	else
		[obj, sol] = obj.forEach(@minHRep);
	end
	return
end

if ~obj.hasHRep
	% convert to Hrep if necessary
	obj.computeHRep();
end
if obj.irredundantHRep
	% nothing to do here
	if nargout>1
		sol.H = obj.H_int;
		sol.He = obj.He_int;
		sol.I = false(size(obj.H_int, 1), 1);
	end
	return
end

if isempty(obj.H_int) && isempty(obj.He_int)
	% still no Hrep available, probably R^n
	sol.H = [];
	sol.He = [];
	sol.I = [];
	return
end

% indices of redundant rows
nold = size(obj.H_int,1); % size of the original system
iredundant_rows = transpose(1:nold);

if obj.irredundantHRep == false
	% Do redundancy elimination

	% Empty (we cannot do much here, only run to numerical difficulties)
	if obj.isEmptySet
		obj.irredundantHRep = true;
		sol.H = obj.H_int;
		sol.He = obj.He_int;
		% logical index set of redundant rows
		sol.I = false(size(obj.H_int,1),1);
		return
	end

	% remove trivially redundant rows with 0*x <= b (note that at this
	% stage we know that the set is not emapy, hence b>=0)
	nA = matNorm(obj.H_int(:,1:end-1));
	zero_rows = nA < MPTOPTIONS.zero_tol;
	if any(zero_rows)
		iredundant_rows(zero_rows) = [];
		obj.H_int(zero_rows, :) = [];
	end
	
	% available heuristics
	% TODO: put these into MPTOPTIONS.modules.geometry.polyhedron.reduce
	bounding_box = true;
	rayshooting = true;
    
    H = obj.H_int;
    He = obj.He_int;
	if bounding_box && ~isempty(H)
		% -1 => Redundant
		%  0 => Unknown
		%  1 => Irredundant
		irr = zeros(size(H,1),1);
		
		A = H(:, 1:end-1);
		b = H(:, end);
		
		% Get the bounds from Polyhedron/outerApprox. This has the added
		% value of automatically storing the bounds in obj.Internal.lb/ub,
		% such that we can reuse them later.
		obj.outerApprox;
		lb = obj.Internal.lb;
		ub = obj.Internal.ub;
		
		% Any inequality that contains all the extreme points of the cube is redundant
		val = (A>0) .* A * (ub - lb) + A * lb - b;
		irr( val < -MPTOPTIONS.rel_tol ) = -1;
		
		% remove detected redundancies
		irem = (irr==-1);
		obj.H_int(irem,:) = []; % original system
		iredundant_rows(irem) = [];
		%irr(irem)     = [];
	end
	
	H = obj.H_int;
	% remove redundances based on uniqueness
	[~,i1] = unique(H,'rows');
	si1 = sort(i1);        % we must sort, otherwise we don't respect order of hyperplanes
	iredundant_rows = iredundant_rows(si1);
	obj.H_int = H(si1,:);
	
	if isempty(He) %%&& ~isempty(x0)
		if rayshooting
			x0 = obj.interiorPoint.x;
		end
		
		H = obj.H_int;
		A = H(:, 1:end-1);
		b = H(:, end);
		
		nc = size(H,1);
		
		if rayshooting
			% try ray-shooting heuristics for full-dimensional polyhedra
			cand=ones(nc,1);
			n = obj.Dim;
			r = ceil(nc/2);
			tol10 = 10*MPTOPTIONS.abs_tol;
			keep = zeros(1,2*r);
			
			for i = 1:r
				d = randn(n,1);
				d = d/norm(d);
				den = A*d;
				if any(norm(den,Inf)<MPTOPTIONS.rel_tol),
					continue
				end
				t = (b-A*x0)./(A*d);
				pos = t>0;
				negpos = ~pos;
				bound_pos = find(pos);
				bound_neg = find(negpos);
				tt = t(bound_pos);
				[t_closest_pos,j1] = min(tt);
				if isempty(j1),
					continue
				end
				if length(find(abs(t-t_closest_pos)<tol10))>1,
					continue
				end
				tt = t(bound_neg);
				[t_closest_neg,j2] = max(tt);
				if isempty(j2),
					continue
				end
				if length(find(abs(t-t_closest_neg)<tol10))>1,
					continue
				end
				if ~isempty(j1),
					keep(2*i-1)= bound_pos(j1);
				else
					t_closest_pos = 0;
				end
				if ~isempty(j2),
					keep(2*i)  = bound_neg(j2);
				else
					t_closest_neg = 0;
				end
				x1 = x0+d*t_closest_pos;
				x2 = x0+d*t_closest_neg;
				x0 = x1+(x2-x1)/2;
			end
			
			keep = unique(keep);
			
			% if "ismembc" is not found on your system, uncomment line below
			cand=(~ismembc(1:length(cand), keep))';
			%cand=(~ismember(1:length(cand), keep))';
		else
			cand = 1:nc;
		end
		
		nonredundant=true(nc,1);
		removerow=[];
		
		k=1;
		lp.f = [];
		lp.A = A;
		lp.b = b;
		lp.Ae = [];
		lp.be = [];
		lp.lb = [];
		lp.ub = [];
		lp.quicklp = true;
		
		while k<=length(cand)
			f_cand = find(cand);
			%if any(k==f_cand)
			if ismembc(k, f_cand)
				b(k) = b(k)+0.1;
				lp.f = -A(k,:);
				lp.A = A(nonredundant, :);
				lp.b = b(nonredundant);
				res = mpt_solve(lp);
				%res = mpt_solve(struct('f',-f1,'A',A(nonredundant,:),'b',b(nonredundant)));
				b(k) = b(k)-0.1;
				nonredundant(k) = false;
				if res.exitflag == MPTOPTIONS.OK
					if -lp.f*res.xopt-b(k)>MPTOPTIONS.abs_tol,
						% non-redundant
						nonredundant(k) = true;
					else
						% redundant
						removerow=[removerow k];
					end
				elseif res.exitflag == MPTOPTIONS.UNBOUNDED
					% non-redundant
					nonredundant(k)= true;
				end
			end
			k=k+1;
		end
		
		iredundant_rows(removerow)=[];
		obj.H_int(removerow,:) = [];
		
	else
		
		
		% for remaining rows test each inequality by solving LP and
		% eliminate rows sequentially
		Ae = obj.He_int(:, 1:end-1);
		be = obj.He_int(:, end);
		
		
		% for remaining rows test each inequality by solving LP and
		% eliminate rows sequentially - this allows better detection of redundant
		% ineqalities as keeping the rows in the next iteration
		
		i=1;
		while ( i<=size(obj.H_int,1) )
			
			A = obj.H_int(:, 1:end-1);
			b = obj.H_int(:, end);
			
			% Solve primal LP
			%     min H(i,1:end-1)*x
			%   s.t.:  H(~i)*[x;-1] <= 0
			%          H(i,1:end-1)*x <= H(i,end)+1
			
			
			H = obj.H_int;
			cost = H(i,1:end-1);
			offset = H(i,end);
			H(i,end) = H(i,end) + 1;
			
			% Setup data for fast call to LP solver
			lpn.f = -cost(:);
			
			lpn.A = H(:,1:end-1);
			lpn.b = H(:,end);
			lpn.Ae = Ae; lpn.be = be;
			
			% try to solve without lb/ub first
			res = mpt_solve(lpn);
			
			% if not feasible add bounds
			if res.exitflag ~=MPTOPTIONS.OK
				lpn.lb = -Inf(obj.Dim,1);
				lpn.ub = Inf(obj.Dim,1);
				res = mpt_solve(lpn);
			end
			
			if res.exitflag ==MPTOPTIONS.OK
				if (-res.obj < offset + MPTOPTIONS.abs_tol)
					% appears to be redundant
					% remove this inequality only if the polyhedron remains feasible
					obj.H_int(i,:) = [];
					iredundant_rows(i) = [];
				else
					% update row counter
					i = i+1;
				end
			else
				% this error should occur only if polyhedron is not feasible
				% (i.e. empty polyhedron-> irrendundant are detected at the beginning)
				% Try to resolve using dual lex LP
				
				% formulate LEX-dual LP
				%  lexmin [[b;be] I]'*y
				% s.t:  [A; Ae]'*y = 0
				%    sum(y) = 1 (all except i-that and equalities Ae corresponds to testing facet )
				%      y_j >= 0 (all y except i and Ae are nonnegative )
				%      y_i, y_Ae are free
				
				%          % this problem formulation was sometimes unbounded
				%          % size of y
				%          ny = size(A,1) + size(Ae,1);
				%          % matrix sum
				%          sy = [ones(1,size(A,1)) zeros(1,size(Ae,1))];
				%          sy(i) = 0;
				%
				%          dp.Ae = [A' Ae'; sy];
				%          dp.be = [zeros(size(dp.Ae,1)-1,1); 1]; %the last row is the sum(y)=1
				%          dp.A = [-eye(size(A,1)) zeros(size(A,1),size(Ae,1))];
				%          dp.A(i,:) = [];
				%          dp.b = zeros(size(dp.A,1),1);
				%          f = [[b;be] eye(ny)];
				
				
				% formulation of lex dual LP with included bound on the objective function
				[n1,n2] = size(A);
				ne = size(Ae,1);
				% all except index i
				ic = setdiff(1:n1,i);
				if isempty(Ae)
					d.Ae = [zeros(1,n2) 1;
						A(ic,:) ones(n1-1,1);
						A(i,:) 0]';
				else
					d.Ae = [zeros(1,n2) 1;
						A(ic,:) ones(n1-1,1);
						A(i,:) 0;
						Ae zeros(ne,1)]';
				end
				d.be = [zeros(n2,1); 1];
				d.A = -eye(n1,n1+ne+1);
				d.b = zeros(n1,1);
				
				% the maximum positive value on t is on the first position
				ff = [[MPTOPTIONS.infbound; b(ic); b(i); be] eye(n1+ne+1)];
				for j=1:size(ff,2)
					%             dp.f = f(:,j); % pick j-th column of f
					d.f = ff(:,j);
					rd = mpt_solve(d);
					%             rp = mpt_solve(dp); % solution without upper bound
					if rd.exitflag ==MPTOPTIONS.OK
						if rd.obj > MPTOPTIONS.lex_tol
							% LP is lex-positive, facet is not redundant
							i = i+1;
							break;
						elseif (rd.obj < -MPTOPTIONS.lex_tol)
							% LP is not lex-positive, facet equation is redundant,
							% remove it
							obj.H_int(i,:) = []; % original system
							iredundant_rows(i) = [];
							break;
						elseif j>=size(ff,2)
							% all possible perturbations failed: this means that
							% facet may be not redundant
							
							% try again with zero tolerance
							for jn=1:size(ff,2)
								% dp.f = f(:,j); % pick j-th column of f
								d.f = ff(:,jn);
								rd = mpt_solve(d);
								if rd.exitflag ==MPTOPTIONS.OK
									if rd.obj > 0
										% LP is lex-positive, facet is not redundant
										i = i+1;
										break;
									elseif (rd.obj < 0)
										% LP is not lex-positive, facet equation is redundant,
										% remove it
										obj.H_int(i,:) = []; % original system
										iredundant_rows(i) = [];
										break;
									elseif j>=size(ff,2)
										% all possible perturbations failed: this means that
										% facet may be not redundant
										i = i+1;
										break;
										%error('Lexicographic approach to redundancy elimination failed.');
									end
								else
									%error('Error solving dual-lex redundancy LP.');
									% if error appears, let it as nonredundant facet
									i = i+1;
									break;
								end
							end
							%error('Lexicographic approach to redundancy elimination failed.');
						end
					else
						%error('Error solving dual-lex redundancy LP.');
						% if error appears, let it as nonredundant facet
						i = i+1;
						break;
					end
				end
				
			end
			
			
		end
	end
	
	obj.irredundantHRep = true;
	% unset obj.optMat since the H-representation might have changed
	obj.optMat = [];
	
end

if isempty(obj.H_int) && isempty(obj.He_int)
	% all rows were removed as redundant to produce R^n
	obj.H_int = [zeros(1, obj.Dim), 1];
	% keep at least one of the original constraints for consistency
	iredundant_rows = 1;
end
sol.H  = obj.H_int;
sol.He = obj.He_int;
sol.I = true(nold,1);
sol.I(iredundant_rows) = false;

end

