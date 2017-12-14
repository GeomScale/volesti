function [ts, iP, iQ] = isNeighbor(P,Q,fP,fQ)
%
%  ISNEIGHBOR: Test if a polyhedron touches another polyhedron along a given facet.
%  ================================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      ts = P.isNeighbor(Q)
%      ts = isNeighbor(P,Q)
%      [ts, iP, iQ] = isNeighbor(P,Q,fP,fQ)
%    
%  
%  DESCRIPTION
%  -----------
%     Return true if the polyhedron P  shares with the polyhedron Q  a part of a
%  common facet. Both polyhedrons must be in H-representation. If they are not, the
%  irredundant H-representation will be computed. The function tests if the
%  intersection of two polyhedra P  and S  in dimension d  is nonempty and is of
%  dimension d-1. If this holds, then the two polyhedra are touching themself along
%  one facet. For closer explanation, see the example below.
%  
%  INPUT
%  -----
%     
%        
%          P  Polyhedron in H-representation           
%             Class: Polyhedron                        
%          Q  Polyhedron in H-representation           
%             Class: Polyhedron                        
%          fP Index of a facet to test from polyhedron 
%             P.                                       
%             Class: double                            
%          fQ Index of a facet to test from polyhedron 
%             Q.                                       
%             Class: double                            
%               
%  
%  
%  OUTPUT
%  ------
%     
%        
%          ts Logical statement if the polyhedron P    
%             touches Q  along the facet.              
%             Class: logical                           
%             Allowed values:                          
%                                                      
%               true                                   
%               false                                  
%                                                      
%          iP Index of a facet from polyhedron P  that 
%             touches the polyhedron Q.                
%             Class: double                            
%          iQ Index of a facet from polyhedron Q  that 
%             touches the polyhedron P.                
%             Class: double                            
%               
%  
%  
%  SEE ALSO
%  --------
%     isAdjacent
%  

%  AUTHOR(s)
%  ---------
%     
%    
%   (c) 2010-2013  Martin Herceg: ETH Zurich
%   mailto:herceg@control.ee.ethz.ch 
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

narginchk(2, 4);

validate_polyhedron(Q);

if numel(P)>1 || numel(Q)>1
    error('P and Q must have the length of 1.');
end

if P.Dim ~= Q.Dim
    error('Polyhedra must be in the same dimension.');
end


% initial values to be returned if adjacency is not proven
ts = false;
iP = [];
iQ = [];

% empty P,Q -> return 0
if isEmptySet(P) || isEmptySet(Q)
    if MPTOPTIONS.verbose>=1 
        disp('P or Q is empty.');
    end
    return
end

T = P.intersect(Q);
if T.isEmptySet
	% non-ontersecting polyhedra cannot be adjacent
	return
end

% make P, Q irredundant
P.minHRep();
Q.minHRep();

% check indices of the facets if provided
if nargin>2
    if isempty(fP)
        fP = 1:size(P.H,1);
    else
        validate_indexset(fP);
        % convert logical index sets to numeric
        if islogical(fP) || any(fP==0)
            fP = find(fP);
        end
        % consider only unique indices
        fP = unique(fP);
        
        n1 = size(P.H,1);
        if numel(fP)>n1 || any(fP>n1)
            error('Facet index set for region P contains indices out of range.');
        end
    end
else
    fP = 1:size(P.H,1);
end
if nargin>3
    if isempty(fQ)
        fQ = 1:size(Q.H,1);
    else
        validate_indexset(fQ);
        % convert logical index sets to numeric
        if islogical(fQ) || any(fQ==0)
            fQ = find(fQ);
        end
        % consider only unique indices
        fQ = unique(fQ);
        n2 = size(Q.H,1);
        if numel(fQ)>n2 || any(fQ>n2)
            error('Facet index set for region Q contains indices out of range.');
        end

    end
else
    fQ = 1:size(Q.H,1);
end



% % test if two polyhedron are adjacent
% 
% H1 = [      0.8905      0.23614           10;
%     -0.055625     0.030184            0;
%     -0.21887     -0.06688   8.8818e-16;
%     0           -1           10;
%     1            0           10;
%     0            1           5;
%     0            1           20];
% 
% H2 = [     0.055625    -0.030184            0;
%     -0.053731    -0.010858   2.2204e-16;
%     0            1           10];
% 
% P1 = Polyhedron('H',H1);
% P2 = Polyhedron('H',H2);
% 
% % first test if i-th row of H1 is a facet (which obviously is except the
% % last row)
% %       max t
% % s.t.: A(ic,:)*th + t <= b(ic,:) , ic is all except i
% %       A(i,:)*th = b(i)
% %
% % if t>0 -> it is a facet
% 
% f = false(size(H1,1),1);
% fn =f;
% 
% for i=1:size(H1,1)
%     
%     Ae = H1(i,1:end-1);
%     be = H1(i,end);
%     ic = setdiff(1:size(H1,1),i);
%     A = H1(ic,1:end-1);
%     b = H1(ic,end);
%     
%     % formulate equality constrained LP
%     lp.Ae = [Ae 0];
%     lp.be = be;
%     lp.A = [A ones(length(ic),1)];
%     lp.b = b;
%     lp.f = zeros(1,3);
%     lp.f(end) = -1;
%     
%     r = mpt_solve(lp);
%     if -r.obj >0
%         f(i) = true;
%     end    
%     
%     % remove the equality and formulate LP
%     [~,p] = max(abs(Ae));
%     gmax = Ae(p);
%     
%     % a1*th1 + a2*th = be    
%     np    = true(2,1);
%     np(p) = false; 
% 
%     % obtain matrices Qnew, qnew
%     lpnew.A = A(:,np) -A(:,p)*Ae(np)'/gmax;
%     lpnew.b = A(:,p)*be/gmax + b;
%     lpnew.f = [0 -1];
%     
%     rn = mpt_solve(lp);
%     if -r.obj >0
%         fn(i) = true;
%     end    
% 
% 
% end

% now check if facet i is adjacent to j
H1 = P.H;
H2 = Q.H;
P_A = H1(:, 1:end-1); P_b = H1(:, end);
Q_A = H2(:, 1:end-1); Q_b = H2(:, end);
P_Ae = P.Ae; P_be = P.be;
Q_Ae = Q.Ae; Q_be = Q.be;

nic = length(P_b)-1;
njc = length(Q_b)-1;

for i=fP
	% formulate equality constrained LP
	% constraints from P1
	Ae1 = P_A(i, :);
	be1 = P_b(i);
	A1 = P_A; A1(i, :) = [];
	b1 = P_b; b1(i) = [];

	for j=fQ
        
        % constraints from P2
		Ae2 = Q_A(j, :);
		be2 = Q_b(j);
		A2 = Q_A; A2(j, :) = [];
		b2 = Q_b; b2(j) = [];

        % formulate problem
        % max t
        %  s.t.: A1*x+t <= b1
        %        A2*x+t <= b2
        %        Ae1_f*x = be1_f  % facet constraint
        %        Ae2_f*x = be2_f  % facet constraint
        
        s.A = [A1 ones(nic,1); A2 ones(njc,1)];
        s.b = [b1; b2];
        % include equalities from P, Q to test low-dim polyhedra
        s.Ae = [Ae1 0; Ae2 0; P_Ae zeros(size(P_Ae,1),1); Q_Ae zeros(size(Q_Ae,1),1)];
        s.be = [be1; be2; P_be; Q_be];
        s.f = zeros(1,P.Dim+1);
        s.f(end) = -1;
		
		% solve the LP quickly
        s.lb = []; s.ub = []; 
		s.quicklp = true;
		
        rs = mpt_solve(s);
        
        if rs.exitflag==MPTOPTIONS.OK
            % if t>0, adjacency has been proven, break and return
            if -rs.obj > MPTOPTIONS.abs_tol
                iP = i;
                iQ = j;
                ts = true;
                break
            end
        end
        
        %         % remove the equality and formulate LP
        %         [~,p] = max(abs(Ae1));
        %         gmax = Ae1(p);
        %
        %         % a1*th1 + a2*th = be
        %         np    = true(2,1);
        %         np(p) = false;
        %
        %         % obtain matrices Qnew, qnew
        %         A1new = A1(:,np) -A1(:,p)*Ae1(np)'/gmax;
        %         b1new = A1(:,p)*be1/gmax + b1;
        %
        %         if norm(Ae2(p),1)<MPTOPTIONS.abs_tol
        %             % if fp is zero, thp can be anything. We choose thp=0
        %             % which means eliminating corresponding column.
        %             A2new = A2;
        %             A2new(:,p) = [];
        %             b2new = b2;
        %         else
        %             A2new = A2(:,np)-A2(:,p)*Ae2(np)'/Ae2(p);
        %             b2new = A2(:,p)*be2/Ae2(p) + b2;
        %         end
        %
        %         sn.A = [A1new ones(length(ic),1); A2new ones(length(jc),1)];
        %         sn.b = [b1new; b2new];
        %         sn.f = [0 -1];
        %
        %         rn = mpt_solve(sn);
        %         if -rn.obj > MPTOPTIONS.abs_tol
        %             adjn{i} = [adjn{i}, j];
        %         end
        
        
        
    end
    if ts
        break
    end
end


% the polyhedra can be neighbors along the facet described as equality
% (e.g. if one of them is affine set)
% In the next part we merge equalities with inequalities and test adjacency
% for the merged constraints. 

if ~ts && ( ~isFullDim(P) || ~isFullDim(Q) )
    % lower-dimensional polyhedra can have double-sided inequalities 
    % mark indices that are equalities
    if ~isempty(P.H_int)        
        [PA,Pb,PAe,Pbe,indeqP]=mpt_ineq2eq(P.A,P.b);
        indP = setdiff(1:size(P.H),indeqP(:));
    else
        PA = []; Pb = []; PAe = []; Pbe = []; indeqP = []; indP = [];
    end
    if ~isempty(Q.H_int)
        [QA,Qb,QAe,Qbe,indeqQ]=mpt_ineq2eq(Q.A,Q.b);
        indQ = setdiff(1:size(Q.H),indeqQ(:));
    else
        QA = []; Qb = []; QAe = []; Qbe = []; indeqQ = []; indQ = [];
    end
        
    % now check if facet i is adjacent to j
    H1 = [PA Pb; PAe Pbe; P.He];
    H2 = [QA Qb; QAe Qbe; Q.He];
    
    for i=1:size(H1,1)
        for j=size(H2,1)
            
            % formulate equality constrained LP
            % constraints from P1
            Ae1 =H1(i,1:end-1);
            be1 = H1(i,end);
            ic = setdiff(1:size(H1,1),i);
            A1 = H1(ic,1:end-1);
            b1 = H1(ic,end);
            
            % constraints from P2
            Ae2 = H2(j,1:end-1);
            be2 = H2(j,end);
            jc = setdiff(1:size(H2,1),j);
            A2 = H2(jc,1:end-1);
            b2 = H2(jc,end);
            
            % formulate problem
            % max t
            %  s.t.: A1*x+t <= b1
            %        A2*x+t <= b2
            %        Ae1_f*x = be1_f  % facet constraint
            %        Ae2_f*x = be2_f  % facet constraint
            
            s.A = [A1 ones(length(ic),1); A2 ones(length(jc),1)];
            s.b = [b1; b2];
            % include equalities from P, Q to test low-dim polyhedra
            s.Ae = [Ae1 0; Ae2 0];
            s.be = [be1; be2];
            s.f = zeros(1,P.Dim+1);
            s.f(end) = -1;
            
            rs = mpt_solve(s);
            
            if rs.exitflag==MPTOPTIONS.OK
                % if t>0, adjacency has been proven, break and return
                if -rs.obj > MPTOPTIONS.abs_tol
                    if i <= numel(indP)
                        iP = indP(i);
                    elseif (i > numel(indP)) && (i<=numel(indP)+size(indeqP,1))
                        % if the adjacent facet belongs do double sided
                        % equalities, return the corresponding couple
                        iP = indeqP( i-numel(indP), : );
                    else
                        iP = [];
                    end
                    if j <= numel(indQ)
                        iQ = indQ(j);
                    elseif (j > numel(indQ)) && (j<=numel(indQ)+size(indeqQ,1))
                        % if the adjacent facet belongs do double sided
                        % equalities, return the corresponding couple
                        iQ = indeqQ( j-numel(indQ), : );
                    else
                        iQ = [];
                    end
                    
                    ts = true;
                    
                    break
                end
            end
            
        end
    end
end


% if the above tests indicate that the polyhedra share a common facet, we
% need to check if their intersection (possibly low-dimensional polyhedron)
% is not empty 
% see test_polyhedron_isAdjacent_13_pass, 14
if ts
    if ~P.isFullDim || ~Q.isFullDim
        % one of the polyhedra is low-dimensional - need to check if T is
        % not empty
        if T.isEmptySet
            % polyhedra are not adjacent
            ts = false;
        end
    else
        % both polyhedra are full-dimensional - check if T is not Full-Dim
        if T.isFullDim || T.isEmptySet
            ts = false;
        end
    end
end


end
 
