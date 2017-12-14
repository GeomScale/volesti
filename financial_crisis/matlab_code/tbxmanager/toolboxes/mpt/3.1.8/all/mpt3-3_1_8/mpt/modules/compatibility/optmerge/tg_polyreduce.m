function [At,Bt,isemptypoly,keptrows]=tg_polyreduce(At,Bt,tolerance,lpsolver,checkemptyrow,candidates,useCdd)
%
% purpose: 
% Given a polyhedron Ax<=B, return an equivalent polyhedron At x<=Bt by
% eliminating redundant constraints.
% Consider only hyperplanes for removal that with indices given in the set 
% 'candidates'.
% If cddUse==1, use cdd('find_interior') to determine if the polyhedron is 
% empty (instead of polyreduce) and use cdd('reduce_h') to determine the
% non-redundant hyperplanes (instead of a loop with lpsolve)
%
% At,Bt        = reduced polyhedron
% isemptypoly  = 1 if (A,B) is empty
% keptrows     = rows of (A,B) kept in (At,Bt)
%

% Copyright is with the following author(s):
%
% (C) 2004 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2003-2004 Tobias Geyer, Automatic Control Laboratory, ETH Zurich,
%          geyer@control.ee.ethz.ch

% ---------------------------------------------------------------------------
% Legal note:
%          This program is free software; you can redistribute it and/or
%          modify it under the terms of the GNU General Public
%          License as published by the Free Software Foundation; either
%          version 2.1 of the License, or (at your option) any later version.
%
%          This program is distributed in the hope that it will be useful,
%          but WITHOUT ANY WARRANTY; without even the implied warranty of
%          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%          General Public License for more details.
% 
%          You should have received a copy of the GNU General Public
%          License along with this library; if not, write to the 
%          Free Software Foundation, Inc., 
%          59 Temple Place, Suite 330, 
%          Boston, MA  02111-1307  USA
%
% ---------------------------------------------------------------------------

global mptOptions

n=size(At,2);
if length(Bt)~=size(At,1),
   error('The A and B matrices must have the same number of rows.')
end

if nargin<7
    % 0 - use reduction through solveLP
    % 1 - use cddmex
    % 2 - use reduce()
    useCdd=1;
end;
if nargin<6
    candidates = 1:length(Bt);
end;
if nargin<5,
    checkemptyrow=0;    % do not check for removal of empty rows
end;
if nargin<4,
    lpsolver=mptOptions.lpsolver;         % use NAG toolbox
    % lpsolver=1;       % Matlab's linprog
end;
if nargin<3,
    tolerance=mptOptions.abs_tol;     % non-negative tolerance to decide if a constraint is redundant
end

isemptypoly=0;
keptrows=(1:length(Bt))';
Aorig = At;
Borig = Bt;


% clear empty rows
% -----------------------------------------------------------------------
if checkemptyrow,

    emptyrows = max(abs(At),[],2) < tolerance;
    if any(emptyrows),
        if any(Bt(emptyrows)<-tolerance),
            isemptypoly=1;
            return;
        else
            % this part is problably not implemented correctly.
            error('tg_polyreduce: debug this part')
            keptrows(emptyrows)=[];
            At(emptyrows,:)=[];
            Bt(emptyrows)=[];
        end
    end
end




% check if polyhedron is empty
% -----------------------------------------------------------------------

if useCdd==1
    
    % use cdd
    Pfull.A = At;
    Pfull.B = Bt;
    Pred = cddmex('find_interior', Pfull);
    R = Pred.objlp;
    if Pred.how~=1 | R<tolerance
        isemptypoly = 1;
        return;
    end
    
elseif useCdd==0
    
    % use polyinnerball
    if ~isemptypoly,
        [xopt,R]=polyinnerball(At,Bt);
        if R<tolerance
            isemptypoly = 1;
            return;
        end
    end
    
end


% remove redundant hyperplanes
% -----------------------------------------------------------------------

if useCdd==1
    
    % use cdd
    try
        [Pred, redrows] = cddmex('reduce_h', Pfull);
    catch
        % CDD failed, try an alternative method (reduce.m)
        Punred = polytope(Aorig, Borig, 0, 2);
        [Pred, keptrows] = reduce(Punred);
        isemptypoly = ~isfulldim(Pred);
        [At, Bt] = double(Pred);
        return
    end
    
    % find the indices of the rows that we want to remove (rows that are
    % both in candidates and redrows)
    remrows = [];
    for i=1:length(redrows)
        if ~isempty(find(candidates==redrows(i)))
            % redrows(i) is in the candidate set
            remrows(end+1) = redrows(i);
        end;
    end;
    
    % remove rows
    keptrows(remrows) = [];
    At = Pred.A;
    Bt = Pred.B;
    
    return
    
elseif useCdd==2
    poly = polytope(At, Bt, 0, 2);
    if ~isfulldim(poly),
        isemptypoly = 1;
        return
    end
    [pred, redrows] = reduce(poly);
    [At, Bt] = double(pred);
    remrows = [];
    for i=1:length(redrows)
        if ~isempty(find(candidates==redrows(i)))
            % redrows(i) is in the candidate set
            remrows(end+1) = redrows(i);
        end;
    end;
    
    % remove rows
    keptrows(remrows) = [];
    
else
    % use loop with lpsolve
    while length(candidates) > 0
        j = candidates(1);
        
        f = At(j,:)';
        
        % is the j-th hyperplane redundant?
        % get indices of all the other hyperplanes 
        ii =[1:j-1, j+1:length(Bt)];
        
        if ~isempty(ii),
            [xopt,dummy,status]=mpt_solveLPi(-f',At(ii,:),Bt(ii),[],[],[],lpsolver);
            obj=f'*xopt-Bt(j);
        else
            status='unbounded';
        end
        
        if strcmp(status,'unbounded') | obj>tolerance,
            % j-th hyperplane is non-redundant
            candidates(1) = [];
        else
            % j-th hyperplane is redundant: remove
            At(j,:)=[];
            Bt(j)=[];
            keptrows(j)=[];
            
            % adapt candidates accordingly
            %candidates = find(candidates>j) - 1;
            candidates(1) = [];
            candidates = candidates - 1; % because we have deleted a row
        end
    end
    
    return
    
end
