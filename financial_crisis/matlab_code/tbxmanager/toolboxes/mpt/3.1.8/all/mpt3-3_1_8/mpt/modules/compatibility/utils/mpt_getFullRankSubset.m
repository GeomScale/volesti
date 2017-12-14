function [Gret, keptrows]=mpt_getFullRankSubset(Gt,AllSubsets)
%MPT_GETFULLRANKSUBSET Removes rows from matrix Gt until it has full row rank
%
% [Gt, keptrows]=mpt_getFullRankSubset(Gt)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% It returns all (or just one) full row rank matrices and the rows which were 
% kept from the original matrix.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% Gt            - matrix
% AllSubsets    - Return all subsets or just one; set to 1 or 0 respectively;
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% Gt        - matrix having full row rank (is cell if all subsets are returned)
% keptrows  - indices of the rows that were NOT eliminated
%

% Copyright is with the following author(s):
%
% (C) 2003 Pascal Grieder, Automatic Control Laboratory, ETH Zurich,
%          grieder@control.ee.ethz.ch

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

narginchk(1, 2);

if(nargin<2 | isempty(AllSubsets))
    AllSubsets=0;
end
rows=size(Gt,1);
if(rank(Gt)==rows)
    %matrix has full row rank 
    keptrows=1:rows;
    Gret=Gt;
    return
end
ctr=0;
if(rank(Gt)==rows-1)
    %one redundant row
    for i=1:rows
        Gtemp=Gt;
        Gtemp(i,:)=[];
        if(rank(Gtemp)<=rows)
            ctr=ctr+1;
            %matrix has full row rank 
            Gret{ctr}=Gtemp;
            keptrows{ctr}=1:rows;
            keptrows{ctr}(i)=[];
        end
    end
    return
end


if(rank(Gt)==rows-2)
    %two redundant rows
    for i=1:rows
        for j=(i+1):rows
            Gtemp=Gt;
            Gtemp([i j],:)=[];
            if(rank(Gtemp)<=rows)
                ctr=ctr+1;
                %matrix has full row rank 
                Gret{ctr}=Gtemp;
                keptrows{ctr}=1:rows;
                keptrows{ctr}([i j])=[];
            end
        end
    end
    return
end
if(rank(Gt)==rows-3)
    %three redundant rows
    for i=1:rows
        for j=(i+1):rows
            for k=(j+1):rows
                Gtemp=Gt;
                Gtemp([i j k],:)=[];
                if(rank(Gtemp)<=rows)
                     ctr=ctr+1;
                    %matrix has full row rank 
                    Gret{ctr}=Gtemp;
                    keptrows{ctr}=1:rows;
                    keptrows{ctr}([i j k])=[];
                end
            end
        end
    end
    return
end
if(rank(Gt)==rows-4)
    %four redundant rows
    for i=1:rows
        for j=(i+1):rows
            for k=(j+1):rows
                for m=(k+1):rows,
                    Gtemp=Gt;
                    Gtemp([i j k m],:)=[];
                    if(rank(Gtemp)<=rows)
                        ctr=ctr+1;
                        %matrix has full row rank 
                        Gret{ctr}=Gtemp;
                        keptrows{ctr}=1:rows;
                        keptrows{ctr}([i j k m])=[];
                    end
                end
            end
        end
    end
    return
end
Gret = {};
ctr = 0;
allcombs = nchoosek(1:rows, rows-rank(Gt));
for ii = 1:size(allcombs,1),
    Gtemp = Gt;
    Gtemp(allcombs(ii, :), :) = [];
    if rank(Gtemp)<=rows,
        ctr = ctr+1;
        Gret{ctr} = Gtemp;
        keptrows{ctr} = 1:rows;
        keptrows{ctr}(allcombs(ii,:)) = [];
        if ctr>2,
            return
        end
    end
end
return
        
disp('Yes, it''s true: I''m lazy. This is why this code was written in such a silly way') 
disp('Feel free to write a proper function to perform this task and send it to me.');
fprintf('\n\n')
disp('Please contact grieder@control.ee.ethz.ch if this occurs')
error('mpt_getFullRankSubset: Primal Degeneracy: Code for four redundant equalities not written yet.')
