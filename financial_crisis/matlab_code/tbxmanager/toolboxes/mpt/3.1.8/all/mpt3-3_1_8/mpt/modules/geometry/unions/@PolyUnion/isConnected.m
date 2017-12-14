function ts = isConnected(obj)
%
%  ISCONNECTED: Test if the union of polyhedra form a connected union. 
%  ====================================================================
%  
%  
%  SYNTAX
%  ------
%     
%      ts = U.isConnected
%      ts = isConnected(U)
%    
%  
%  DESCRIPTION
%  -----------
%     Return true if the union U of polyhedra is connected and false otherwise.
%  Once this method has been called, the information about the connectivity can be
%  retrieved from U.Internal.Connected property. Note tha if the union U is convex,
%  it implies that the union is connected. Note that this function is
%  computationally demanding and is suitable for unions with small number of
%  polyhedra.
%  
%  INPUT
%  -----
%     
%        
%          U Union of polyhedra in the same           
%            dimension.                               
%            Class: PolyUnion                         
%              
%  
%  
%  OUTPUT
%  ------
%     
%        
%          ts True if union of polyhedra is connected  
%             and false otherwise.                     
%             Class: logical                           
%             Allowed values:                          
%                                                      
%               true                                   
%               false                                  
%                                                      
%               
%  
%  
%  SEE ALSO
%  --------
%     isConvex,  isOverlapping,  isFullDim,  isBounded
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

% deal with arrays
if numel(obj)>1
    ts = -ones(size(obj));
    for i=1:numel(obj)
        ts(i) = obj(i).isConnected;
    end
    return;
end

% empty obj
if obj.Num==0
    ts = false;
    return;
elseif obj.Num==1
	% single region is connected
	ts = true;
	return
end

% check connectivity
if isempty(obj.Internal.Connected)        
    % if there is info about the convexity, use this info to determine
    % connectivity
    if ~isempty(obj.Internal.Convex)        
        if obj.Internal.Convex
            % if the union is convex, then it is connected as well
            % otherwise we have to check for connectivity
            ts = true;
            obj.Internal.Connected = ts;
            return
        end
    end
    
    % Each polyhedron inside the array must contain a common point
    % with any of the remaining sets.
    % The easiest would be to check connectivity between neighbors,
    % but since we don't know what are the neighbors, we run a
    % consecutive test for each polyhedron in the array reducing of
    % those polyhedra that have been checked.
    
    % get all combinations with
    c = cbn(1:obj.Num);
    
    % run a test consecutively for each combination of polyhedra in the set
    t = false(1,obj.Num);
    while size(c,1)>0
        if MPTOPTIONS.verbose >= 1
            fprintf('Regions: %d-%d \n',c(1,1),c(1,2));
        end
        
        % take always the first combination
        reg = c(1,:);
        IC = intersect(obj.Set(reg(1)),obj.Set(reg(2)));
        if ~IC.isEmptySet
            % region c(i,1)=reg(1)  is connected with c(i,2)=reg(2)
            t(reg(1)) = true;
            t(reg(2)) = true;
            % kick out region c(i,1) and c(i,2) from the list
            if ~isempty(c)
                c(c(:,1)==reg(1),:) = [];
            end
            if ~isempty(c)
                c(c(:,1)==reg(2),:) = [];
            end
        else
            % this combination has been checked, remove it from the list
            c(1,:)=[];
        end
        
    end
    % if all are connected, return true
    ts = all(t);
    obj.Internal.Connected = ts;
    
else
   ts = obj.Internal.Connected; 
end



end

function y = cbn(v)
%
% create two combinations of sets to be checked
% v - vector of indices
%

N = numel(v);

y = zeros(N*(N-1),2);
k=1;
for i=1:N
    for j=setdiff(1:N,i)
        y(k,:) = [v(i),v(j)];
        k = k+1;
    end
end



end
