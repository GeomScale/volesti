function tf = isEmptySet(P)
%
%  ISEMPTYSET: Test if a polyhedron has a non-empty interior. 
%  ===========================================================
%  
%  
%  SYNTAX
%  ------
%     
%      tf = P.isEmptySet
%      tf = isEmptySet(P)
%    
%  
%  DESCRIPTION
%  -----------
%     Return true if the polyhedron P  is empty or false otherwise. A polyhedron is
%  considered as empty if the diameter of the inscribed ball is less than
%  region_tol or if there exist no feasible point inside the set.
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
%          tf True if the polyhedron P  is empty and   
%             false otherwise.                         
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
%     isEmptySet,  isBounded,  isFullDim
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

% Deal with arrays of polyhedra
tf = false(size(P));
for i=1:length(P)
    % compute interior point only if P(i).Internal.Empty property has not been set
	if P(i).Dim==0 || (isempty(P(i).H_int) && isempty(P(i).He_int) && ...
			isempty(P(i).V_int) && isempty(P(i).R_int))
		P(i).Internal.Empty = true;
	elseif isempty(P(i).Internal.Empty)
        ip = P(i).interiorPoint;
        if ~isempty(ip.r) && 2*ip.r<MPTOPTIONS.region_tol
            % note that low-dimensional polyhedra have radius=0
            % we need to check how far is the interior point from the
            % inequalities
            if ip.isStrict
                % full-dim
                P(i).Internal.Empty = true;
            else          
                % low-dim
                if P(i).hasVRep
                    % calculate the maximum distance to vertices
                    hn = P(i).V-repmat(x',2,1);
                    d = sqrt(sum(hn.*hn,2));
                    dmax = max(d);
                else
                    % calculate the maximum violation of inequalities
                    hn = P(i).H*[ip.x;-1];
                    dmax = max(hn);                    
                end
                if dmax>MPTOPTIONS.abs_tol
                    % maximum allowable tolerance is violated, region is empty
                    P(i).Internal.Empty = true;
                else
                    % no constraints violated, region is not empty
                    P(i).Internal.Empty = false;

                end
			end
		elseif any(isnan(ip.x))
            % NaN in interior point means empty polyhedron
            P(i).Internal.Empty = true;
		else
			P(i).Internal.Empty = isempty(ip.x);
        end
    end
    tf(i) = P(i).Internal.Empty;
end
end
