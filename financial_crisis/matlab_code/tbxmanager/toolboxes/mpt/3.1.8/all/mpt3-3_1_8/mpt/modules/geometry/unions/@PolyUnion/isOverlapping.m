function ts = isOverlapping(obj)
%
%  ISOVERLAPPING: Test if the union of polyhedra contains overlaps. 
%  =================================================================
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
%     Return true if the union U of polyhedra contains overlaps and false
%  otherwise. Once this method has been called, the information about the overlaps
%  can be retrieved from U.Internal.Overlaps property. This function considers
%  following two cases to detect overlaps: 
%    
%     1. If two full-dimensional polyhedra overlap, then the intersection of these
%     polyhedra must be full-dimensional. 
%     2. If low-dimensional and full-dimensional polyhedra overlap, then the
%     intersection of these polyhedra must not be empty. 
%     Note that this function is computationally demanding and is suitable for
%  unions with small number of polyhedra.
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
%          ts True if union of polyhedra has overlaps  
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
%     isConvex,  isConnected,  isFullDim,  isBounded
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
        ts(i) = obj(i).isOverlapping;
    end
    return;
end

% empty obj
if obj.Num <= 1
	% no regions = no overlaps
	% single region = no overlaps
    ts = false;
    return;
end

if isempty(obj.Internal.Overlaps)
    ts = false;
    % get all combinations without repetitions
    % the function "combnk" comes from statistical toolbox
    %c = combnk(1:obj.Num,2);
    c = nchoosek(1:obj.Num,2);
    
    % run a test consecutively for each combination of polyhedra in the set    
    for i=1:size(c,1)
        if MPTOPTIONS.verbose >= 1
            fprintf('Regions: %d-%d \n',c(i,1),c(i,2));
        end

        IC = intersect(obj.Set(c(i,1)),obj.Set(c(i,2)));
        if ~obj.Set(c(i,1)).isFullDim || ~obj.Set(c(i,2)).isFullDim
            if ~IC.isEmptySet
                ts = true;
                break;
            end
        else
            if IC.isFullDim
                ts = true;
                break;
            end
        end
    end
    obj.Internal.Overlaps = ts;
else
    ts = obj.Internal.Overlaps;
end



end
