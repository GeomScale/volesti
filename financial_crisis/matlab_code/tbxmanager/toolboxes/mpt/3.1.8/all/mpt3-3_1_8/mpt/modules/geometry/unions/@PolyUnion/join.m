function out = join(obj)
%
%  JOIN: Merges arrays of polyunions to one polyunion object. 
%  ===========================================================
%  
%  
%  SYNTAX
%  ------
%     
%      Un = U.join
%      Un = join(U)
%    
%  
%  DESCRIPTION
%  -----------
%     Merging of array of PolyUnion objecs in the same dimension to single
%  PolyUnion object. If the union has Bounded and FullDim properties set as true,
%  then all polyhedra must be bounded and full-dimensional.
%  
%  INPUT
%  -----
%     
%        
%          U Array of PolyUnion objects in the same   
%            dimension.                               
%            Class: PolyUnion                         
%              
%  
%  
%  OUTPUT
%  ------
%     
%        
%          Un Single PolyUnion object.                 
%             Class: PolyUnion                         
%               
%  
%  
%  SEE ALSO
%  --------
%     merge,  isFullDim,  isBounded
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
 
 
dm = cell(size(obj));
[dm{:}] = obj.Dim;
obj(cellfun('isempty',dm)) = [];
obj([dm{:}]==0);

no = numel(obj);
if no<=1
    % nothing to do
    out = obj;
else    
    dims = zeros(no, 1);
    for i=1:no
        dims(i) = obj(i).Dim;
    end
    
    if any(diff(dims))
        error('Only polyunions of the identical dimensions can be joined.');
    end

    % do not check convexity, overlaps, and connectivity
    % check boundedness
    bnd = cell(no, 1);
    for i=1:no
        bnd{i} = obj(i).Internal.Bounded;
    end
    ie = cellfun('isempty',bnd);
    if any(~ie)
        % one of unions has Bounded property set
        if any([bnd{:}]==1)
            % check for bounded sets
            b = obj.isBounded;
            if ~all(b)
                error('All unions must be bounded.')
            end
        end
        indexb = find(~ie);
        isb = bnd{indexb(1)};
    else
        isb = [];
    end
    
    % check dimensionality
    d = cell(no, 1);
    for i=1:no
        d{i} = obj(i).Internal.FullDim;
    end
    id = cellfun('isempty',d);
    if any(~id)
        % one of unions has FullDim property set
        if any([d{:}]==1)
            % check for fulldimensionality
            df = obj.isFullDim;
            if ~all(df)
                error('All unions must be full-dimensional.')
            end
        end
        indexfd = find(~id);
        isfd = d{indexfd(1)};
    else
        isfd = [];
    end
    
    out = PolyUnion('Set', obj(1).Set);
    if numel(obj)>1
        for i = 2:length(obj)
            out.add(obj(i).Set);
        end
    end
    
    % set internal properties
    out.Internal.Bounded = isb;
    out.Internal.FullDim = isfd;
    
end


end
