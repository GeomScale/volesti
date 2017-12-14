function closest = closestRegion(P,x0)
%CLOSESTREGION Returns index of a region which is closest to the point x0
%
% closest = closestRegion(P, x0)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Returns index of a region which is closest to the point x0
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P                 - Polytopes
% x0                - given point, if not specified: Options.clickx0=1
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------

% closest   - if the point does not belong to P, index of the polytope which is
%             closest to the point x is returned
%
% see also POLYHEDRON/ISIN
%

% Copyright is with the following author(s):
%
%     2010, Revised by Martin Herceg, Automatic Control Lab, ETH Zurich
% (C) 2006 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch

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

lenP = length(P);
if lenP < 2,
    closest = 1;
else
    d = zeros(1,lenP);
    for i=1:lenP,
        ret=P(i).distance(x0);
        d(i) = ret.dist;
    end
    [~,closest] = min(d);
    closest = closest(1);
end
