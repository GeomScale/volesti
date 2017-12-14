function [isequal] = mpt_iscombequal(vector1,vector2)
%MPT_ISCOMBEQUAL Are two vectors combinatorially equal
%
%[isequal] = mpt_iscombequal(vector1,vector2)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Checks if two vectors are combinatorially equal
% e.g. [1 2 3] and [2 3 1] would return true, since all elements are identical
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% vector1, vector2   - vectors
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% isequal  - 1 if vectors are combinatorially equal, 0 otherwise
%

% Copyright is with the following author(s):
%
% (C) 2003 Mato Baotic, Automatic Control Laboratory, ETH Zurich,
%          baotic@control.ee.ethz.ch
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

isequal = (all(ismember(vector1,vector2))) & (all(ismember(vector2,vector1)));