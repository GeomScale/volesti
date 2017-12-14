function Options = mpt_defaultOptions(Options, varargin)
%MPT_DEFAULTOPTIONS Sets default values of undefined fields of Options structure
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Sets default values of fields which are not defined in the Options structure.
%
% Options = mpt_defaultOptions(Options, 'f1', v1, 'f2', v2)
%
%   will set Options.f1=v1 and Options.f2=v2 if Options.f1 and/or Options.f2 do
%   not exist in Options.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
%

% Copyright is with the following author(s):
%
% (C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich
%     kvasnica@control.ee.ethz.ch

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

narginchk(2, Inf);

nargs = length(varargin);
if rem(nargs, 2)~=0,
    error('mpt_defaultOptions: wrong number of input arguments!');
end
for ia = 1:2:nargs,
    field = varargin{ia};
    default = varargin{ia+1};
    if ~isfield(Options, field),
        Options = setfield(Options, field, default);
    end
end
